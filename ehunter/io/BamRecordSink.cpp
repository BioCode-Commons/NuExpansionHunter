//
// NuExpansionHunter
// Copyright 2026 NuExpansionHunter contributors
//
// This file is part of NuExpansionHunter, a modified version of
// Illumina ExpansionHunter v5.0.0 (Copyright 2016-2019 Illumina, Inc.).
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#include "io/BamRecordSink.hh"

#include <condition_variable>
#include <atomic>
#include <deque>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>

namespace ehunter
{

constexpr std::size_t BamRecordSink::defaultBudgetBytes;

struct BamRecordSink::State
{
    struct Charge
    {
        explicit Charge(size_t bytes) : remaining(bytes), requested(bytes) {}
        size_t remaining; // Protected by State::mutex.
        std::atomic<size_t> requested;
        std::atomic<bool> active {false};
    };

    // Called under mutex after admission is published. A callback may have
    // finished synchronously inside prepareAccounted, before the charge exists.
    void applyCapacity(Charge& charge)
    {
        const size_t requested = charge.requested.load();
        if (charge.active.load() && requested < charge.remaining)
        {
            resident -= charge.remaining - requested;
            charge.remaining = requested;
        }
    }

    struct Envelope
    {
        Records records;
        Bytes data;
        std::future<Bytes> prepared;
        size_t count;
        std::size_t bytes;
        std::shared_ptr<Charge> compressedCharge;
    };

    State(WriteRecord callback, std::size_t budget, WriteBytes encodedWriter, PrepareAccounted accounted)
        : writeRecord(std::move(callback))
        , writeBytes(std::move(encodedWriter))
        , prepareAccounted(std::move(accounted))
        , maxBytes(budget)
    {
        if (!writeRecord)
            throw std::invalid_argument("BAM record sink requires a writer");
        consumer = std::thread([this] { consume(); });
    }

    void consume() noexcept
    {
        for (;;)
        {
            Envelope current;
            {
                std::unique_lock<std::mutex> lock(mutex);
                ready.wait(lock, [&] { return closing || !batches.empty(); });
                if (batches.empty())
                    return;
                current = std::move(batches.front());
                batches.pop_front();
            }

            try
            {
                if (current.prepared.valid())
                    current.data = current.prepared.get();
                if (!current.data.empty())
                    writeBytes(current.data.data(), current.data.size());
                for (auto& record : current.records)
                {
                    writeRecord(record.get());
                    record.reset();
                }
            }
            catch (...)
            {
                // Cancel and release accepted work without allocating another
                // container on the failure path. Waiting producers then wake
                // and release their own unsubmitted records.
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    error = std::current_exception();
                    resident = 0;
                }
                credits.notify_all();
                ready.notify_all();
                Bytes().swap(current.data);
                current.records = Records();
                // Dropping a future can join the encoder. An encoder returns its
                // memory credits under mutex, so destroy cancelled futures outside it.
                for (;;)
                {
                    Envelope cancelled;
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        if (batches.empty()) break;
                        cancelled = std::move(batches.front());
                        batches.pop_front();
                    }
                }
                return;
            }

            Bytes().swap(current.data);
            current.records = Records();
            bool wake;
            {
                std::lock_guard<std::mutex> lock(mutex);
                if (current.compressedCharge) current.compressedCharge->active.store(false);
                resident -= current.compressedCharge ? current.compressedCharge->remaining : current.bytes;
                // Wake waiting producers once half the budget is free (an empty
                // sink included), not after every batch: with many producers of
                // single-record batches, a wake per batch keeps them contending
                // for the mutex. Producers wait only while batches are resident,
                // so the drain always reaches this point.
                wake = waiting && resident <= maxBytes / 2;
            }
            if (wake) credits.notify_all();
        }
    }

    WriteRecord writeRecord;
    WriteBytes writeBytes;
    PrepareAccounted prepareAccounted;
    const std::size_t maxBytes;
    mutable std::mutex mutex;
    std::mutex finishing;
    std::condition_variable ready, credits;
    std::deque<Envelope> batches;
    std::size_t resident = 0;
    std::size_t waiting = 0; // Producers waiting for credits.
    bool closing = false; // Explicit finish only; error separately cancels admission.
    std::exception_ptr error;
    std::thread consumer;
};

BamRecordSink::BamRecordSink(WriteRecord writeRecord, std::size_t maxBytes, WriteBytes writeBytes, PrepareAccounted prepareAccounted)
    : state_(new State(std::move(writeRecord), maxBytes, std::move(writeBytes), std::move(prepareAccounted)))
{
}

BamRecordSink::~BamRecordSink()
{
    try
    {
        finish();
    }
    catch (...)
    {
    }
}

std::size_t BamRecordSink::recordBytes(const bam1_t& record)
{
    return sizeof(bam1_t) + static_cast<std::size_t>(record.m_data);
}

void BamRecordSink::submit(Records records, WriteErrorPolicy policy)
{
    if (records.empty())
        return;

    std::size_t bytes = records.capacity() * sizeof(Record);
    for (const auto& record : records)
    {
        if (!record)
            throw std::invalid_argument("Cannot submit a null BAM record");
        const auto weight = recordBytes(*record);
        if (weight > std::numeric_limits<std::size_t>::max() - bytes)
            throw std::length_error("BAM record batch is too large");
        bytes += weight;
    }

    enqueue(std::move(records), {}, 0, bytes, policy);
}

void BamRecordSink::submitEncoded(Bytes data, size_t count, size_t workspaceBytes)
{
    if (!state_->writeBytes || data.empty() || !count)
        throw std::invalid_argument("Invalid encoded BAM batch");
    if (workspaceBytes > std::numeric_limits<size_t>::max() - data.capacity())
        throw std::length_error("BAM encoding workspace is too large");
    const size_t bytes = data.capacity() + workspaceBytes;
    enqueue({}, std::move(data), count, bytes, WriteErrorPolicy::kThrowImmediately);
}

void BamRecordSink::enqueue(Records records, Bytes data, size_t count, size_t bytes, WriteErrorPolicy policy)
{
    auto& state = *state_;
    {
        std::unique_lock<std::mutex> lock(state.mutex);
        const auto available = [&] {
            return state.error || state.closing || !state.maxBytes || !state.resident
                || (bytes <= state.maxBytes && state.resident <= state.maxBytes - bytes);
        };
        if (!available())
        {
            ++state.waiting;
            state.credits.wait(lock, available);
            --state.waiting;
        }
        if (state.error)
        {
            if (policy == WriteErrorPolicy::kDeferUntilFinish && !state.closing)
                return;
            std::rethrow_exception(state.error);
        }
        if (state.closing)
            throw std::logic_error("Cannot submit BAM records after finish");
        if (bytes > std::numeric_limits<std::size_t>::max() - state.resident)
            throw std::length_error("Resident BAM batches are too large");
        state.batches.push_back({ std::move(records), std::move(data), {}, count, bytes, {} });
        auto& entry = state.batches.back();
        if (!entry.data.empty() && state.prepareAccounted)
        {
            try
            {
                auto charge = std::make_shared<State::Charge>(bytes);
                entry.compressedCharge = charge;
                std::weak_ptr<State> owner = state_;
                entry.prepared = state.prepareAccounted(std::move(entry.data), [owner, charge](size_t remaining) {
                    charge->requested.store(remaining);
                    if (!charge->active.load()) return;
                    if (auto state = owner.lock())
                    {
                        {
                            std::lock_guard<std::mutex> lock(state->mutex);
                            if (state->error) return;
                            state->applyCapacity(*charge);
                        }
                        state->credits.notify_all();
                    }
                });
            }
            catch (...) { state.batches.pop_back(); throw; }
        }
        state.resident += bytes;
        if (entry.compressedCharge)
        {
            entry.compressedCharge->active.store(true);
            state.applyCapacity(*entry.compressedCharge);
        }
    }
    state.ready.notify_one();
}

void BamRecordSink::finish()
{
    auto& state = *state_;
    // Serialize join() as well as closing, including concurrent finish callers.
    std::lock_guard<std::mutex> finishLock(state.finishing);
    {
        std::lock_guard<std::mutex> lock(state.mutex);
        state.closing = true;
    }
    state.ready.notify_all();
    state.credits.notify_all();
    if (state.consumer.joinable())
        state.consumer.join();
    std::lock_guard<std::mutex> lock(state.mutex);
    if (state.error)
        std::rethrow_exception(state.error);
}

std::size_t BamRecordSink::residentBytes() const
{
    std::lock_guard<std::mutex> lock(state_->mutex);
    return state_->resident;
}

}
