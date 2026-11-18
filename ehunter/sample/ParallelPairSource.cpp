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
#include "sample/ParallelPairSource.hh"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <future>
#include <limits>
#include <mutex>
#include <thread>

#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "absl/strings/string_view.h"
#include "core/HtsHelpers.hh"
#include "core/ThreadPool.hh"
#include "sample/GenomeMask.hh"
#include "sample/HtsFileStreamer.hh"
#include "spdlog/spdlog.h"

namespace ehunter
{
namespace
{
struct RoutedRead
{
    explicit RoutedRead(const htshelpers::HtsFileStreamer::StreamingRecord& record)
        : packed(record.name, record.flag, record.length, record.sequence, record.quality, true),
        nameHash(absl::Hash<absl::string_view>()(packed.fragmentId())) {}
    PackedRead packed;
    size_t nameHash;
};
struct NameHash
{
    size_t operator()(const RoutedRead& read) const { return read.nameHash; }
};
struct NameEqual
{
    bool operator()(const RoutedRead& a, const RoutedRead& b) const
    {
        return a.packed.fragmentId() == b.packed.fragmentId();
    }
};
using Pending = absl::flat_hash_set<RoutedRead, NameHash, NameEqual>;

struct Record
{
    int32_t tid, mtid;
    int64_t position, matePosition;
    uint32_t ordinal;
    RoutedRead read;
};

struct Batch
{
    explicit Batch(size_t laneCount) : lanes(laneCount), completed(laneCount) {}
    std::vector<std::vector<Record>> lanes;
    struct OrderedPair
    {
        uint32_t ordinal;
        StreamingInputPair pair;
    };
    std::vector<std::vector<OrderedPair>> completed;
    uint32_t admitted = 0;
    std::vector<StreamingInputPair*> results;
    bool ended = false;

    void fill(htshelpers::HtsFileStreamer& input, const GenomeMask& mask, const std::function<void()>& check,
              int64_t rangeEnd = -1)
    {
        ended = false;
        size_t payloadBytes = 0;
        admitted = 0;
        for (auto& records : lanes)
            records.clear();
        unsigned scanned = 0;
        while (admitted < kStreamingBatchRecords && payloadBytes < 8 * 1024 * 1024)
        {
            if (rangeEnd >= 0 && input.virtualOffset() >= rangeEnd)
                break;
            htshelpers::HtsFileStreamer::StreamingRecord record;
            if (!input.nextStreamingRecord(record))
            {
                ended = true;
                break;
            }
            if (record.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
                continue;
            if (record.tid == -1)
            {
                ended = true;
                break;
            }
            if (++scanned % 1024 == 0)
                check();
            if (!(mask.query(record.tid, record.pos) || mask.query(record.mtid, record.mpos))
                || !(record.flag & BAM_FPAIRED))
                continue;
            RoutedRead read(record);
            const size_t lane = read.nameHash % lanes.size();
            payloadBytes += read.packed.retainedBytes() + sizeof(read.nameHash);
            lanes[lane].push_back({record.tid, record.mtid, record.pos, record.mpos,
                                   admitted++, std::move(read)});
        }
    }

    void pairLane(size_t lane, Pending& pending, const AnalyzerFinder& finder)
    {
        auto& pairs = completed[lane];
        pairs.clear();
        // Each admitted record can complete at most one pair. Reserving once
        // keeps published pointers stable until the consumer releases the slot.
        pairs.reserve(lanes[lane].size());
        for (auto& record : lanes[lane])
        {
            bool retained = false;
            const auto found = pending.lazy_emplace(record.read, [&](const Pending::constructor& insert) {
                insert(std::move(record.read));
                retained = true;
            });
            if (retained) continue;
            auto& read = record.read.packed;
            auto targets = finder.queryPair(record.tid, record.position, record.position + static_cast<int64_t>(read.length()),
                record.mtid, record.matePosition, record.matePosition + static_cast<int64_t>(found->packed.length()));
            if (targets.empty())
            {
                pending.erase(found);
                continue;
            }
            auto first = pending.extract(found);
            pairs.push_back({record.ordinal, StreamingInputPair{
                std::move(read),
                std::move(first.value().packed), std::move(targets)}});
        }
    }

    void gather()
    {
        // One thread writes this ordinal array after all lane jobs are ready.
        // The mate tables and per-lane completed sequences have independent owners.
        results.assign(admitted, nullptr);
        for (auto& lane : completed)
            for (auto& entry : lane)
                results[entry.ordinal] = &entry.pair;
    }
};

struct PairingBatch
{
    explicit PairingBatch(size_t lanes) : batch(lanes), errors(lanes) {}
    Batch batch;
    std::vector<std::exception_ptr> errors;
    std::atomic<size_t> remaining { 0 };
    std::promise<void> done;
};

// Every name lane visits the same bounded ring in FIFO order. Publish one
// slot and complete one future per batch, independent of the lane count.
class PairingLanes
{
public:
    PairingLanes(size_t lanes, size_t window, const AnalyzerFinder& finder)
        : ring_(window), finder_(finder)
    {
        threads_.reserve(lanes);
        try
        {
            for (size_t lane = 0; lane < lanes; ++lane)
                threads_.emplace_back(&PairingLanes::run, this, lane);
        }
        catch (...)
        {
            finish();
            throw;
        }
    }

    ~PairingLanes() { finish(); }

    std::future<void> publish(PairingBatch& job)
    {
        job.done = std::promise<void>();
        auto completion = job.done.get_future();
        job.remaining.store(threads_.size(), std::memory_order_relaxed);
        {
            std::lock_guard<std::mutex> lock(mutex_);
            // The owner reuses a slot only after its published results drain.
            ring_[published_ % ring_.size()] = &job;
            ++published_;
        }
        ready_.notify_all();
        return completion;
    }

private:
    void finish()
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        ready_.notify_all();
        for (auto& thread : threads_) thread.join();
    }

    void run(size_t lane)
    {
        Pending pending;
        for (size_t sequence = 0;; ++sequence)
        {
            PairingBatch* job;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                ready_.wait(lock, [&] { return sequence < published_ || stopping_; });
                if (sequence == published_) return;
                job = ring_[sequence % ring_.size()];
            }
            try { job->batch.pairLane(lane, pending, finder_); }
            catch (...) { job->errors[lane] = std::current_exception(); }
            if (job->remaining.fetch_sub(1, std::memory_order_acq_rel) == 1)
            {
                auto done = std::move(job->done);
                auto error = std::find_if(job->errors.begin(), job->errors.end(),
                    [](const std::exception_ptr& e) { return static_cast<bool>(e); });
                if (error != job->errors.end())
                    done.set_exception(*error);
                else
                    done.set_value();
            }
        }
    }

    std::vector<PairingBatch*> ring_;
    const AnalyzerFinder& finder_;
    std::vector<std::thread> threads_;
    std::mutex mutex_;
    std::condition_variable ready_;
    size_t published_ = 0;
    bool stopping_ = false;
};

// A bounded window of independent HTS handles scans physical file ranges.
// Results are swapped back into input order before consulting any mate table.
// A partial range can yield many batches without buffering the entire range.
class RangeReader
{
    struct Slot
    {
        Slot(size_t lanes, std::unique_ptr<htshelpers::HtsFileStreamer> reader)
            : batch(lanes), input(std::move(reader)) {}
        Batch batch;
        std::unique_ptr<htshelpers::HtsFileStreamer> input;
        int64_t end = 0;
        std::future<void> job;
    };

public:
    RangeReader(htshelpers::HtsFileStreamer& input, const GenomeMask& mask,
                const std::function<void()>& check, size_t lanes, unsigned workers)
        : main_(input), mask_(mask), check_(check),
          starts_(workers > 1 ? input.bamRangeStarts(workers) : std::vector<int64_t>())
    {
        const size_t count = std::min(size_t(workers), starts_.size());
        if (!count)
            return;
        pool_.reset(new ctpl::thread_pool(static_cast<int>(count)));
        // Keep further ranges in flight when an early range takes longer.
        // Handle count controls lookahead; the decode pool still has count threads.
        const size_t lookahead = std::min(starts_.size(), 4 * count);
        for (size_t i = 0; i < lookahead; ++i)
        {
            std::unique_ptr<htshelpers::HtsFileStreamer> reader;
            try { reader = input.independentReader(); }
            catch (const std::exception&)
            {
                if (i < count)
                    throw;
                break; // Optional lookahead must fit the process file-descriptor limit.
            }
            slots_.emplace_back(new Slot(lanes, std::move(reader)));
            scheduleNext(*slots_.back());
        }
        spdlog::info("Indexed BAM scan: ranges={} readers={}", starts_.size(), count);
    }

    void fill(Batch& output)
    {
        if (pending_.empty())
        {
            output.fill(main_, mask_, check_);
            return;
        }
        auto& slot = *pending_.front();
        pending_.pop_front();
        slot.job.get();
        std::swap(slot.batch, output);
        if (output.ended)
            return;
        const int64_t offset = slot.input->virtualOffset();
        if (offset < slot.end)
        {
            launch(slot, offset);
            pending_.push_front(&slot);
        }
        else if (offset == slot.end)
            scheduleNext(slot);
        else
        {
            // A stale index supplied a cut inside the last record just read.
            // All records through offset came from a verified start. Discard
            // speculative ranges and continue at that exact record boundary.
            pending_.clear();
            main_.seekVirtualOffset(offset);
        }
    }
private:
    void scheduleNext(Slot& slot)
    {
        if (next_ == starts_.size())
            return;
        const auto start = starts_[next_++];
        slot.end = next_ < starts_.size() ? starts_[next_] : std::numeric_limits<int64_t>::max();
        launch(slot, start);
        pending_.push_back(&slot);
    }
    void launch(Slot& slot, int64_t start)
    {
        slot.job = pool_->push([this, &slot, start](int) {
            if (slot.input->virtualOffset() != start)
                slot.input->seekVirtualOffset(start);
            slot.batch.fill(*slot.input, mask_, check_, slot.end);
        });
    }

    htshelpers::HtsFileStreamer& main_;
    const GenomeMask& mask_;
    const std::function<void()>& check_;
    std::vector<int64_t> starts_;
    size_t next_ = 0;
    std::vector<std::unique_ptr<Slot>> slots_;
    std::deque<Slot*> pending_;
    // Join outstanding jobs before releasing their slots, including on error.
    std::unique_ptr<ctpl::thread_pool> pool_;
};
}

size_t streamingWindowSlots(size_t laneCount) { return std::max(size_t(4), std::min(size_t(96), 4 * laneCount)); }

// Contract [K-34]
// Ensures:  for every record sequence, the pairing lanes emit the same pairs in the same order as StreamingPairSource
//           (v5.0.0's encounter-order pairing).
// Reference: StreamingPairSource::next.
void streamParallelPairs(htshelpers::HtsFileStreamer& input, const GenomeMask& mask, const AnalyzerFinder& finder,
    size_t laneCount, const std::function<void()>& check, unsigned decodeWorkers, const PrepareStreamingBatch& prepare)
{
    struct Slot : PairingBatch
    {
        using PairingBatch::PairingBatch;
        std::future<void> paired;
        std::future<std::function<void()>> prepared;
    };
    const size_t window = streamingWindowSlots(laneCount);
    std::deque<Slot> slots;
    for (size_t i = 0; i < window; ++i)
        slots.emplace_back(laneCount);
    RangeReader ranges(input, mask, check, laneCount, decodeWorkers);
    PairingLanes pairing(laneCount, window, finder);
    std::unique_ptr<ctpl::thread_pool> staging(
        new ctpl::thread_pool(static_cast<int>(std::max(size_t(1), std::min(laneCount, size_t(decodeWorkers))))));
    std::mutex ringMutex;
    std::condition_variable canFill, canPublish;
    size_t written = 0, consumed = 0;
    size_t refillAt = 0;
    bool publisherWaiting = false;
    bool stopping = false;
    std::exception_ptr readError;
    // One reader owns the input order for the entire scan. A consumed slot
    // returns directly to it without allocating and dispatching another task.
    std::thread reader([&] {
        try
        {
            for (size_t sequence = 0;; ++sequence)
            {
                {
                    std::unique_lock<std::mutex> lock(ringMutex);
                    if (sequence - consumed == window)
                    {
                        // Returning several slots together amortizes wakeups
                        // while the publisher still has the rest of the window.
                        refillAt = consumed + std::max(size_t(1), window / 4);
                        canFill.wait(lock, [&] { return stopping || consumed >= refillAt; });
                        refillAt = 0;
                    }
                    if (stopping) return;
                }
                Slot& slot = slots[sequence % window];
                ranges.fill(slot.batch);
                const bool ended = slot.batch.ended;
                slot.paired = pairing.publish(slot);
                {
                    auto* ptr = &slot;
                    slot.prepared = staging->push([&, ptr](int) {
                        ptr->paired.get();
                        ptr->batch.gather();
                        return prepare(ptr->batch.results);
                    });
                }
                {
                    std::lock_guard<std::mutex> lock(ringMutex);
                    ++written;
                    if (publisherWaiting) canPublish.notify_one();
                }
                if (ended) return;
            }
        }
        catch (...)
        {
            std::lock_guard<std::mutex> lock(ringMutex);
            readError = std::current_exception();
            canPublish.notify_one();
        }
    });
    auto joinReader = [&] {
        {
            std::lock_guard<std::mutex> lock(ringMutex);
            stopping = true;
        }
        canFill.notify_one();
        reader.join();
    };
    try
    {
        for (size_t sequence = 0;; ++sequence)
        {
            auto& slot = slots[sequence % window];
            {
                std::unique_lock<std::mutex> lock(ringMutex);
                if (sequence == written)
                {
                    publisherWaiting = true;
                    canPublish.wait(lock, [&] { return sequence < written || readError; });
                    publisherWaiting = false;
                }
                if (sequence == written) std::rethrow_exception(readError);
            }
            std::function<void()> commit = slot.prepared.get();
            commit();
            commit = {};
            if (slot.batch.ended) break;
            {
                std::lock_guard<std::mutex> lock(ringMutex);
                ++consumed;
                if (consumed == refillAt) canFill.notify_one();
            }
        }
    }
    catch (...)
    {
        joinReader();
        throw;
    }
    joinReader();
}
}
