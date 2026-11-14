#include "io/BamRecordSink.hh"

#include <condition_variable>
#include <chrono>
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
    struct Envelope
    {
        Records records;
        std::size_t bytes;
    };

    State(WriteRecord callback, std::size_t budget)
        : writeRecord(std::move(callback))
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

            uint64_t written = 0;
            try
            {
                for (auto& record : current.records)
                {
                    writeRecord(record.get());
                    ++written;
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
                    writtenRecords += written;
                    batches.clear();
                }
                credits.notify_all();
                ready.notify_all();
                current.records = Records();
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    resident = 0;
                }
                return;
            }

            current.records = Records();
            {
                std::lock_guard<std::mutex> lock(mutex);
                resident -= current.bytes;
                writtenRecords += written;
            }
            credits.notify_all();
        }
    }

    WriteRecord writeRecord;
    const std::size_t maxBytes;
    mutable std::mutex mutex;
    std::mutex finishing;
    std::condition_variable ready, credits;
    std::deque<Envelope> batches;
    std::size_t resident = 0;
    std::size_t peakResident = 0;
    uint64_t writtenRecords = 0;
    double submitWaitSeconds = 0, drainSeconds = 0;
    bool closing = false; // Explicit finish only; error separately cancels admission.
    std::exception_ptr error;
    std::thread consumer;
};

BamRecordSink::BamRecordSink(WriteRecord writeRecord, std::size_t maxBytes)
    : state_(new State(std::move(writeRecord), maxBytes))
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

    auto& state = *state_;
    {
        std::unique_lock<std::mutex> lock(state.mutex);
        const auto available = [&] {
            return state.error || state.closing || !state.maxBytes || !state.resident
                || (bytes <= state.maxBytes && state.resident <= state.maxBytes - bytes);
        };
        if (!available())
        {
            const auto started = std::chrono::steady_clock::now();
            state.credits.wait(lock, available);
            state.submitWaitSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
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
        state.batches.push_back({ std::move(records), bytes });
        state.resident += bytes;
        if (state.resident > state.peakResident)
            state.peakResident = state.resident;
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
    const bool joining = state.consumer.joinable();
    const auto started = std::chrono::steady_clock::now();
    if (joining)
        state.consumer.join();
    std::lock_guard<std::mutex> lock(state.mutex);
    if (joining)
        state.drainSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    if (state.error)
        std::rethrow_exception(state.error);
}

std::size_t BamRecordSink::residentBytes() const
{
    std::lock_guard<std::mutex> lock(state_->mutex);
    return state_->resident;
}

BamRecordSink::Statistics BamRecordSink::statistics() const
{
    std::lock_guard<std::mutex> lock(state_->mutex);
    return { state_->peakResident, state_->writtenRecords, state_->submitWaitSeconds, state_->drainSeconds };
}

}
