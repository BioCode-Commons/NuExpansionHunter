// Ordered, bounded streaming read preprocessing. HTS access stays on the caller.
#pragma once

#include <condition_variable>
#include <atomic>
#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>
#include <boost/optional.hpp>

#include "core/Read.hh"
#include "sample/GenomeMask.hh"
#include "sample/HtsFileStreamer.hh"

namespace ehunter
{
struct DecodedStreamRecord
{
    int32_t contig, mateContig;
    int64_t position, matePosition;
    Read read;
};

class DecodedReadStream
{
public:
    // Decimal MB. This limits new per-worker batch storage, not process RSS,
    // existing mate/evidence queues, HTS buffers, allocator overhead or stacks.
    static constexpr std::size_t workerCacheLimit = 20000000;
    DecodedReadStream(htshelpers::HtsFileStreamer& input, const GenomeMask& mask, unsigned workers);
    ~DecodedReadStream();
    DecodedReadStream(const DecodedReadStream&) = delete;
    DecodedReadStream& operator=(const DecodedReadStream&) = delete;
    // Borrowed until next() is called again or the stream is destroyed. The
    // caller may move the owned Read out, but must not retain this pointer.
    DecodedStreamRecord* next();
    std::size_t peakWorkerCacheBytes() const;
    std::size_t synchronousLargeRecords() const { return largeRecords_; }

private:
    static constexpr std::size_t maxRecords = 1024;
    static constexpr std::size_t maxRawBytes = 512 * 1024;
    struct RawRecord
    {
        int32_t contig, mateContig, length;
        int64_t position, matePosition;
        uint32_t flags;
        std::size_t name, sequence, quality;
    };
    struct Slot
    {
        std::vector<uint8_t> raw;
        std::vector<RawRecord> records;
        std::vector<DecodedStreamRecord> decoded;
        std::mutex mutex;
        std::condition_variable cv;
        bool done = false;
        std::exception_ptr error;
        std::atomic<std::size_t> highWater { 0 };
        std::size_t owner = 0;
    };
    struct Worker
    {
        std::mutex mutex;
        std::condition_variable cv;
        std::array<Slot*, 2> jobs {{ nullptr, nullptr }};
        std::size_t jobHead = 0, jobCount = 0;
        bool stop = false;
        std::thread thread;
    };
    bool seekSelected();
    bool fill(Slot& slot);
    void decode(Worker& worker);
    void shutdown() noexcept;
    void account(Slot& slot);
    DecodedStreamRecord* decodeCurrent();

    htshelpers::HtsFileStreamer& input_;
    const GenomeMask& mask_;
    std::vector<std::unique_ptr<Slot>> slots_;
    std::vector<std::unique_ptr<Worker>> workers_;
    std::deque<std::size_t> pending_;
    Slot* consuming_ = nullptr;
    std::size_t consumed_ = 0, largeRecords_ = 0;
    bool started_ = false, hasRecord_ = false, finished_ = false, largePending_ = false;
    std::exception_ptr inputError_;
    boost::optional<DecodedStreamRecord> directRecord_;
};
}
