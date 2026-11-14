// Owned BAM batches with a dedicated output consumer and byte backpressure.
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

#include "htslib/sam.h"

namespace ehunter
{

class BamRecordSink
{
public:
    struct DestroyRecord
    {
        void operator()(bam1_t* record) const noexcept { bam_destroy1(record); }
    };
    using Record = std::unique_ptr<bam1_t, DestroyRecord>;
    using Records = std::vector<Record>;
    using WriteRecord = std::function<void(const bam1_t*)>;
    enum class WriteErrorPolicy
    {
        kThrowImmediately,
        kDeferUntilFinish
    };
    static constexpr std::size_t defaultBudgetBytes = 8 * 1024 * 1024;
    struct Statistics
    {
        std::size_t peakResidentBytes;
        uint64_t writtenRecords;
        double submitWaitSeconds;
        double drainSeconds;
    };

    // Callback runs on one dedicated consumer thread. The callback and any
    // borrowed file/header must remain valid until finish() returns or throws.
    // Zero explicitly selects an unbounded sink. A batch larger than a nonzero
    // budget is allowed only when no other batch is resident.
    explicit BamRecordSink(WriteRecord writeRecord, std::size_t maxBytes = defaultBudgetBytes);
    ~BamRecordSink();
    BamRecordSink(const BamRecordSink&) = delete;
    BamRecordSink& operator=(const BamRecordSink&) = delete;

    // Ownership is consumed on success and failure. Producers can block for
    // credits, but the consumer never requires a producer's compute executor.
    // Deferred mode preserves seeking compatibility: after a consumer write
    // failure, release submissions without throwing until explicit finish().
    // Validation/allocation errors and submissions after finish still throw.
    void submit(Records records, WriteErrorPolicy policy = WriteErrorPolicy::kThrowImmediately);
    // Close admission, drain accepted batches, join, and propagate write errors.
    // Concurrent/repeated finish calls are supported. Producers should normally
    // be joined first; a producer waiting during finish is awakened and rejected.
    void finish();

    static std::size_t recordBytes(const bam1_t& record);
    // Includes queued and consumer-owned batches, excluding producer-local
    // records, HTS compression buffers, and allocator/container overhead.
    std::size_t residentBytes() const;
    Statistics statistics() const;

private:
    struct State;
    std::unique_ptr<State> state_;
};

}
