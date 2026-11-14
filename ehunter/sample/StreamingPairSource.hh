// Sequential BAM admission, encounter-order pairing and locus routing.
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

extern "C"
{
#include "htslib/sam.h"
}

#include "core/Read.hh"
#include "sample/AnalyzerFinder.hh"

namespace ehunter
{

class GenomeMask;
namespace htshelpers
{
class HtsFileStreamer;
}

// Own only the BAM fields needed by the existing sequence decoder. Copying this
// payload never retains the streamer's reusable BAM allocation. A raw payload
// contains BAM bases/qualities; a cold payload contains exact decoded-case bases.
class PackedRead
{
public:
    explicit PackedRead(const bam1_t& alignment, bool decodedCold = false);
    Read decode() const;
    // Consume the owned sequence/name. Cold storage moves directly into Read;
    // raw storage is decoded and released. Afterwards only destruction or
    // assignment is supported; no input buffers survive into the returned Read.
    Read releaseRead();
    size_t length() const { return static_cast<size_t>(length_); }
    const std::string& fragmentId() const { return name_; }
    bool isDecoded() const { return (flags_ & decodedFlag) != 0; }
    // Charge retained string capacity, including its terminator, conservatively
    // counting even inline string storage. Cold/consumed records charge zero.
    size_t rawPayloadBytes() const { return isDecoded() || payload_.empty() ? 0 : payload_.capacity() + 1; }
    // Conservative retained storage, excluding allocator bookkeeping.
    size_t retainedBytes() const { return sizeof(*this) + name_.capacity() + 1 + payload_.capacity() + 1; }

private:
    friend class StreamingPairSource;
    void makeDecoded();
    // SAM flags occupy 16 bits, leaving this private storage tag disjoint.
    static constexpr uint32_t decodedFlag = uint32_t(1) << 31;
    std::string name_;
    std::string payload_;
    int32_t length_;
    uint32_t flags_;
};

struct StreamingInputPair
{
    // The arriving record is read; the earlier admitted record is mate.
    PackedRead read;
    PackedRead mate;
    std::vector<AnalyzerBundle> bundles;
};

// A single caller owns input admission, the pending-QNAME table and the mutable
// coordinate query. No input buffers or pending-table references escape next().
// Input, mask and finder must outlive the source. Pending unmatched records are
// intentionally retained until a matching encounter or EOF; they are not
// bounded by genomic distance or by a transient worker-batch budget.
class StreamingPairSource
{
public:
    static constexpr std::size_t defaultPendingRawBytes = 2000000;
    StreamingPairSource(htshelpers::HtsFileStreamer& input, const GenomeMask& mask, AnalyzerFinder& finder,
        std::size_t pendingRawBudget = defaultPendingRawBytes,
        std::function<void()> cancellationCheck = {});
    ~StreamingPairSource();
    StreamingPairSource(const StreamingPairSource&) = delete;
    StreamingPairSource& operator=(const StreamingPairSource&) = delete;

    // Return only pairs with nonempty routing, in completion order. Empty means
    // permanent EOF; final unmatched records have been released.
    // The optional cancellation check runs every 1024 scanned primary records,
    // including filtered records, with the interval preserved across calls.
    // It cannot interrupt blocking HTS I/O or the streamer's internal scan
    // past secondary/supplementary records. A thrown check aborts admission;
    // do not resume next() after that exception.
    boost::optional<StreamingInputPair> next();
    std::size_t pendingReadCount() const;
    std::size_t pendingRawBytes() const { return pendingRawBytes_; }
    std::size_t pendingColdReadCount() const { return pendingColdReads_; }
    std::size_t peakPendingReadCount() const { return peakPendingReads_; }
    std::size_t peakPendingRawBytes() const { return peakPendingRawBytes_; }
    std::size_t peakColdPendingReadCount() const { return peakColdPendingReads_; }
    std::size_t coldReadCount() const { return coldReads_; }

private:
    struct PendingReads;
    htshelpers::HtsFileStreamer& input_;
    const GenomeMask& mask_;
    AnalyzerFinder& finder_;
    std::unique_ptr<PendingReads> pending_;
    const std::size_t pendingRawBudget_;
    std::function<void()> cancellationCheck_;
    unsigned recordsUntilCancellationCheck_ = 1024;
    std::size_t pendingRawBytes_ = 0, pendingColdReads_ = 0;
    std::size_t peakPendingReads_ = 0, peakPendingRawBytes_ = 0, peakColdPendingReads_ = 0, coldReads_ = 0;
    bool finished_ = false;
};

}
