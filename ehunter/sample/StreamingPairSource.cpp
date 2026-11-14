#include "sample/StreamingPairSource.hh"

#include <algorithm>
#include <utility>

#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "absl/strings/string_view.h"

#include "core/HtsHelpers.hh"
#include "sample/GenomeMask.hh"
#include "sample/HtsFileStreamer.hh"

namespace ehunter
{
namespace
{

struct ReadNameHash
{
    using is_transparent = void;
    std::size_t operator()(absl::string_view name) const { return absl::Hash<absl::string_view>()(name); }
    std::size_t operator()(const PackedRead& read) const { return (*this)(absl::string_view(read.fragmentId())); }
};

struct ReadNameEqual
{
    using is_transparent = void;
    bool operator()(const PackedRead& left, const PackedRead& right) const { return left.fragmentId() == right.fragmentId(); }
    bool operator()(const PackedRead& left, absl::string_view right) const { return left.fragmentId() == right; }
    bool operator()(absl::string_view left, const PackedRead& right) const { return left == right.fragmentId(); }
};

}

constexpr uint32_t PackedRead::decodedFlag;
constexpr std::size_t StreamingPairSource::defaultPendingRawBytes;

PackedRead::PackedRead(const bam1_t& alignment, bool decodedCold)
    : name_(bam_get_qname(&alignment))
    , length_(alignment.core.l_qseq)
    , flags_(alignment.core.flag)
{
    if (name_.empty() || length_ <= 0)
    {
        // Keep ReadId/Read validation and its error text without decoding any
        // valid sequence. Invalid names/lengths throw before reading the arrays.
        (void)htshelpers::decodeReadParts(name_.c_str(), flags_, length_, nullptr, nullptr);
    }
    if (decodedCold)
    {
        auto read = htshelpers::decodeReadParts(
            name_.c_str(), flags_, length_, bam_get_seq(&alignment), bam_get_qual(&alignment));
        payload_ = read.releaseSequence();
        flags_ |= decodedFlag;
    }
    else
    {
        const size_t packedBytes = (length() + 1) / 2;
        payload_.assign(reinterpret_cast<const char*>(bam_get_seq(&alignment)), packedBytes + length());
    }
}

Read PackedRead::decode() const
{
    if (isDecoded())
        return Read(ReadId(name_, flags_ & BAM_FREAD1 ? MateNumber::kFirstMate : MateNumber::kSecondMate),
            payload_, (flags_ & BAM_FREVERSE) != 0);
    const auto* data = reinterpret_cast<const uint8_t*>(payload_.data());
    return htshelpers::decodeReadParts(
        name_.c_str(), flags_, length_, data, data + (length() + 1) / 2);
}

void PackedRead::makeDecoded()
{
    if (isDecoded())
        return;
    auto read = decode();
    auto sequence = read.releaseSequence();
    payload_.swap(sequence); // Release the previous raw capacity, even for SSO.
    flags_ |= decodedFlag;
}

Read PackedRead::releaseRead()
{
    if (isDecoded())
    {
        Read read(ReadId(std::move(name_), flags_ & BAM_FREAD1 ? MateNumber::kFirstMate : MateNumber::kSecondMate),
            std::move(payload_), (flags_ & BAM_FREVERSE) != 0);
        length_ = 0;
        return read;
    }
    auto read = decode();
    std::string().swap(name_);
    std::string().swap(payload_);
    length_ = 0;
    return read;
}

struct StreamingPairSource::PendingReads
{
    absl::flat_hash_set<PackedRead, ReadNameHash, ReadNameEqual> reads;
};

StreamingPairSource::StreamingPairSource(
    htshelpers::HtsFileStreamer& input, const GenomeMask& mask, AnalyzerFinder& finder, std::size_t pendingRawBudget,
    std::function<void()> cancellationCheck)
    : input_(input)
    , mask_(mask)
    , finder_(finder)
    , pending_(new PendingReads())
    , pendingRawBudget_(pendingRawBudget)
    , cancellationCheck_(std::move(cancellationCheck))
{
    pending_->reads.reserve(1024);
}

StreamingPairSource::~StreamingPairSource() = default;

std::size_t StreamingPairSource::pendingReadCount() const
{
    return pending_ ? pending_->reads.size() : 0;
}

boost::optional<StreamingInputPair> StreamingPairSource::next()
{
    if (finished_)
        return boost::none;

    while (input_.trySeekingToNextPrimaryAlignment() && input_.isStreamingAlignedReads())
    {
        // Poll before filtering or retaining another pending mate. The
        // counter belongs to the source, not to an individual next() call.
        if (--recordsUntilCancellationCheck_ == 0)
        {
            recordsUntilCancellationCheck_ = 1024;
            if (cancellationCheck_)
                cancellationCheck_();
        }
        // Evaluation order is part of admission compatibility: read-mask hit
        // short-circuits the mate mask, and both precede paired-flag filtering.
        if (!(mask_.query(input_.currentReadContigId(), input_.currentReadPosition())
                || mask_.query(input_.currentMateContigId(), input_.currentMatePosition()))
            || !input_.currentIsPaired())
            continue;

        const bam1_t* raw = input_.currentAlignment();
        const char* name = bam_get_qname(raw);
        // An irrelevant completed pair must still reject an empty admitted
        // QNAME/query. Delegate exception text and validation order to Read.
        if (!*name || raw->core.l_qseq <= 0)
            (void)htshelpers::decodeReadParts(name, raw->core.flag, raw->core.l_qseq, nullptr, nullptr);

        const auto found = pending_->reads.find(absl::string_view(name));
        if (found == pending_->reads.end())
        {
            const auto length = static_cast<std::size_t>(raw->core.l_qseq);
            const auto available = pendingRawBudget_ - pendingRawBytes_;
            PackedRead first(*raw, (length + 1) / 2 + length > available);
            // String capacity can exceed the requested byte length. Account its
            // actual retained capacity before admitting a raw pending payload.
            if (first.rawPayloadBytes() > available)
                first.makeDecoded();
            const auto rawCharge = first.rawPayloadBytes();
            const bool cold = first.isDecoded();
            pending_->reads.emplace(std::move(first));
            pendingRawBytes_ += rawCharge;
            if (cold)
            {
                ++pendingColdReads_;
                ++coldReads_;
            }
            peakPendingReads_ = std::max(peakPendingReads_, pending_->reads.size());
            peakPendingRawBytes_ = std::max(peakPendingRawBytes_, pendingRawBytes_);
            peakColdPendingReads_ = std::max(peakColdPendingReads_, pendingColdReads_);
            continue;
        }

        const int64_t position = input_.currentReadPosition();
        const int64_t matePosition = input_.currentMatePosition();
        auto bundles = finder_.queryStreaming(
            input_.currentReadContigId(), position, position + raw->core.l_qseq,
            input_.currentMateContigId(), matePosition,
            matePosition + static_cast<int64_t>(found->length()));
        if (bundles.empty())
        {
            // Consume both encounters, even when neither routes to a locus.
            // A later repetition of this QNAME starts an independent pair.
            pendingRawBytes_ -= found->rawPayloadBytes();
            pendingColdReads_ -= found->isDecoded() ? 1 : 0;
            pending_->reads.erase(found);
            continue;
        }

        PackedRead arriving(*raw);
        const auto rawCharge = found->rawPayloadBytes();
        const bool cold = found->isDecoded();
        auto earlier = pending_->reads.extract(found);
        pendingRawBytes_ -= rawCharge;
        pendingColdReads_ -= cold ? 1 : 0;
        return StreamingInputPair { std::move(arriving), std::move(earlier.value()), std::move(bundles) };
    }

    finished_ = true;
    pending_.reset();
    pendingRawBytes_ = 0;
    pendingColdReads_ = 0;
    return boost::none;
}

}
