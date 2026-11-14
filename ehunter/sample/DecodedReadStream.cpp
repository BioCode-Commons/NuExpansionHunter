#include "sample/DecodedReadStream.hh"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <stdexcept>
#include "core/HtsHelpers.hh"

namespace ehunter
{
constexpr std::size_t DecodedReadStream::workerCacheLimit;
constexpr std::size_t DecodedReadStream::maxRecords;
constexpr std::size_t DecodedReadStream::maxRawBytes;

DecodedReadStream::DecodedReadStream(
    htshelpers::HtsFileStreamer& input, const GenomeMask& mask, unsigned workers)
    : input_(input), mask_(mask)
{
    // Allocate all control blocks before starting any thread. Partial thread
    // creation failure must join the threads already started.
    workers_.reserve(workers);
    for (unsigned i = 0; i < workers; ++i)
        workers_.emplace_back(new Worker);
    slots_.reserve(std::size_t(workers) * 2);
    for (std::size_t i = 0; i < std::size_t(workers) * 2; ++i)
    {
        slots_.emplace_back(new Slot);
        slots_.back()->owner = i % workers;
    }
    try
    {
        for (auto& worker : workers_)
            worker->thread = std::thread(&DecodedReadStream::decode, this, std::ref(*worker));
    }
    catch (...)
    {
        shutdown();
        throw;
    }
}

DecodedReadStream::~DecodedReadStream() { shutdown(); }

void DecodedReadStream::shutdown() noexcept
{
    for (auto& worker : workers_)
    {
        { std::lock_guard<std::mutex> lock(worker->mutex); worker->stop = true; }
        worker->cv.notify_one();
    }
    for (auto& worker : workers_)
        if (worker->thread.joinable())
            worker->thread.join();
}

bool DecodedReadStream::seekSelected()
{
    if (hasRecord_)
        return true;
    if (finished_)
        return false;
    // A run of rejected records often stays in the same two mask bins. Keep
    // these tiny memos only during this seek: the caller can change the mask
    // between next() calls, and such changes must be observed immediately.
    struct BinResult
    {
        int32_t contig = -1;
        uint64_t bin = 0;
        bool result = false;
        bool valid = false;
    };
    BinResult readBin, mateBin;
    const auto queryBin = [&](BinResult& cached, int32_t contig, int64_t position) {
        if (contig < 0 || position < 0)
            return false;
        const auto bin = GenomeMask::binIndex(static_cast<uint64_t>(position));
        if (!cached.valid || cached.contig != contig || cached.bin != bin)
        {
            const bool result = mask_.query(contig, position);
            cached.contig = contig;
            cached.bin = bin;
            cached.result = result;
            cached.valid = true;
        }
        return cached.result;
    };
    try
    {
        while (input_.trySeekingToNextPrimaryAlignment() && input_.isStreamingAlignedReads())
        {
            if ((queryBin(readBin, input_.currentReadContigId(), input_.currentReadPosition())
                || queryBin(mateBin, input_.currentMateContigId(), input_.currentMatePosition()))
                && input_.currentIsPaired())
            {
                hasRecord_ = true;
                return true;
            }
        }
    }
    catch (...)
    {
        if (slots_.empty())
            throw;
        // Earlier batches and their errors must be observed before this later
        // input error, even though the reader has looked ahead.
        inputError_ = std::current_exception();
    }
    finished_ = true;
    return false;
}

DecodedStreamRecord* DecodedReadStream::decodeCurrent()
{
    directRecord_.emplace(DecodedStreamRecord { input_.currentReadContigId(), input_.currentMateContigId(),
        input_.currentReadPosition(), input_.currentMatePosition(), input_.decodeRead() });
    hasRecord_ = false;
    return &*directRecord_;
}

const bam1_t* DecodedReadStream::nextUndecoded()
{
    if (!slots_.empty())
        throw std::logic_error("Undecoded streaming requires synchronous read preprocessing");
    directRecord_ = boost::none;
    if (!seekSelected())
        return nullptr;
    const auto* bam = input_.currentAlignment();
    // Routing does not inspect bases, but skipping an eager Read must not hide
    // its validation errors. Delegate unusual invalid records to that decoder.
    if (bam->core.l_qseq <= 0 || !*bam_get_qname(bam))
        input_.decodeRead();
    hasRecord_ = false;
    return bam;
}

void DecodedReadStream::account(Slot& slot)
{
    std::size_t bytes = sizeof(Slot) + slot.raw.capacity()
        + slot.records.capacity() * sizeof(RawRecord)
        + slot.decoded.capacity() * sizeof(DecodedStreamRecord);
    // Counting SSO storage again is conservative. Moved records have left this
    // cache and belong to the original routing/mate queues.
    for (const auto& record : slot.decoded)
        bytes += record.read.fragmentId().capacity() + 1 + record.read.sequence().capacity() + 1;
    slot.highWater.store(std::max(slot.highWater.load(), bytes));
    // Two fixed slots belong to each worker. Sum their individual high-water
    // marks conservatively, including slots waiting for ordered consumption.
    const auto owner = slot.owner;
    const auto combined = slots_[owner]->highWater.load()
        + slots_[owner + workers_.size()]->highWater.load() + sizeof(Worker);
    if (combined > workerCacheLimit)
        throw std::runtime_error("Read preprocessing cache exceeded its 20 MB worker limit");
}

bool DecodedReadStream::fill(Slot& slot)
{
    if (largePending_ || finished_)
        return false;
    slot.raw.clear();
    slot.records.clear();
    slot.decoded.clear();
    // Reuse a modest batch, not the full authorized 20 MB per worker.
    slot.raw.reserve(maxRawBytes);
    slot.records.reserve(maxRecords);
    slot.decoded.reserve(maxRecords);
    while (slot.records.size() < maxRecords && seekSelected())
    {
        const bam1_t* bam = input_.currentAlignment();
        const std::size_t nameBytes = std::strlen(bam_get_qname(bam)) + 1;
        const std::size_t length = bam->core.l_qseq;
        const std::size_t sequenceBytes = (length + 1) / 2;
        const std::size_t needed = nameBytes + sequenceBytes + length;
        if (needed > maxRawBytes)
        {
            // Keep this record in the reader until all earlier batches drain;
            // then decode synchronously without copying its raw payload.
            if (slot.records.empty())
                largePending_ = true;
            break;
        }
        if (needed > maxRawBytes - slot.raw.size())
            break;
        const std::size_t offset = slot.raw.size();
        slot.raw.resize(offset + needed);
        std::memcpy(slot.raw.data() + offset, bam_get_qname(bam), nameBytes);
        if (sequenceBytes)
            std::memcpy(slot.raw.data() + offset + nameBytes, bam_get_seq(bam), sequenceBytes);
        if (length)
            std::memcpy(slot.raw.data() + offset + nameBytes + sequenceBytes, bam_get_qual(bam), length);
        slot.records.push_back({ bam->core.tid, bam->core.mtid, bam->core.l_qseq,
            bam->core.pos, bam->core.mpos, bam->core.flag, offset,
            offset + nameBytes, offset + nameBytes + sequenceBytes });
        hasRecord_ = false;
    }
    account(slot);
    if (slot.records.empty())
        return false;
    {
        std::lock_guard<std::mutex> lock(slot.mutex);
        slot.error = nullptr;
        slot.done = false;
    }
    auto& worker = *workers_[slot.owner];
    {
        std::lock_guard<std::mutex> lock(worker.mutex);
        assert(worker.jobCount < worker.jobs.size());
        worker.jobs[(worker.jobHead + worker.jobCount) % worker.jobs.size()] = &slot;
        ++worker.jobCount;
    }
    worker.cv.notify_one();
    return true;
}

void DecodedReadStream::decode(Worker& worker)
{
    for (;;)
    {
        Slot* current;
        {
            std::unique_lock<std::mutex> lock(worker.mutex);
            worker.cv.wait(lock, [&] { return worker.stop || worker.jobCount; });
            if (worker.stop)
                return;
            current = worker.jobs[worker.jobHead];
            worker.jobHead = (worker.jobHead + 1) % worker.jobs.size();
            --worker.jobCount;
        }
        auto& slot = *current;
        std::exception_ptr error;
        try
        {
            for (const auto& raw : slot.records)
                slot.decoded.push_back({ raw.contig, raw.mateContig, raw.position, raw.matePosition,
                    htshelpers::decodeReadParts(reinterpret_cast<const char*>(slot.raw.data() + raw.name),
                        raw.flags, raw.length, slot.raw.data() + raw.sequence, slot.raw.data() + raw.quality) });
            account(slot);
        }
        catch (...)
        {
            error = std::current_exception();
        }
        {
            std::lock_guard<std::mutex> lock(slot.mutex);
            slot.error = error;
            slot.done = true;
        }
        slot.cv.notify_one();
    }
}

DecodedStreamRecord* DecodedReadStream::next()
{
    // Release the previous synchronous/large result before looking ahead or
    // allocating its successor. Large records must not become hidden cache.
    directRecord_ = boost::none;
    if (slots_.empty())
        return seekSelected() ? decodeCurrent() : nullptr;

    if (consuming_ && consumed_ < consuming_->decoded.size())
        return &consuming_->decoded[consumed_++];
    if (consuming_)
    {
        if (consuming_->error)
            std::rethrow_exception(consuming_->error);
        const auto released = pending_.front();
        pending_.pop_front();
        consuming_->decoded.clear();
        consuming_ = nullptr;
        if (fill(*slots_[released]))
            pending_.push_back(released);
    }
    if (!started_ || pending_.empty())
    {
        started_ = true;
        if (largePending_)
        {
            largePending_ = false;
            ++largeRecords_;
            return decodeCurrent();
        }
        for (std::size_t i = 0; i < slots_.size(); ++i)
        {
            if (fill(*slots_[i]))
                pending_.push_back(i);
            else
                break;
        }
    }
    if (pending_.empty())
    {
        if (largePending_)
        {
            largePending_ = false;
            ++largeRecords_;
            return decodeCurrent();
        }
        if (inputError_)
            std::rethrow_exception(inputError_);
        return nullptr;
    }
    auto& slot = *slots_[pending_.front()];
    {
        std::unique_lock<std::mutex> lock(slot.mutex);
        slot.cv.wait(lock, [&] { return slot.done; });
        if (slot.error && slot.decoded.empty())
            std::rethrow_exception(slot.error);
    }
    consuming_ = &slot;
    consumed_ = 1;
    return &slot.decoded.front();
}

std::size_t DecodedReadStream::peakWorkerCacheBytes() const
{
    std::size_t peak = 0;
    for (std::size_t i = 0; i < workers_.size(); ++i)
        peak = std::max(peak, slots_[i]->highWater.load()
            + slots_[i + workers_.size()]->highWater.load() + sizeof(Worker));
    return peak;
}
}
