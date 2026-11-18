//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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
//
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: streaming exposes physical BAM ranges and borrowed record fields with htslib fallback.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "sample/HtsFileStreamer.hh"

#include "core/HtsHelpers.hh"
#include "core/OrderedBlocks.hh"
#include "htslib/bgzf.h"
#include "htslib/hts_endian.h"
#include <array>
#include <fstream>
#include <sys/stat.h>

using std::string;

namespace ehunter
{

namespace htshelpers
{

void HtsFileStreamer::openHtsFile(const unsigned decompressionThreads)
{
    htsFilePtr_ = sam_open(htsFilePath_.c_str(), "r");

    if (!htsFilePtr_)
    {
        throw std::runtime_error("Failed to read BAM file " + htsFilePath_);
    }

    // Required step for parsing of some CRAMs
    if (hts_set_fai_filename(htsFilePtr_, htsReferencePath_.c_str()) != 0)
    {
        throw std::runtime_error("Failed to set index of: " + htsReferencePath_);
    }

    /// Create thread pool for bgzf block decompression. The behavior of htslib seems to be to use this pool instead of
    /// (not in addition to) the calling thread, therefore there is no point in creating a decompression thread-pool
    /// with less than 2 threads.
    if (decompressionThreads > 1)
    {
        htsThreadPool_.pool = hts_tpool_init(decompressionThreads);
        if (not htsThreadPool_.pool)
        {
            throw std::runtime_error(
                "HtsFileStreamer: Error creating htslib threadpool with " + std::to_string(decompressionThreads)
                + " threads.");
        }
        hts_set_opt(htsFilePtr_, HTS_OPT_THREAD_POOL, &htsThreadPool_);
    }
}

void HtsFileStreamer::loadHeader()
{
    htsHeaderPtr_ = sam_hdr_read(htsFilePtr_);

    if (!htsHeaderPtr_)
    {
        throw std::runtime_error("Failed to read header of " + htsFilePath_);
    }

    contigInfo_ = htshelpers::decodeContigInfo(htsHeaderPtr_);
}

void HtsFileStreamer::prepareForStreamingAlignments() { htsAlignmentPtr_ = bam_init1(); }

bool HtsFileStreamer::trySeekingToNextPrimaryAlignment()
{
    if (status_ != Status::kStreamingReads)
    {
        return false;
    }

    while (nextPhysicalRecord())
        if (currentIsPrimary())
            return true;
    return false;
}

bool HtsFileStreamer::nextPhysicalRecord()
{
    const auto result = sam_read1(htsFilePtr_, htsHeaderPtr_, htsAlignmentPtr_);
    if (result >= 0)
        return true;
    if (result < -1)
        throw std::runtime_error("Failed to extract a record from " + htsFilePath_);
    status_ = Status::kFinishedStreaming;
    return false;
}

// A common short-read record has one M operation and fits inside the current
// verified BGZF block. Its sequence and quality can be consumed in place. Hand
// every boundary, complex CIGAR, repairable name and malformed record to HTSlib.
static bool readContiguousMatch(BGZF& block, int targetCount, HtsFileStreamer::StreamingRecord& result)
{
    const int available = block.block_length - block.block_offset;
    if (available <= 36) return false;
    const auto* bytes = static_cast<const uint8_t*>(block.uncompressed_block) + block.block_offset;
    const auto size = le_to_u32(bytes);
    if (size < 32 || size >= static_cast<uint32_t>(available - 4)) return false;
    const uint8_t* core = bytes + 4;
    const uint32_t nameLength = core[8];
    const int32_t length = le_to_i32(core + 16);
    if (!nameLength || length <= 0 || le_to_u16(core + 12) != 1) return false;
    const uint64_t sequenceBytes = (static_cast<uint64_t>(length) + 1) / 2;
    if (32 + nameLength + 4 + sequenceBytes + length > size) return false;
    const uint8_t* name = core + 32;
    if (name[nameLength - 1] || le_to_u32(name + nameLength) != (uint64_t(length) << BAM_CIGAR_SHIFT))
        return false;
    const int32_t tid = le_to_i32(core), mtid = le_to_i32(core + 20);
    if (tid < -1 || mtid < -1 || tid >= targetCount || mtid >= targetCount) return false;
    const uint8_t* sequence = name + nameLength + 4;
    result = { tid, mtid, le_to_i32(core + 4), le_to_i32(core + 24), le_to_u16(core + 14), length,
        reinterpret_cast<const char*>(name), sequence, sequence + sequenceBytes };
    // The strict bound above leaves a byte in the current block. HTSlib owns
    // all block transitions, EOF detection, CRC checking and virtual seeks.
    block.block_offset += size + 4;
    block.uncompressed_address += size + 4;
    return true;
}

// Contract [K-30]
// Ensures:  for a well-formed record wholly inside one BGZF block the in-place view's fields equal htslib's
//           sam_read1; otherwise the sam_read1 fallback runs; no byte outside the block is read.
// Reference: sam_read1 (the fallback path), with htslib 1.10.2 as pinned by the superbuild.
bool HtsFileStreamer::nextStreamingRecord(StreamingRecord& result)
{
    if (htsFilePtr_->format.format == bam && readContiguousMatch(*htsFilePtr_->fp.bgzf, htsHeaderPtr_->n_targets, result))
        return true;
    if (!nextPhysicalRecord()) return false;
    const auto& c = htsAlignmentPtr_->core;
    result = { c.tid, c.mtid, c.pos, c.mpos, c.flag, c.l_qseq, bam_get_qname(htsAlignmentPtr_),
        bam_get_seq(htsAlignmentPtr_), bam_get_qual(htsAlignmentPtr_) };
    return true;
}

bool HtsFileStreamer::currentIsPrimary() const { return isPrimaryAlignment(htsAlignmentPtr_); }

int64_t HtsFileStreamer::virtualOffset() const { return bgzf_tell(htsFilePtr_->fp.bgzf); }

void HtsFileStreamer::seekVirtualOffset(int64_t offset)
{
    if (bgzf_seek(htsFilePtr_->fp.bgzf, offset, SEEK_SET) < 0)
        throw std::runtime_error("Cannot seek in " + htsFilePath_);
    status_ = Status::kStreamingReads;
}

std::unique_ptr<HtsFileStreamer> HtsFileStreamer::independentReader() const
{
    return std::unique_ptr<HtsFileStreamer>(new HtsFileStreamer(htsFilePath_, htsReferencePath_));
}

// Recognize two complete mapped records at the start of a probed BGZF block.
// This is a hint only. RangeReader accepts this cut only when the preceding
// physical scan reaches its exact virtual offset at a BAM record boundary.
static size_t candidateRecordSize(const uint8_t* data, size_t available, int targets)
{
    if (available < 36) return 0;
    const size_t bytes = le_to_u32(data);
    if (bytes < 32 || bytes > available - 4) return 0;
    const auto* core = data + 4;
    const int32_t tid = le_to_i32(core), mateTid = le_to_i32(core + 20), length = le_to_i32(core + 16);
    const size_t name = core[8], cigar = le_to_u16(core + 12);
    if (tid < 0 || tid >= targets || mateTid < -1 || mateTid >= targets || le_to_i32(core + 4) < 0
        || length <= 0 || name < 2 || !cigar
        || 32 + name + 4 * cigar + (uint64_t(length) + 1) / 2 + length > bytes)
        return 0;
    if (core[32 + name - 1] != 0) return 0;
    for (size_t j = 0; j + 1 < name; ++j)
        if (core[32 + j] < 33 || core[32 + j] > 126) return 0;
    return bytes + 4;
}

static void splitLargeBamRanges(const string& path, int targets, unsigned workers, std::vector<int64_t>& offsets)
{
    constexpr int64_t stride = 1024 * 1024;
    std::vector<int64_t> probes;
    for (size_t i = 1; i < offsets.size(); ++i)
    {
        const auto begin = offsets[i - 1] >> 16, end = offsets[i] >> 16;
        if (end - begin >= 2 * stride)
            for (auto address = begin + stride; address < end - stride / 2; address += stride)
                probes.push_back(address);
    }
    if (probes.empty()) return;
    auto scan = [&](int, size_t begin, size_t end) {
        std::vector<int64_t> cuts;
        std::ifstream raw(path, std::ios::binary);
        std::unique_ptr<BGZF, decltype(&bgzf_close)> block(bgzf_open(path.c_str(), "r"), bgzf_close);
        if (!raw || !block) return cuts;
        // Common canonical BGZF header. Files written with another valid gzip
        // header keep their index-derived cuts.
        const uint8_t header[] = {31,139,8,4,0,0,0,0,0,255,6,0,'B','C',2,0};
        std::array<uint8_t, BGZF_MAX_BLOCK_SIZE + sizeof(header)> buffer;
        for (size_t i = begin; i < end; ++i)
        {
            raw.clear();
            raw.seekg(probes[i]);
            raw.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
            const auto last = buffer.begin() + raw.gcount();
            const auto hit = std::search(buffer.begin(), last, std::begin(header), std::end(header));
            if (hit == last) continue;
            const int64_t offset = (probes[i] + (hit - buffer.begin())) << 16;
            if (bgzf_seek(block.get(), offset, SEEK_SET) < 0 || bgzf_read_block(block.get()) < 0) continue;
            const auto* bytes = static_cast<const uint8_t*>(block->uncompressed_block);
            const size_t size = block->block_length;
            const size_t first = candidateRecordSize(bytes, size, targets);
            if (first && candidateRecordSize(bytes + first, size - first, targets))
                cuts.push_back(offset);
        }
        return cuts;
    };
    orderedBlocks(probes.size(), 32, workers, scan, [&](const std::vector<int64_t>& cuts) {
        offsets.insert(offsets.end(), cuts.begin(), cuts.end());
    });
    std::sort(offsets.begin(), offsets.end());
}

std::vector<int64_t> HtsFileStreamer::bamRangeStarts(unsigned workers) const
{
    struct stat file;
    if (htsFilePtr_->format.format != bam || stat(htsFilePath_.c_str(), &file) || !S_ISREG(file.st_mode))
        return {};
    std::unique_ptr<hts_idx_t, decltype(&hts_idx_destroy)> index(
        sam_index_load3(htsFilePtr_, htsFilePath_.c_str(), nullptr, HTS_IDX_SILENT_FAIL), hts_idx_destroy);
    if (!index)
        return {};
    std::vector<int64_t> offsets {virtualOffset()};
    for (int tid = 0; tid < htsHeaderPtr_->n_targets; ++tid)
    {
        for (int64_t pos = 0; pos < htsHeaderPtr_->target_len[tid]; pos += 32768)
        {
            std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> it(
                sam_itr_queryi(index.get(), tid, pos, pos + 1), hts_itr_destroy);
            if (!it || !it->n_off)
                continue;
            const int64_t offset = it->off[0].u;
            if ((offset >> 16) - (offsets.back() >> 16) >= 256 * 1024 && (offset >> 16) < file.st_size)
                offsets.push_back(offset);
        }
    }
    if (offsets.size() == 1)
        offsets.clear();
    else if (workers > 1)
        splitLargeBamRanges(htsFilePath_, htsHeaderPtr_->n_targets, workers, offsets);
    return offsets;
}

bool HtsFileStreamer::isStreamingAlignedReads() const
{
    return status_ != Status::kFinishedStreaming && currentReadContigId() != -1;
}

Read HtsFileStreamer::decodeRead() const { return htshelpers::decodeRead(htsAlignmentPtr_); }

HtsFileStreamer::~HtsFileStreamer()
{
    bam_destroy1(htsAlignmentPtr_);
    htsAlignmentPtr_ = nullptr;

    bam_hdr_destroy(htsHeaderPtr_);
    htsHeaderPtr_ = nullptr;

    sam_close(htsFilePtr_);
    htsFilePtr_ = nullptr;

    if (htsThreadPool_.pool)
    {
        hts_tpool_destroy(htsThreadPool_.pool);
    }
}

}

}
