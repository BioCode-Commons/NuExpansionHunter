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

#pragma once

#include <memory>
#include <string>
#include <vector>

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"

namespace ehunter
{

namespace htshelpers
{

class HtsFileStreamer
{
public:
    /// \param[in] decompressionThreads Total size of thread pool used for bgzip decompression of the hts file. Has no
    /// effect if the file is uncompressed. When set to one or less the calling thread handles all decompression and no
    /// thread pool is used.
    ///
    HtsFileStreamer(
        const std::string& htsFilePath, const std::string& htsReferencePath, const unsigned decompressionThreads = 1)
        : htsFilePath_(htsFilePath)
        , htsReferencePath_(htsReferencePath)
        , contigInfo_({})
    {
        openHtsFile(decompressionThreads);
        loadHeader();
        prepareForStreamingAlignments();
    }
    ~HtsFileStreamer();

    bool trySeekingToNextPrimaryAlignment();

    // Physical BAM ranges include all records, even records discarded by the
    // primary filter. Index offsets are hints until the preceding scan reaches
    // the exact virtual offset at a record boundary.
    std::vector<int64_t> bamRangeStarts(unsigned workers = 1) const;
    std::unique_ptr<HtsFileStreamer> independentReader() const;
    int64_t virtualOffset() const;
    void seekVirtualOffset(int64_t offset);
    bool nextPhysicalRecord();
    struct StreamingRecord
    {
        int32_t tid, mtid;
        hts_pos_t pos, mpos;
        uint32_t flag;
        int32_t length;
        const char* name;
        const uint8_t* sequence;
        const uint8_t* quality;
    };
    // The returned fields borrow the current block or HTS record until the next
    // stream operation. Use this view directly instead of currentAlignment().
    bool nextStreamingRecord(StreamingRecord& record);
    bool currentIsPrimary() const;

    int32_t currentReadContigId() const { return htsAlignmentPtr_->core.tid; }
    hts_pos_t currentReadPosition() const { return htsAlignmentPtr_->core.pos; }
    int32_t currentReadLength() const;
    int32_t currentMateContigId() const { return htsAlignmentPtr_->core.mtid; }
    hts_pos_t currentMatePosition() const { return htsAlignmentPtr_->core.mpos; }

    bool currentIsPaired() const { return (htsAlignmentPtr_->core.flag & BAM_FPAIRED); }

    bool isStreamingAlignedReads() const;

    Read decodeRead() const;
    // Borrowed until the next seek; asynchronous users must copy needed fields.
    const bam1_t* currentAlignment() const { return htsAlignmentPtr_; }

private:
    enum class Status
    {
        kStreamingReads,
        kFinishedStreaming
    };

    void openHtsFile(unsigned decompressionThreads);
    void loadHeader();
    void prepareForStreamingAlignments();

    const std::string htsFilePath_;
    const std::string htsReferencePath_;
    ReferenceContigInfo contigInfo_;
    Status status_ = Status::kStreamingReads;

    htsFile* htsFilePtr_ = nullptr;
    bam1_t* htsAlignmentPtr_ = nullptr;
    bam_hdr_t* htsHeaderPtr_ = nullptr;
    htsThreadPool htsThreadPool_ = { nullptr, 0 };
};

}
}
