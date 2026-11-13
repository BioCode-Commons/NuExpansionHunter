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

#include "core/HtsHelpers.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>

#include "spdlog/spdlog.h"

using std::pair;
using std::string;
using std::vector;

namespace spd = spdlog;

namespace ehunter
{
namespace htshelpers
{

string decodeQuals(bam1_t* htsAlignPtr)
{
    string quals;
    uint8_t* htsQualsPtr = bam_get_qual(htsAlignPtr);
    const int readLength = htsAlignPtr->core.l_qseq;
    quals.resize(readLength);

    for (int index = 0; index < readLength; ++index)
    {
        quals[index] = static_cast<char>(33 + htsQualsPtr[index]);
    }

    return quals;
}

string decodeBases(bam1_t* htsAlignPtr)
{
    string bases;
    uint8_t* htsSeqPtr = bam_get_seq(htsAlignPtr);
    const int32_t readLength = htsAlignPtr->core.l_qseq;
    bases.resize(readLength);

    for (int32_t index = 0; index < readLength; ++index)
    {
        bases[index] = seq_nt16_str[bam_seqi(htsSeqPtr, index)];
    }

    return bases;
}

LinearAlignmentStats decodeAlignmentStats(bam1_t* htsAlignPtr)
{
    LinearAlignmentStats alignmentStats;
    alignmentStats.chromId = htsAlignPtr->core.tid;
    alignmentStats.pos = htsAlignPtr->core.pos;
    alignmentStats.mapq = htsAlignPtr->core.qual;
    alignmentStats.mateChromId = htsAlignPtr->core.mtid;
    alignmentStats.matePos = htsAlignPtr->core.mpos;

    uint32_t samFlag = htsAlignPtr->core.flag;
    alignmentStats.isPaired = samFlag & BAM_FPAIRED;
    alignmentStats.isMapped = !(samFlag & BAM_FUNMAP);
    alignmentStats.isMateMapped = !(samFlag & BAM_FMUNMAP);

    return alignmentStats;
}

bool isPrimaryAlignment(bam1_t* htsAlignPtr)
{
    return !((htsAlignPtr->core.flag & BAM_FSECONDARY) || (htsAlignPtr->core.flag & BAM_FSUPPLEMENTARY));
}

/// Lower-case version of htslib seq_nt16_str table;
class Cache_seq_nt16_str_lc
{
public:
    Cache_seq_nt16_str_lc()
    {
        for (unsigned i(0); i < 16; ++i)
        {
            data[i] = std::tolower(seq_nt16_str[i]);
        }
    }

    char data[16];
};

static Cache_seq_nt16_str_lc seq_nt16_str_lc;

Read decodeRead(bam1_t* htsAlignPtr)
{
    return decodeReadParts(bam_get_qname(htsAlignPtr), htsAlignPtr->core.flag, htsAlignPtr->core.l_qseq,
        bam_get_seq(htsAlignPtr), bam_get_qual(htsAlignPtr));
}

Read decodeReadParts(const char* qname, uint32_t samFlag, int32_t readLength,
    const uint8_t* htsSeqPtr, const uint8_t* htsQualPtr)
{
    const bool isFirstMate = samFlag & BAM_FREAD1;
    const bool isReversed = samFlag & BAM_FREVERSE;

    MateNumber mateNumber = isFirstMate ? MateNumber::kFirstMate : MateNumber::kSecondMate;
    ReadId readId(qname, mateNumber);

    string bases;

    {
        // Decode bases and convert low-quality bases to lowercase:
        static const uint8_t lowBaseQualityCutoff(20);

        bases.resize(readLength);

        // Each BAM byte encodes two bases. Decode both nibbles together so the
        // hot loop does not reload the byte or compute a variable shift for
        // each base. Quality casing remains independent for the two bases.
        int32_t index = 0;
        for (; index < readLength - 1; index += 2)
        {
            const uint8_t packed = htsSeqPtr[index / 2];
            const char* first((htsQualPtr[index] <= lowBaseQualityCutoff) ? seq_nt16_str_lc.data : seq_nt16_str);
            const char* second((htsQualPtr[index + 1] <= lowBaseQualityCutoff) ? seq_nt16_str_lc.data : seq_nt16_str);
            bases[index] = first[packed >> 4];
            bases[index + 1] = second[packed & 15];
        }
        if (index < readLength)
        {
            const char* converter((htsQualPtr[index] <= lowBaseQualityCutoff) ? seq_nt16_str_lc.data : seq_nt16_str);
            bases[index] = converter[htsSeqPtr[index / 2] >> 4];
        }
    }

    return { std::move(readId), std::move(bases), isReversed };
}

ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr)
{
    vector<pair<string, int64_t>> contigNamesAndSizes;
    const int32_t numContigs = htsHeaderPtr->n_targets;
    contigNamesAndSizes.reserve(numContigs);

    for (int32_t contigIndex = 0; contigIndex != numContigs; ++contigIndex)
    {
        const string contig = htsHeaderPtr->target_name[contigIndex];
        int64_t size = htsHeaderPtr->target_len[contigIndex];
        contigNamesAndSizes.push_back(std::make_pair(contig, size));
    }

    return ReferenceContigInfo(contigNamesAndSizes);
}

} // namespace htshelpers
}
