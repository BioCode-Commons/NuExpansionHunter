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
// This file differs from ExpansionHunter v5.0.0: read decoding exposes packed-array inputs and SIMD decoding with a scalar reference.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "core/HtsHelpers.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#endif

#include "core/ReadDecodeKernels.hh"
#include "core/VerifyPolicy.hh"
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

namespace
{
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
bool readDecodeAvx2Supported()
{
    if (verify::kNoAvx2)
        return false;
    static const bool supported = __builtin_cpu_supports("avx2");
    return supported;
}
#endif
}


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

    return { std::move(readId), decodeReadBases(readLength, htsSeqPtr, htsQualPtr), isReversed };
}

// Contract [K-22]
// Requires: readLength >= 0; htsSeqPtr holds (readLength + 1) / 2 packed bytes and htsQualPtr readLength qualities.
// Ensures:  the result equals decodeReadBasesReference(readLength, htsSeqPtr, htsQualPtr); no byte beyond
//           htsSeqPtr[(readLength + 1) / 2) or htsQualPtr[readLength) is read.
// Reference: decodeReadBasesReference, ExpansionHunter v5.0.0's per-base loop in decodeRead.
string decodeReadBasesReference(int32_t readLength, const uint8_t* htsSeqPtr, const uint8_t* htsQualPtr)
{
    string bases;
    bases.resize(readLength);
    rd_decode_reference(readLength, htsSeqPtr, htsQualPtr, &bases[0], seq_nt16_str, seq_nt16_str_lc.data);
    return bases;
}

static string decodeReadBasesOptimized(int32_t readLength, const uint8_t* htsSeqPtr, const uint8_t* htsQualPtr)
{
    string bases;
    bases.resize(readLength);
    int32_t index = 0;
    if (readLength >= 32 && !verify::kScalarOnly)
    {
#if defined(__aarch64__) && defined(__ARM_NEON)
        index = rd_decode_neon_blocks(readLength, htsSeqPtr, htsQualPtr, &bases[0], seq_nt16_str);
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
        if (readDecodeAvx2Supported())
            index = rd_decode_avx2_blocks(readLength, htsSeqPtr, htsQualPtr, &bases[0], seq_nt16_str);
#endif
    }
    // Complete vector blocks never read beyond the actual BAM arrays.
    rd_decode_tail(index, readLength, htsSeqPtr, htsQualPtr, &bases[0], seq_nt16_str, seq_nt16_str_lc.data);
    return bases;
}

string decodeReadBases(int32_t readLength, const uint8_t* htsSeqPtr, const uint8_t* htsQualPtr)
{
    string bases = decodeReadBasesOptimized(readLength, htsSeqPtr, htsQualPtr);
    if (verify::kShadow && bases != decodeReadBasesReference(readLength, htsSeqPtr, htsQualPtr))
        throw std::logic_error("Read decoding differs from the per-base reference");
    return bases;
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
