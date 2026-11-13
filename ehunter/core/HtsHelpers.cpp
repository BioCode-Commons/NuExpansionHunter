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

#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#endif

#if defined(EH_READ_DECODE_TEST_FORCE_SCALAR) && defined(EH_READ_DECODE_TEST_FORCE_AVX2)
#error "Conflicting read decoder test dispatch overrides"
#endif
#if defined(EH_READ_DECODE_TEST_FORCE_AVX2) && !(defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)))
#error "The AVX2 read decoder test override requires an x86 compiler with target attributes"
#endif

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
#ifdef EH_READ_DECODE_TEST_INSTRUMENT_SIMD
thread_local uint64_t decodedSimdBlocks = 0;
#endif

#if defined(__aarch64__) && defined(__ARM_NEON) && !defined(EH_READ_DECODE_TEST_FORCE_SCALAR)
int32_t decodeReadNeonBlocks(int32_t length, const uint8_t* packed, const uint8_t* qualities, char* bases)
{
    const uint8x16_t alphabet = vld1q_u8(reinterpret_cast<const uint8_t*>(seq_nt16_str));
    const uint8x16_t nibbleMask = vdupq_n_u8(15);
    const uint8x16_t cutoff = vdupq_n_u8(20);
    const uint8x16_t caseBit = vdupq_n_u8(0x20);
    int32_t index = 0;
    for (; index <= length - 32; index += 32)
    {
        const uint8x16_t encoded = vld1q_u8(packed + index / 2);
        const uint8x16_t high = vshrq_n_u8(encoded, 4);
        const uint8x16_t low = vandq_u8(encoded, nibbleMask);
        const uint8x16_t first = vqtbl1q_u8(alphabet, vzip1q_u8(high, low));
        const uint8x16_t second = vqtbl1q_u8(alphabet, vzip2q_u8(high, low));
        // Every BAM symbol is an ASCII uppercase letter or '='. ORing the
        // lowercase bit preserves '=' and exactly matches the scalar table.
        const uint8x16_t firstCase = vandq_u8(vcleq_u8(vld1q_u8(qualities + index), cutoff), caseBit);
        const uint8x16_t secondCase = vandq_u8(vcleq_u8(vld1q_u8(qualities + index + 16), cutoff), caseBit);
        vst1q_u8(reinterpret_cast<uint8_t*>(bases + index), vorrq_u8(first, firstCase));
        vst1q_u8(reinterpret_cast<uint8_t*>(bases + index + 16), vorrq_u8(second, secondCase));
#ifdef EH_READ_DECODE_TEST_INSTRUMENT_SIMD
        ++decodedSimdBlocks;
#endif
    }
    return index;
}
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(EH_READ_DECODE_TEST_FORCE_SCALAR)
bool readDecodeAvx2Supported()
{
#ifdef EH_READ_DECODE_TEST_FORCE_AVX2
    // Test-only override for an explicitly validated ISA/emulation target.
    return true;
#else
    return __builtin_cpu_supports("avx2");
#endif
}

// Only scalar arguments cross the baseline-ISA/AVX2 boundary. No global
// -mavx2 flag or AVX2 instructions are required on unsupported x86 machines.
int32_t __attribute__((target("avx2"), noinline)) decodeReadAvx2Blocks(
    int32_t length, const uint8_t* packed, const uint8_t* qualities, char* bases)
{
    const __m128i alphabet128 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(seq_nt16_str));
    const __m256i alphabet = _mm256_broadcastsi128_si256(alphabet128);
    const __m128i nibbleMask = _mm_set1_epi8(15);
    const __m256i cutoff = _mm256_set1_epi8(20);
    const __m256i caseBit = _mm256_set1_epi8(0x20);
    const __m256i zero = _mm256_setzero_si256();
    int32_t index = 0;
    for (; index <= length - 32; index += 32)
    {
        const __m128i encoded = _mm_loadu_si128(reinterpret_cast<const __m128i*>(packed + index / 2));
        const __m128i high = _mm_and_si128(_mm_srli_epi16(encoded, 4), nibbleMask);
        const __m128i low = _mm_and_si128(encoded, nibbleMask);
        const __m128i first = _mm_unpacklo_epi8(high, low);
        const __m128i second = _mm_unpackhi_epi8(high, low);
        const __m256i nibbles = _mm256_inserti128_si256(_mm256_castsi128_si256(first), second, 1);
        const __m256i decoded = _mm256_shuffle_epi8(alphabet, nibbles);
        const __m256i quality = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(qualities + index));
        // Saturating unsigned subtraction handles all qualities 0..255,
        // including the BAM absent-quality value 255, without signed compares.
        const __m256i lower = _mm256_cmpeq_epi8(_mm256_subs_epu8(quality, cutoff), zero);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(bases + index),
            _mm256_or_si256(decoded, _mm256_and_si256(lower, caseBit)));
#ifdef EH_READ_DECODE_TEST_INSTRUMENT_SIMD
        ++decodedSimdBlocks;
#endif
    }
    return index;
}
#endif
}

#ifdef EH_READ_DECODE_TEST_INSTRUMENT_SIMD
namespace detail
{
uint64_t readDecodeSimdBlocksForTesting() { return decodedSimdBlocks; }
bool readDecodeSimdSupportedForTesting()
{
#if defined(EH_READ_DECODE_TEST_FORCE_SCALAR)
    return false;
#elif defined(__aarch64__) && defined(__ARM_NEON)
    return true;
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    return readDecodeAvx2Supported();
#else
    return false;
#endif
}
}
#endif

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

        int32_t index = 0;
        if (readLength >= 32)
        {
#if defined(__aarch64__) && defined(__ARM_NEON) && !defined(EH_READ_DECODE_TEST_FORCE_SCALAR)
            index = decodeReadNeonBlocks(readLength, htsSeqPtr, htsQualPtr, &bases[0]);
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(EH_READ_DECODE_TEST_FORCE_SCALAR)
            if (readDecodeAvx2Supported())
                index = decodeReadAvx2Blocks(readLength, htsSeqPtr, htsQualPtr, &bases[0]);
#endif
        }
        // Complete vector blocks never read beyond the actual BAM arrays.
        // Each remaining BAM byte encodes two bases. Decode both nibbles so the
        // hot loop does not reload the byte or compute a variable shift for
        // each base. Quality casing remains independent for the two bases.
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
