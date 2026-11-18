//
// NuExpansionHunter
// Copyright 2026 NuExpansionHunter contributors
//
// This file is part of NuExpansionHunter, a modified version of
// Illumina ExpansionHunter v5.0.0 (Copyright 2016-2019 Illumina, Inc.).
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
// BAM sequence and quality packing (kernel K-31), as pure functions over plain arrays.
//
// This header is valid C and C++. It declares no intrinsics itself: the file
// that includes it provides them (<immintrin.h> in the product, the lane-wise
// models of verification/simd in a proof harness).
#pragma once

#include <ctype.h>
#include <stddef.h>
#include <stdint.h>

// Packs query[index, size) two bases per byte with htslib's nucleotide table nt16 (seq_nt16_table) and writes
// quality 40 for an upper-case base and 0 otherwise. For an odd size the low nibble of the last byte is 0.
static inline void be_encode_tail(
    size_t index, const char* query, size_t size, uint8_t* sequence, uint8_t* quality, const unsigned char* nt16)
{
    for (; index + 1 < size; index += 2)
    {
        const unsigned char first = (unsigned char)query[index];
        const unsigned char second = (unsigned char)query[index + 1];
        sequence[index / 2] = (uint8_t)((nt16[first] << 4) | nt16[second]);
        quality[index] = isupper(first) ? 40 : 0;
        quality[index + 1] = isupper(second) ? 40 : 0;
    }
    if (index < size)
    {
        sequence[index / 2] = (uint8_t)(nt16[(unsigned char)query[index]] << 4);
        quality[index] = isupper((unsigned char)query[index]) ? 40 : 0;
    }
}

// Contract [K-31] reference: the scalar packer from the first base.
static inline void be_encode_reference(
    const char* query, size_t size, uint8_t* sequence, uint8_t* quality, const unsigned char* nt16)
{
    be_encode_tail(0, query, size, sequence, quality, nt16);
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
// Contract [K-31] AVX2 packer.
// Requires: nt16 is htslib's seq_nt16_table, which gives both cases of each letter the same code; the C
//           locale (isupper is true exactly for 'A'..'Z').
// Ensures:  returns an even index <= size; sequence[0, index / 2) and quality[0, index) equal those of
//           be_encode_reference; be_encode_tail from index completes the reference result. No byte outside
//           query[0, size), sequence[0, (size + 1) / 2) or quality[0, size) is touched.
// Reference: be_encode_reference.
static inline size_t __attribute__((target("avx2"))) be_encode_avx2_blocks(
    const char* query, size_t size, uint8_t* sequence, uint8_t* quality, const unsigned char* nt16)
{
    const size_t evenLength = size & ~(size_t)1;
    if (evenLength < 32)
        return 0;
    const __m256i first = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i*)(nt16 + 64)));
    const __m256i second = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i*)(nt16 + 80)));
    size_t index = 0;
    while (index < evenLength)
    {
        // Read the final full vector at an even base offset. Its overlap with
        // the preceding vector writes identical sequence bytes and qualities.
        const size_t start = index < evenLength - 32 ? index : evenLength - 32;
        const __m256i bases = _mm256_loadu_si256((const __m256i*)(query + start));
        const __m256i upper = _mm256_and_si256(bases, _mm256_set1_epi8((char)0xdf));
        const __m256i letters = _mm256_sub_epi8(upper, _mm256_set1_epi8('A'));
        const __m256i valid = _mm256_cmpeq_epi8(letters, _mm256_min_epu8(letters, _mm256_set1_epi8(25)));
        if (_mm256_movemask_epi8(valid) != -1)
            break;
        // htslib assigns the same nucleotide code to either case of each letter.
        const __m256i codes = _mm256_blendv_epi8(_mm256_shuffle_epi8(first, upper),
            _mm256_shuffle_epi8(second, upper), _mm256_cmpgt_epi8(upper, _mm256_set1_epi8('O')));
        const __m256i packedWords = _mm256_maddubs_epi16(codes, _mm256_set1_epi16(0x0110));
        const __m256i packed = _mm256_packus_epi16(packedWords, _mm256_setzero_si256());
        const __m128i bytes = _mm_unpacklo_epi64(_mm256_castsi256_si128(packed), _mm256_extracti128_si256(packed, 1));
        _mm_storeu_si128((__m128i*)(sequence + start / 2), bytes);
        _mm256_storeu_si256(
            (__m256i*)(quality + start), _mm256_and_si256(_mm256_cmpeq_epi8(bases, upper), _mm256_set1_epi8(40)));
        index = start + 32;
    }
    return index;
}
#endif

#if defined(__aarch64__) && defined(__ARM_NEON)
// Contract [K-31] NEON packer: the same contract as be_encode_avx2_blocks, in blocks of 16 bases.
// Reference: be_encode_reference.
static inline size_t be_encode_neon_blocks(
    const char* query, size_t size, uint8_t* sequence, uint8_t* quality, const unsigned char* nt16)
{
    // nt16 codes of '@'..'_' (upper-case letters at 65..90), indexed by upper - 64.
    const uint8x16x2_t table = { { vld1q_u8(nt16 + 64), vld1q_u8(nt16 + 80) } };
    size_t index = 0;
    for (; index + 16 <= size; index += 16)
    {
        const uint8x16_t bases = vld1q_u8((const uint8_t*)(query + index));
        const uint8x16_t upper = vandq_u8(bases, vdupq_n_u8(0xdf));
        const uint8x16_t valid = vcleq_u8(vsubq_u8(upper, vdupq_n_u8('A')), vdupq_n_u8(25));
        if (vminvq_u8(valid) != 255)
            break;
        // htslib assigns the same nucleotide code to either case of each letter.
        const uint8x16_t codes = vqtbl2q_u8(table, vsubq_u8(upper, vdupq_n_u8(64)));
        const uint8x8_t even = vget_low_u8(vuzp1q_u8(codes, codes));
        const uint8x8_t odd = vget_low_u8(vuzp2q_u8(codes, codes));
        vst1_u8(sequence + index / 2, vorr_u8(vshl_n_u8(even, 4), odd));
        vst1q_u8(quality + index, vandq_u8(vceqq_u8(bases, upper), vdupq_n_u8(40)));
    }
    return index;
}
#endif
