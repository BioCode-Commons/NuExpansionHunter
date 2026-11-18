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
// Matching base run (kernel K-19), as a pure function over plain arrays.
//
// This header is valid C and C++. It declares no intrinsics itself: the file
// that includes it provides them (<emmintrin.h> in the product, the lane-wise
// models of ehunter/verification/simd in a proof harness).
#pragma once

#include <stddef.h>
#include <stdint.h>

// codes::kReferenceQueryCodeMatchLookupTable[reference code][query code] of graphutils/BaseMatching.hh;
// a unit test checks that the two tables agree.
static const uint8_t bm_code_match[16][5] = {
    { 1, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0 }, { 0, 0, 1, 0, 0 }, { 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0 }, { 0, 1, 1, 1, 0 },
    { 1, 0, 1, 1, 0 }, { 1, 1, 0, 1, 0 }, { 0, 0, 1, 1, 0 }, { 1, 1, 0, 0, 0 }, { 1, 1, 1, 1, 0 }, { 1, 0, 1, 0, 0 },
    { 0, 1, 1, 0, 0 }, { 1, 1, 1, 0, 0 }, { 1, 0, 0, 1, 0 }, { 0, 1, 0, 1, 0 }
};

// Contract [K-19]
// Requires: reference[0, length) and query[0, length) readable; referenceCodes and queryCodes are the
//           256-entry base encoding tables of graphutils/BaseMatching.hh.
// Ensures:  the result equals bm_matching_base_run_reference(...) with the same arguments: the length of the
//           longest prefix (reverse = 0) or suffix (reverse = 1) of matching base pairs.
// Reference: bm_matching_base_run_reference, the per-base loop of v5.0.0's path extension.
static inline size_t bm_matching_base_run_reference(const char* reference, const char* query, size_t length,
    int reverse, const uint8_t* referenceCodes, const uint8_t* queryCodes)
{
    size_t matched = 0;
    while (matched < length)
    {
        const size_t pos = reverse ? length - matched - 1 : matched;
        if (!bm_code_match[referenceCodes[(uint8_t)reference[pos]]][queryCodes[(uint8_t)query[pos]]])
            break;
        ++matched;
    }
    return matched;
}

// Equal canonical bases (either case) have an unambiguous result in the scalar
// lookup table. Every other byte pair, including IUPAC reference symbols, uses that table.
static inline size_t bm_matching_base_run(const char* reference, const char* query, size_t length, int reverse,
    int vector, const uint8_t* referenceCodes, const uint8_t* queryCodes)
{
    size_t matched = 0;
    while (matched < length)
    {
#if defined(__SSE2__) && (defined(__GNUC__) || defined(__clang__))
        if (vector && length - matched >= 16)
        {
            const size_t begin = reverse ? length - matched - 16 : matched;
            const __m128i fold = _mm_set1_epi8(32);
            const __m128i r = _mm_or_si128(fold, _mm_loadu_si128((const __m128i*)(reference + begin)));
            const __m128i q = _mm_or_si128(fold, _mm_loadu_si128((const __m128i*)(query + begin)));
            const __m128i ac = _mm_or_si128(_mm_cmpeq_epi8(q, _mm_set1_epi8('a')), _mm_cmpeq_epi8(q, _mm_set1_epi8('c')));
            const __m128i gt = _mm_or_si128(_mm_cmpeq_epi8(q, _mm_set1_epi8('g')), _mm_cmpeq_epi8(q, _mm_set1_epi8('t')));
            const __m128i plain = _mm_and_si128(_mm_cmpeq_epi8(r, q), _mm_or_si128(ac, gt));
            const unsigned other = (unsigned)_mm_movemask_epi8(plain) ^ 0xffffU;
            if (!other)
            {
                matched += 16;
                continue;
            }
            matched += reverse ? (size_t)(__builtin_clz(other) - 16) : (size_t)__builtin_ctz(other);
        }
#elif defined(__aarch64__) && defined(__ARM_NEON)
        if (vector && length - matched >= 16)
        {
            const size_t begin = reverse ? length - matched - 16 : matched;
            const uint8x16_t fold = vdupq_n_u8(32);
            const uint8x16_t r = vorrq_u8(fold, vld1q_u8((const uint8_t*)(reference + begin)));
            const uint8x16_t q = vorrq_u8(fold, vld1q_u8((const uint8_t*)(query + begin)));
            const uint8x16_t ac = vorrq_u8(vceqq_u8(q, vdupq_n_u8('a')), vceqq_u8(q, vdupq_n_u8('c')));
            const uint8x16_t gt = vorrq_u8(vceqq_u8(q, vdupq_n_u8('g')), vceqq_u8(q, vdupq_n_u8('t')));
            const uint8x16_t plain = vandq_u8(vceqq_u8(r, q), vorrq_u8(ac, gt));
            // Four mask bits per lane, lane 0 in the lowest nibble.
            const uint64_t other
                = ~vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(plain), 4)), 0);
            if (!other)
            {
                matched += 16;
                continue;
            }
            matched += reverse ? (size_t)(__builtin_clzll(other) / 4) : (size_t)(__builtin_ctzll(other) / 4);
        }
#else
        (void)vector;
#endif
        const size_t pos = reverse ? length - matched - 1 : matched;
        if (!bm_code_match[referenceCodes[(uint8_t)reference[pos]]][queryCodes[(uint8_t)query[pos]]])
            break;
        ++matched;
    }
    return matched;
}
