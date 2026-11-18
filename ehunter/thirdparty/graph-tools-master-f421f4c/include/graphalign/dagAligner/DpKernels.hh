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
// Vector kernels of the affine DP, as pure functions over plain arrays.
//
// This header is valid C and C++. It declares no intrinsics itself: the file
// that includes it provides them (<immintrin.h> in the product, the lane-wise
// models of ehunter/verification/simd in a proof harness), so that the proof
// harnesses check exactly the code the product runs. Each kernel states its
// contract next to the reference it must equal.
#pragma once

#include <stdint.h>

#include "graphutils/BaseMatchingKernels.hh"

#define DP_SCORE_MIN ((int16_t)-32768)

static inline int16_t dp_max16(int16_t a, int16_t b) { return a < b ? b : a; }

// Contract [K-01] reference: one cell of the insertion recurrence of ExpansionHunter v5.0.0
// (recomputeForInsertion), with v5's int16 narrowing of every candidate. Returns the trace bit
// (v == f, the insertion predicate when the gap-opening sum is representable).
static inline int dp_insertion_cell(int16_t* v, int16_t* previousF, int16_t* previousV, int16_t gapOpen,
                                    int16_t gapExt)
{
    const int16_t extension = (int16_t)(*previousF + gapExt);
    const int16_t opening = (int16_t)(*previousV + gapOpen + gapExt);
    const int16_t currentF = dp_max16(extension, opening);
    const int16_t currentV = dp_max16(*v, currentF);
    *v = currentV;
    *previousF = currentF;
    *previousV = currentV;
    return currentV == currentF;
}

// Contract [K-03]
// Requires: nothing (all int16 inputs, including sums that wrap in 16 bits; any gap scores).
// Ensures:  for a row with exactly one predecessor p, dp_fused_lane gives the same E and V cell as the three
//           general passes dp_deletion_lane, dp_align_lane and dp_consolidate_lane applied to cells that
//           start at DP_SCORE_MIN (E, G and V of the current row).
// Reference: dp_deletion_lane + dp_align_lane + dp_consolidate_lane (the general row code of v5.0.0).
static inline void dp_fused_lane(int16_t ep, int16_t vp, int16_t vpPrevious, int16_t penalty, int16_t gapOpen,
                                 int16_t gapExt, int16_t* e, int16_t* v)
{
    // Narrow each candidate before taking max, just as in the general
    // recurrence (including scores near the 16-bit limits).
    const int16_t extension = (int16_t)(ep + gapExt);
    const int16_t opening = (int16_t)(vp + gapOpen + gapExt);
    const int16_t deletion = dp_max16(extension, opening);
    const int16_t match = (int16_t)(vpPrevious + penalty);
    *e = deletion;
    *v = dp_max16(deletion, match);
}

// One predecessor's contribution to the E cell (v5.0.0 recomputeForDeletion).
static inline void dp_deletion_lane(int16_t ep, int16_t vp, int16_t gapOpen, int16_t gapExt, int16_t* e)
{
    const int16_t fromE = (int16_t)(ep + gapExt);
    const int16_t fromV = (int16_t)(vp + gapOpen + gapExt);
    const int16_t candidate = fromE > fromV ? fromE : fromV;
    *e = *e > candidate ? *e : candidate;
}

// One predecessor's contribution to the G cell (v5.0.0 recomputeForAlign).
static inline void dp_align_lane(int16_t vpPrevious, int16_t penalty, int16_t* g)
{
    const int16_t candidate = (int16_t)(vpPrevious + penalty);
    *g = *g > candidate ? *g : candidate;
}

// V = max(V, G, E) (v5.0.0 consolidate).
static inline void dp_consolidate_lane(int16_t g, int16_t e, int16_t* v)
{
    const int16_t candidate = g > e ? g : e;
    *v = *v > candidate ? *v : candidate;
}

// Row passes over n lanes (a multiple of 16), built from the lane functions above. vp is the predecessor
// row at query offset 0; lane i also reads vp[i - 1], the cell to its left (the row's boundary cell for i = 0).
//
// Contract [K-03] rows
// Requires: n >= 0 is a multiple of 16; the arrays written (e, v or g, n lanes each) do not overlap the arrays
//           read, which holds for a current row t and a predecessor row p != t; penalties has n values.
// Ensures:  every SIMD row kernel below (SSE2, AVX2, NEON) writes, lane by lane, exactly what the
//           matching reference row writes, for all int16 inputs and gap scores.
// Reference: dp_fused_row_reference, dp_deletion_row_reference, dp_align_row_reference and
//            dp_consolidate_row_reference, i.e. the lane functions applied to each lane.
static inline void dp_fused_row_reference(const int16_t* ep, const int16_t* vp, const int16_t* penalties,
                                          int16_t gapOpen, int16_t gapExt, int16_t* e, int16_t* v, int n)
{
    for (int i = 0; i < n; ++i)
        dp_fused_lane(ep[i], vp[i], vp[i - 1], penalties[i], gapOpen, gapExt, &e[i], &v[i]);
}

static inline void dp_deletion_row_reference(const int16_t* ep, const int16_t* vp, int16_t gapOpen, int16_t gapExt,
                                             int16_t* e, int n)
{
    for (int i = 0; i < n; ++i)
        dp_deletion_lane(ep[i], vp[i], gapOpen, gapExt, &e[i]);
}

static inline void dp_align_row_reference(const int16_t* vp, const int16_t* penalties, int16_t* g, int n)
{
    for (int i = 0; i < n; ++i)
        dp_align_lane(vp[i - 1], penalties[i], &g[i]);
}

static inline void dp_consolidate_row_reference(const int16_t* g, const int16_t* e, int16_t* v, int n)
{
    for (int i = 0; i < n; ++i)
        dp_consolidate_lane(g[i], e[i], &v[i]);
}

// Contract [K-09] reference: position of the only mismatch of query[0, n) against target[0, n), n if there is
// none, -1 if there are two or more. Codes: target 0..15, query 0..4.
static inline int dp_linear_first_mismatch_reference(const uint8_t* query, const uint8_t* target, int n)
{
    int firstMismatch = n;
    for (int i = 0; i < n; ++i)
        if (!bm_code_match[target[i]][query[i]])
        {
            if (firstMismatch != n)
                return -1;
            firstMismatch = i;
        }
    return firstMismatch;
}

// K-03 row kernels. Every lane operation is a wrapping 16-bit add or a signed 16-bit max, which is what the
// reference does after its int16 narrowing; (int16_t)(vp + gapOpen + gapExt) equals vp plus the wrapped sum
// (int16_t)(gapOpen + gapExt) in 16-bit arithmetic.
#define DP_ROW_KERNELS(SUFFIX, REST, ATTRIBUTE, VEC, LOAD, STORE, ADD, MAX, SET1, WIDTH)                                  \
    static inline void ATTRIBUTE dp_fused_row_##SUFFIX(const int16_t* ep, const int16_t* vp,                         \
        const int16_t* penalties, int16_t gapOpen, int16_t gapExt, int16_t* e, int16_t* v, int n)                    \
    {                                                                                                                  \
        const VEC extend = SET1(gapExt), opening = SET1((int16_t)(gapOpen + gapExt));                                    \
        int i = 0;                                                                                                     \
        for (; n - i >= WIDTH; i += WIDTH)                                                                             \
        {                                                                                                              \
            const VEC deletion = MAX(ADD(LOAD(ep + i), extend), ADD(LOAD(vp + i), opening));                            \
            const VEC match = ADD(LOAD(vp + i - 1), LOAD(penalties + i));                                             \
            STORE(e + i, deletion);                                                                                    \
            STORE(v + i, MAX(deletion, match));                                                                        \
        }                                                                                                              \
        dp_fused_row_##REST(ep + i, vp + i, penalties + i, gapOpen, gapExt, e + i, v + i, n - i);                 \
    }                                                                                                                  \
    static inline void ATTRIBUTE dp_deletion_row_##SUFFIX(const int16_t* ep, const int16_t* vp, int16_t gapOpen,     \
        int16_t gapExt, int16_t* e, int n)                                                                            \
    {                                                                                                                  \
        const VEC extend = SET1(gapExt), opening = SET1((int16_t)(gapOpen + gapExt));                                    \
        int i = 0;                                                                                                     \
        for (; n - i >= WIDTH; i += WIDTH)                                                                             \
            STORE(e + i, MAX(LOAD(e + i), MAX(ADD(LOAD(ep + i), extend), ADD(LOAD(vp + i), opening))));                 \
        dp_deletion_row_##REST(ep + i, vp + i, gapOpen, gapExt, e + i, n - i);                                     \
    }                                                                                                                  \
    static inline void ATTRIBUTE dp_align_row_##SUFFIX(const int16_t* vp, const int16_t* penalties, int16_t* g,     \
        int n)                                                                                                         \
    {                                                                                                                  \
        int i = 0;                                                                                                     \
        for (; n - i >= WIDTH; i += WIDTH)                                                                             \
            STORE(g + i, MAX(LOAD(g + i), ADD(LOAD(vp + i - 1), LOAD(penalties + i))));                              \
        dp_align_row_##REST(vp + i, penalties + i, g + i, n - i);                                                  \
    }                                                                                                                  \
    static inline void ATTRIBUTE dp_consolidate_row_##SUFFIX(const int16_t* g, const int16_t* e, int16_t* v, int n) \
    {                                                                                                                  \
        int i = 0;                                                                                                     \
        for (; n - i >= WIDTH; i += WIDTH)                                                                             \
            STORE(v + i, MAX(LOAD(v + i), MAX(LOAD(g + i), LOAD(e + i))));                                           \
        dp_consolidate_row_##REST(g + i, e + i, v + i, n - i);                                                     \
    }

#define DP_SSE2_LOAD(p) _mm_loadu_si128((const __m128i*)(p))
#define DP_SSE2_STORE(p, x) _mm_storeu_si128((__m128i*)(p), x)
#define DP_AVX2_LOAD(p) _mm256_loadu_si256((const __m256i*)(p))
#define DP_AVX2_STORE(p, x) _mm256_storeu_si256((__m256i*)(p), x)
#define DP_NO_ATTRIBUTE

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))

// Contract [K-01] AVX2 block, 16 cells.
// Requires: gap open -8, gap extension -2; on entry to a block, previousF, previousV and every v[i]
//           are >= DP_SCORE_MIN + 40 (otherwise the block is not taken and the caller continues with
//           the scalar cells).
// Ensures:  for every block taken, v[i..i+16), previousF, previousV and trace bits (i - start) .. +15 equal
//           16 applications of dp_insertion_cell(.., -8, -2).
// Reference: dp_insertion_cell.
static inline int __attribute__((target("avx2"))) dp_insertion_avx2(
    int16_t* v, int i, int end, int16_t* previousF, int16_t* previousV, uint64_t* trace, int start)
{
    // Across sixteen cells the earliest opening can lose 10 + 15*2
    // points. This guard also protects every dominated scalar opening,
    // so int16 narrowing cannot change which candidate wins.
    const int16_t safeMinimum = DP_SCORE_MIN + 40;
    const __m256i safeMinimumVector = _mm256_set1_epi16(safeMinimum);
    const __m256i bottom = _mm256_set1_epi16(DP_SCORE_MIN);
    for (; end - i >= 16; i += 16)
    {
        const __m256i base = _mm256_loadu_si256((const __m256i*)(v + i));
        if (*previousF < safeMinimum || *previousV < safeMinimum
            || _mm256_movemask_epi8(_mm256_cmpgt_epi16(safeMinimumVector, base)))
        {
            break;
        }
        // permute + alignr shifts across the 128-bit lane boundary.
        // The vacated positions contain negative infinity, not zero.
        __m256i insertion = _mm256_add_epi16(
            _mm256_insert_epi16(
                _mm256_alignr_epi8(base, _mm256_permute2x128_si256(bottom, base, 0x20), 14), *previousV, 0),
            _mm256_set1_epi16(-10));
        insertion = _mm256_insert_epi16(
            insertion, dp_max16((int16_t)(*previousF - 2), (int16_t)(*previousV - 10)), 0);

        __m256i shifted = _mm256_adds_epi16(insertion, _mm256_set1_epi16(-2));
        insertion = _mm256_max_epi16(insertion,
            _mm256_alignr_epi8(shifted, _mm256_permute2x128_si256(bottom, shifted, 0x20), 14));
        shifted = _mm256_adds_epi16(insertion, _mm256_set1_epi16(-4));
        insertion = _mm256_max_epi16(insertion,
            _mm256_alignr_epi8(shifted, _mm256_permute2x128_si256(bottom, shifted, 0x20), 12));
        shifted = _mm256_adds_epi16(insertion, _mm256_set1_epi16(-8));
        insertion = _mm256_max_epi16(insertion,
            _mm256_alignr_epi8(shifted, _mm256_permute2x128_si256(bottom, shifted, 0x20), 8));
        shifted = _mm256_adds_epi16(insertion, _mm256_set1_epi16(-16));
        insertion = _mm256_max_epi16(insertion, _mm256_permute2x128_si256(bottom, shifted, 0x20));

        const __m256i current = _mm256_max_epi16(base, insertion);
        _mm256_storeu_si256((__m256i*)(v + i), current);
        const __m256i equal = _mm256_cmpeq_epi16(current, insertion);
        const __m128i traceBytes
            = _mm_packs_epi16(_mm256_castsi256_si128(equal), _mm256_extracti128_si256(equal, 1));
        const uint64_t sixteenBits = (uint64_t)(unsigned)_mm_movemask_epi8(traceBytes);
        *trace |= sixteenBits << (i - start);
        *previousF = (int16_t)_mm256_extract_epi16(insertion, 15);
        *previousV = (int16_t)_mm256_extract_epi16(current, 15);
    }
    return i;
}

static inline __m512i __attribute__((target("avx512f,avx512bw"), always_inline))
dp_extend_insertions_avx512(__m512i insertion, __m512i indices, int distance)
{
    const __m512i earlier
        = _mm512_permutexvar_epi16(_mm512_sub_epi16(indices, _mm512_set1_epi16((short)distance)), insertion);
    const __m512i extended = _mm512_subs_epi16(earlier, _mm512_set1_epi16((short)(2 * distance)));
    return _mm512_mask_max_epi16(insertion, (__mmask32)(0xffffffffu << distance), insertion, extended);
}

// Contract [K-01] AVX-512BW block, 32 cells.
// Requires: gap open -8, gap extension -2; on entry to a block, previousF, previousV and every v[i]
//           are >= DP_SCORE_MIN + 72 (otherwise the block is not taken).
// Ensures:  for every block taken, v[i..i+32), previousF, previousV and trace bits (i - start) .. +31 equal
//           32 applications of dp_insertion_cell(.., -8, -2).
// Reference: dp_insertion_cell.
static inline int __attribute__((target("avx512f,avx512bw"))) dp_insertion_avx512(
    int16_t* v, int i, int end, int16_t* previousF, int16_t* previousV, uint64_t* trace, int start)
{
    const __m512i indices = _mm512_set_epi16(31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    // At most 10 opening points and 31 two-point extensions are lost.
    // Leave rows close to the scalar narrowing boundary to the fallback.
    const int16_t safeMinimum = DP_SCORE_MIN + 72;
    for (; end - i >= 32; i += 32)
    {
        const __m512i base = _mm512_loadu_si512(v + i);
        if (*previousF < safeMinimum || *previousV < safeMinimum
            || _mm512_cmplt_epi16_mask(base, _mm512_set1_epi16(safeMinimum)))
            break;
        const __m512i preceding = _mm512_permutexvar_epi16(_mm512_sub_epi16(indices, _mm512_set1_epi16(1)), base);
        __m512i insertion = _mm512_sub_epi16(preceding, _mm512_set1_epi16(10));
        const int16_t first = dp_max16((int16_t)(*previousF - 2), (int16_t)(*previousV - 10));
        insertion = _mm512_mask_mov_epi16(insertion, 1, _mm512_set1_epi16(first));
        insertion = dp_extend_insertions_avx512(insertion, indices, 1);
        insertion = dp_extend_insertions_avx512(insertion, indices, 2);
        insertion = dp_extend_insertions_avx512(insertion, indices, 4);
        insertion = dp_extend_insertions_avx512(insertion, indices, 8);
        insertion = dp_extend_insertions_avx512(insertion, indices, 16);
        const __m512i current = _mm512_max_epi16(base, insertion);
        _mm512_storeu_si512(v + i, current);
        *trace |= (uint64_t)_mm512_cmpeq_epi16_mask(current, insertion) << (i - start);
        *previousF = (int16_t)_mm256_extract_epi16(_mm512_extracti64x4_epi64(insertion, 1), 15);
        *previousV = (int16_t)_mm256_extract_epi16(_mm512_extracti64x4_epi64(current, 1), 15);
    }
    return i;
}



// SSE2 is part of x86-64, so this tier needs no run-time check. There is no AVX-512 tier: with rows of a few
// hundred cells, 512-bit row passes were slower than 256-bit ones on Zen 4.
DP_ROW_KERNELS(sse2, reference, DP_NO_ATTRIBUTE, __m128i, DP_SSE2_LOAD, DP_SSE2_STORE, _mm_add_epi16, _mm_max_epi16,
               _mm_set1_epi16, 8)
DP_ROW_KERNELS(avx2, reference, __attribute__((target("avx2"))), __m256i, DP_AVX2_LOAD, DP_AVX2_STORE, _mm256_add_epi16,
               _mm256_max_epi16, _mm256_set1_epi16, 16)

// Contract [K-09] AVX2 scan.
// Requires: target codes 0..15, query codes 0..4, n >= 0.
// Ensures:  the result equals dp_linear_first_mismatch_reference(query, target, n).
// Reference: dp_linear_first_mismatch_reference.
static inline int __attribute__((target("avx2"))) dp_linear_first_mismatch_avx2(
    const uint8_t* query, const uint8_t* target, int queryLength)
{
    // Bit sets encode exactly the reference/query IUPAC table.
    const uint8_t targetMasks[16] = { 1, 2, 4, 8, 0, 14, 13, 11, 12, 3, 15, 5, 6, 7, 9, 10 };
    const uint8_t queryMasks[16] = { 1, 2, 4, 8, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const __m256i targetTable = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i*)targetMasks));
    const __m256i queryTable = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i*)queryMasks));
    int firstMismatch = queryLength;
    int i = 0;
    for (; queryLength - i >= 32; i += 32)
    {
        const __m256i targetBits = _mm256_shuffle_epi8(targetTable, _mm256_loadu_si256((const __m256i*)(target + i)));
        const __m256i queryBits = _mm256_shuffle_epi8(queryTable, _mm256_loadu_si256((const __m256i*)(query + i)));
        const uint32_t mismatch = (uint32_t)_mm256_movemask_epi8(
            _mm256_cmpeq_epi8(_mm256_and_si256(targetBits, queryBits), _mm256_setzero_si256()));
        if (mismatch)
        {
            if ((mismatch & (mismatch - 1)) || firstMismatch != queryLength)
                return -1;
            firstMismatch = i + __builtin_ctz(mismatch);
        }
    }
    for (; i < queryLength; ++i)
        if (!bm_code_match[target[i]][query[i]])
        {
            if (firstMismatch != queryLength)
                return -1;
            firstMismatch = i;
        }
    return firstMismatch;
}

#endif

#if defined(__aarch64__) && defined(__ARM_NEON)
// Contract [K-01] NEON block, 8 cells.
// Requires: gap open -8, gap extension -2; on entry to a block, previousF, previousV and every v[i]
//           are >= DP_SCORE_MIN + 24 (otherwise the block is not taken and the caller continues with
//           the scalar cells).
// Ensures:  for every block taken, v[i..i+8), previousF, previousV and trace bits (i - start) .. +7 equal
//           8 applications of dp_insertion_cell(.., -8, -2).
// Reference: dp_insertion_cell.
static inline int dp_insertion_neon(int16_t* v, int i, int end, int16_t* previousF, int16_t* previousV,
                                    uint64_t* trace, int start)
{
    // With nonpositive opening scores and no narrowing overflow, opening from the preceding F cannot beat
    // extending it. The eight-cell recurrence is therefore a max-plus prefix scan over the unmodified V
    // cells. A 24-point margin covers every scalar candidate in this block, including the dominated gap
    // openings.
    const int16_t safeMinimum = DP_SCORE_MIN + 24;
    const int16x8_t bottom = vdupq_n_s16(DP_SCORE_MIN);
    const uint8x8_t traceWeights = { 1, 2, 4, 8, 16, 32, 64, 128 };
    for (; end - i >= 8; i += 8)
    {
        const int16x8_t base = vld1q_s16(v + i);
        // Preserve v5.0.0's modular int16 behavior when the prefix-scan algebra is not valid.
        if (*previousF < safeMinimum || *previousV < safeMinimum || vminvq_s16(base) < safeMinimum)
            break;
        int16x8_t insertion = vaddq_s16(vextq_s16(vdupq_n_s16(*previousV), base, 7), vdupq_n_s16(-10));
        insertion = vsetq_lane_s16(dp_max16((int16_t)(*previousF - 2), (int16_t)(*previousV - 10)), insertion, 0);
        // Saturation applies only to negative-infinity padding lanes; the range guard protects scores.
        insertion = vmaxq_s16(insertion, vextq_s16(bottom, vqaddq_s16(insertion, vdupq_n_s16(-2)), 7));
        insertion = vmaxq_s16(insertion, vextq_s16(bottom, vqaddq_s16(insertion, vdupq_n_s16(-4)), 6));
        insertion = vmaxq_s16(insertion, vextq_s16(bottom, vqaddq_s16(insertion, vdupq_n_s16(-8)), 4));
        const int16x8_t current = vmaxq_s16(base, insertion);
        vst1q_s16(v + i, current);
        const uint8x8_t traceLanes = vmovn_u16(vceqq_s16(current, insertion));
        const uint64_t eightBits = vaddv_u8(vand_u8(traceLanes, traceWeights));
        *trace |= eightBits << (i - start);
        *previousF = vgetq_lane_s16(insertion, 7);
        *previousV = vgetq_lane_s16(current, 7);
    }
    return i;
}

// Contract [K-09] NEON scan.
// Requires: target codes 0..15, query codes 0..4, n >= 0.
// Ensures:  the result equals dp_linear_first_mismatch_reference(query, target, n).
// Reference: dp_linear_first_mismatch_reference.
static inline int dp_linear_first_mismatch_neon(const uint8_t* query, const uint8_t* target, int queryLength)
{
    // Bit sets encode exactly the reference/query IUPAC table.
    const uint8_t targetMasks[16] = { 1, 2, 4, 8, 0, 14, 13, 11, 12, 3, 15, 5, 6, 7, 9, 10 };
    const uint8_t queryMasks[16] = { 1, 2, 4, 8, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const uint8_t lanes[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    const uint8x16_t targetTable = vld1q_u8(targetMasks), queryTable = vld1q_u8(queryMasks);
    int firstMismatch = queryLength;
    int i = 0;
    for (; queryLength - i >= 16; i += 16)
    {
        const uint8x16_t targetBits = vqtbl1q_u8(targetTable, vld1q_u8(target + i));
        const uint8x16_t queryBits = vqtbl1q_u8(queryTable, vld1q_u8(query + i));
        const uint8x16_t mismatch = vceqq_u8(vandq_u8(targetBits, queryBits), vdupq_n_u8(0));
        const unsigned count = vaddvq_u8(vshrq_n_u8(mismatch, 7));
        if (count)
        {
            if (count > 1 || firstMismatch != queryLength)
                return -1;
            firstMismatch = i + vminvq_u8(vbslq_u8(mismatch, vld1q_u8(lanes), vdupq_n_u8(255)));
        }
    }
    for (; i < queryLength; ++i)
        if (!bm_code_match[target[i]][query[i]])
        {
            if (firstMismatch != queryLength)
                return -1;
            firstMismatch = i;
        }
    return firstMismatch;
}

#define DP_NEON_LOAD(p) vld1q_s16(p)
#define DP_NEON_STORE(p, x) vst1q_s16(p, x)
// NEON is part of AArch64, so this tier is selected at compile time.
DP_ROW_KERNELS(neon, reference, DP_NO_ATTRIBUTE, int16x8_t, DP_NEON_LOAD, DP_NEON_STORE, vaddq_s16, vmaxq_s16,
               vdupq_n_s16, 8)
#endif
