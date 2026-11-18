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
// Read decoding (kernel K-22), as pure functions over plain arrays.
//
// This header is valid C and C++. It declares no intrinsics itself: the file
// that includes it provides them (<immintrin.h> in the product, the lane-wise
// models of verification/simd in a proof harness).
#pragma once

#include <stdint.h>

// Contract [K-22] reference: ExpansionHunter v5.0.0's per-base loop. upper is htslib's seq_nt16_str and lower
// its lower-case copy; a base whose quality is at most 20 is written in lower case.
static inline void rd_decode_reference(
    int32_t length, const uint8_t* packed, const uint8_t* qualities, char* bases, const char* upper, const char* lower)
{
    for (int32_t index = 0; index < length; ++index)
    {
        const char* converter = qualities[index] <= 20 ? lower : upper;
        bases[index] = converter[packed[index >> 1] >> ((~index & 1) << 2) & 0xf];
    }
}

// Decodes bases [index, length) two nibbles at a time; the scalar part of the optimized decoder.
static inline void rd_decode_tail(int32_t index, int32_t length, const uint8_t* packed, const uint8_t* qualities,
    char* bases, const char* upper, const char* lower)
{
    // Each BAM byte encodes two bases. Decode both nibbles so the loop does not
    // reload the byte or compute a variable shift for each base. Quality casing
    // remains independent for the two bases.
    for (; index < length - 1; index += 2)
    {
        const uint8_t pair = packed[index / 2];
        const char* first = qualities[index] <= 20 ? lower : upper;
        const char* second = qualities[index + 1] <= 20 ? lower : upper;
        bases[index] = first[pair >> 4];
        bases[index + 1] = second[pair & 15];
    }
    if (index < length)
    {
        const char* converter = qualities[index] <= 20 ? lower : upper;
        bases[index] = converter[packed[index / 2] >> 4];
    }
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
// Contract [K-22] AVX2 blocks of 32 bases.
// Requires: length >= 0; packed[0, (length + 1) / 2) and qualities[0, length) readable; upper is htslib's
//           seq_nt16_str (every symbol an upper-case ASCII letter or '=').
// Ensures:  returns index = 32 * floor(length / 32); bases[0, index) equal rd_decode_reference's; no byte
//           outside packed[0, index / 2) or qualities[0, index) is read. Followed by rd_decode_tail from
//           index, the whole result equals rd_decode_reference.
// Reference: rd_decode_reference.
static inline int32_t __attribute__((target("avx2"))) rd_decode_avx2_blocks(
    int32_t length, const uint8_t* packed, const uint8_t* qualities, char* bases, const char* upper)
{
    const __m128i alphabet128 = _mm_loadu_si128((const __m128i*)upper);
    const __m256i alphabet = _mm256_broadcastsi128_si256(alphabet128);
    const __m128i nibbleMask = _mm_set1_epi8(15);
    const __m256i cutoff = _mm256_set1_epi8(20);
    const __m256i caseBit = _mm256_set1_epi8(0x20);
    const __m256i zero = _mm256_setzero_si256();
    int32_t index = 0;
    for (; index <= length - 32; index += 32)
    {
        const __m128i encoded = _mm_loadu_si128((const __m128i*)(packed + index / 2));
        const __m128i high = _mm_and_si128(_mm_srli_epi16(encoded, 4), nibbleMask);
        const __m128i low = _mm_and_si128(encoded, nibbleMask);
        const __m128i first = _mm_unpacklo_epi8(high, low);
        const __m128i second = _mm_unpackhi_epi8(high, low);
        const __m256i nibbles = _mm256_inserti128_si256(_mm256_castsi128_si256(first), second, 1);
        const __m256i decoded = _mm256_shuffle_epi8(alphabet, nibbles);
        const __m256i quality = _mm256_loadu_si256((const __m256i*)(qualities + index));
        // Saturating unsigned subtraction handles all qualities 0..255,
        // including the BAM absent-quality value 255, without signed compares.
        // Every BAM symbol is an ASCII upper-case letter or '='; setting the
        // lower-case bit preserves '=' and matches the lower-case table.
        const __m256i lower = _mm256_cmpeq_epi8(_mm256_subs_epu8(quality, cutoff), zero);
        _mm256_storeu_si256((__m256i*)(bases + index), _mm256_or_si256(decoded, _mm256_and_si256(lower, caseBit)));
    }
    return index;
}
#endif

#if defined(__aarch64__) && defined(__ARM_NEON)
// Contract [K-22] NEON blocks of 32 bases: the same contract as rd_decode_avx2_blocks.
// Reference: rd_decode_reference.
static inline int32_t rd_decode_neon_blocks(int32_t length, const uint8_t* packed, const uint8_t* qualities,
                                            char* bases, const char* upper)
{
    const uint8x16_t alphabet = vld1q_u8((const uint8_t*)upper);
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
        // Every BAM symbol is an ASCII uppercase letter or '='. ORing the lowercase bit preserves '=' and
        // exactly matches the scalar table.
        const uint8x16_t firstCase = vandq_u8(vcleq_u8(vld1q_u8(qualities + index), cutoff), caseBit);
        const uint8x16_t secondCase = vandq_u8(vcleq_u8(vld1q_u8(qualities + index + 16), cutoff), caseBit);
        vst1q_u8((uint8_t*)(bases + index), vorrq_u8(first, firstCase));
        vst1q_u8((uint8_t*)(bases + index + 16), vorrq_u8(second, secondCase));
    }
    return index;
}
#endif
