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
// Lane-wise C models of the x86 SSE2/SSSE3/AVX2/AVX-512BW intrinsics that the proved kernels use.
//
// A proof harness includes this header in place of <immintrin.h>, then the product's kernel header, so the
// model checker analyses exactly the product's kernel code. Every model is checked against the hardware
// intrinsic by simd/crosscheck_x86.cpp (all models, random and boundary inputs); a proof is only as good as
// that check. Vectors are little-endian byte arrays, as in the hardware registers.
#pragma once

#include <stdint.h>
#include <string.h>

typedef struct { uint8_t b[16]; } __m128i;
typedef struct { uint8_t b[32]; } __m256i;
typedef struct { uint8_t b[64]; } __m512i;
typedef uint32_t __mmask32;

#define EH_SIMD_MODEL static inline

// Element access
EH_SIMD_MODEL int16_t m_get16(const uint8_t* b, int lane) { return (int16_t)(uint16_t)(b[2 * lane] | (b[2 * lane + 1] << 8)); }
EH_SIMD_MODEL void m_set16(uint8_t* b, int lane, int16_t value)
{
    b[2 * lane] = (uint8_t)((uint16_t)value & 0xff);
    b[2 * lane + 1] = (uint8_t)((uint16_t)value >> 8);
}
EH_SIMD_MODEL int16_t m_sat16(int32_t v) { return (int16_t)(v > 32767 ? 32767 : v < -32768 ? -32768 : v); }
EH_SIMD_MODEL uint8_t m_satu8_from16(int16_t v) { return (uint8_t)(v < 0 ? 0 : v > 255 ? 255 : v); }
EH_SIMD_MODEL int8_t m_sat8_from16(int16_t v) { return (int8_t)(v < -128 ? -128 : v > 127 ? 127 : v); }

// ---------------------------------------------------------------- 128-bit
EH_SIMD_MODEL __m128i _mm_loadu_si128(const __m128i* p) { __m128i r; memcpy(r.b, p, 16); return r; }
EH_SIMD_MODEL void _mm_storeu_si128(__m128i* p, __m128i a) { memcpy(p, a.b, 16); }
EH_SIMD_MODEL __m128i _mm_set1_epi8(char v) { __m128i r; for (int i = 0; i < 16; ++i) r.b[i] = (uint8_t)v; return r; }
EH_SIMD_MODEL __m128i _mm_set1_epi16(short v) { __m128i r; for (int i = 0; i < 8; ++i) m_set16(r.b, i, v); return r; }
EH_SIMD_MODEL __m128i _mm_add_epi16(__m128i a, __m128i c) { __m128i r; for (int i = 0; i < 8; ++i) m_set16(r.b, i, (int16_t)(m_get16(a.b, i) + m_get16(c.b, i))); return r; }
EH_SIMD_MODEL __m128i _mm_max_epi16(__m128i a, __m128i c) { __m128i r; for (int i = 0; i < 8; ++i) { int16_t x = m_get16(a.b, i), y = m_get16(c.b, i); m_set16(r.b, i, x > y ? x : y); } return r; }
EH_SIMD_MODEL __m128i _mm_or_si128(__m128i a, __m128i c) { __m128i r; for (int i = 0; i < 16; ++i) r.b[i] = a.b[i] | c.b[i]; return r; }
EH_SIMD_MODEL __m128i _mm_and_si128(__m128i a, __m128i c) { __m128i r; for (int i = 0; i < 16; ++i) r.b[i] = a.b[i] & c.b[i]; return r; }
EH_SIMD_MODEL __m128i _mm_cmpeq_epi8(__m128i a, __m128i c) { __m128i r; for (int i = 0; i < 16; ++i) r.b[i] = a.b[i] == c.b[i] ? 0xff : 0; return r; }
EH_SIMD_MODEL int _mm_movemask_epi8(__m128i a) { int m = 0; for (int i = 0; i < 16; ++i) m |= (a.b[i] >> 7) << i; return m; }
EH_SIMD_MODEL __m128i _mm_srli_epi16(__m128i a, int count)
{
    __m128i r;
    for (int i = 0; i < 8; ++i) m_set16(r.b, i, count > 15 ? 0 : (int16_t)((uint16_t)m_get16(a.b, i) >> count));
    return r;
}
EH_SIMD_MODEL __m128i _mm_unpacklo_epi8(__m128i a, __m128i c) { __m128i r; for (int i = 0; i < 8; ++i) { r.b[2 * i] = a.b[i]; r.b[2 * i + 1] = c.b[i]; } return r; }
EH_SIMD_MODEL __m128i _mm_unpackhi_epi8(__m128i a, __m128i c) { __m128i r; for (int i = 0; i < 8; ++i) { r.b[2 * i] = a.b[8 + i]; r.b[2 * i + 1] = c.b[8 + i]; } return r; }
EH_SIMD_MODEL __m128i _mm_unpacklo_epi64(__m128i a, __m128i c) { __m128i r; memcpy(r.b, a.b, 8); memcpy(r.b + 8, c.b, 8); return r; }
EH_SIMD_MODEL __m128i _mm_packs_epi16(__m128i a, __m128i c)
{
    __m128i r;
    for (int i = 0; i < 8; ++i) { r.b[i] = (uint8_t)m_sat8_from16(m_get16(a.b, i)); r.b[8 + i] = (uint8_t)m_sat8_from16(m_get16(c.b, i)); }
    return r;
}

// ---------------------------------------------------------------- 256-bit
EH_SIMD_MODEL __m256i _mm256_loadu_si256(const __m256i* p) { __m256i r; memcpy(r.b, p, 32); return r; }
EH_SIMD_MODEL void _mm256_storeu_si256(__m256i* p, __m256i a) { memcpy(p, a.b, 32); }
EH_SIMD_MODEL __m256i _mm256_setzero_si256(void) { __m256i r; memset(r.b, 0, 32); return r; }
EH_SIMD_MODEL __m256i _mm256_set1_epi8(char v) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = (uint8_t)v; return r; }
EH_SIMD_MODEL __m256i _mm256_set1_epi16(short v) { __m256i r; for (int i = 0; i < 16; ++i) m_set16(r.b, i, v); return r; }
EH_SIMD_MODEL __m256i _mm256_broadcastsi128_si256(__m128i a) { __m256i r; memcpy(r.b, a.b, 16); memcpy(r.b + 16, a.b, 16); return r; }
EH_SIMD_MODEL __m128i _mm256_castsi256_si128(__m256i a) { __m128i r; memcpy(r.b, a.b, 16); return r; }
// The upper half is undefined in hardware; the kernels overwrite it (inserti128) before use. Model it as zero.
EH_SIMD_MODEL __m256i _mm256_castsi128_si256(__m128i a) { __m256i r; memcpy(r.b, a.b, 16); memset(r.b + 16, 0, 16); return r; }
EH_SIMD_MODEL __m128i _mm256_extracti128_si256(__m256i a, int imm) { __m128i r; memcpy(r.b, a.b + 16 * (imm & 1), 16); return r; }
EH_SIMD_MODEL __m256i _mm256_inserti128_si256(__m256i a, __m128i c, int imm) { __m256i r = a; memcpy(r.b + 16 * (imm & 1), c.b, 16); return r; }
EH_SIMD_MODEL __m256i _mm256_and_si256(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = a.b[i] & c.b[i]; return r; }
EH_SIMD_MODEL __m256i _mm256_or_si256(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = a.b[i] | c.b[i]; return r; }
EH_SIMD_MODEL __m256i _mm256_cmpeq_epi8(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = a.b[i] == c.b[i] ? 0xff : 0; return r; }
EH_SIMD_MODEL __m256i _mm256_cmpgt_epi8(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = (int8_t)a.b[i] > (int8_t)c.b[i] ? 0xff : 0; return r; }
EH_SIMD_MODEL __m256i _mm256_sub_epi8(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = (uint8_t)(a.b[i] - c.b[i]); return r; }
EH_SIMD_MODEL __m256i _mm256_min_epu8(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = a.b[i] < c.b[i] ? a.b[i] : c.b[i]; return r; }
EH_SIMD_MODEL __m256i _mm256_subs_epu8(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = a.b[i] > c.b[i] ? (uint8_t)(a.b[i] - c.b[i]) : 0; return r; }
EH_SIMD_MODEL int _mm256_movemask_epi8(__m256i a) { uint32_t m = 0; for (int i = 0; i < 32; ++i) m |= (uint32_t)(a.b[i] >> 7) << i; return (int)m; }
// Per 128-bit lane: index bit 7 set gives 0, otherwise the byte of the same lane at index & 15.
EH_SIMD_MODEL __m256i _mm256_shuffle_epi8(__m256i a, __m256i c)
{
    __m256i r;
    for (int i = 0; i < 32; ++i) r.b[i] = (c.b[i] & 0x80) ? 0 : a.b[(i & 16) + (c.b[i] & 15)];
    return r;
}
EH_SIMD_MODEL __m256i _mm256_blendv_epi8(__m256i a, __m256i c, __m256i mask) { __m256i r; for (int i = 0; i < 32; ++i) r.b[i] = (mask.b[i] & 0x80) ? c.b[i] : a.b[i]; return r; }
// Unsigned bytes of a times signed bytes of c, adjacent products added with signed saturation.
EH_SIMD_MODEL __m256i _mm256_maddubs_epi16(__m256i a, __m256i c)
{
    __m256i r;
    for (int i = 0; i < 16; ++i)
        m_set16(r.b, i, m_sat16((int32_t)a.b[2 * i] * (int8_t)c.b[2 * i] + (int32_t)a.b[2 * i + 1] * (int8_t)c.b[2 * i + 1]));
    return r;
}
// Per 128-bit lane: 8 words of a then 8 words of c, unsigned saturation.
EH_SIMD_MODEL __m256i _mm256_packus_epi16(__m256i a, __m256i c)
{
    __m256i r;
    for (int lane = 0; lane < 2; ++lane)
        for (int i = 0; i < 8; ++i)
        {
            r.b[16 * lane + i] = m_satu8_from16(m_get16(a.b, 8 * lane + i));
            r.b[16 * lane + 8 + i] = m_satu8_from16(m_get16(c.b, 8 * lane + i));
        }
    return r;
}
EH_SIMD_MODEL __m256i _mm256_add_epi16(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 16; ++i) m_set16(r.b, i, (int16_t)(m_get16(a.b, i) + m_get16(c.b, i))); return r; }
EH_SIMD_MODEL __m256i _mm256_adds_epi16(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 16; ++i) m_set16(r.b, i, m_sat16((int32_t)m_get16(a.b, i) + m_get16(c.b, i))); return r; }
EH_SIMD_MODEL __m256i _mm256_max_epi16(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 16; ++i) { int16_t x = m_get16(a.b, i), y = m_get16(c.b, i); m_set16(r.b, i, x > y ? x : y); } return r; }
EH_SIMD_MODEL __m256i _mm256_cmpgt_epi16(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 16; ++i) m_set16(r.b, i, m_get16(a.b, i) > m_get16(c.b, i) ? -1 : 0); return r; }
EH_SIMD_MODEL __m256i _mm256_cmpeq_epi16(__m256i a, __m256i c) { __m256i r; for (int i = 0; i < 16; ++i) m_set16(r.b, i, m_get16(a.b, i) == m_get16(c.b, i) ? -1 : 0); return r; }
EH_SIMD_MODEL __m256i _mm256_insert_epi16(__m256i a, int value, int index) { __m256i r = a; m_set16(r.b, index & 15, (int16_t)value); return r; }
EH_SIMD_MODEL int _mm256_extract_epi16(__m256i a, int index) { return (int)(uint16_t)m_get16(a.b, index & 15); }
// imm 0x20: low 128 bits of a, then low 128 bits of c (the only control the kernels use; others are rejected).
EH_SIMD_MODEL __m256i _mm256_permute2x128_si256(__m256i a, __m256i c, int imm)
{
    __m256i r;
    const uint8_t* sources[4] = { a.b, a.b + 16, c.b, c.b + 16 };
    for (int half = 0; half < 2; ++half)
    {
        const int control = (imm >> (4 * half)) & 0xf;
        if (control & 8) memset(r.b + 16 * half, 0, 16);
        else memcpy(r.b + 16 * half, sources[control & 3], 16);
    }
    return r;
}
// Per 128-bit lane: the 32-byte concatenation (a lane : c lane) shifted right by imm bytes, low 16 bytes kept.
EH_SIMD_MODEL __m256i _mm256_alignr_epi8(__m256i a, __m256i c, int imm)
{
    __m256i r;
    for (int lane = 0; lane < 2; ++lane)
        for (int i = 0; i < 16; ++i)
        {
            const int k = i + imm;
            r.b[16 * lane + i] = k < 16 ? c.b[16 * lane + k] : k < 32 ? a.b[16 * lane + k - 16] : 0;
        }
    return r;
}

// ---------------------------------------------------------------- 512-bit (AVX-512BW)
EH_SIMD_MODEL __m512i _mm512_loadu_si512(const void* p) { __m512i r; memcpy(r.b, p, 64); return r; }
EH_SIMD_MODEL void _mm512_storeu_si512(void* p, __m512i a) { memcpy(p, a.b, 64); }
EH_SIMD_MODEL __m512i _mm512_set1_epi16(short v) { __m512i r; for (int i = 0; i < 32; ++i) m_set16(r.b, i, v); return r; }
// Arguments from the highest lane (e31) to the lowest (e0), as in the hardware intrinsic.
EH_SIMD_MODEL __m512i _mm512_set_epi16(short e31, short e30, short e29, short e28, short e27, short e26, short e25,
    short e24, short e23, short e22, short e21, short e20, short e19, short e18, short e17, short e16, short e15,
    short e14, short e13, short e12, short e11, short e10, short e9, short e8, short e7, short e6, short e5, short e4,
    short e3, short e2, short e1, short e0)
{
    const short v[32] = { e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19,
        e20, e21, e22, e23, e24, e25, e26, e27, e28, e29, e30, e31 };
    __m512i r;
    for (int i = 0; i < 32; ++i) m_set16(r.b, i, v[i]);
    return r;
}
EH_SIMD_MODEL __m512i _mm512_sub_epi16(__m512i a, __m512i c) { __m512i r; for (int i = 0; i < 32; ++i) m_set16(r.b, i, (int16_t)(m_get16(a.b, i) - m_get16(c.b, i))); return r; }
EH_SIMD_MODEL __m512i _mm512_subs_epi16(__m512i a, __m512i c) { __m512i r; for (int i = 0; i < 32; ++i) m_set16(r.b, i, m_sat16((int32_t)m_get16(a.b, i) - m_get16(c.b, i))); return r; }
EH_SIMD_MODEL __m512i _mm512_max_epi16(__m512i a, __m512i c) { __m512i r; for (int i = 0; i < 32; ++i) { int16_t x = m_get16(a.b, i), y = m_get16(c.b, i); m_set16(r.b, i, x > y ? x : y); } return r; }
// Result lane i takes lane (index[i] & 31) of a.
EH_SIMD_MODEL __m512i _mm512_permutexvar_epi16(__m512i index, __m512i a) { __m512i r; for (int i = 0; i < 32; ++i) m_set16(r.b, i, m_get16(a.b, m_get16(index.b, i) & 31)); return r; }
EH_SIMD_MODEL __m512i _mm512_mask_mov_epi16(__m512i src, __mmask32 k, __m512i a) { __m512i r; for (int i = 0; i < 32; ++i) m_set16(r.b, i, ((k >> i) & 1) ? m_get16(a.b, i) : m_get16(src.b, i)); return r; }
EH_SIMD_MODEL __m512i _mm512_mask_max_epi16(__m512i src, __mmask32 k, __m512i a, __m512i c)
{
    __m512i r;
    for (int i = 0; i < 32; ++i) { int16_t x = m_get16(a.b, i), y = m_get16(c.b, i); m_set16(r.b, i, ((k >> i) & 1) ? (x > y ? x : y) : m_get16(src.b, i)); }
    return r;
}
EH_SIMD_MODEL __mmask32 _mm512_cmplt_epi16_mask(__m512i a, __m512i c) { __mmask32 m = 0; for (int i = 0; i < 32; ++i) m |= (uint32_t)(m_get16(a.b, i) < m_get16(c.b, i)) << i; return m; }
EH_SIMD_MODEL __mmask32 _mm512_cmpeq_epi16_mask(__m512i a, __m512i c) { __mmask32 m = 0; for (int i = 0; i < 32; ++i) m |= (uint32_t)(m_get16(a.b, i) == m_get16(c.b, i)) << i; return m; }
EH_SIMD_MODEL __m256i _mm512_extracti64x4_epi64(__m512i a, int imm) { __m256i r; memcpy(r.b, a.b + 32 * (imm & 1), 32); return r; }
