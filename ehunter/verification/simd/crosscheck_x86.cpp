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
// Checks every model of x86_models.h against the hardware intrinsic on this CPU.
// Build: g++ -O1 -mavx2 -mavx512f -mavx512bw crosscheck_x86.cpp; run without arguments. Exit status 0 when
// every model agrees on all inputs (boundary values, then random vectors).
#include <immintrin.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <stdint.h>
#include <string.h>

namespace model
{
#include "x86_models.h"
}

namespace
{
std::mt19937_64 rng(20261117);
int failures = 0;
long checks = 0;

const int16_t edges16[] = { 0, 1, -1, 2, -2, 127, 128, -128, -129, 255, 256, 32767, -32768, 32766, -32767, 32727,
    -32728, -32696, -32697, 10, -10, 40, -40 };
const uint8_t edges8[] = { 0, 1, 0x7f, 0x80, 0x81, 0xff, 0xfe, 0x20, 0xdf, 'A', 'Z', 'a', 'z', 'O', 'P', '@', '[', 20,
    21, 19, 32, 0x0f, 0x10 };

template <int N> void fill(uint8_t (&b)[N], int mode)
{
    for (int i = 0; i < N; ++i)
    {
        if (mode == 0) b[i] = uint8_t(rng());
        else if (mode == 1) b[i] = edges8[rng() % sizeof(edges8)];
        else
        {
            const int16_t v = edges16[rng() % (sizeof(edges16) / 2)];
            b[i] = (i & 1) ? uint8_t(uint16_t(v) >> 8) : uint8_t(uint16_t(v) & 0xff);
        }
    }
}

template <int N> bool same(const void* a, const void* b, const char* name)
{
    ++checks;
    if (std::memcmp(a, b, N) == 0) return true;
    if (failures++ < 20) std::printf("MISMATCH %s\n", name);
    return false;
}

#define V128(n) uint8_t n##_b[16]; fill(n##_b, mode); __m128i n; std::memcpy(&n, n##_b, 16); model::__m128i n##_m; std::memcpy(&n##_m, n##_b, 16)
#define V256(n) uint8_t n##_b[32]; fill(n##_b, mode); __m256i n; std::memcpy(&n, n##_b, 32); model::__m256i n##_m; std::memcpy(&n##_m, n##_b, 32)
#define V512(n) uint8_t n##_b[64]; fill(n##_b, mode); __m512i n; std::memcpy(&n, n##_b, 64); model::__m512i n##_m; std::memcpy(&n##_m, n##_b, 64)
#define CHECK(N, hw, md, name) do { auto h_ = (hw); auto m_ = (md); same<N>(&h_, &m_, name); } while (0)

__attribute__((target("avx2,avx512f,avx512bw"))) void round(int mode)
{
    V128(a); V128(c); V256(x); V256(y); V256(z); V512(p); V512(q); V512(r);
    const char s8 = char(edges8[rng() % sizeof(edges8)] ^ (mode == 0 ? uint8_t(rng()) : 0));
    const short s16 = short(edges16[rng() % (sizeof(edges16) / 2)] + (mode == 0 ? short(rng()) : 0));
    const uint32_t k = uint32_t(rng());
    // 128-bit
    CHECK(16, _mm_loadu_si128(reinterpret_cast<const __m128i*>(a_b)), model::_mm_loadu_si128(reinterpret_cast<const model::__m128i*>(a_b)), "_mm_loadu_si128");
    { uint8_t h[16], m[16]; _mm_storeu_si128(reinterpret_cast<__m128i*>(h), a); model::_mm_storeu_si128(reinterpret_cast<model::__m128i*>(m), a_m); same<16>(h, m, "_mm_storeu_si128"); }
    CHECK(16, _mm_set1_epi8(s8), model::_mm_set1_epi8(s8), "_mm_set1_epi8");
    CHECK(16, _mm_set1_epi16(s16), model::_mm_set1_epi16(s16), "_mm_set1_epi16");
    CHECK(16, _mm_add_epi16(a, c), model::_mm_add_epi16(a_m, c_m), "_mm_add_epi16");
    CHECK(16, _mm_max_epi16(a, c), model::_mm_max_epi16(a_m, c_m), "_mm_max_epi16");
    CHECK(16, _mm_or_si128(a, c), model::_mm_or_si128(a_m, c_m), "_mm_or_si128");
    CHECK(16, _mm_and_si128(a, c), model::_mm_and_si128(a_m, c_m), "_mm_and_si128");
    CHECK(16, _mm_cmpeq_epi8(a, c), model::_mm_cmpeq_epi8(a_m, c_m), "_mm_cmpeq_epi8");
    CHECK(4, _mm_movemask_epi8(a), model::_mm_movemask_epi8(a_m), "_mm_movemask_epi8");
    CHECK(16, _mm_srli_epi16(a, 4), model::_mm_srli_epi16(a_m, 4), "_mm_srli_epi16");
    CHECK(16, _mm_unpacklo_epi8(a, c), model::_mm_unpacklo_epi8(a_m, c_m), "_mm_unpacklo_epi8");
    CHECK(16, _mm_unpackhi_epi8(a, c), model::_mm_unpackhi_epi8(a_m, c_m), "_mm_unpackhi_epi8");
    CHECK(16, _mm_unpacklo_epi64(a, c), model::_mm_unpacklo_epi64(a_m, c_m), "_mm_unpacklo_epi64");
    CHECK(16, _mm_packs_epi16(a, c), model::_mm_packs_epi16(a_m, c_m), "_mm_packs_epi16");
    // 256-bit
    CHECK(32, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x_b)), model::_mm256_loadu_si256(reinterpret_cast<const model::__m256i*>(x_b)), "_mm256_loadu_si256");
    { uint8_t h[32], m[32]; _mm256_storeu_si256(reinterpret_cast<__m256i*>(h), x); model::_mm256_storeu_si256(reinterpret_cast<model::__m256i*>(m), x_m); same<32>(h, m, "_mm256_storeu_si256"); }
    CHECK(32, _mm256_setzero_si256(), model::_mm256_setzero_si256(), "_mm256_setzero_si256");
    CHECK(32, _mm256_set1_epi8(s8), model::_mm256_set1_epi8(s8), "_mm256_set1_epi8");
    CHECK(32, _mm256_set1_epi16(s16), model::_mm256_set1_epi16(s16), "_mm256_set1_epi16");
    CHECK(32, _mm256_broadcastsi128_si256(a), model::_mm256_broadcastsi128_si256(a_m), "_mm256_broadcastsi128_si256");
    CHECK(16, _mm256_castsi256_si128(x), model::_mm256_castsi256_si128(x_m), "_mm256_castsi256_si128");
    CHECK(32, _mm256_inserti128_si256(_mm256_castsi128_si256(a), c, 1), model::_mm256_inserti128_si256(model::_mm256_castsi128_si256(a_m), c_m, 1), "_mm256_castsi128_si256+inserti128");
    CHECK(16, _mm256_extracti128_si256(x, 1), model::_mm256_extracti128_si256(x_m, 1), "_mm256_extracti128_si256(1)");
    CHECK(16, _mm256_extracti128_si256(x, 0), model::_mm256_extracti128_si256(x_m, 0), "_mm256_extracti128_si256(0)");
    CHECK(32, _mm256_inserti128_si256(x, a, 0), model::_mm256_inserti128_si256(x_m, a_m, 0), "_mm256_inserti128_si256(0)");
    CHECK(32, _mm256_inserti128_si256(x, a, 1), model::_mm256_inserti128_si256(x_m, a_m, 1), "_mm256_inserti128_si256(1)");
    CHECK(32, _mm256_and_si256(x, y), model::_mm256_and_si256(x_m, y_m), "_mm256_and_si256");
    CHECK(32, _mm256_or_si256(x, y), model::_mm256_or_si256(x_m, y_m), "_mm256_or_si256");
    CHECK(32, _mm256_cmpeq_epi8(x, y), model::_mm256_cmpeq_epi8(x_m, y_m), "_mm256_cmpeq_epi8");
    CHECK(32, _mm256_cmpgt_epi8(x, y), model::_mm256_cmpgt_epi8(x_m, y_m), "_mm256_cmpgt_epi8");
    CHECK(32, _mm256_sub_epi8(x, y), model::_mm256_sub_epi8(x_m, y_m), "_mm256_sub_epi8");
    CHECK(32, _mm256_min_epu8(x, y), model::_mm256_min_epu8(x_m, y_m), "_mm256_min_epu8");
    CHECK(32, _mm256_subs_epu8(x, y), model::_mm256_subs_epu8(x_m, y_m), "_mm256_subs_epu8");
    CHECK(4, _mm256_movemask_epi8(x), model::_mm256_movemask_epi8(x_m), "_mm256_movemask_epi8");
    CHECK(32, _mm256_shuffle_epi8(x, y), model::_mm256_shuffle_epi8(x_m, y_m), "_mm256_shuffle_epi8");
    CHECK(32, _mm256_blendv_epi8(x, y, z), model::_mm256_blendv_epi8(x_m, y_m, z_m), "_mm256_blendv_epi8");
    CHECK(32, _mm256_maddubs_epi16(x, y), model::_mm256_maddubs_epi16(x_m, y_m), "_mm256_maddubs_epi16");
    CHECK(32, _mm256_packus_epi16(x, y), model::_mm256_packus_epi16(x_m, y_m), "_mm256_packus_epi16");
    CHECK(32, _mm256_add_epi16(x, y), model::_mm256_add_epi16(x_m, y_m), "_mm256_add_epi16");
    CHECK(32, _mm256_adds_epi16(x, y), model::_mm256_adds_epi16(x_m, y_m), "_mm256_adds_epi16");
    CHECK(32, _mm256_max_epi16(x, y), model::_mm256_max_epi16(x_m, y_m), "_mm256_max_epi16");
    CHECK(32, _mm256_cmpgt_epi16(x, y), model::_mm256_cmpgt_epi16(x_m, y_m), "_mm256_cmpgt_epi16");
    CHECK(32, _mm256_cmpeq_epi16(x, y), model::_mm256_cmpeq_epi16(x_m, y_m), "_mm256_cmpeq_epi16");
    CHECK(32, _mm256_insert_epi16(x, s16, 0), model::_mm256_insert_epi16(x_m, s16, 0), "_mm256_insert_epi16(0)");
    CHECK(4, _mm256_extract_epi16(x, 15), model::_mm256_extract_epi16(x_m, 15), "_mm256_extract_epi16(15)");
    CHECK(4, _mm256_extract_epi16(x, 0), model::_mm256_extract_epi16(x_m, 0), "_mm256_extract_epi16(0)");
    CHECK(32, _mm256_permute2x128_si256(x, y, 0x20), model::_mm256_permute2x128_si256(x_m, y_m, 0x20), "_mm256_permute2x128_si256(0x20)");
    CHECK(32, _mm256_alignr_epi8(x, y, 14), model::_mm256_alignr_epi8(x_m, y_m, 14), "_mm256_alignr_epi8(14)");
    CHECK(32, _mm256_alignr_epi8(x, y, 12), model::_mm256_alignr_epi8(x_m, y_m, 12), "_mm256_alignr_epi8(12)");
    CHECK(32, _mm256_alignr_epi8(x, y, 8), model::_mm256_alignr_epi8(x_m, y_m, 8), "_mm256_alignr_epi8(8)");
    // 512-bit
    CHECK(64, _mm512_loadu_si512(p_b), model::_mm512_loadu_si512(p_b), "_mm512_loadu_si512");
    { uint8_t h[64], m[64]; _mm512_storeu_si512(h, p); model::_mm512_storeu_si512(m, p_m); same<64>(h, m, "_mm512_storeu_si512"); }
    CHECK(64, _mm512_set1_epi16(s16), model::_mm512_set1_epi16(s16), "_mm512_set1_epi16");
    CHECK(64, _mm512_set_epi16(31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0),
          model::_mm512_set_epi16(31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0), "_mm512_set_epi16");
    CHECK(64, _mm512_sub_epi16(p, q), model::_mm512_sub_epi16(p_m, q_m), "_mm512_sub_epi16");
    CHECK(64, _mm512_subs_epi16(p, q), model::_mm512_subs_epi16(p_m, q_m), "_mm512_subs_epi16");
    CHECK(64, _mm512_max_epi16(p, q), model::_mm512_max_epi16(p_m, q_m), "_mm512_max_epi16");
    CHECK(64, _mm512_permutexvar_epi16(p, q), model::_mm512_permutexvar_epi16(p_m, q_m), "_mm512_permutexvar_epi16");
    CHECK(64, _mm512_mask_mov_epi16(p, k, q), model::_mm512_mask_mov_epi16(p_m, k, q_m), "_mm512_mask_mov_epi16");
    CHECK(64, _mm512_mask_max_epi16(p, k, q, r), model::_mm512_mask_max_epi16(p_m, k, q_m, r_m), "_mm512_mask_max_epi16");
    CHECK(4, _mm512_cmplt_epi16_mask(p, q), model::_mm512_cmplt_epi16_mask(p_m, q_m), "_mm512_cmplt_epi16_mask");
    CHECK(4, _mm512_cmpeq_epi16_mask(p, q), model::_mm512_cmpeq_epi16_mask(p_m, q_m), "_mm512_cmpeq_epi16_mask");
    CHECK(32, _mm512_extracti64x4_epi64(p, 1), model::_mm512_extracti64x4_epi64(p_m, 1), "_mm512_extracti64x4_epi64(1)");
    CHECK(32, _mm512_extracti64x4_epi64(p, 0), model::_mm512_extracti64x4_epi64(p_m, 0), "_mm512_extracti64x4_epi64(0)");
}
}

int main()
{
    if (!__builtin_cpu_supports("avx2") || !__builtin_cpu_supports("avx512bw"))
    {
        std::printf("crosscheck needs AVX2 and AVX-512BW\n");
        return 2;
    }
    for (int i = 0; i < 300000; ++i)
        round(i % 3);
    std::printf("simd model crosscheck: %ld comparisons, %d mismatches\n", checks, failures);
    return failures ? 1 : 0;
}
