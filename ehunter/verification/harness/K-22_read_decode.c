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
// Obligation [K-22]: for all packed bases and qualities, the AVX2 decoder followed by the two-nibble tail
// writes the bases of v5.0.0's per-base loop, and reads no byte outside packed[0, (n + 1) / 2) or
// qualities[0, n) (checked by the model checker's bounds checks on exactly sized arrays). Complete for one
// 32-base block (48 symbolic bytes) plus TAIL bases (default 1); further blocks repeat the same state.
#include "harness.h"
#include "tables.h"
#include "core/ReadDecodeKernels.hh"

#if defined(__aarch64__)
#define K22_KERNEL rd_decode_neon_blocks
#else
#define K22_KERNEL rd_decode_avx2_blocks
#endif

#ifndef TAIL
#define TAIL 1
#endif
#define N (32 + TAIL)

static void check_once(void)
{
    uint8_t packed[(N + 1) / 2], qualities[N];
    for (int i = 0; i < (N + 1) / 2; ++i)
        packed[i] = IN_U8();
    for (int i = 0; i < N; ++i)
        qualities[i] = IN_U8();
    const char* upper = (const char*)hts_nt16_upper;
    const char* lower = (const char*)hts_nt16_lower;
    char bases[N], expected[N];
    const int32_t index = K22_KERNEL(N, packed, qualities, bases, upper);
    OBLIGATION(index == 32);
    rd_decode_tail(index, N, packed, qualities, bases, upper, lower);
    rd_decode_reference(N, packed, qualities, expected, upper, lower);
    for (int i = 0; i < N; ++i)
        OBLIGATION(bases[i] == expected[i]);
}

int main(void)
{
#ifndef __CPROVER__
#if !defined(__aarch64__)
    if (!__builtin_cpu_supports("avx2"))
        return 77;
#endif
#endif
    for (long round = 0; round < HARNESS_ROUNDS; ++round)
        check_once();
    return 0;
}
