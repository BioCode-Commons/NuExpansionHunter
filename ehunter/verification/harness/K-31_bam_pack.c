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
// Obligation [K-31]: for all query strings, the AVX2 packer followed by the scalar tail writes the sequence
// and quality bytes of the scalar reference, and touches no byte outside the exactly sized arrays. Complete
// for one 32-base vector plus TAIL bases (default 1, odd so that the overlapping final vector and the odd
// nibble are both exercised).
#include "harness.h"
#include "tables.h"
#include "io/BamEncodingKernels.hh"

#if defined(__aarch64__)
#define K31_KERNEL be_encode_neon_blocks
#else
#define K31_KERNEL be_encode_avx2_blocks
#endif

#ifndef TAIL
#define TAIL 1
#endif
#define N (32 + TAIL)

static void check_once(void)
{
    char query[N];
    for (int i = 0; i < N; ++i)
        query[i] = (char)IN_U8();
    uint8_t sequence[(N + 1) / 2], quality[N], expectedSequence[(N + 1) / 2], expectedQuality[N];
    const size_t index = K31_KERNEL(query, N, sequence, quality, hts_nt16_table);
    be_encode_tail(index, query, N, sequence, quality, hts_nt16_table);
    be_encode_reference(query, N, expectedSequence, expectedQuality, hts_nt16_table);
    for (int i = 0; i < (N + 1) / 2; ++i)
        OBLIGATION(sequence[i] == expectedSequence[i]);
    for (int i = 0; i < N; ++i)
        OBLIGATION(quality[i] == expectedQuality[i]);
}

int main(void)
{
#ifndef __CPROVER__
#if !defined(__aarch64__)
    if (!__builtin_cpu_supports("avx2"))
        return 77;
#endif
    for (int c = 0; c < 256; ++c)
        OBLIGATION(!isupper(c) == !(c >= 'A' && c <= 'Z'));
#endif
    for (long round = 0; round < HARNESS_ROUNDS; ++round)
        check_once();
    return 0;
}
