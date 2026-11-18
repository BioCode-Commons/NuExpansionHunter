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
// Obligation [K-01]: for all 16-bit inputs, a vector block of the insertion scan equals the same number of
// steps of the one-cell recurrence (gap open -8, extension -2), including the carried F and V and the trace
// bits; a block the guard rejects changes nothing. Complete per block: every input is symbolic.
// Build define K01_AVX512 selects the 32-lane AVX-512BW block; the default is the 16-lane AVX2 block. On
// aarch64 the 8-lane NEON block is checked (executable mode; the NEON intrinsics have no model).
#include "harness.h"
#include "graphalign/dagAligner/DpKernels.hh"

#if defined(__aarch64__)
#define LANES 8
#define KERNEL dp_insertion_neon
#elif defined(K01_AVX512)
#define LANES 32
#define KERNEL dp_insertion_avx512
#else
#define LANES 16
#define KERNEL dp_insertion_avx2
#endif

static void check_once(void)
{
    int16_t v[LANES], expected[LANES];
    for (int i = 0; i < LANES; ++i)
        expected[i] = v[i] = IN_I16();
    int16_t f = IN_I16(), h = IN_I16();
    int16_t expectedF = f, expectedH = h;
    // Bits already in the trace word are kept; the block's bits go to (i - start) .. + LANES - 1.
    const int start = 0;
    uint64_t trace = 0, expectedTrace = 0;

    const int end = KERNEL(v, 0, LANES, &f, &h, &trace, start);

    if (end == 0)
    {
        for (int i = 0; i < LANES; ++i)
            OBLIGATION(v[i] == expected[i]);
        OBLIGATION(f == expectedF && h == expectedH && trace == 0);
        return;
    }
    OBLIGATION(end == LANES);
    for (int i = 0; i < LANES; ++i)
        expectedTrace |= (uint64_t)dp_insertion_cell(&expected[i], &expectedF, &expectedH, -8, -2) << i;
    for (int i = 0; i < LANES; ++i)
        OBLIGATION(v[i] == expected[i]);
    OBLIGATION(f == expectedF);
    OBLIGATION(h == expectedH);
    OBLIGATION(trace == expectedTrace);
}

int main(void)
{
#if !defined(__CPROVER__) && !defined(__aarch64__)
    if (!__builtin_cpu_supports("avx2") || !__builtin_cpu_supports("avx512bw"))
        return 77;
#endif
    for (long round = 0; round < HARNESS_ROUNDS; ++round)
        check_once();
    return 0;
}
