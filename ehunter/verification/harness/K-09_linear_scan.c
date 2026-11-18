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
// Obligation [K-09]: for all code arrays (target codes 0..15, query codes 0..4), the AVX2 linear mismatch scan
// returns the reference result. Complete for one vector block plus a scalar tail: n = 32 + TAIL symbolic codes
// (TAIL defaults to 1); the carried state of the block loop is the first-mismatch position, which a second
// block reads exactly as the tail does.
#include "harness.h"
#include "graphalign/dagAligner/DpKernels.hh"

#if defined(__aarch64__)
#define K09_KERNEL dp_linear_first_mismatch_neon
#else
#define K09_KERNEL dp_linear_first_mismatch_avx2
#endif

#ifndef TAIL
#define TAIL 1
#endif
#define N (32 + TAIL)

static void check_once(void)
{
    uint8_t query[N], target[N];
    for (int i = 0; i < N; ++i)
    {
        query[i] = IN_U8() % 5;
        target[i] = IN_U8() % 16;
#ifdef __CPROVER__
        ASSUME(query[i] <= 4 && target[i] <= 15);
#endif
    }
    for (int n = 32; n <= N; ++n)
        OBLIGATION(K09_KERNEL(query, target, n) == dp_linear_first_mismatch_reference(query, target, n));
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
