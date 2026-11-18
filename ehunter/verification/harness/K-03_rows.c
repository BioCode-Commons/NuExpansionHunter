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

// Obligation [K-03 rows]: for all 16-bit inputs and gap scores, each SIMD row pass (fused, deletion, align,
// consolidate) writes exactly what its reference row writes. Complete per block: every input is symbolic.
// Build define K03_TIER selects the kernels: 1 = SSE2 (two 8-lane vectors), 2 = AVX2 (one 16-lane vector). On
// aarch64 the NEON rows (two 8-lane vectors) are checked in executable mode.
#include "harness.h"
#include "graphalign/dagAligner/DpKernels.hh"

#if defined(__aarch64__)
#define LANES 16
#define TIER neon
#elif K03_TIER == 2
#define LANES 16
#define TIER avx2
#else
#define LANES 16
#define TIER sse2
#endif
#define KERNEL_NAME(PASS, SUFFIX) dp_##PASS##_row_##SUFFIX
#define KERNEL(PASS, SUFFIX) KERNEL_NAME(PASS, SUFFIX)

static void check_once(void)
{
    // vpRow[0] is the predecessor row's boundary cell, read as vp[-1] by lane 0.
    int16_t epRow[LANES], vpRow[LANES + 1], penalties[LANES];
    for (int i = 0; i < LANES; ++i)
    {
        epRow[i] = IN_I16();
        penalties[i] = IN_I16();
    }
    for (int i = 0; i < LANES + 1; ++i)
        vpRow[i] = IN_I16();
    const int16_t* vp = vpRow + 1;
    const int16_t gapOpen = IN_I16(), gapExt = IN_I16();

    int16_t e[LANES], v[LANES], eRef[LANES], vRef[LANES];
    KERNEL(fused, TIER)(epRow, vp, penalties, gapOpen, gapExt, e, v, LANES);
    dp_fused_row_reference(epRow, vp, penalties, gapOpen, gapExt, eRef, vRef, LANES);
    for (int i = 0; i < LANES; ++i)
        OBLIGATION(e[i] == eRef[i] && v[i] == vRef[i]);

    int16_t g[LANES], gRef[LANES];
    for (int i = 0; i < LANES; ++i)
    {
        e[i] = eRef[i] = IN_I16();
        g[i] = gRef[i] = IN_I16();
        v[i] = vRef[i] = IN_I16();
    }
    KERNEL(deletion, TIER)(epRow, vp, gapOpen, gapExt, e, LANES);
    dp_deletion_row_reference(epRow, vp, gapOpen, gapExt, eRef, LANES);
    KERNEL(align, TIER)(vp, penalties, g, LANES);
    dp_align_row_reference(vp, penalties, gRef, LANES);
    KERNEL(consolidate, TIER)(g, e, v, LANES);
    dp_consolidate_row_reference(gRef, eRef, vRef, LANES);
    for (int i = 0; i < LANES; ++i)
        OBLIGATION(e[i] == eRef[i] && g[i] == gRef[i] && v[i] == vRef[i]);
}

int main(void)
{
#if !defined(__CPROVER__) && !defined(__aarch64__)
#if K03_TIER == 2
    if (!__builtin_cpu_supports("avx2"))
        return 77;
#endif
#endif
    for (long round = 0; round < HARNESS_ROUNDS; ++round)
        check_once();
    return 0;
}
