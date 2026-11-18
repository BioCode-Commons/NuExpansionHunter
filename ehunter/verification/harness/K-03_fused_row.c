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
// Obligation [K-03]: for every predecessor cell (E, V, V of the previous query position), penalty and pair of
// gap scores, the fused single-predecessor lane writes the same E and V as the three general passes started
// from DP_SCORE_MIN. Complete: all six 16-bit inputs are symbolic. Lanes of a row are independent, so the
// row follows lane by lane.
#include "harness.h"
#include "graphalign/dagAligner/DpKernels.hh"

static void check_once(void)
{
    const int16_t ep = IN_I16(), vp = IN_I16(), vpPrevious = IN_I16(), penalty = IN_I16();
    const int16_t gapOpen = IN_I16(), gapExt = IN_I16();

    int16_t fusedE, fusedV;
    dp_fused_lane(ep, vp, vpPrevious, penalty, gapOpen, gapExt, &fusedE, &fusedV);

    int16_t e = DP_SCORE_MIN, g = DP_SCORE_MIN, v = DP_SCORE_MIN;
    dp_deletion_lane(ep, vp, gapOpen, gapExt, &e);
    dp_align_lane(vpPrevious, penalty, &g);
    dp_consolidate_lane(g, e, &v);

    OBLIGATION(fusedE == e);
    OBLIGATION(fusedV == v);
}

int main(void)
{
    for (long round = 0; round < HARNESS_ROUNDS; ++round)
        check_once();
    return 0;
}
