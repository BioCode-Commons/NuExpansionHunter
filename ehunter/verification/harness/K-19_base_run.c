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
// Obligation [K-19]: for all byte strings, the SSE2 matching-base run equals the per-base reference, in both
// directions. Complete for one vector block plus a scalar tail: 16 + TAIL symbolic bytes of reference and of
// query (TAIL defaults to 2); a block either accepts 16 bases, after which the loop state is the same as
// before the block with matched + 16, or hands its first rejected position to the scalar step.
#include "harness.h"
#include "tables.h"
#include "graphutils/BaseMatchingKernels.hh"

#ifndef TAIL
#define TAIL 2
#endif
#define N (16 + TAIL)

static void check_once(void)
{
    char reference[N], query[N];
    for (int i = 0; i < N; ++i)
    {
        reference[i] = (char)IN_U8();
        query[i] = (char)IN_U8();
    }
    for (int reverse = 0; reverse < 2; ++reverse)
        for (size_t length = 16; length <= N; ++length)
            OBLIGATION(bm_matching_base_run(reference, query, length, reverse, 1, gt_reference_codes, gt_query_codes)
                == bm_matching_base_run_reference(reference, query, length, reverse, gt_reference_codes, gt_query_codes));
}

int main(void)
{
    // The kernel's own copy of the match table is the product's table.
    for (int r = 0; r < 16; ++r)
        for (int q = 0; q < 5; ++q)
            OBLIGATION(bm_code_match[r][q] == gt_code_match[5 * r + q]);
    for (long round = 0; round < HARNESS_ROUNDS; ++round)
        check_once();
    return 0;
}
