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
// Verification switches.
//
// Each optimized kernel that has a reference (a scalar path, the dense matrix,
// an unskipped loop) can be pinned to one side for the whole binary, so that
// "release output = reference output" can be checked on real data at both
// extremes of every policy. A release build defines none of the EH_VERIFY_*
// macros. They are set only through the CMake option EH_VERIFY, for example
// -DEH_VERIFY="SIMD=SCALAR;BUDGET_SEARCH=NEVER", which defines
// EH_VERIFY_<NAME>=EH_VERIFY_<VALUE>. Kernel code reads the constants below and
// never tests the macros directly. ehunter/verification/README.md lists the
// switches with the kernels each one covers.
#pragma once

#include <string>

#define EH_VERIFY_DEFAULT 0
#define EH_VERIFY_OFF 1
#define EH_VERIFY_ON 2
#define EH_VERIFY_SCALAR 3
#define EH_VERIFY_AVX2 4
#define EH_VERIFY_NEVER 5
#define EH_VERIFY_ALWAYS 6
#define EH_VERIFY_REFERENCE 7
#define EH_VERIFY_SSE42 8

#ifndef EH_VERIFY_SIMD
#define EH_VERIFY_SIMD EH_VERIFY_DEFAULT
#endif
#ifndef EH_VERIFY_FAST_PATHS
#define EH_VERIFY_FAST_PATHS EH_VERIFY_DEFAULT
#endif
#ifndef EH_VERIFY_BUDGET_SEARCH
#define EH_VERIFY_BUDGET_SEARCH EH_VERIFY_DEFAULT
#endif
#ifndef EH_VERIFY_ROW_STOP
#define EH_VERIFY_ROW_STOP EH_VERIFY_DEFAULT
#endif
#ifndef EH_VERIFY_FUSED_ROWS
#define EH_VERIFY_FUSED_ROWS EH_VERIFY_DEFAULT
#endif
#ifndef EH_VERIFY_GRAPH_CACHE
#define EH_VERIFY_GRAPH_CACHE EH_VERIFY_DEFAULT
#endif
#ifndef EH_VERIFY_FOLD_SKIP
#define EH_VERIFY_FOLD_SKIP EH_VERIFY_DEFAULT
#endif
#ifndef EH_VERIFY_LAYOUT
#define EH_VERIFY_LAYOUT EH_VERIFY_DEFAULT
#endif
#ifndef EH_VERIFY_SHADOW
#define EH_VERIFY_SHADOW EH_VERIFY_DEFAULT
#endif

#if EH_VERIFY_SIMD != EH_VERIFY_DEFAULT && EH_VERIFY_SIMD != EH_VERIFY_SCALAR && EH_VERIFY_SIMD != EH_VERIFY_SSE42     \
    && EH_VERIFY_SIMD != EH_VERIFY_AVX2
#error "EH_VERIFY: SIMD must be SCALAR, SSE42 or AVX2"
#endif
#if EH_VERIFY_FAST_PATHS != EH_VERIFY_DEFAULT && EH_VERIFY_FAST_PATHS != EH_VERIFY_OFF
#error "EH_VERIFY: FAST_PATHS must be OFF"
#endif
#if EH_VERIFY_BUDGET_SEARCH != EH_VERIFY_DEFAULT && EH_VERIFY_BUDGET_SEARCH != EH_VERIFY_NEVER                       \
    && EH_VERIFY_BUDGET_SEARCH != EH_VERIFY_ALWAYS
#error "EH_VERIFY: BUDGET_SEARCH must be NEVER or ALWAYS"
#endif
#if EH_VERIFY_ROW_STOP != EH_VERIFY_DEFAULT && EH_VERIFY_ROW_STOP != EH_VERIFY_OFF
#error "EH_VERIFY: ROW_STOP must be OFF"
#endif
#if EH_VERIFY_FUSED_ROWS != EH_VERIFY_DEFAULT && EH_VERIFY_FUSED_ROWS != EH_VERIFY_OFF
#error "EH_VERIFY: FUSED_ROWS must be OFF"
#endif
#if EH_VERIFY_GRAPH_CACHE != EH_VERIFY_DEFAULT && EH_VERIFY_GRAPH_CACHE != EH_VERIFY_OFF
#error "EH_VERIFY: GRAPH_CACHE must be OFF"
#endif
#if EH_VERIFY_FOLD_SKIP != EH_VERIFY_DEFAULT && EH_VERIFY_FOLD_SKIP != EH_VERIFY_OFF
#error "EH_VERIFY: FOLD_SKIP must be OFF"
#endif
#if EH_VERIFY_LAYOUT != EH_VERIFY_DEFAULT && EH_VERIFY_LAYOUT != EH_VERIFY_REFERENCE
#error "EH_VERIFY: LAYOUT must be REFERENCE"
#endif
#if EH_VERIFY_SHADOW != EH_VERIFY_DEFAULT && EH_VERIFY_SHADOW != EH_VERIFY_ON
#error "EH_VERIFY: SHADOW must be ON"
#endif

namespace graphtools
{
namespace verify
{
// SIMD=SCALAR: every dispatched kernel takes its scalar path (K-01, K-03 rows, K-09, K-19, K-22, K-23, K-31).
constexpr bool kScalarOnly = EH_VERIFY_SIMD == EH_VERIFY_SCALAR;
// SIMD=SSE42: the x86-64-v2 tier; no AVX2 or AVX-512 path (the SSE kernels of K-03 rows, K-19, K-23 and the
// k-mer probe stay; K-01, K-09, K-22 and K-31 run scalar).
constexpr bool kNoAvx2 = EH_VERIFY_SIMD == EH_VERIFY_SCALAR || EH_VERIFY_SIMD == EH_VERIFY_SSE42;
// SIMD=AVX2: no AVX-512 path; AVX2 paths still require the CPU to support AVX2.
constexpr bool kNoAvx512 = EH_VERIFY_SIMD != EH_VERIFY_DEFAULT;
// FAST_PATHS=OFF: every extension goes to the DP with no proven minimum cost (K-09, K-11, K-12, K-17).
constexpr bool kFastPaths = EH_VERIFY_FAST_PATHS == EH_VERIFY_DEFAULT;
// BUDGET_SEARCH=NEVER: every eligible call uses the dense matrix (K-07, K-10).
constexpr bool kBudgetSearchNever = EH_VERIFY_BUDGET_SEARCH == EH_VERIFY_NEVER;
// BUDGET_SEARCH=ALWAYS: every eligible call searches, with no skip, no measured cap and no cost hint (K-10).
constexpr bool kBudgetSearchAlways = EH_VERIFY_BUDGET_SEARCH == EH_VERIFY_ALWAYS;
// ROW_STOP=OFF: the dense fill computes every row (K-05).
constexpr bool kRowStop = EH_VERIFY_ROW_STOP == EH_VERIFY_DEFAULT;
// FUSED_ROWS=OFF: rows with one predecessor use the general row code (K-03, K-08).
constexpr bool kFusedRows = EH_VERIFY_FUSED_ROWS == EH_VERIFY_DEFAULT;
// GRAPH_CACHE=OFF: the target graph is prepared afresh for every call (K-13).
constexpr bool kGraphCache = EH_VERIFY_GRAPH_CACHE == EH_VERIFY_DEFAULT;
// FOLD_SKIP=OFF: the canonical fold runs every step (K-21).
constexpr bool kFoldSkip = EH_VERIFY_FOLD_SKIP == EH_VERIFY_DEFAULT;
// LAYOUT=REFERENCE: hash layout for every k-mer index, owned catalog in streaming (K-15, K-28).
constexpr bool kReferenceLayout = EH_VERIFY_LAYOUT == EH_VERIFY_REFERENCE;
// SHADOW=ON: kernels with an in-tree reference also run it and abort on a difference (K-05, K-21, K-22, K-31).
constexpr bool kShadow = EH_VERIFY_SHADOW == EH_VERIFY_ON;

constexpr bool kAnySet = EH_VERIFY_SIMD != EH_VERIFY_DEFAULT || EH_VERIFY_FAST_PATHS != EH_VERIFY_DEFAULT
    || EH_VERIFY_BUDGET_SEARCH != EH_VERIFY_DEFAULT || EH_VERIFY_ROW_STOP != EH_VERIFY_DEFAULT
    || EH_VERIFY_FUSED_ROWS != EH_VERIFY_DEFAULT || EH_VERIFY_GRAPH_CACHE != EH_VERIFY_DEFAULT
    || EH_VERIFY_FOLD_SKIP != EH_VERIFY_DEFAULT || EH_VERIFY_LAYOUT != EH_VERIFY_DEFAULT
    || EH_VERIFY_SHADOW != EH_VERIFY_DEFAULT;

// Settings of this build, for the log; empty in a release build.
inline std::string describeSettings()
{
    std::string text;
    const auto add = [&text](const char* item) { text += (text.empty() ? "" : " ") + std::string(item); };
    if (kScalarOnly) add("SIMD=SCALAR");
    else if (kNoAvx2) add("SIMD=SSE42");
    else if (kNoAvx512) add("SIMD=AVX2");
    if (!kFastPaths) add("FAST_PATHS=OFF");
    if (kBudgetSearchNever) add("BUDGET_SEARCH=NEVER");
    if (kBudgetSearchAlways) add("BUDGET_SEARCH=ALWAYS");
    if (!kRowStop) add("ROW_STOP=OFF");
    if (!kFusedRows) add("FUSED_ROWS=OFF");
    if (!kGraphCache) add("GRAPH_CACHE=OFF");
    if (!kFoldSkip) add("FOLD_SKIP=OFF");
    if (kReferenceLayout) add("LAYOUT=REFERENCE");
    if (kShadow) add("SHADOW=ON");
    return text;
}
}
}
