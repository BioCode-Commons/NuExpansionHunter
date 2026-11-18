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
// Two modes from one driver source.
//  - Proof (__CPROVER__ defined, run by CBMC): inputs are nondeterministic, preconditions are assumed, the
//    obligation is asserted; intrinsics are the lane-wise models of simd/x86_models.h.
//  - Executable (built by a C compiler with the hardware intrinsics): the same body runs over boundary values
//    and a seeded generator, and the program exits non-zero on the first violated obligation.
#pragma once

#include <stdint.h>
#include <stddef.h>

#ifdef __CPROVER__
#include <assert.h>
#include "../simd/x86_models.h"
#include <ctype.h>
// The C locale's isupper (the kernels' contract requires the C locale). Executable mode checks that the
// library's isupper agrees with it on every byte.
#undef isupper
int isupper(int c) { return c >= 'A' && c <= 'Z'; }
uint8_t nondet_u8(void);
uint16_t nondet_u16(void);
int nondet_int(void);
#define IN_U8() nondet_u8()
#define IN_I16() ((int16_t)nondet_u16())
#define IN_INT() nondet_int()
#define ASSUME(c) __CPROVER_assume(c)
#define OBLIGATION(c) assert(c)
#define HARNESS_ROUNDS 1
#else
#if defined(__aarch64__)
#include <arm_neon.h>
#else
#include <immintrin.h>
#endif
#include <stdio.h>
#include <stdlib.h>
static uint64_t harness_state = 0x9e3779b97f4a7c15ull;
static uint64_t harness_next(void)
{
    harness_state ^= harness_state << 13;
    harness_state ^= harness_state >> 7;
    harness_state ^= harness_state << 17;
    return harness_state;
}
// One value in four is a boundary value; the rest are uniform.
static __attribute__((unused)) uint8_t harness_u8(void)
{
    static const uint8_t edges[] = { 0, 1, 19, 20, 21, 32, 0x7f, 0x80, 0xff, 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', '=' };
    const uint64_t r = harness_next();
    return (r & 3) ? (uint8_t)(r >> 8) : edges[(r >> 8) % sizeof(edges)];
}
static __attribute__((unused)) int16_t harness_i16(void)
{
    static const int16_t edges[] = { 0, 1, -1, 32767, -32768, -32728, -32729, -32696, -32697, -32727, 100, -100 };
    const uint64_t r = harness_next();
    return (r & 3) ? (int16_t)(r >> 16) : edges[(r >> 16) % (sizeof(edges) / sizeof(edges[0]))];
}
#define IN_U8() harness_u8()
#define IN_I16() harness_i16()
#define IN_INT() ((int)harness_next())
#define ASSUME(c) do { if (!(c)) return; } while (0)
#define OBLIGATION(c) do { if (!(c)) { fprintf(stderr, "%s:%d: obligation violated: %s\n", __FILE__, __LINE__, #c); exit(1); } } while (0)
#ifndef HARNESS_ROUNDS
#define HARNESS_ROUNDS 2000000
#endif
#endif
