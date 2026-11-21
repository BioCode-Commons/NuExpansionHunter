<!--
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
-->

# Verification

For inputs with defined v5.0.0 behavior within the documented package limits,
the claim is identical VCF and JSON bytes, and identical BAM headers and multisets
of decoded records. BAM compression bytes and record order are outside the claim;
the unused low nibble of an odd-length sequence is immaterial to decoded records.
Evidence rows must fit a 16-bit best-cell index and fewer than 2^32 cells per locus.
The argument composes alignment equivalence (K-01–K-14), read equivalence
(K-15–K-23), locus equivalence (K-24–K-29, K-34–K-35), and output equivalence
(K-30–K-33, K-36). See [OBLIGATIONS.md](OBLIGATIONS.md) for individual contracts.

## Trusted base

The compiler and flags (`-ffp-contract=off`, without `-ffast-math`), C/C++ libraries
including `log` and `lgamma`, htslib 1.10.2 and its three package patches, zlib,
libdeflate, ISA-L, Boost, Abseil and spdlog are trusted. Compression is trusted only
through htslib's decoding contract. SHA-1 supplies catalog-cache integrity checks.
SIMD proofs trust the lane-wise intrinsic models in `simd/`; hardware cross-checks
connect those models to the instruction implementations. Testing a bounded case
set does not establish an unbounded equivalence claim.

## References and drivers

| Kind | Definition | Location |
| --- | --- | --- |
| R1 | Scalar pure function or restated v5 loop | Beside the optimized kernel |
| R2 | Executable v5 computation retained by the product | Dense matrix, general edge builder, owned catalog, serial paths |
| R3 | v5 text with namespace, include and necessary type adaptors | `reference/` |

The R3 matrix uses the R3 `StrAlign` calculator and filter, so it does not share the
optimized summary calculation. `JsonWriter_v5.hh` extracts the whole-document
`JsonWriter::write` body; its callable adaptor returns a string to the unit test.
The variant visitors remain shared. `check_reference.sh` records every upstream
path and this extraction explicitly, and prints the unfiltered unified diffs.
Copyright notices are preserved; namespace, includes, callable and borrowed-view
adaptors are the intended differences. The pinned wrapper uses the in-tree dense
matrix, whose separate obligations remain part of the trusted composition.

The C harnesses have two modes from one assertion body: proof mode assumes the
contract preconditions on symbolic inputs, and executable mode enumerates or
samples inputs. `harness/harness.h` selects the mode. The K-25 C++ driver is an
executable differential test: a seeded generator applies public add, remove,
filter and in-repeat padding operations and compares every row, best alignment,
and saturated allele lookup after each mutation. A symbolic K-25 driver remains
an open obligation. `probes/` supplies additional deterministic case sets and
in-process oracles. Internal second-best scores and nonoptimal sparse cells are
outside the alignment contract; canonical traceback output includes optimal paths
in their returned order.

## Build and run

After the parent project enables the opt-in subdirectory, configure its normal
package build with:

```sh
cmake -S ehunter -B build -DEH_BUILD_VERIFICATION=ON
cmake --build build
ctest --test-dir build -R '^verification_' --output-on-failure
```

Use the same dependency and toolchain options as the package build. A verification
build uses, for example, `-DEH_VERIFY="SIMD=SCALAR;BUDGET_SEARCH=NEVER"`.
`EH_VERIFY` is empty in a release build, and nonempty builds refuse installation.
Compile the libraries and probes with matching switches; changing only a probe's
flags does not switch a kernel compiled into a library.

```sh
bash ehunter/verification/scripts/run_proofs.sh
bash ehunter/verification/scripts/check_reference.sh
python3 ehunter/verification/scripts/run_probes.py --help
python3 ehunter/verification/scripts/run_probes.py \
  --package-build build --v5-build build-v5 --out probe-output \
  --verification-build 'SIMD=SCALAR=build-scalar' \
  --verification-build 'BUDGET_SEARCH=NEVER;FAST_PATHS=OFF=build-dense' \
  --verification-build 'FAST_PATHS=OFF;SIMD=SCALAR=build-fastpaths' \
  --verification-build 'BUDGET_SEARCH=ALWAYS;FUSED_ROWS=OFF=build-general' \
  --verification-build 'GRAPH_CACHE=OFF;FAST_PATHS=OFF=build-uncached' \
  --verification-build 'LAYOUT=REFERENCE=build-layout'
```

Prepare each supplied policy build with exactly the corresponding `EH_VERIFY`
value; all library paths, dependency includes, compiler and linker flags are
arguments or environment variables documented in the runner. The runner builds
only probes against these supplied libraries. It creates and removes a detached
v5.0.0 worktree for probes using v5 APIs. Package-only APIs use in-process oracles
and compare default with the specified policy build. One line per probe reports
`IDENTICAL`, `DIFFERS`, or `BUILD FAILED`; any failure gives a nonzero exit.
Commands and complete outputs remain under `--out`. Run builds and executable
checks under the shared build/run lock and the designated CPU affinity.

`dense_dp_cases` compares v5.0.0 only on inputs that meet the K-04/K-05
precondition: ordered DAG rows, where every predecessor row is lower than its
row. Its invalid-order family (self and future edges, every scoring setting) is
outside that precondition; v5.0.0's output on those inputs depends on prior
matrix storage, so neither implementation is a function of its input there. The
probe runs that family after the compared output, writes only a case count to
stderr, and checks only that it terminates; run it under AddressSanitizer for
out-of-bounds accesses (`dense_dp_cases invalid-order` runs the family alone).
`catalog_row_order <reads.bam> <reference.fa> <catalog.json>...` shows that
ExpansionHunter's own graphs meet the precondition: for every locus as the
program loads it, it runs the edge discovery that `PinnedDagAligner` uses,
forward and reversed, from every node at several start offsets and extensions,
and fails on any edge list without ordered rows. It is built with the
verification targets and is not registered with ctest because it needs data.

`run_neon_drivers.sh` takes `CC` (an aarch64 compiler), `RUNNER` (for example
`qemu-aarch64` on an x86-64 host) and `ROUNDS`.
`run_proofs.sh` accepts kernel names and the `CBMC`, `CC`, `CXX`, `TASKSET`,
`HTSLIB_INCLUDE`, and `HTSLIB_LIBS` environment variables. It writes proof and
hardware-model results under `results/`. The CMake targets register all ten
probes plus K-25 with ctest. Fixture generators in `fixtures/` accept their output
and converter paths as arguments; `generate_complex_fixtures.py` includes late
mismatches, complex indels, skewed loads and a heavy third locus.
`scripts/compare_outputs.py` compares writer outputs, while
`scripts/bam_multiset_compare.cpp` compares decoded BAM multisets. Build that
standalone helper with C++17 and zlib as shown in its header; it uses filesystem
and pooled string-view storage. The package probes and K-25 use C++11.

## EH_VERIFY switches

Use the name without `EH_VERIFY_` in the CMake option.

| Name | Values | Kernels and effect |
| --- | --- | --- |
| SIMD | SCALAR, SSE42, AVX2 | K-01, K-03 rows, K-09, K-19, K-22, K-23, K-31: pin dispatch. SCALAR: every kernel scalar. SSE42: the x86-64-v2 tier (SSE row kernels, K-19, K-23 and the k-mer probe; K-01, K-09, K-22 and K-31 scalar). AVX2: no AVX-512. AVX2 and AVX-512 still check CPU support |
| FAST_PATHS | OFF | K-09, K-11, K-12, K-17: general DP without a fast-path minimum-cost hint |
| BUDGET_SEARCH | NEVER, ALWAYS | K-07, K-10: dense fallback, or search without skip, measured cap or cost hint |
| ROW_STOP | OFF | K-05: fill all dense rows |
| FUSED_ROWS | OFF | K-03, K-08: general dense and sparse row code |
| GRAPH_CACHE | OFF | K-13: prepare the target afresh |
| FOLD_SKIP | OFF | K-21: execute every canonical fold step |
| LAYOUT | REFERENCE | K-15, K-28: hash k-mer layout and owned streaming catalog |
| SHADOW | ON | K-05, K-21, K-22, K-31: compare supported in-tree references and poison dropped dense rows |
| EXTENSION_CACHE | ON | Experimental cross-locus reuse of pinned single-node extension traces |

With `EXTENSION_CACHE=ON`, the streaming pipeline stores a bounded, in-process
table of optimal scores and ordered raw CIGARs for target alignments of shared
inputs. Hits require the same input pair, physical mate, read orientation,
extension direction, complete oriented query bytes and certified target prefix.
Donor and recipient independently pass the default-score, single-node and
full-target score-width guards. The receiving extension still remaps paths to
its graph. A busy bucket, a missing donor or a trace that exceeds the table's
capacity uses the ordinary extension alignment. A hit skips sparse/dense
search measurements and uses the optimum to update the ordinary cost hint.

For query length Q and donor score S, a linear pinned endpoint consuming t
target bases costs at least 2(t − Q) relative to the all-match score 5Q.
The cached prefix of length L = Q + floor((5Q − S)/2) + 1 therefore contains
every endpoint and predecessor on an optimal traceback. All later endpoints
also strictly exceed the optimum cost on the recipient. The cache preserves
the entire returned CIGAR vector and its order, including ties.

`EH_VERIFY="EXTENSION_CACHE=ON;SHADOW=ON"` also runs the uncached alignment
on every hit and checks the optimum score and the complete ordered raw vector.
`verification_extension_cache_cases` exercises both extension directions,
remapping across graphs, tied traces, byte guards and lock contention. The
application logs only aggregate table counts when this switch is enabled.

## K-mer index boundary profiles

`kmer_index_oracle` runs every profile by default; `audit PROFILE` selects one.
It compares independent path enumeration and the R3 index, checks mapped/copy/move
lifetime, and queries bucket boundaries and 100,000 seeded random keys per profile.
The reference-layout build expects hash storage throughout.

| Profiles | Boundary | Expected package representation |
| --- | --- | --- |
| `start2047`, `start2048` | Maximum single-node start position 2,047 vs. 2,048 | Prefix / fallback |
| `node3`, `node4` | Single-node path uses node ID 3 vs. 4 | Prefix / fallback |
| `dense14`, `dense10` | Homopolymer graph with many repeated path occurrences | Fallback |
| `skew14`, `skew10` | More than 32 keys share a prefix bucket | Fallback |
| `entries65535`, `entries65536` | Complete eight-base de Bruijn sequence, trimmed by one key or untrimmed | Fallback; 16-bit / 32-bit hash buckets |
| `complex14`, `complex10` | Four-node branching graph, two repeat loops, one ambiguous flank base | Prefix |
| `ambiguous14`, `ambiguous10` | Ambiguous repeat motif and flank bases; bounded expansion | Prefix |
| `tail14`, `tail10` | Last prefix bucket contains suffix 1; absent suffix 0 must not match padding | Prefix |
| `empty14` | No complete 14-mer exists | Fallback |
| `maximum16` | All-zero and all-one 32-bit keys, repeated path multiplicity | Fallback |

The largest sequence is 65,543 bases with at most 65,536 distinct keys.
Ambiguity is bounded to keep fixture construction finite.

## Proof status

`scripts/run_proofs.sh` runs each proof with CBMC 6.11 and the CaDiCaL solver (`--bounds-check --pointer-check --div-by-zero-check
--undefined-shift-check --unwinding-assertions`), then the same driver in executable mode on the host CPU, and
checks that one seeded fault in a copy of the kernel header makes the proof fail. The SIMD models were checked
against the hardware (AMD Zen 4, AVX2 and AVX-512BW) on 19.5 million comparisons with no mismatch. Logs are in
`results/`. OBLIGATIONS.md separates complete kernel proofs from case testing.

| Obligation | Tool | Bound | Result |
| --- | --- | --- | --- |
| K-01 AVX2 block | CBMC + SIMD models | 16 cells, all 18 inputs symbolic | proved (CaDiCaL, 6 min), seeded fault detected |
| K-01 AVX-512BW block | CBMC + SIMD models | 32 cells, all 34 inputs symbolic | proved (CaDiCaL, 110 min), seeded fault detected |
| K-03 fused row | CBMC | one lane, all 6 inputs symbolic | proved, seeded fault detected |
| K-03 row kernels, SSE2 | CBMC + SIMD models | 16 lanes, all inputs symbolic | proved (20 s), seeded fault detected |
| K-03 row kernels, AVX2 | CBMC + SIMD models | 16 lanes, all inputs symbolic | proved (21 s), seeded fault detected |
| K-09 linear scan | CBMC + SIMD models | 32-byte block + 1-base tail | proved, seeded fault detected |
| K-19 base run | CBMC + SIMD models | 16-byte block + 2-byte tail, both directions | proved, seeded fault detected |
| K-22 read decoding | CBMC + SIMD models | 32-base block + 1-base tail, exact-size arrays | proved, seeded fault detected |
| K-31 BAM packing | CBMC + SIMD models | 32-base vector + 1-base tail, exact-size arrays | proved, seeded fault detected |

NEON tested under emulation; intrinsic models not cross-checked on ARM hardware. The NEON kernels (K-01, K-03
rows, K-09, K-19, K-22, K-23 and K-31) live next to their scalar references in the same headers.
`scripts/run_neon_drivers.sh` builds the harnesses for aarch64 in executable mode, where each selects its NEON
kernel, and compares it with the reference on random and boundary inputs; K-23 is covered by
`OrientationPredictorTest` and `orientation_oracle`. On this project's x86-64 host they run under qemu-user, with
the aarch64 unit tests and an identity subset; no NEON proof is claimed. A block proof covers a whole row by
induction on the carried state: the state one block passes to the next is a state the harness starts from.
