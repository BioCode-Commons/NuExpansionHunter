#!/bin/bash
#
# NuExpansionHunter
# Copyright 2026 NuExpansionHunter contributors
#
# This file is part of NuExpansionHunter, a modified version of
# Illumina ExpansionHunter v5.0.0 (Copyright 2016-2019 Illumina, Inc.).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# usage: run_proofs.sh [kernel ...]   (default: all bit-precise proofs)
# Runs the complete bit-precise proofs with CBMC, then each harness in executable mode on this CPU, and the
# SIMD model cross-check. Results go to verification/results/<kernel>.log; the script prints one line per
# proof and exits non-zero if any proof or check fails.
# Environment: CBMC (default: cbmc on PATH), CXX/CC, HTSLIB_INCLUDE and HTSLIB_LIBS for the table generator
# (optional: the committed harness/tables.h is used when they are not set), TASKSET (optional CPU pinning).
set -u
V=$(cd "$(dirname "$0")/.." && pwd); E=$(dirname "$V"); G=$E/thirdparty/graph-tools-master-f421f4c/include
CBMC=${CBMC:-cbmc}; CC=${CC:-gcc}; CXX=${CXX:-g++}; PIN=${TASKSET:-}
R=$V/results; mkdir -p "$R"; W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
# Loops whose trip count depends on symbolic data are unwound to a bound, and the unwinding assertions prove
# that the bound is never reached, so the proofs stay complete.
FLAGS="--bounds-check --pointer-check --div-by-zero-check --undefined-shift-check --unwinding-assertions --object-bits 16 --sat-solver cadical"
declare -A PROOFS=(
  [K-01-avx2]="K-01_insertion_scan.c --unwind 70"
  [K-01-avx512]="K-01_insertion_scan.c -DK01_AVX512 --unwind 70"
  [K-03]="K-03_fused_row.c --unwind 70"
  [K-03-rows-sse2]="K-03_rows.c -DK03_TIER=1 --unwind 50"
  [K-03-rows-avx2]="K-03_rows.c -DK03_TIER=2 --unwind 50"
  [K-09]="K-09_linear_scan.c --unwind 70"
  [K-19]="K-19_base_run.c --unwind 20"
  [K-22]="K-22_read_decode.c --unwind 70"
  [K-31]="K-31_bam_pack.c --unwind 70"
)
# One seeded fault per proof, in a copy of the product kernel header: the proof must then fail. This checks
# that each proof constrains the code it names (a vacuous proof would still pass).
declare -A MUTANTS=(
  [K-01-avx2]="thirdparty/graph-tools-master-f421f4c/include/graphalign/dagAligner/DpKernels.hh|s/_mm256_set1_epi16(-10))/_mm256_set1_epi16(-9))/"
  [K-01-avx512]="thirdparty/graph-tools-master-f421f4c/include/graphalign/dagAligner/DpKernels.hh|s/_mm512_sub_epi16(preceding, _mm512_set1_epi16(10))/_mm512_sub_epi16(preceding, _mm512_set1_epi16(9))/"
  [K-03]="thirdparty/graph-tools-master-f421f4c/include/graphalign/dagAligner/DpKernels.hh|s/const int16_t match = (int16_t)(vpPrevious + penalty);/const int16_t match = (int16_t)(vpPrevious + penalty + 1);/"
  [K-03-rows-sse2]="thirdparty/graph-tools-master-f421f4c/include/graphalign/dagAligner/DpKernels.hh|s/_mm_add_epi16, _mm_max_epi16,/_mm_add_epi16, _mm_add_epi16,/"
  [K-03-rows-avx2]="thirdparty/graph-tools-master-f421f4c/include/graphalign/dagAligner/DpKernels.hh|s/DP_AVX2_STORE, _mm256_add_epi16,/DP_AVX2_STORE, _mm256_max_epi16,/"
  [K-09]="thirdparty/graph-tools-master-f421f4c/include/graphalign/dagAligner/DpKernels.hh|s/14, 13, 11, 12, 3, 15/14, 13, 10, 12, 3, 15/"
  [K-19]="thirdparty/graph-tools-master-f421f4c/include/graphutils/BaseMatchingKernels.hh|s/_mm_set1_epi8('g')/_mm_set1_epi8('n')/"
  [K-22]="core/ReadDecodeKernels.hh|s/_mm256_set1_epi8(20)/_mm256_set1_epi8(21)/"
  [K-31]="io/BamEncodingKernels.hh|s/_mm256_set1_epi8(40)/_mm256_set1_epi8(41)/"
)
# Logs name files relative to the repository, so that no host path reaches the package.
sanitize() { sed -i "s|$(dirname "$E")/||g; s|$W/[^ ]*/include/|mutant:|g; s|$W/|tmp/|g" "$1"; }
names=("$@"); [ ${#names[@]} -eq 0 ] && names=(K-01-avx2 K-01-avx512 K-03 K-03-rows-sse2 K-03-rows-avx2 K-09 K-19 K-22 K-31)
status=0
if [ -n "${HTSLIB_INCLUDE:-}" ]; then
  $CXX -O1 -I"$HTSLIB_INCLUDE" -I"$G" -o "$W/gen" "$V/scripts/gen_tables.cpp" ${HTSLIB_LIBS:-} && "$W/gen" > "$W/tables.h" \
    && cmp -s "$W/tables.h" "$V/harness/tables.h" && echo "tables.h: matches the product and htslib tables" \
    || { echo "tables.h: DIFFERS from the product or htslib tables"; status=1; }
fi
$CXX -O1 -mavx2 -mavx512f -mavx512bw -o "$W/crosscheck" "$V/simd/crosscheck_x86.cpp" && $PIN "$W/crosscheck" > "$R/simd-crosscheck.log" 2>&1
case $? in 0) echo "simd models: $(tail -1 "$R/simd-crosscheck.log")";; *) echo "simd models: FAILED (see results/simd-crosscheck.log)"; status=1;; esac
for name in "${names[@]}"; do
  set -- ${PROOFS[$name]}; src=$1; shift
  start=$(date +%s)
  ( cd "$V/harness" && $PIN $CBMC -I"$G" -I"$E" "$@" $FLAGS "$src" ) > "$R/$name.log" 2>&1; rc=$?
  secs=$(( $(date +%s) - start ))
  sanitize "$R/$name.log"
  summary=$(grep -E "^\*\* [0-9]+ of [0-9]+ failed|VERIFICATION (SUCCESSFUL|FAILED)" "$R/$name.log" | tail -1)
  if [ $rc -eq 0 ] && grep -q "VERIFICATION SUCCESSFUL" "$R/$name.log"; then echo "$name: proof PASSED in ${secs}s ($summary)"
  else echo "$name: proof FAILED (rc $rc, ${secs}s, $summary)"; status=1; fi
  IFS='|' read -r header expr <<< "${MUTANTS[$name]}"
  M=$W/mutant-$name; mkdir -p "$M/$(dirname "$header")" "$M/include/$(dirname "${header#*include/}")"
  target=$M/$header; [[ $header == thirdparty/* ]] && target=$M/include/${header#*include/}
  sed "$expr" "$E/$header" > "$target"
  if cmp -s "$E/$header" "$target"; then echo "$name: mutant NOT APPLIED"; status=1
  else
    ( cd "$V/harness" && $PIN $CBMC -I"$M/include" -I"$M" -I"$G" -I"$E" "$@" $FLAGS "$src" ) > "$R/$name.mutant.log" 2>&1
    sanitize "$R/$name.mutant.log"
    if grep -q "VERIFICATION FAILED" "$R/$name.mutant.log"; then echo "$name: seeded fault detected"
    else echo "$name: seeded fault NOT detected"; status=1; fi
  fi
  cdefs=(); for a in "$@"; do [[ $a == -D* ]] && cdefs+=("$a"); done
  $CC -O2 -I"$G" -I"$E" "${cdefs[@]}" -o "$W/$name" "$V/harness/$src" && $PIN "$W/$name"; rc=$?
  case $rc in 0) echo "$name: executable mode passed";; 77) echo "$name: executable mode skipped (CPU lacks the ISA)";;
    *) echo "$name: executable mode FAILED"; status=1;; esac
done
exit $status
