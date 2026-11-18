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
# usage: run_neon_drivers.sh
# Builds the kernel harnesses for aarch64 in executable mode, where they select the NEON kernels (K-01, K-03
# rows, K-09, K-19, K-22, K-31), and runs each against its scalar reference. On an x86-64 host set CC to an
# aarch64 cross compiler and RUNNER to a user-mode emulator (for example qemu-aarch64); on an aarch64 host the
# defaults (cc, no runner) run natively. ROUNDS sets the number of random rounds per driver.
set -u
V=$(cd "$(dirname "$0")/.." && pwd); E=$(dirname "$V"); G=$E/thirdparty/graph-tools-master-f421f4c/include
CC=${CC:-cc}; RUNNER=${RUNNER:-}; ROUNDS=${ROUNDS:-300000}
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
status=0
for h in K-01_insertion_scan K-03_rows K-09_linear_scan K-19_base_run K-22_read_decode K-31_bam_pack; do
  if ! $CC -O2 -Wall -Wextra -Werror -DHARNESS_ROUNDS=$ROUNDS -I"$G" -I"$E" -o "$W/$h" "$V/harness/$h.c" > "$W/$h.log" 2>&1; then
    echo "$h: BUILD FAILED"; cat "$W/$h.log"; status=1; continue
  fi
  if ! file "$W/$h" 2>/dev/null | grep -q aarch64; then echo "$h: not an aarch64 binary (set CC)"; status=1; continue; fi
  ( cd "$V/harness" && $RUNNER "$W/$h" ); rc=$?
  case $rc in 0) echo "$h: NEON driver passed ($ROUNDS rounds)";; *) echo "$h: NEON driver FAILED (rc $rc)"; status=1;; esac
done
exit $status
