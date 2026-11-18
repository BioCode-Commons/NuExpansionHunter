#!/bin/bash
#
# Modified by NuExpansionHunter contributors, 2026.
# This file differs from ExpansionHunter v5.0.0: the binary and the output directory are arguments, and the result is compared with output/.
# Modifications are licensed under the Apache License, Version 2.0.
#
# Runs ExpansionHunter on the example data and compares the result with output/.
# output/ holds the output of ExpansionHunter v5.0.0 (the Linux x86-64 release
# binary) for the command below, run in this directory with
# --output-prefix output/repeats.
#
# usage: example/run.sh [binary] [output directory]
#   binary            default: $EXPANSIONHUNTER, else build/install/bin/ExpansionHunter
#   output directory  default: a new temporary directory
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
binary=${1:-${EXPANSIONHUNTER:-$here/../build/install/bin/ExpansionHunter}}
out=${2:-$(mktemp -d)}
if [ ! -x "$binary" ]; then
  echo "ExpansionHunter not found at $binary; pass the binary as the first argument or set EXPANSIONHUNTER" >&2
  exit 1
fi
mkdir -p "$out"
out=$(cd "$out" && pwd)
cd "$here"

"$binary" \
  --reads input/variants.bam \
  --reference input/reference.fa \
  --variant-catalog input/variants.json \
  --output-prefix "$out/repeats"

status=0
for file in repeats.vcf repeats.json; do
  if cmp -s "$out/$file" "output/$file"; then
    echo "$file: identical to output/$file"
  else
    echo "$file: differs from output/$file"
    status=1
  fi
done
# BAM files are compared by header and records, not by compressed bytes.
if command -v samtools > /dev/null; then
  bam_digest() { { samtools view -H "$1" | grep -v '^@PG'; samtools view "$1" | LC_ALL=C sort; } | cksum; }
  if [ "$(bam_digest "$out/repeats_realigned.bam")" = "$(bam_digest output/repeats_realigned.bam)" ]; then
    echo "repeats_realigned.bam: same header and records as output/repeats_realigned.bam"
  else
    echo "repeats_realigned.bam: differs from output/repeats_realigned.bam"
    status=1
  fi
else
  echo "repeats_realigned.bam: not compared (samtools not found)"
fi
echo "Output written to $out"
exit $status
