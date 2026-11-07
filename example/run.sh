#!/bin/bash
# Modified for NuExpansionHunter (2026): run the example from any working directory.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p output

"${EXPANSIONHUNTER:-ExpansionHunter}" \
  --reads input/variants.bam \
  --reference input/reference.fa \
  --variant-catalog input/variants.json \
  --output-prefix output/repeats
