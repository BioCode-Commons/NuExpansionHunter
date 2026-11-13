# Four-active-thread streaming regression benchmark

The comparison is between the former DP PR (`66200f6`), Happy-Friday (`d9173a3`), and a candidate. Use identical Release toolchain/dependencies and separate build directories. Every command is recorded with binary/input hashes. These are synthetic 150 bp diagnostics, not a real-WGS or cross-hardware speed guarantee.

## Inputs

Compile `tools/streaming_validation/sam_to_bam.c` against the same HTSlib as the application. Create the base fixture using:

```sh
python3 tools/frontend_validation/generate_fixture.py --kind mixed \
  --out /path/to/base-variable --bam-helper /path/to/sam_to_bam \
  --loci 1000 --covered-loci 1000 --target-pairs 50 \
  --background-pairs 450 --near-background-pairs 50 --variable-qualities
python3 tools/round4_validation/generate_complex_fixtures.py \
  --source-fixture /path/to/base-variable --out /path/to/complex-inputs \
  --sam-to-bam /path/to/sam_to_bam
python3 tools/round4_validation/generate_heavy_third.py \
  --source-fixture /path/to/base-variable --out /path/to/complex-inputs \
  --sam-to-bam /path/to/sam_to_bam
```

The generated cases are late high-quality mismatches, three-repeat multi-indel graphs, skewed heavy/light loci, and a heavy third locus that stresses reader-help scheduling. Every case contains exactly 10% coordinate-target reads, 90% background, and only 150 bp reads. Their manifests record independent SAM, pairing, coordinate and CIGAR checks. The heavy-third case is an artificial scheduling stress, not typical coverage.

For the earlier mixed-indel control, generate the same base command without
`--variable-qualities` into a different directory, then run
`generate_mixed_indels.py --source-fixture /path/to/base-fixed --out /path/to/mixed-indels --sam-to-bam /path/to/sam_to_bam`.
Pass that directory as another `--input`. The original ordinary-variable base can
also be included to check that gains on the difficult inputs preserve its speed.

## Measurements

```sh
python3 tools/round4_validation/benchmark_three.py \
  --pr /path/to/pr/ExpansionHunter --happy /path/to/happy/ExpansionHunter \
  --candidate /path/to/candidate/ExpansionHunter \
  --pr-source /path/to/pr-source --happy-source /path/to/happy-source \
  --candidate-source /path/to/candidate-source \
  --input /path/to/complex-inputs/late-mismatches \
  --input /path/to/complex-inputs/complex-indels \
  --input /path/to/complex-inputs/skewed-heavy-light \
  --input /path/to/complex-inputs/heavy-third \
  --out /path/to/new-results --cycles 9 --target-speedup 1.2
```

Run serially with no concurrent builds, generators, probes or benchmarks. The runner excludes one predetermined warmup triple, balances variant order, and reports paired median speed ratios plus bootstrap intervals. There is no page-cache clearing. It separately prepares compatible catalog caches for each binary before measurement.

All variants receive `--analysis-mode streaming --threads 2 --read-decode-threads 0 --decode-threads 0 --encode-threads 0 --region-extension-length 1000`. This budgets two compute workers, one reader and one BAM writer, including reader assistance. It does not enforce physical-core affinity on macOS.

The primary phase includes input, routing, graph alignment/DP and worker drain. Locus preparation is excluded. Read plus output time and the sum of all post-preparation phases are also reported. Peak RSS includes the entire process; both the median and maximum peak-RSS ratios must be at most 1.05 against each baseline.

After every triple, the runner checks VCF/JSON byte identity and the complete BAM header and raw-record multiset, preserving every field/tag and duplicate count. BAM record order and compressed bytes may differ. Comparisons and hashes are outside timed processes. A failed speed, output or RAM gate must be reported; a speed average cannot substitute for passing each input against both baselines.

## Exact algorithm and bounded storage

The candidate preserves the original scoring, seed rules, per-locus read order,
ordered traceback and canonical alignment. Its shortcuts are guarded proofs;
uncertainty, ties or a traversal cap return to the established DP.

For the default pinned alignment, rewrite score as `5 * query_length - cost`.
Matches cost 0, mismatches 9, ordinary deletion opening/extension 10/2 and
insertion opening/extension 15/7. Existing virtual-row boundary rules are kept.
An exhaustive capped ungapped search can return a unique cost-0/9 result, or
prove that no cost-at-most-9 result exists. Only the latter proof skips the first
budget attempt. A capped/incomplete search never supplies that proof.

A separate, single-node shortcut accepts exactly two strictly internal
mismatches only after ruling out every cost-at-most-18 single deletion of length
1..5 and every single insertion of length 1, using matching prefix/suffix lengths.
Five spare target bases keep all such competitors inside that node. Legacy
boundary I1+D1/D1+I1 paths cost 17 but skip only the first query/target base,
leaving both internal mismatches; they cannot compete. All remaining ambiguous,
boundary or out-of-range cases fall back.

Sparse DP fuses initialization, relaxation and support/endpoint collection for
single-predecessor rows (gapped widths at most 64). It computes the same H/E/F
values and insertion-equality flags. The previous operation-budget accounting
and dense fallback remain. A separate performance hint starts at budget 36
only after three qualifying DP calls have actual optimum cost above 18; a
qualifying cost at most 18 resets it. The hint is not a lower-bound assertion:
it searches a wider complete range containing every optimum and every tied
path. It stores one small counter and never memoizes read alignments.

Graph preparation is an owned four-slot LRU keyed by graph revision, start node,
start offset, extension and direction. Every graph mutation invalidates its
identity; copies initially share equivalent content. Each slot retains at most
64 KiB of actual container capacity, including EdgeMap buffers. Oversized
preparations remain local to that call. This bounds retained preparation data
at 256 KiB per alignment workspace, about 768 KiB for the configured three
processing roles; allocator overhead and transient allocations still require
RSS measurement.

The reader helps at most 32 read pairs before publishing the unique continuation
at the queue tail. Per-locus FIFO/sole ownership and the existing active-queue
limit stay intact. Workers still drain their locus to empty. No additional read
pairing cache or worker is introduced.

Additional changes eliminate temporary index lists and Path/Cigar copies,
project graph alignments in one operation pass, defer decoding of provably
unrouted second reads, and use equivalent common-base encoding paths. The dense
matrix now keeps its SIMD padding within vector size while separately retaining
the actual query length, including the previous padding-value history. This
fixes a pre-existing sanitizer-detected container overread.
