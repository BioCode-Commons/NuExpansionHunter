# Streaming frontend diagnostics

These tools measure the complete application in streaming mode while preserving
all VCF, JSON and realigned BAM outputs. They reuse the existing deterministic
catalog generator, isolated `posix_spawn`/`wait4` measurement, and exact BAM
header/full-record multiset comparator from `tools/streaming_validation`.

The fixtures are local simulations, **not WGS or expanded-allele truth sets**.
On a 10-logical-CPU Mac, `--threads 16 32` deliberately oversubscribes compute
workers; it does not simulate 16 or 32 physical cores or x86 AVX2 throughput.

Generate inputs serially, outside all timing and compilation:

```sh
python3 tools/frontend_validation/generate_fixture.py --kind frontend \
  --out /new/fixtures/frontend --bam-helper /path/sam_to_bam
python3 tools/frontend_validation/generate_fixture.py --kind target \
  --out /new/fixtures/target --bam-helper /path/sam_to_bam
python3 tools/frontend_validation/generate_fixture.py --kind outside \
  --out /new/fixtures/outside --bam-helper /path/sam_to_bam
python3 tools/frontend_validation/generate_fixture.py --kind candidate-negative \
  --out /new/fixtures/candidate-negative --bam-helper /path/sam_to_bam
```

Defaults are 10,000 independent loci, 1,000 covered loci, 50 target pairs per
covered locus (100,000 target records), 150 bp reads, mixed qualities, one
substitution per target read, and periodic missing qualities. `frontend` adds
2,000,000 records inside the coarse mask but outside exact extraction regions;
`outside` adds the same number with both mates outside the coarse mask. `target`
adds no background. All reference repeats are 10 CAGs. Extension length must
remain 1000. Use `--loci 100 --covered-loci 10 --target-pairs 10
--background-pairs 10` for lightweight correctness calibration.

`candidate-negative` starts with the same target layout and no extra background.
For pairs where `pair_index % 10 != 0`, it replaces both sequences with independent
deterministic ACGT strings while retaining coordinates, flags, CIGAR and qualities.
For a pair count divisible by ten, 90% of pairs are replaced and 10% retain the
original positive-control sequences. The manifest records actual selected pairs,
changed records and fractions. These reads deliberately disagree with their mapped
coordinates and may still pass orientation/seed tests by chance; this is a rejection
path diagnostic, not a claim of 90% measured rejection or realistic WGS error rates.

Run all processes serially:

```sh
python3 tools/frontend_validation/run_matrix.py \
  --baseline /path/baseline/ExpansionHunter --candidate /path/candidate/ExpansionHunter \
  --baseline-source /path/baseline-source --candidate-source /path/candidate-source \
  --input /new/fixtures/frontend --input /new/fixtures/target --input /new/fixtures/outside \
  --out /new/results --threads 16 32 --preprocess 0 --decode 2 --encode 1 --pairs 5
```

The default creates separate compatible prepared caches for each binary and
fixture, outside timed runs. Cache loading remains inside each process timing.
To reuse existing compatible archives, supply both `--baseline-cache PATH` and
`--candidate-cache PATH`. These paths are reused for every specified fixture;
the application validates compatibility on every run. The harness records each
archive's path, size and hash without interpreting its header. These options
cannot be combined with `--uncached`.
Each configuration has one excluded warmup followed by alternating AB/BA
pairs. Every pair compares all VCF/JSON bytes and complete BAM header/record
bytes, ignoring record order only and retaining duplicates. Warmups also check
equality across worker configurations. Raw outputs, logs, process wall/user/sys
time, peak RSS, phase timers, hashes and environment are retained. The summary
uses paired ratios and descriptive bootstrap intervals without an automatic
acceptance rule. Five pairs provide local evidence, not a general speed claim.

For a small correctness matrix, add `--preprocess 0 1`; repeat in a new output
directory with `--uncached` to exercise owned indexes. Run the bundled example
as another input only after naming its files `reference.fa`, `reference.fa.fai`,
`catalog.json` and `reads.bam`. The matrix does not drop BAM output or reduce
evidence thresholds. Do not run it concurrently with builds or other benchmarks.

`read_processing_seconds` includes HTS reading/decompression, decoding, pairing,
query/dispatch and asynchronous alignment completion. Its phase time is not a
DP-specific timer. Peak RSS is a process observation, not a lifetime memory
bound. Cold-cache performance is not measured because the OS cache is not reset.

When catalog initialization is amortized across samples, use
`summary.read_processing` as the primary comparison. It excludes catalog setup,
genotyping and the final output phase, including the final BAM writer drain;
retain total wall and output-phase measurements as secondary observations.

The expanded local run uses `--background-pairs 500` for `frontend` and `outside`
(10,100,000 records each), and `--covered-loci 10000` for `target` and
`candidate-negative` (1,000,000 records each). All retain 10,000 loci. Pass all
four input directories to the matrix with `--pairs 9 --threads 16 32
--preprocess 0 --decode 2 --encode 1`. Keep the initial smaller run separately;
do not pool its measurements with the expanded workloads.
