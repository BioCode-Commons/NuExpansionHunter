# Streaming benchmark results

The candidate reduces read-processing time while preserving the specified VCF/JSON and BAM outputs. Five difficult synthetic cases show paired median speedups of **1.243–1.734× over the earlier DP PR** and **1.127–1.279× over Happy-Friday**. The ordinary mixed-read control shows 2.227× and 1.742×, respectively. **The original requirement of at least 1.2× against both baselines on every difficult input is not fully met.**

The baselines are PR `66200f6` and Happy-Friday `d9173a3`. Correctness uses unmodified official ExpansionHunter v5.0.0 (`9f5db1d`) as an additional oracle. The measured candidate binary SHA256 is `aafe0d09728513ad86345b696304dd2a5580ab478780c14088c387a1883b0015`.

## Timing and memory

Only 150 bp reads, with 10% coordinate-target reads and 90% background. Streaming uses two compute workers, one reader and one BAM writer; additional read/HTS decode/encode threads are zero. This is a four-active-thread budget, without physical-core affinity on macOS ARM64. Read-processing includes graph alignment/DP and worker drain; loci preprocessing is excluded. Peak RSS includes the entire process.

Each fixture has one excluded warmup triple and nine serial, balanced-order triples. Builds, generators, probes and other benchmark jobs did not run concurrently. Page caches were not cleared. Host timing varied substantially; the following are paired medians, not cross-hardware guarantees. Speedup is the median of the per-cycle baseline/candidate time ratios, which need not equal a ratio of medians.

| Fixture | vs earlier PR | vs Happy-Friday |
|---|---:|---:|
| Mixed indels | 1.734× | 1.264× |
| Late mismatches | 1.581× | 1.178× |
| Multi-repeat complex indels | 1.491× | 1.127× |
| Skewed heavy/light loci | 1.417× | 1.279× |
| Heavy third locus | 1.243× | 1.153× |
| Ordinary mixed-read control | 2.227× | 1.742× |

All 54 measured pairs are faster than the earlier PR; 53/54 are faster than Happy-Friday. One heavy-third run is 1.87% slower than Happy-Friday. The first three-cycle screen had Happy-Friday about 9.2% slower than the old PR on that artificial scheduling stress; later nine-cycle measurements did not reproduce a median regression for Happy-Friday itself. That observation should not be represented as a deterministic causal guarantee.

For every fixture and both baselines, the median and maximum process peak RSS pass the 1.05× gate. The largest median increase is about 1.14%; the largest maximum-to-maximum increase is about 1.05%. All per-cycle phase times, peak RSS, output hashes, paired statistics and bootstrap intervals are in [benchmark_summary.json](benchmark_summary.json).

## Correctness and limitations

- VCF/JSON are byte-identical. The complete BAM header and raw-record multiset are identical, including all fields, tags and duplicate counts. Only record order and compression bytes may differ; the default compression level was not changed.
- All 73 final pipeline comparisons directly against official v5.0.0 pass (68 existing cases and five complex fixtures). All 120 comparisons against the two baselines in the formal runs, including warmups, pass.
- All 4,320 final orientation→seed→gapped cases (11,224 ordered alignments) match official v5 and d917. The two-mismatch shortcut has 2,500 directed comparisons with active hits and forced fallback coverage.
- Dense padding has 912 state comparisons against official v5/d917, including varying query lengths, selected default/nondefault scores, full padding values and empty-input state. Actual changed-source ASan/UBSan probes pass. This is not a whole-application sanitizer claim; some linked dependencies were not instrumented.
- GraphTools passes 30/30 executable suites. EH UnitTests passes 167/168. The sole failure, `PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime`, also occurs on unmodified d917 and concerns edge-label JSON ordering. The newly added padding regression was separately rebuilt and passes.
- Extreme-score, nonproduction local traceback has a separately identified pre-existing limitation; no universal equivalence claim is made for that API domain.
- Synthetic read-level perturbations are not consistent diploid truth, real WGS or clinical validation. Heavy-third is an artificial scheduling stress. No 250 bp/long-read pipeline performance or native x86 AVX2 performance was measured in this round.

See [README.md](README.md) for the exact algorithm, bounded storage and reproduction commands. Parameterized input generators reproduce the original five difficult fixtures: 30 reference/FAI/catalog/SAM/BAM/BAI hashes agree. Later proposals to relax the short-query threshold or directly construct single-gap results are not included.
