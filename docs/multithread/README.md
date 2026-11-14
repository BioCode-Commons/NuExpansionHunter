# Streaming parallel execution: final validation report

**The frozen v6 benchmark campaign is complete.** On the primary workload of
100,000 repeat loci and 20 million BAM records, v6 reduced median paired
end-to-end time by 51.03% and logged read-processing time by 53.30%.
Complete outputs matched. Peak RSS increased in two of three primary pairs;
the maximum observed peak rose by 303.516 MiB. The measurements do not
establish a RAM reduction or absence of runtime regressions: six selected
configurations had median paired end-to-end slowdowns.

All 13 original stress cases and the independent T8 timing recheck completed.
The main statistics cover 27 configurations, 81 formal pairs, and 27 warmup
pairs. Including the retained original T8 runs, 84 formal pairs and 28 warmup
pairs were executed: 224 sample-analysis processes, all in pairs that passed
complete-output comparison. Further optimization and profiling stopped at the
user's request. Historical component experiments below remain separate and
are not pooled with this final comparison.

[Complete measurement overview](measurement-overview.md) ·
[Every executed formal analysis](individual-formal-runs.md) ·
[Detailed measurements and raw-data links](measurements/report.md) ·
[Machine-readable matrix](measurement-matrix.json) ·
[Complete regression audit](regression-audit.md)

## Versions and scope

| Role | Identity | Executable SHA-256 |
|---|---|---|
| Fixed RAMOPT baseline | `b00650d297754f68d78778a441e2fcbf5096250b` | `2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b` |
| Frozen v6 candidate | `Jiapeng-multithread`, changes based on `b00650d` | `26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b` |

The baseline executable was frozen before the RAMOPT commit was made; its
recorded product-source files match that commit. The candidate source-file
hashes identify the benchmarked working tree. Freezing v6 for measurement does
not establish that its memory change has succeeded. The original experimental
manifest remains historical provenance. [Final provenance manifest](manifest.json)

The benchmark scope is **coordinate-sorted BAM, streaming mode, and cached
sample analysis**. Each binary uses a matching prepared catalog archive.
Reported RAM is the **whole sample-analysis process peak RSS**, including
loading and using the prepared archive. Read-processing time is also reported
separately. Catalog/index preparation time and RAM are excluded from the
sample-analysis comparison. Hashing, output comparison, and fixture/output cleanup are outside the
measured application process; application-owned input cleanup remains timed.

The original extra-memory reference was 20,000,000 bytes per compute worker;
the user later allowed a modest excess if disclosed accurately. This is a
comparison reference, not a process RAM ceiling or an established outcome.
The original 50% memory-reduction target was not achieved in these measurements.
The primary-workload speedup also does not satisfy the no-runtime-regression
requirement across all configurations; all observed slowdowns remain visible
in the stress results. MB denotes 1,000,000 bytes; MiB denotes 1,048,576 bytes.

## Execution design and compatibility

One input owner preserves primary-record admission, consecutive admitted-QNAME
pairing, mate roles, routing coordinates, and quality-dependent base case.
Pending raw sequence/quality storage is capped at 2,000,000 bytes; additional
unmatched reads retain exact decoded sequences. Missing and distant mates are
not discarded to enforce an artificial memory limit. The pending-mate table
still grows with unmatched names: the raw-payload cap is not a bound on its
total storage.

A shared executor performs initialization, sequence decoding, alignment, and
genotyping. Alignment chunks can run independently; each locus consumes its
results in original encounter order. Scheduling storage belongs to active
chunks rather than every catalog locus. Per-worker alignment workspaces are
constructed independently. A separate publisher workspace supports assistance
when admission would otherwise wait.

The v6 change releases a consumed item's alignments and oriented reads before
processing the next item. It preserves the retained evidence and serialized
BAM records. Its effect on process RSS is not established by the screening
experiment below.

Reduction batches transfer owned BAM records to a dedicated consumer with a
default 8 MiB streaming queue budget. Write failure wakes blocked producers;
non-batched seeking writes retain deferred error reporting through `finish()`.
Borrowed adjacency accessors avoid archive-owner reference-count traffic
inside synchronous traversals, while public owning accessors retain their
lifetime guarantee. Streaming selectors optionally capture one locale per
analysis; other callers retain fresh-global-locale behavior.

`--threads` counts compute workers. Codec pools and coordination threads are
additional. The fixed-codec cases use two decode and one encode worker;
automatic codec settings are a separate comparison. Queue/cache charges and
per-worker budgets are accounting estimates, not allocator-exact RSS limits.
The BAM budget excludes producer-local batches and HTS buffers; one oversized
batch can run alone so that large records are preserved.
[Design and diagnosis](diagnosis-and-design.md)

## Measurement method and limits

The local host is macOS ARM64 with 10 reported logical CPUs. CPU affinity was
not enforced, and exclusive host use or stable OS memory pressure was not
established. Benchmark application processes run serially; owned builds and
fixture-generation jobs are stopped during measurement. This controls our
jobs without establishing that the rest of the host is idle.

The runner counterbalances version order, records warmups separately, and
reports wall time, user-plus-system CPU time, read-processing time, and peak
RSS. A paired ratio means candidate divided by baseline within one pair;
the reported median paired ratio is not necessarily the ratio of the two
independent medians. Warmups are excluded from formal summaries. Short or
empty-input cases are primarily compatibility checks, not reliable speed
ratios. Few-pair intervals are descriptive rather than hardware-wide guarantees.

**Read-phase boundaries differ.** RAMOPT logs the end of read processing
before destroying its unmatched-read storage and BAM reader. The candidate
logs that phase after closing the input scope, including that cleanup. Raw
phase times remain reported without adjustment, but they do not isolate the
same work. In particular, the observed roughly 20–25% longer missing-mate
read phases, with similar end-to-end time, must not be interpreted entirely
as slower alignment computation. They also must not be removed from the
reported measurements. End-to-end time includes cleanup in both versions.

**A timing anomaly was retained and independently rechecked.** In the original
`genotyped10k-fixed` T8 configuration, zero-based round 2, the baseline's logged
read phase was 410.054 seconds while the external `perf_counter` interval was
15.1607 seconds; the UTC span was 411.266 seconds. Its cause has not been
established. The final integrity audit checked all 224 sample-analysis runs
and 24 separate preparation runs and flagged only this run. The complete
original T8 configuration remains in the evidence, including all three formal
pairs and its warmup. Main statistics use all three formal pairs and the
warmup from a separate T8 recheck, rather than selecting individual favorable
repetitions. No recorded measurements were changed. The recheck passed exact
output comparison and had no timing-integrity flags.
[Selection policy and source files](evidence/selection-provenance.json) ·
[Timing-integrity audit](evidence/timing-integrity-audit.json) ·
[Original fixed-codec results](evidence/results/genotyped10k-fixed.json) ·
[Independent T8 recheck](evidence/results/genotyped10k-t8-recheck.json)

Large RSS variation was observed for the same frozen binary and input. In
the v5/v6 screening, v5 peak RSS was 3,765.234 MiB in warmup and 3,028.141 MiB
in the measured run. A concurrent host snapshot reported about 5.66 GiB
occupied by the memory compressor; this is host-wide context, not memory
attributable to either application. Per-process anonymous/file-backed/
compressed-memory decomposition was not performed. The optional sampler
remains uncompiled and unrun. These limitations neither prove nor dismiss a
candidate RAM regression; the observed process peaks remain reported.
[Screening host snapshot](evidence/screening-host-memory.json)

Local results cannot establish 16-to-32-core scaling on the HG002 server or
prove that the published HappyFriday regression is removed. Most fixtures
are synthetic. Their agreement establishes preservation relative to the
baseline, not genotype accuracy on real WGS.

## Exact output policy

Each completed pair must match all VCF bytes and all JSON bytes. BAM comparison
requires identical complete header bytes and an identical multiset of complete
record bytes, including auxiliary fields and duplicate multiplicity. Only
BAM record order may differ. The comparator performs an exact external sort;
diagnostic hashes are not the sole equality check. Successful exit or matching
record counts alone are insufficient. No record cutoff, sampling, or dropped
mate is used to make comparisons pass.

## Final v6 primary comparison — complete

The primary fixture contains 100,000 covered repeat loci and 20,000,000 BAM
records, with 100 read pairs per locus and 150-base reads. The completed
configuration used eight compute workers, fixed codecs 2/1, one warmup pair,
and three measured pairs. The following summary excludes warmup.

| Final primary metric | RAMOPT baseline | Frozen v6 |
|---|---:|---:|
| Median end-to-end wall time, seconds | 184.0256 | 90.1127 |
| Median read-processing time, seconds | 169.5318 | 79.7250 |
| Median CPU time, seconds | 1,487.8427 | 742.9978 |
| Median whole cached-analysis peak RSS, MiB | 3,776.234375 | 3,720.890625 |
| Maximum observed whole cached-analysis peak RSS, MiB | 3,839.406250 | 4,142.921875 |
| Range of whole cached-analysis peak RSS, MiB | 3,339.187500–3,839.406250 | 3,453.765625–4,142.921875 |

All three measured pairs were faster in both end-to-end and read-processing
time. The median within-pair candidate/baseline ratios were **0.489675 for
wall time** and **0.467038 for read-processing time**, corresponding to 51.03%
and 53.30% paired reductions. The ratio of the independent read-processing
medians is 0.470266, or a 52.97% reduction; this is a different summary and
must not be labeled the paired reduction. The median paired CPU ratio was
0.522131, distinct from the ratio of independent CPU medians.

| Measured pair | RAMOPT peak RSS, MiB | v6 peak RSS, MiB | v6 − RAMOPT, MiB |
|---|---:|---:|---:|
| 0 | 3,339.187500 | 3,720.890625 | +381.703125 |
| 1 | 3,839.406250 | 3,453.765625 | −385.640625 |
| 2 | 3,776.234375 | 4,142.921875 | +366.687500 |

**Peak RSS increased in two of three measured pairs.** Both increases exceed
the original 160,000,000-byte reference allowance at eight compute workers.
Although the candidate's independent median RSS is lower, its maximum
observed RSS is 303.515625 MiB higher and its median paired RSS delta is
+366.687500 MiB. The observed variation does not establish a RAM reduction or
justify treating the original memory reference as met in every run.

All three measured pairs and the separate warmup passed the complete-output
policy: identical JSON/VCF bytes, identical complete BAM headers, and equal
multisets of all **20,000,000 BAM records** per application run. Compressed BAM
file bytes differed, as allowed by the record-order policy.

A separate genotype/evidence check on measured pair 0 passed for both
versions: all 100,000 loci were called, all 100,000 VCF repeat records were
called, and every locus reported repeat evidence. Coverage ranged from
15.97297 to 16.21622, above the validator's minimum of 10. This confirms
nonempty genotype computation on the synthetic fixture; it is not an
independent biological-accuracy assessment.

[Final primary measurements](evidence/final-v6-primary.json) ·
[Primary genotype/evidence validation](evidence/final-v6-primary-genotype-validation.json)

## Final stress campaign — complete

All 13 original cases and 26 original thread configurations completed, with
one warmup pair and three formal pairs per configuration. The additional T8
recheck contributes another three formal pairs and one warmup. The main
statistics substitute that entire recheck for the original anomalous T8
configuration. Together with the primary case, this yields 14 fixture/settings
cases and 27 selected configurations; fixed and automatic I/O are separate
settings of the same 10,000-locus fixture. Total executed counts retain both
T8 attempts. Cache-preparation processes are excluded from those analysis
counts and from the time/RAM comparison.

Every original case and the recheck passed complete-output comparison for
all formal and warmup pairs. All inputs were coordinate-sorted. Fixed-codec
cases used two decode and one encode worker; `genotyped10k-auto-io` used each
version's automatic codec settings. The table below contains every selected
stress configuration. B is RAMOPT; C is frozen v6. Time changes are medians
of paired C/B ratios minus one, not ratios of independent medians. RSS columns
show the maximum whole-process peak observed among the three formal runs.
Read-phase cleanup boundaries differ between versions, as detailed above;
end-to-end time covers complete process cleanup in both.

| Stress case/settings | Compute workers | E2E median B / C, s | Paired E2E change | Read median B / C, s | Paired read change | Peak RSS max B / C, MiB |
|---|---:|---:|---:|---:|---:|---:|
| `genotyped10k-fixed` | 1 | 14.269 / 13.013 | -8.80% | 12.056 / 10.730 | -11.00% | 433.172 / 435.469 |
| `genotyped10k-fixed` | 4 | 11.879 / 9.112 | -22.22% | 10.699 / 7.927 | -23.49% | 435.984 / 445.500 |
| `genotyped10k-auto-io` | 8 | 21.835 / 10.900 | -51.23% | 20.458 / 9.418 | -55.05% | 442.750 / 449.516 |
| `repository-example` | 1 | 0.150 / 0.111 | -26.35% | 0.098 / 0.059 | -40.03% | 16.344 / 17.984 |
| `repository-example` | 8 | 0.151 / 0.074 | -51.07% | 0.099 / 0.021 | -78.33% | 16.438 / 23.000 |
| `longread1000` | 1 | 0.791 / 0.735 | +0.61% | 0.682 / 0.625 | -0.62% | 24.328 / 26.031 |
| `longread1000` | 8 | 0.910 / 0.454 | -48.98% | 0.334 / 0.327 | -1.81% | 61.641 / 44.641 |
| `higherror-indels150` | 1 | 2.145 / 2.124 | +0.14% | 2.036 / 2.013 | +0.23% | 23.891 / 25.000 |
| `higherror-indels150` | 8 | 0.851 / 0.913 | +7.20% | 0.775 / 0.825 | +6.63% | 29.641 / 33.344 |
| `iupac-missingquality150` | 1 | 0.632 / 0.580 | -8.80% | 0.513 / 0.467 | -9.17% | 24.406 / 25.266 |
| `iupac-missingquality150` | 8 | 0.281 / 0.256 | -9.94% | 0.195 / 0.173 | -11.25% | 28.141 / 31.328 |
| `expandedrepeat150-genotyped` | 1 | 5.195 / 4.809 | -7.42% | 5.007 / 4.592 | -8.28% | 46.125 / 49.469 |
| `expandedrepeat150-genotyped` | 8 | 4.496 / 2.253 | -49.01% | 4.409 / 2.150 | -50.36% | 50.938 / 66.734 |
| `dense-fanout-large` | 1 | 10.882 / 8.877 | -19.27% | 10.010 / 8.105 | -20.20% | 136.609 / 142.078 |
| `dense-fanout-large` | 8 | 10.766 / 6.296 | -42.12% | 10.512 / 5.999 | -43.49% | 167.203 / 185.078 |
| `hot-locus` | 1 | 0.480 / 0.325 | -33.00% | 0.386 / 0.228 | -41.68% | 35.328 / 32.062 |
| `hot-locus` | 8 | 0.455 / 0.237 | -50.40% | 0.374 / 0.157 | -63.42% | 36.391 / 43.016 |
| `empty-input` | 1 | 0.051 / 0.050 | -2.59% | 0.000105 / 0.000179 | +38.52% | 15.859 / 15.781 |
| `empty-input` | 8 | 0.041 / 0.039 | -3.86% | 0.000158 / 0.000171 | +9.87% | 15.922 / 15.844 |
| `pending-million-missing-mates` | 1 | 0.622 / 0.616 | +0.83% | 0.444 / 0.550 | +25.07% | 858.609 / 795.219 |
| `pending-million-missing-mates` | 8 | 0.595 / 0.572 | -2.41% | 0.438 / 0.521 | +20.02% | 858.922 / 795.422 |
| `pending-large-far-apart` | 1 | 84.627 / 85.454 | +0.98% | 82.568 / 82.913 | +0.42% | 117.078 / 116.484 |
| `pending-large-far-apart` | 8 | 34.038 / 36.142 | +6.18% | 33.118 / 35.077 | +5.92% | 211.859 / 237.422 |
| `mixed100k` | 1 | 17.658 / 16.652 | -5.70% | 1.936 / 1.663 | -14.08% | 2,381.750 / 2,384.062 |
| `mixed100k` | 8 | 10.998 / 9.844 | -10.88% | 1.992 / 0.901 | -55.05% | 2,381.906 / 2,397.094 |
| `genotyped10k-t8-recheck` | 8 | 21.084 / 9.772 | -53.65% | 19.717 / 8.559 | -56.59% | 441.812 / 449.250 |

**Six configurations had positive median paired end-to-end changes:**
`longread1000` T1 (+0.61%), `higherror-indels150` T1 (+0.14%) and T8 (+7.20%),
`pending-million-missing-mates` T1 (+0.83%), and `pending-large-far-apart` T1
(+0.98%) and T8 (+6.18%). These observations are retained even where absolute
times are short or the independent medians move in the opposite direction.
The far-apart T8 median end-to-end time increased from 34.038 to 36.142 seconds.

Eight configurations had positive median paired read-phase changes:
both high-error configurations, both empty-input configurations, both
missing-mate configurations, and both far-apart configurations. Missing-mate
read phases increased 25.07% at T1 and 20.02% at T8; the cleanup-boundary
difference described above may contribute, but the entire delta is not
attributed to cleanup. Empty-input phases were roughly 0.1–0.2 milliseconds,
so their percentages do not support a useful compute-performance conclusion.
End-to-end time includes cleanup in both versions.

All selected stress configurations stayed within the original observed
20,000,000-byte-per-worker extra-memory reference; this does not mean their
RSS decreased. The primary case exceeded that reference, as shown above.
Every per-pair increase, RSS range, and formal maximum remains available in
the complete matrix and detailed report. Output equality and completed
measurements establish neither universal performance acceptance nor a
hardware-wide absence of regression.
[Stress plan and fixture sizes](evidence/final-v6-stress-plan.json) ·
[Completed original campaign](evidence/final-v6-stress-results.json) ·
[Complete time/RSS matrix](measurement-overview.md) ·
[Individual formal runs](individual-formal-runs.md) ·
[Detailed report with raw evidence](measurements/report.md)

## Completed correctness evidence

The latest application suite ran 239 tests: **238 passed and one pre-existing
test failed**. `PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime`
compares graph JSON containing the same edge labels in different orders
(`a,b` versus `b,a`). The same failure appears in the earlier baseline and
RAMOPT evidence. It is disclosed rather than counted as a full-suite pass.
[Latest unit log](evidence/unit-v6.log) ·
[Earlier baseline failure](evidence/preexisting-baseline-test-failure.log)

All 30 GraphTools test programs passed in the latest GraphTools run. That run
precedes the v6 item-release change, which changes the streaming reducer and
does not change GraphTools. [GraphTools log](evidence/graph-tests.log)

| Component validation | Instrumentation and scope | Completed result |
|---|---|---|
| Executor and publisher assistance | Standalone executor/test translation units; ASan and TSan; additional UBSan scope unverified because exact compile flags were not retained | 19 tests × 20 repetitions per configuration; 760 total case executions |
| BAM sink, including deferred seeking errors | Standalone sink/test translation units; ASan+UBSan and TSan | 11 tests × 20 repetitions per configuration; 440 total case executions |
| Locale snapshot lifetime and semantics | Current aligner/test translation units; ASan+UBSan | 6 focused cases passed |
| Borrowed adjacency and owning-range lifetime | Focused Graph tests under ASan | 22 cases passed |

These component runs reported no sanitizer diagnostics. They are not whole-
application sanitizer coverage, and the final v6 item-release edit was not
included in those earlier sanitizer binaries. Available commands and manifests
accompany the evidence; provenance gaps are listed per component. Exact compile commands and validation-time source hashes
were not retained for the executor and borrowed-adjacency runs. The locale
command manifest records instrumented sources but not their validation-time
hashes. These missing records are not reconstructed as exact provenance.
[Component sanitizer evidence](evidence/component-sanitizers.json)

A seeking-mode smoke comparison of RAMOPT against v5 used the repository
example. JSON and VCF bytes matched; the full BAM header and all 1,098 records
matched as a multiset. This is a completed v5 compatibility smoke test, not a
seeking performance campaign or a new v6 measurement. The v6 edit is confined
to streaming reduction. [Seeking smoke comparison](evidence/seeking-v5.json)

## Appendix A: historical RAMOPT versus v4

This completed experiment used the same 100,000-locus/20-million-record primary
fixture and eight compute workers with fixed codecs 2/1. Its candidate was v4,
SHA-256 `fac4e0369a4df15b7e876bfa5490dda495a35ebe53cf5d1c1a9d6ee54b0520ab`.
All three measured pairs and the separate warmup matched the full output
policy. These measurements must not be pooled with the final v6 comparison.

| Historical v4 metric | RAMOPT median | v4 median | Median within-pair ratio |
|---|---:|---:|---:|
| Wall time, seconds | 153.828 | 81.212 | 0.519202 |
| Read-processing time, seconds | 143.520 | 70.672 | 0.487697 |
| CPU time, seconds | 1,309.810 | 651.139 | 0.495389 |
| Whole cached-analysis peak RSS, MiB | 3,431.422 | 3,909.750 | 1.141590 |

The three paired peak-RSS increases were **484.922, 143.797, and 512.734 MiB**.
All measured v4 runs were faster, and all used more peak RSS than their paired
RAMOPT runs. The median paired wall-time reduction was 48.08%; this is a v4
result, separate from the completed v6 primary comparison above.
[Historical v4 measurements](evidence/historical-v4-primary.json)

End-of-read logs reported peak charged pipeline storage of approximately
19.8 MB, at most 32 active chunks, pending raw payload of 23,200 bytes, and an
8,388,576-byte peak BAM sink charge. These are internal accounting counters,
not RSS at the end of reading. No per-phase OS RSS was measured. Sink wait
seconds sum blocked producer intervals; sink drain time excludes the later
HTS close/compression flush.

## Appendix B: v5 versus v6 item-release screening

This is a separate component experiment. Its baseline is **v5**, SHA-256
`f304d6f09e2d4b4a9b50541b551abae4307641bcb597b0a6507c65700b18b6b9`, rather
than RAMOPT. It has one measured pair and one separate warmup on the primary
fixture, eight compute workers, and fixed codecs 2/1. Both pairs matched the
complete output policy.

| Screening run | v5 wall, s | v6 wall, s | v5 peak RSS, MiB | v6 peak RSS, MiB | v6 − v5 RSS, MiB |
|---|---:|---:|---:|---:|---:|
| Warmup, v6 then v5 | 100.000 | 70.749 | 3,765.234 | 3,549.359 | −215.875 |
| Measured, v5 then v6 | 106.842 | 93.291 | 3,028.141 | 3,684.188 | +656.047 |

In the measured pair, read-processing time was 93.441 versus 82.010 seconds;
CPU time was 783.103 versus 784.310 seconds. Wall time decreased 12.68%, but
the RSS direction reversed between warmup and measurement. **This screening
does not establish a RAM reduction from early item release.** It also cannot
be combined with Appendix A to estimate v6 versus RAMOPT: the baselines,
sampling counts, and run conditions differ.
[v5/v6 screening measurements](evidence/v5-v6-screening.json)

## Measurement outcome

The frozen v6 comparison is complete, with exact full outputs across the
primary workload, all original stress cases, and the independent T8 recheck.
The large primary workload became substantially faster, but its peak RSS
regressed in two of three formal pairs. Six selected configurations also had
median paired end-to-end slowdowns. The current evidence therefore does not
meet the original combination of 50% lower memory and no runtime regression.
It also does not establish that the original extra-memory reference is met
in every run.

These conclusions apply to the recorded local cached-analysis measurements.
They retain the known unit-test failure, sanitizer provenance gaps, clock
anomaly and complete recheck, read-phase boundary difference, and host-memory
variability described above. Real HG002 testing on the intended 32-core
allocation remains necessary before claiming that the published scaling
regression has been removed. No unperformed memory diagnostics are used to
explain away the observed RSS increases.
