# Stress timing regression audit

Captured: 2026-11-14T20:45:12.740484+00:00

This audit covers 13 completed stress cases and 26 worker configurations: 78 formal pairs and 26 excluded warmup pairs. Every completed pair is verified under the exact-output policy. All 13 stress cases are complete; the clock-qualified original T8 configuration remains disclosed. This is an observational audit, not an unconditional performance acceptance statement.

B = frozen RAMOPT `b00650d`; C = frozen v6. Positive changes mean that C is slower. Every comparison below was recomputed from individual formal pairs; stored summary ratios and confidence-interval inputs agree. Warmups are excluded from regression summaries.

## Metrics and comparability

- **End-to-end** is the stored process `wall_seconds`, measured with Python `time.perf_counter()` around process creation and `wait4`; it includes archive loading and process teardown. Preparation, output comparison and hashing are outside this measurement.
- **Read phase** is the application log field `read_processing_seconds`, computed with C++ `std::chrono::steady_clock`. It is elapsed wall time, not summed worker time or CPU time.
- **Paired change** is `100 × (median(C/B) − 1)`. B/C medians are independent medians and must not be used to recreate the paired change. **Paired delta** is `median(C − B)` in seconds. Three-pair bootstrap intervals are descriptive and do not prove an invariant performance difference.
- The read-phase cleanup boundary differs: RAMOPT declares `unpairedReads` and `readStreamer` at function scope and logs the phase before their destruction; v6 destroys its scoped input/source/pipeline before logging the phase. This is a confirmed scope difference, particularly relevant with many unmatched reads. It does not quantify how much of any observed difference is cleanup. End-to-end covers both versions' complete process lifetime.
- Short and empty-input ratios can amplify microsecond or millisecond differences. All positive observations remain disclosed; they are not dismissed as proof of noise.

## All positive median paired changes

| Case | T | Metric | B / C median (s) | Paired change | Paired delta (s) | 95% descriptive interval | Slower pairs |
| --- | ---: | --- | ---: | ---: | ---: | --- | ---: |
| `empty-input` | 1 | Read phase | 0.000105 / 0.000179 | +38.52% | +0.000050 | [+35.65%, +73.75%] | 3/3 |
| `empty-input` | 8 | Read phase | 0.000158 / 0.000171 | +9.87% | +0.000016 | [-12.58%, +15.80%] | 2/3 |
| `higherror-indels150` | 1 | End-to-end | 2.144888 / 2.123863 | +0.14% | +0.002949 | [-3.10%, +1.24%] | 2/3 |
| `higherror-indels150` | 1 | Read phase | 2.035567 / 2.012876 | +0.23% | +0.004824 | [-3.19%, +1.40%] | 2/3 |
| `higherror-indels150` | 8 | End-to-end | 0.851498 / 0.912797 | +7.20% | +0.061298 | [-7.50%, +11.48%] | 2/3 |
| `higherror-indels150` | 8 | Read phase | 0.774738 / 0.825073 | +6.63% | +0.051346 | [-7.84%, +10.83%] | 2/3 |
| `longread1000` | 1 | End-to-end | 0.790517 / 0.734940 | +0.61% | +0.004899 | [-20.20%, +3.40%] | 2/3 |
| `pending-large-far-apart` | 1 | End-to-end | 84.626849 / 85.454341 | +0.98% | +0.827491 | [-2.61%, +7.33%] | 2/3 |
| `pending-large-far-apart` | 1 | Read phase | 82.567784 / 82.912817 | +0.42% | +0.345033 | [-2.88%, +6.55%] | 2/3 |
| `pending-large-far-apart` | 8 | End-to-end | 34.038496 / 36.141901 | +6.18% | +2.103405 | [+3.49%, +6.76%] | 3/3 |
| `pending-large-far-apart` | 8 | Read phase | 33.117766 / 35.077498 | +5.92% | +1.959732 | [+3.37%, +6.59%] | 3/3 |
| `pending-million-missing-mates` | 1 | End-to-end | 0.621567 / 0.615825 | +0.83% | +0.005159 | [-24.95%, +1.17%] | 2/3 |
| `pending-million-missing-mates` | 1 | Read phase | 0.444087 / 0.550379 | +25.07% | +0.110314 | [-5.39%, +26.89%] | 2/3 |
| `pending-million-missing-mates` | 8 | Read phase | 0.438014 / 0.521046 | +20.02% | +0.086107 | [+18.90%, +24.13%] | 3/3 |

There are 6 positive end-to-end median comparisons and 8 positive read-phase median comparisons in the completed original campaign.

For `pending-million-missing-mates`, T1 is +0.83% end-to-end (+0.005159 s paired delta) but +25.07% read phase (+0.110314 s). T8 is −2.41% end-to-end (−0.014001 s) but +20.02% read phase (+0.086107 s), with all three T8 read-phase observations slower. Preserve both measurements and disclose the cleanup-boundary difference; do not claim that the entire +20–25% is additional alignment computation. This fixture contains no completed read pairs.

## Positive individual pairs despite nonpositive median change

| Case | T | Metric | Median paired change | Individual slower formal pairs (change; delta seconds) |
| --- | ---: | --- | ---: | --- |
| `empty-input` | 8 | End-to-end | -3.86% | 2: +18.96%; +0.007128 s |
| `longread1000` | 1 | Read phase | -0.62% | 3: +3.02%; +0.018321 s |
| `longread1000` | 8 | Read phase | -1.81% | 3: +10.66%; +0.035644 s |
| `mixed100k` | 1 | End-to-end | -5.70% | 1: +7.41%; +1.165560 s |
| `mixed100k` | 8 | End-to-end | -10.88% | 3: +0.68%; +0.066629 s |
| `pending-million-missing-mates` | 8 | End-to-end | -2.41% | 3: +0.85%; +0.005055 s |

## Timing-integrity exception

`genotyped10k-fixed`, T8, formal pair 3 (stored round 2), baseline has a raw read phase of **410.054144792 s** but process wall time of **15.160689833 s**. The UTC wrapper interval is **411.266445 s**, a **396.105755167 s** difference from process wall. Its phase sum is **411.083287625 s**. This is the only flagged run in the final integrity audit.

The original 1,291-byte stdout log itself contains `Phase read_processing_seconds=410.054144792`; stderr is empty. Logged timestamps span 12:01:52 to 12:08:43. The runner extracts the exact numeric literal with a regular expression, so this is not a parsing error. The duration code uses steady-clock interfaces in both paths, not subtraction of UTC/system-clock timestamps; a current Python clock query identifies `mach_absolute_time()` for `perf_counter`. The original runtime clock query was not captured. A host suspension or different elapsed-clock behavior is plausible, but the precise OS event is not established.

The affected stored read-phase paired changes are −50.97%, −49.72%, and **−97.80%**. The stored median is −50.97%, but the third sample is inconsistent and must not support a clean read-phase claim. The associated end-to-end/RSS sample remains original evidence but should be labeled as the same clock-affected run, rather than silently accepted because its numeric value looks plausible.

Publication treatment: all original JSON/logs and equality records are preserved, and the original T8 configuration is timing-qualified. The completed supplemental configuration below is a separate experiment; it is not pooled with selected original repetitions. The final integrity audit includes every completed run. No other configuration is flagged for clock inconsistency.

## Supplemental T8 recheck — separate experiment

Completed: three formal pairs plus one excluded warmup, all exact-output comparisons verified. These observations are separate from the 26 original stress configurations and are not pooled into their medians.

| Case | T | End-to-end change | Read-phase change | Evidence |
| --- | ---: | ---: | ---: | --- |
| `genotyped10k-t8-recheck-v6` | 8 | -53.65% | -56.59% | [JSON](evidence/results/genotyped10k-t8-recheck.json) |

| Metric | B / C median (s) | Paired change | Paired delta (s) | 95% descriptive interval | Slower pairs |
| --- | ---: | ---: | ---: | --- | ---: |
| End-to-end | 21.083803 / 9.772231 | -53.65% | -11.311572 | [-53.89%, -51.30%] | 0/3 |
| Read phase | 19.716700 / 8.559377 | -56.59% | -11.157323 | [-57.43%, -54.70%] | 0/3 |

Supplemental timing-integrity flags: 0.

## Integrity coverage and tolerances

[Machine-readable integrity snapshot](evidence/timing-integrity-audit.json) inspected 248 completed application records: 224 sample-analysis runs and 24 preparation runs, including the completed primary comparison, all original stress cases, and the completed supplemental T8 recheck. It flags exactly one preserved original run; the supplemental run has no flags.

- Each phase and their sum may exceed process wall only within `max(0.05 s, 2% × wall)` tolerance.
- UTC wrapper duration may differ from process wall within `max(2 s, 5% × wall)` because it also includes driver/metadata overhead.
- These are integrity diagnostics, not performance regression allowances. No flags does not establish exclusive host use or stable power/cache/memory conditions.
- The [audit script](evidence/audit_timing_integrity.py) reads small JSON files only and launches no benchmark, build or test.

## Original completed configuration overview

Values are median paired percentage changes. Positive means slower. The original clock-qualified configuration is explicitly marked.

| Case | T | End-to-end change | Read-phase change | Evidence |
| --- | ---: | ---: | ---: | --- |
| `dense-fanout-large` | 1 | -19.27% | -20.20% | [JSON](evidence/results/dense-fanout-large.json) |
| `dense-fanout-large` | 8 | -42.12% | -43.49% | [JSON](evidence/results/dense-fanout-large.json) |
| `empty-input` | 1 | -2.59% | +38.52% | [JSON](evidence/results/empty-input.json) |
| `empty-input` | 8 | -3.86% | +9.87% | [JSON](evidence/results/empty-input.json) |
| `expandedrepeat150-genotyped` | 1 | -7.42% | -8.28% | [JSON](evidence/results/expandedrepeat150-genotyped.json) |
| `expandedrepeat150-genotyped` | 8 | -49.01% | -50.36% | [JSON](evidence/results/expandedrepeat150-genotyped.json) |
| `genotyped10k-auto-io` | 8 | -51.23% | -55.05% | [JSON](evidence/results/genotyped10k-auto-io.json) |
| `genotyped10k-fixed` | 1 | -8.80% | -11.00% | [JSON](evidence/results/genotyped10k-fixed.json) |
| `genotyped10k-fixed` | 4 | -22.22% | -23.49% | [JSON](evidence/results/genotyped10k-fixed.json) |
| `genotyped10k-fixed` | 8 | -46.65% (clock-qualified) | -50.97% (clock-qualified) | [JSON](evidence/results/genotyped10k-fixed.json) |
| `higherror-indels150` | 1 | +0.14% | +0.23% | [JSON](evidence/results/higherror-indels150.json) |
| `higherror-indels150` | 8 | +7.20% | +6.63% | [JSON](evidence/results/higherror-indels150.json) |
| `hot-locus` | 1 | -33.00% | -41.68% | [JSON](evidence/results/hot-locus.json) |
| `hot-locus` | 8 | -50.40% | -63.42% | [JSON](evidence/results/hot-locus.json) |
| `iupac-missingquality150` | 1 | -8.80% | -9.17% | [JSON](evidence/results/iupac-missingquality150.json) |
| `iupac-missingquality150` | 8 | -9.94% | -11.25% | [JSON](evidence/results/iupac-missingquality150.json) |
| `longread1000` | 1 | +0.61% | -0.62% | [JSON](evidence/results/longread1000.json) |
| `longread1000` | 8 | -48.98% | -1.81% | [JSON](evidence/results/longread1000.json) |
| `mixed100k` | 1 | -5.70% | -14.08% | [JSON](evidence/results/mixed100k.json) |
| `mixed100k` | 8 | -10.88% | -55.05% | [JSON](evidence/results/mixed100k.json) |
| `pending-large-far-apart` | 1 | +0.98% | +0.42% | [JSON](evidence/results/pending-large-far-apart.json) |
| `pending-large-far-apart` | 8 | +6.18% | +5.92% | [JSON](evidence/results/pending-large-far-apart.json) |
| `pending-million-missing-mates` | 1 | +0.83% | +25.07% | [JSON](evidence/results/pending-million-missing-mates.json) |
| `pending-million-missing-mates` | 8 | -2.41% | +20.02% | [JSON](evidence/results/pending-million-missing-mates.json) |
| `repository-example` | 1 | -26.35% | -40.03% | [JSON](evidence/results/repository-example.json) |
| `repository-example` | 8 | -51.07% | -78.33% | [JSON](evidence/results/repository-example.json) |

## Source references

- Frozen RAMOPT `b00650d:ehunter/sample/HtsStreamingSampleAnalysis.cpp`: lines 304–309 create the read timer and function-scope pending/read-stream storage; lines 401–404 log the phase before function exit.
- Frozen-v6 working source `ehunter/sample/HtsStreamingSampleAnalysis.cpp`: lines 396–412 enclose input/source/pipeline cleanup before logging the phase.
- `tools/streaming_validation/measure.py`: `measure()` uses `time.perf_counter()` around `posix_spawn`/`wait4`.
- `tools/parallel_validation/run.py`: lines 150–155 parse phase literals from retained stdout/stderr; lines 293–318 store UTC wrapper timestamps around the measurement driver.
- [Clock-affected raw stdout](evidence/timing-anomalies/01-measure.stdout).

The audit used only small JSON/code/log reads. No product, runner, plan, result, or benchmark output was modified. Published evidence paths are relative to this report.
