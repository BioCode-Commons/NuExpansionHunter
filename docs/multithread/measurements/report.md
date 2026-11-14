# Streaming parallel comparison report

Generated: 2026-11-14T20:45:15.350060+00:00

This is an observational report, not an acceptance declaration. B = baseline; C = candidate; T = compute workers, not a process-wide CPU cap. Positive paired wall changes are actual observed median slowdowns. One warmup pair per configuration is excluded from timing/RSS statistics but must be verified. Failed application runs are excluded from timing statistics; unverified completed pairs can appear and are identified separately. Comparisons with different baselines remain separate.

Seconds are end-to-end process wall time. Changes are the median of paired C/B ratios, not the ratio of medians. RSS is process peak memory; MiB = 1,048,576 bytes. Each memory allowance is baseline peak RSS + T × the recorded decimal-MB allowance (default 20,000,000 bytes per worker); both paired and median/maximum differences are checked. A partial result cannot establish a gate. Exact outputs require every VCF/JSON byte and complete BAM header/full-record multiset, preserving duplicates; only BAM record order is ignored. These checks are read from the runner, not repeated here.

Confidence intervals are copied only when the runner's stored paired ratios match the displayed pair set. Missing or stale intervals are shown as unavailable. Few-pair bootstrap intervals are descriptive and do not establish a hardware-general absence of regression.

## primary100k

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary100k | 20,000,000 | 100,000 | 8 | 3/3 | Complete | 184.03 / 90.11 | -51.03% | \[-56.67%, -46.21%\] | 3,776.23 / 3,839.41 | 3,720.89 / 4,142.92 | +366.69 / +381.70 | Allowance exceeded (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/01-primary100k.json); [Full raw metadata](data/01-primary100k.raw.json); [Fixture manifest](data/01-primary100k.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic, independent loci and distinct deterministic flanks; not a WGS performance claim

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 3 | 1,487.84 / 743.00 | -47.79% | 3 | 169.53 / 79.73 | -53.30% |

## genotyped10k-fixed

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| genotyped10k-fixed | 2,000,000 | 10,000 | 1 | 3/3 | Complete | 14.27 / 13.01 | -8.80% | \[-13.04%, -5.29%\] | 432.95 / 433.17 | 434.80 / 435.47 | +1.70 / +3.98 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| genotyped10k-fixed | 2,000,000 | 10,000 | 4 | 3/3 | Complete | 11.88 / 9.11 | -22.22% | \[-23.29%, -17.16%\] | 435.42 / 435.98 | 442.95 / 445.50 | +7.53 / +9.52 | Within observed allowance (+76.29 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/02-genotyped10k-fixed.json); [Full raw metadata](data/02-genotyped10k-fixed.raw.json); [Fixture manifest](data/02-genotyped10k-fixed.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic, independent loci and distinct deterministic flanks; not a WGS performance claim

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 33.49 / 28.30 | -15.50% | 3 | 12.06 / 10.73 | -11.00% |
| 4 | 3 | 65.88 / 48.88 | -25.11% | 3 | 10.70 / 7.93 | -23.49% |

## genotyped10k-auto-io

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| genotyped10k-auto-io | 2,000,000 | 10,000 | 8 | 3/3 | Complete | 21.83 / 10.90 | -51.23% | \[-52.58%, -49.67%\] | 440.53 / 442.75 | 446.30 / 449.52 | +5.78 / +6.77 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/03-genotyped10k-auto-io.json); [Full raw metadata](data/03-genotyped10k-auto-io.raw.json); [Fixture manifest](data/03-genotyped10k-auto-io.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=auto; encode=auto; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic, independent loci and distinct deterministic flanks; not a WGS performance claim

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 3 | 190.34 / 87.71 | -54.61% | 3 | 20.46 / 9.42 | -55.05% |

## repository-example

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| repository-example | Not recorded | 1 | 1 | 3/3 | Complete | 0.15 / 0.11 | -26.35% | \[-26.61%, -26.32%\] | 16.28 / 16.34 | 17.92 / 17.98 | +1.64 / +1.64 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| repository-example | Not recorded | 1 | 8 | 3/3 | Complete | 0.15 / 0.07 | -51.07% | \[-51.17%, -50.50%\] | 16.41 / 16.44 | 22.89 / 23.00 | +6.48 / +6.67 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/04-repository-example.json); [Full raw metadata](data/04-repository-example.raw.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Counts unavailable; BAM/catalog files were not scanned; loci supplied with --loci

Fixture scope: Not recorded

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 0.13 / 0.13 | +0.45% | 3 | 0.10 / 0.06 | -40.03% |
| 8 | 3 | 0.13 / 0.16 | +17.00% | 3 | 0.10 / 0.02 | -78.33% |

## longread1000

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| longread1000 | 25,600 | 200 | 1 | 3/3 | Complete | 0.79 / 0.73 | +0.61% slower | \[-20.20%, +3.40%\] | 24.25 / 24.33 | 26.02 / 26.03 | +1.73 / +1.78 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| longread1000 | 25,600 | 200 | 8 | 3/3 | Complete | 0.91 / 0.45 | -48.98% | \[-53.61%, -43.38%\] | 61.52 / 61.64 | 44.45 / 44.64 | -17.00 / -16.97 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/05-longread1000.json); [Full raw metadata](data/05-longread1000.raw.json); [Fixture manifest](data/05-longread1000.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic streaming input extremes; not representative WGS or a guaranteed alignment fast-path/fallback witness

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 2.09 / 1.75 | -9.39% | 3 | 0.68 / 0.63 | -0.62% |
| 8 | 3 | 3.48 / 2.01 | -42.11% | 3 | 0.33 / 0.33 | -1.81% |

## higherror-indels150

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| higherror-indels150 | 25,600 | 200 | 1 | 3/3 | Complete | 2.14 / 2.12 | +0.14% slower | \[-3.10%, +1.24%\] | 23.81 / 23.89 | 24.94 / 25.00 | +1.19 / +1.22 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| higherror-indels150 | 25,600 | 200 | 8 | 3/3 | Complete | 0.85 / 0.91 | +7.20% slower | \[-7.50%, +11.48%\] | 29.44 / 29.64 | 33.33 / 33.34 | +3.91 / +3.95 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/06-higherror-indels150.json); [Full raw metadata](data/06-higherror-indels150.raw.json); [Fixture manifest](data/06-higherror-indels150.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic streaming input extremes; not representative WGS or a guaranteed alignment fast-path/fallback witness

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 4.26 / 4.19 | +0.38% | 3 | 2.04 / 2.01 | +0.23% |
| 8 | 3 | 6.51 / 6.72 | +3.59% | 3 | 0.77 / 0.83 | +6.63% |

## iupac-missingquality150

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| iupac-missingquality150 | 25,600 | 200 | 1 | 3/3 | Complete | 0.63 / 0.58 | -8.80% | \[-12.13%, -1.99%\] | 24.31 / 24.41 | 25.20 / 25.27 | +0.95 / +0.98 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| iupac-missingquality150 | 25,600 | 200 | 8 | 3/3 | Complete | 0.28 / 0.26 | -9.94% | \[-10.77%, -3.60%\] | 28.00 / 28.14 | 30.67 / 31.33 | +2.64 / +3.56 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/07-iupac-missingquality150.json); [Full raw metadata](data/07-iupac-missingquality150.raw.json); [Fixture manifest](data/07-iupac-missingquality150.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic streaming input extremes; not representative WGS or a guaranteed alignment fast-path/fallback witness

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 1.23 / 1.09 | -12.17% | 3 | 0.51 / 0.47 | -9.17% |
| 8 | 3 | 1.78 / 1.59 | -11.25% | 3 | 0.19 / 0.17 | -11.25% |

## expandedrepeat150-genotyped

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| expandedrepeat150-genotyped | 105,600 | 200 | 1 | 3/3 | Complete | 5.19 / 4.81 | -7.42% | \[-13.71%, -4.97%\] | 46.06 / 46.12 | 49.31 / 49.47 | +3.22 / +3.58 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| expandedrepeat150-genotyped | 105,600 | 200 | 8 | 3/3 | Complete | 4.50 / 2.25 | -49.01% | \[-51.77%, -49.01%\] | 50.77 / 50.94 | 66.38 / 66.73 | +15.61 / +15.80 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/08-expandedrepeat150-genotyped.json); [Full raw metadata](data/08-expandedrepeat150-genotyped.raw.json); [Fixture manifest](data/08-expandedrepeat150-genotyped.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic streaming input extremes; not representative WGS or a guaranteed alignment fast-path/fallback witness

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 10.63 / 9.71 | -8.69% | 3 | 5.01 / 4.59 | -8.28% |
| 8 | 3 | 38.40 / 18.02 | -53.06% | 3 | 4.41 / 2.15 | -50.36% |

## dense-fanout-large

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dense-fanout-large | 50,000 | 24 | 1 | 3/3 | Complete | 10.88 / 8.88 | -19.27% | \[-19.33%, -17.79%\] | 136.22 / 136.61 | 141.09 / 142.08 | +5.08 / +5.47 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| dense-fanout-large | 50,000 | 24 | 8 | 3/3 | Complete | 10.77 / 6.30 | -42.12% | \[-42.41%, -40.24%\] | 164.67 / 167.20 | 177.05 / 185.08 | +12.67 / +20.41 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/09-dense-fanout-large.json); [Full raw metadata](data/09-dense-fanout-large.raw.json); [Fixture manifest](data/09-dense-fanout-large.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Bounded synthetic streaming extreme diagnostics; not real WGS or genotype truth.

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 25.86 / 19.78 | -24.48% | 3 | 10.01 / 8.10 | -20.20% |
| 8 | 3 | 88.53 / 55.22 | -37.35% | 3 | 10.51 / 6.00 | -43.49% |

## hot-locus

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| hot-locus | 40,000 | 128 | 1 | 3/3 | Complete | 0.48 / 0.33 | -33.00% | \[-33.68%, -32.19%\] | 35.30 / 35.33 | 31.80 / 32.06 | -3.53 / -3.23 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| hot-locus | 40,000 | 128 | 8 | 3/3 | Complete | 0.45 / 0.24 | -50.40% | \[-56.97%, -47.82%\] | 35.45 / 36.39 | 42.70 / 43.02 | +6.62 / +7.25 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/10-hot-locus.json); [Full raw metadata](data/10-hot-locus.raw.json); [Fixture manifest](data/10-hot-locus.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Bounded synthetic streaming extreme diagnostics; not real WGS or genotype truth.

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 0.59 / 0.60 | -1.30% | 3 | 0.39 / 0.23 | -41.68% |
| 8 | 3 | 0.58 / 1.06 | +69.85% | 3 | 0.37 / 0.16 | -63.42% |

## empty-input

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| empty-input | 0 | 128 | 1 | 3/3 | Complete | 0.05 / 0.05 | -2.59% | \[-7.05%, -0.80%\] | 15.83 / 15.86 | 15.78 / 15.78 | -0.05 / -0.03 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| empty-input | 0 | 128 | 8 | 3/3 | Complete | 0.04 / 0.04 | -3.86% | \[-7.24%, +18.96%\] | 15.92 / 15.92 | 15.81 / 15.84 | -0.11 / -0.03 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/11-empty-input.json); [Full raw metadata](data/11-empty-input.raw.json); [Fixture manifest](data/11-empty-input.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Bounded synthetic streaming extreme diagnostics; not real WGS or genotype truth.

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 0.03 / 0.03 | -1.64% | 3 | 0.00 / 0.00 | +38.52% |
| 8 | 3 | 0.04 / 0.04 | +4.11% | 3 | 0.00 / 0.00 | +9.87% |

## pending-million-missing-mates

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pending-million-missing-mates | 2,097,152 | 128 | 1 | 3/3 | Complete | 0.62 / 0.62 | +0.83% slower | \[-24.95%, +1.17%\] | 858.59 / 858.61 | 795.22 / 795.22 | -63.39 / -63.36 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| pending-million-missing-mates | 2,097,152 | 128 | 8 | 3/3 | Complete | 0.60 / 0.57 | -2.41% | \[-3.95%, +0.85%\] | 858.77 / 858.92 | 795.42 / 795.42 | -63.36 / -63.34 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/12-pending-million-missing-mates.json); [Full raw metadata](data/12-pending-million-missing-mates.raw.json); [Fixture manifest](data/12-pending-million-missing-mates.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Bounded synthetic streaming extreme diagnostics; not real WGS or genotype truth.

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 0.99 / 0.97 | -0.13% | 3 | 0.44 / 0.55 | +25.07% |
| 8 | 3 | 0.97 / 0.94 | -1.78% | 3 | 0.44 / 0.52 | +20.02% |

## pending-large-far-apart

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pending-large-far-apart | 262,144 | 128 | 1 | 3/3 | Complete | 84.63 / 85.45 | +0.98% slower | \[-2.61%, +7.33%\] | 100.83 / 117.08 | 102.84 / 116.48 | +2.45 / +15.66 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| pending-large-far-apart | 262,144 | 128 | 8 | 3/3 | Complete | 34.04 / 36.14 | +6.18% slower | \[+3.49%, +6.76%\] | 209.64 / 211.86 | 234.67 / 237.42 | +25.03 / +25.56 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/13-pending-large-far-apart.json); [Full raw metadata](data/13-pending-large-far-apart.raw.json); [Fixture manifest](data/13-pending-large-far-apart.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Bounded synthetic streaming extreme diagnostics; not real WGS or genotype truth.

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 164.87 / 155.63 | -5.60% | 3 | 82.57 / 82.91 | +0.42% |
| 8 | 3 | 262.74 / 257.30 | -2.07% | 3 | 33.12 / 35.08 | +5.92% |

## mixed100k

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mixed100k | 2,200,000 | 100,000 | 1 | 3/3 | Complete | 17.66 / 16.65 | -5.70% | \[-9.75%, +7.41%\] | 2,381.73 / 2,381.75 | 2,380.80 / 2,384.06 | -0.94 / +2.31 | Within observed allowance (+19.07 MiB) | 3/3 verified; warmup verified |
| mixed100k | 2,200,000 | 100,000 | 8 | 3/3 | Complete | 11.00 / 9.84 | -10.88% | \[-17.72%, +0.68%\] | 2,381.89 / 2,381.91 | 2,395.44 / 2,397.09 | +15.19 / +19.48 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/14-mixed100k.json); [Full raw metadata](data/14-mixed100k.raw.json); [Fixture manifest](data/14-mixed100k.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic, independent loci and distinct deterministic flanks; not a WGS performance claim

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 20.54 / 18.14 | -9.70% | 3 | 1.94 / 1.66 | -14.08% |
| 8 | 3 | 44.21 / 34.62 | -22.74% | 3 | 1.99 / 0.90 | -55.05% |

## genotyped10k-t8-recheck

Input run status: complete; declared complete: true.

| Fixture | Records | Loci | T | Pairs measured/requested | State | Wall B / C (s) | Paired wall change | 95% CI of change | RSS B median / max (MiB) | RSS C median / max (MiB) | Paired RSS delta median / largest (MiB) | RAM allowance | Exact outputs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| genotyped10k-t8-recheck | 2,000,000 | 10,000 | 8 | 3/3 | Complete | 21.08 / 9.77 | -53.65% | \[-53.89%, -51.30%\] | 441.19 / 441.81 | 449.05 / 449.25 | +8.06 / +9.25 | Within observed allowance (+152.59 MiB) | 3/3 verified; warmup verified |

[Compact case JSON](data/15-genotyped10k-t8-recheck.json); [Full raw metadata](data/15-genotyped10k-t8-recheck.raw.json); [Fixture manifest](data/15-genotyped10k-t8-recheck.manifest.json).


Baseline SHA-256: 2a74137bfa9e123a1fbb7e4f3374b9b296c676588b53839649fd45d3ca11188b.

Candidate SHA-256: 26f3f58964e42c0d9b332f7c0ab2fb77fe95bc3e0db67a3e8c7992a0a5b9095b.

Timing policy: Cache preparation, hashes, warmup, comparison and cleanup are excluded. Archive loading remains timed. No OS page-cache clearing. All application processes run serially on the same inherited CPU mask.

Requested I/O settings: decode=2; encode=1; BAM queue override=version defaults.

Fixture counts: Fixture manifest; records supplied with --records; loci supplied with --loci

Fixture scope: Synthetic, independent loci and distinct deterministic flanks; not a WGS performance claim

### CPU and read-phase details

Read-phase boundaries differ: b00650d logs the phase before unmatched-read and input-stream destruction; v6 logs it after the input scope closes. These diagnostic elapsed times are not strictly equivalent intervals. Pending-mate cleanup may contribute to an observed difference; no entire delta is attributed to cleanup. End-to-end process wall time is the primary comparison.

CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker alignment/reduction elapsed counters are not CPU time and are not substituted here.

| T | CPU pairs | CPU B / C (s) | Paired CPU change | Read-phase pairs | Read phase B / C (s) | Paired read-phase change |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 3 | 179.73 / 74.93 | -58.80% | 3 | 19.72 / 8.56 | -56.59% |
