# Streaming memory optimization: Jiapeng-RAMOPT

Baseline: `f3b1256286d5620a4dd8310786793a9c1913d841`. Candidate: the frozen v1 product implementation on this branch. Product source is unchanged from the benchmarked executable; publication adds documentation and portable paths.

This change compresses eligible exact k-mer indexes, replaces per-locus inline queues with pooled active queues and bounded shared mutexes, and releases streaming scratch/evidence earlier. Archive format is **v4**: rebuild existing catalog caches with the new executable.

Primary scope: **coordinate-sorted BAM, streaming, prepared catalog cache**. The data are synthetic diagnostics plus the repository example, measured on macOS ARM64. Uncached and unsorted rows are supplementary. The approximately 54% peak-RAM reduction applies to sparse100k; genotyped10k saves about 38–39%. There is no universal 50% memory or zero-runtime-regression guarantee.

All **17 datasets, 53 configurations and 185 measured pairs** preserve exact VCF/JSON bytes and the complete BAM header/full-record multiset, including duplicate multiplicity. Each configuration also has one excluded and output-checked warmup pair. Genotype-positive checks cover every formal/warmup output of the supplemental genotyping fixtures.

EH unit tests: baseline 167/168, candidate 172/173. Both have the same pre-existing `PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime` edge-label-order failure. All 30 GraphTools test executables pass. Queue TSan and index ASan/UBSan scope, extreme microbenchmark regressions and rejected variants are retained in the evidence.

- [Detailed HTML table, intervals, maximum peak RSS and cache preparation](benchmark-comparison.html) (download and open locally).
- [Full measurement and comparison records, JSON gzip](benchmark-data.json.gz).
- [Validation findings and limitations](validation-report.txt).
- [Reproduction instructions](reproduction.txt).
- [Evidence logs and component measurements](evidence/).

Public copies replace local paths with `<WORKSPACE>` / `<TMPDIR>`; measurement values and recorded input/output/binary hashes are unchanged. The compressed JSON contains records and hashes; generated BAM/FASTA and every raw output file are reproduced using the tools rather than checked into Git. Publication file checksums are in [manifest.json](manifest.json).

## End-to-end comparisons

Time is process start to exit; cache preparation is separate. RAM is the median of individual process peak RSS values, in MiB. Changes use the median of paired candidate/baseline ratios; **negative means lower** and may differ from dividing the two marginal medians. Confidence intervals and full ranges are in the HTML/JSON. Runs with intervals crossing parity do not establish a speed difference. Thread counts are compute/decode/encode workers, not physical-core limits. Additional sequence-preprocessing workers are zero.

| Dataset | Cache | BAM order | C/D/E | Pairs | Baseline s | Candidate s | Time change | Baseline MiB | Candidate MiB | RAM change |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| covered10k | mapped | sorted | 8/2/1 | 5 | 7.074 | 6.930 | -2.04% | 629.12 | 348.55 | -44.71% |
| dense-fanout | mapped | sorted | 1/2/1 | 3 | 1.147 | 1.160 | +0.46% | 43.08 | 42.58 | -2.59% |
| dense-fanout | mapped | sorted | 8/2/1 | 3 | 1.241 | 1.251 | +1.15% | 52.83 | 51.84 | -1.59% |
| dense-fanout-large | mapped | sorted | 1/2/1 | 9 | 5.848 | 6.005 | -5.28% | 135.34 | 137.11 | +1.18% |
| dense-fanout-large | mapped | sorted | 8/2/1 | 9 | 9.395 | 9.335 | +0.53% | 171.41 | 165.34 | -3.70% |
| empty-input | mapped | sorted | 1/2/1 | 3 | 0.048 | 0.038 | -22.83% | 19.44 | 15.88 | -18.33% |
| empty-input | mapped | sorted | 8/2/1 | 3 | 0.034 | 0.032 | -7.81% | 19.44 | 15.91 | -18.17% |
| expandedrepeat150 | mapped | sorted | 1/2/1 | 3 | 0.738 | 0.717 | -2.91% | 32.97 | 27.38 | -16.37% |
| expandedrepeat150 | mapped | sorted | 8/2/1 | 3 | 0.730 | 0.708 | +0.18% | 36.41 | 31.12 | -14.35% |
| expandedrepeat150-genotyped | mapped | sorted | 8/2/1 | 5 | 4.011 | 4.065 | +0.56% | 56.33 | 50.92 | -9.67% |
| far-apart | mapped | sorted | 1/2/1 | 3 | 34.590 | 33.970 | -3.52% | 72.42 | 66.86 | -7.68% |
| far-apart | mapped | sorted | 8/2/1 | 3 | 14.859 | 14.469 | -4.72% | 137.33 | 127.14 | -7.26% |
| genotyped10k | mapped | sorted | 8/2/1 | 5 | 19.965 | 19.403 | -4.77% | 724.17 | 441.39 | -39.02% |
| genotyped10k | mapped | sorted | 8/auto/auto | 5 | 19.720 | 19.234 | -4.14% | 729.36 | 447.64 | -38.43% |
| higherror-indels150 | mapped | sorted | 1/2/1 | 3 | 1.263 | 1.253 | -0.21% | 29.73 | 23.98 | -19.34% |
| higherror-indels150 | mapped | sorted | 8/2/1 | 3 | 0.519 | 0.514 | -0.56% | 35.31 | 29.31 | -16.84% |
| hot-locus | mapped | sorted | 1/2/1 | 3 | 0.421 | 0.411 | -1.32% | 38.86 | 35.31 | -9.13% |
| hot-locus | mapped | sorted | 8/2/1 | 3 | 0.366 | 0.358 | -0.49% | 39.11 | 35.44 | -9.39% |
| iupac-missingquality150 | mapped | sorted | 1/2/1 | 3 | 0.388 | 0.375 | -3.34% | 30.23 | 24.34 | -19.48% |
| iupac-missingquality150 | mapped | sorted | 8/2/1 | 3 | 0.190 | 0.190 | -0.46% | 33.80 | 27.88 | -17.03% |
| longread1000 | mapped | sorted | 1/2/1 | 3 | 0.455 | 0.436 | -3.49% | 30.61 | 24.36 | -20.42% |
| longread1000 | mapped | sorted | 8/2/1 | 3 | 0.585 | 0.577 | -1.76% | 67.28 | 61.59 | -7.72% |
| missing-mates | mapped | sorted | 1/2/1 | 3 | 0.091 | 0.077 | -15.05% | 72.28 | 68.67 | -5.01% |
| missing-mates | mapped | sorted | 8/2/1 | 3 | 0.077 | 0.074 | -4.01% | 72.48 | 68.91 | -4.93% |
| sparse100k | mapped | sorted | 8/2/1 | 5 | 11.899 | 7.858 | -35.26% | 5128.95 | 2362.11 | -54.08% |
| shuffled | mapped | unsorted | 1/2/1 | 3 | 1.240 | 1.256 | +3.55% | 60.59 | 58.28 | -3.82% |
| shuffled | mapped | unsorted | 8/2/1 | 3 | 1.191 | 1.267 | +9.79% | 71.00 | 63.92 | -5.80% |
| stripe-collisions | mapped | unsorted | 2/2/1 | 3 | 4.185 | 2.460 | -45.42% | 1496.83 | 687.30 | -54.08% |
| stripe-collisions | mapped | unsorted | 8/2/1 | 3 | 2.551 | 1.783 | -31.67% | 1502.84 | 697.22 | -53.60% |
| covered10k | owned | sorted | 8/2/1 | 5 | 5.948 | 6.118 | +1.72% | 860.48 | 612.61 | -28.06% |
| dense-fanout | owned | sorted | 1/2/1 | 3 | 1.198 | 1.361 | +15.39% | 43.53 | 43.84 | +1.85% |
| dense-fanout | owned | sorted | 8/2/1 | 3 | 0.614 | 0.614 | -1.23% | 54.62 | 52.91 | -1.86% |
| empty-input | owned | sorted | 1/2/1 | 3 | 0.147 | 0.143 | -2.91% | 19.94 | 16.98 | -14.81% |
| empty-input | owned | sorted | 8/2/1 | 3 | 0.070 | 0.071 | -0.79% | 20.94 | 18.06 | -13.73% |
| example | owned | sorted | 8/2/1 | 5 | 0.100 | 0.101 | +0.15% | 16.50 | 16.61 | +0.66% |
| expandedrepeat150 | owned | sorted | 1/2/1 | 3 | 0.885 | 0.879 | -0.59% | 32.52 | 28.03 | -13.79% |
| expandedrepeat150 | owned | sorted | 8/2/1 | 3 | 0.342 | 0.359 | +3.39% | 37.09 | 32.69 | -11.43% |
| far-apart | owned | sorted | 1/2/1 | 3 | 34.530 | 34.951 | +0.11% | 79.34 | 73.19 | -7.89% |
| far-apart | owned | sorted | 8/2/1 | 3 | 15.812 | 16.093 | +1.78% | 139.14 | 125.61 | -9.72% |
| higherror-indels150 | owned | sorted | 1/2/1 | 3 | 1.438 | 1.449 | +0.34% | 29.08 | 24.52 | -15.95% |
| higherror-indels150 | owned | sorted | 8/2/1 | 3 | 0.578 | 0.566 | -2.13% | 34.97 | 30.31 | -13.32% |
| hot-locus | owned | sorted | 1/2/1 | 3 | 0.511 | 0.500 | +2.88% | 39.11 | 36.02 | -7.95% |
| hot-locus | owned | sorted | 8/2/1 | 3 | 0.470 | 0.476 | -0.58% | 40.03 | 37.09 | -7.19% |
| iupac-missingquality150 | owned | sorted | 1/2/1 | 3 | 0.545 | 0.546 | +0.24% | 30.84 | 25.17 | -18.39% |
| iupac-missingquality150 | owned | sorted | 8/2/1 | 3 | 0.241 | 0.249 | +3.48% | 35.11 | 30.12 | -15.88% |
| longread1000 | owned | sorted | 1/2/1 | 3 | 0.599 | 0.591 | -0.06% | 29.64 | 24.86 | -16.04% |
| longread1000 | owned | sorted | 8/2/1 | 3 | 0.597 | 0.603 | +1.15% | 73.56 | 70.25 | -6.63% |
| missing-mates | owned | sorted | 1/2/1 | 3 | 0.250 | 0.266 | +0.28% | 73.05 | 70.12 | -4.00% |
| missing-mates | owned | sorted | 8/2/1 | 3 | 0.148 | 0.151 | +5.52% | 74.25 | 71.41 | -3.83% |
| shuffled | owned | unsorted | 1/2/1 | 3 | 1.360 | 1.491 | +10.73% | 62.12 | 60.55 | -2.37% |
| shuffled | owned | unsorted | 8/2/1 | 3 | 1.150 | 1.222 | +2.62% | 103.02 | 102.53 | +0.35% |
| stripe-collisions | owned | unsorted | 2/2/1 | 3 | 26.452 | 29.652 | +16.19% | 1407.19 | 793.33 | -45.30% |
| stripe-collisions | owned | unsorted | 8/2/1 | 3 | 14.470 | 14.472 | +0.02% | 1510.31 | 811.81 | -46.18% |

At compute=8, auto resolves to decode=8 and encode=8. `covered10k` is LowDepth (~8x); `genotyped10k` actually calls all 10,000 loci (~16x). `expandedrepeat150-genotyped` includes fragment-length calibration and calls all 200 loci; this checks execution, not biological genotype truth.
