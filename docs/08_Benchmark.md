# NuExpansionHunter benchmark

NuExpansionHunter was **7.7–13.9× faster on average** and used **87–91% less peak memory** than ExpansionHunter v5.0.0 in these completed WGS comparisons. **All 22 comparisons passed VCF, JSON and BAM-record identity checks.**

Both tools used the same number of threads and the same CPU allocation for each comparison. Tests used streaming mode and the whole-genome GangSTR hg38 v13 catalog (832,152 loci) on eight AMD EPYC 9634 machines. Results are from the November 19, 2026 snapshot.

## Summary

Only comparisons where **both tools finished successfully** are included. Each row averages the same matched samples for both tools.

| Threads | Samples | v5 time (min) | NuEH time (min) | Speed-up | v5 peak RAM (GiB) | NuEH peak RAM (GiB) | RAM saved |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 9 | 204.8 | 26.7 | 7.69× | 73.9 | 9.7 | 86.9% |
| 16 | 9 | 114.2 | 13.5 | 8.50× | 91.2 | 11.5 | 87.1% |
| 32 | 4 | 100.6 | 7.2 | 13.94× | 157.8 | 14.3 | 90.9% |

Speed-up and memory saving are averaged per sample. The 32-thread row uses a smaller sample set and should not be read as a scaling comparison with the other rows. Eight additional v5 runs at 32 threads ran out of memory under the 172 GiB cap; these are excluded from the performance averages.

## Per-sample results

Times are in minutes; memory is peak RAM in GiB. Every row below passed all three output identity checks.

| Sample | Threads | v5 time (min) | NuEH time (min) | Speed-up | v5 peak RAM (GiB) | NuEH peak RAM (GiB) | RAM saved |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ERR1955394 | 8 | 176.97 | 23.82 | 7.43× | 72.02 | 8.72 | 87.9% |
| ERR1955394 | 16 | 96.78 | 12.05 | 8.03× | 84.05 | 10.55 | 87.5% |
| ERR1955394 | 32 | 85.35 | 6.07 | 14.07× | 155.39 | 14.25 | 90.8% |
| ERR1955395 | 8 | 220.63 | 29.22 | 7.55× | 72.21 | 8.73 | 87.9% |
| ERR1955395 | 16 | 125.70 | 14.75 | 8.52× | 82.20 | 10.58 | 87.1% |
| ERR1955395 | 32 | 107.12 | 7.42 | 14.43× | 164.21 | 14.26 | 91.3% |
| ERR1955396 | 8 | 205.10 | 27.54 | 7.45× | 72.09 | 8.74 | 87.9% |
| ERR1955396 | 16 | 109.23 | 13.88 | 7.87× | 72.51 | 10.54 | 85.5% |
| ERR1955396 | 32 | 94.58 | 6.95 | 13.61× | 148.10 | 14.15 | 90.4% |
| ERR1955398 | 8 | 240.88 | 33.46 | 7.20× | 72.10 | 8.91 | 87.6% |
| ERR1955398 | 16 | 128.87 | 16.87 | 7.64× | 74.74 | 10.75 | 85.6% |
| ERR1955398 | 32 | 115.33 | 8.45 | 13.64× | 163.64 | 14.40 | 91.2% |
| SRR21753324 | 8 | 188.32 | 23.60 | 7.98× | 73.54 | 9.59 | 87.0% |
| SRR21753324 | 16 | 108.97 | 11.99 | 9.09× | 107.77 | 11.40 | 89.4% |
| SRR21753325 | 8 | 173.67 | 22.35 | 7.77× | 72.76 | 9.47 | 87.0% |
| SRR21753325 | 16 | 103.62 | 11.33 | 9.15× | 105.85 | 11.27 | 89.4% |
| SRR21753326 | 8 | 177.63 | 22.40 | 7.93× | 72.65 | 9.07 | 87.5% |
| SRR21753326 | 16 | 101.43 | 11.30 | 8.98× | 101.25 | 10.90 | 89.2% |
| SRR21753327 | 8 | 206.30 | 25.72 | 8.02× | 80.88 | 12.62 | 84.4% |
| SRR21753327 | 16 | 110.98 | 13.06 | 8.50× | 82.34 | 14.42 | 82.5% |
| SRR21753328 | 8 | 253.42 | 32.19 | 7.87× | 76.89 | 11.25 | 85.4% |
| SRR21753328 | 16 | 141.83 | 16.29 | 8.71× | 110.27 | 13.05 | 88.2% |

Only GangSTR had completed matched comparisons in this snapshot. BAM identity refers to decoded records, not compressed file bytes. These results cover the completed comparisons, not the entire planned sample cohort.

## HG002 chromosome subsets

Each row compares one ExpansionHunter v5.0.0 run with one NuExpansionHunter
run. All runs used the same HG002 35× PCR-free BAM restricted to chromosomes
21 and 22 and the corresponding GRCh38 no-alt reference. The catalog subset
determines the loci evaluated.

| Catalog | Loci | v5 time (s) | NuEH time (s) | Speedup | v5 peak RSS (GiB) | NuEH peak RSS (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GangSTR hg38 v13, chr21–22 | 22,519 | 387.01 | 16.96 | 22.82× | 6.74 | 1.42 |
| TRExplorer v2.1, chr22 | 86,149 | 1,126.22 | 127.94 | 8.80× | 17.86 | 4.34 |
| TRExplorer v2.1, chr21–22 | 164,507 | 2,578.22 | 292.71 | 8.81× | 26.60 | 7.95 |

Both programs ran with `--analysis-mode streaming --sex male --threads 32`
and the same BAM, reference, and catalog within each pair. Runs were pinned
to logical CPUs `14-29,98-113` on an AMD EPYC 9634, spanning 16 physical cores
with simultaneous multithreading. The NuExpansionHunter release build used
GCC 15.2 with `-march=znver4`, `-ffp-contract=off`, libdeflate, ISA-L, and jemalloc.

GNU time measured wall time and process peak resident set size (RSS) for
each invocation, including catalog loading and output writing. Speedup
divides v5 wall time by NuEH wall time. All three pairs matched on VCF and
JSON byte digests, BAM header digests, and order-independent digests of
decoded BAM records.
