# NuExpansionHunter v1.1.0 versus ExpansionHunter v5.0.0

Across 22 completed whole-genome comparisons, mean speedups for
NuExpansionHunter v1.1.0 ranged from **10.09× to 18.79×** by thread count.
Mean peak RSS reductions ranged from **86.9% to 90.9%**.
All 22 comparisons passed the four output identity checks.

## Summary

Each row uses the same samples for both versions. Wall time and peak RSS are
arithmetic means. Speedup (v5 wall time divided by v1.1.0 wall time) and RSS
reduction (one minus v1.1.0 peak RSS divided by v5 peak RSS) are calculated for
each pair and then averaged within that thread count.

| Threads | Samples | v5 wall (min) | v1.1.0 wall (min) | Mean speedup | v5 peak RSS (GiB) | v1.1.0 peak RSS (GiB) | Mean RSS reduction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 9 | 204.8 | 20.3 | 10.09× | 73.9 | 9.7 | 86.9% |
| 16 | 9 | 114.2 | 10.3 | 11.06× | 91.2 | 11.6 | 87.1% |
| 32 | 4 | 100.6 | 5.3 | 18.79× | 157.8 | 14.3 | 90.9% |

The 32-thread group covers four samples, so comparisons across thread counts
also reflect the change in sample mix. Nine of the 21 B1 ExpansionHunter
v5.0.0 runs at 32 threads hit the 172 GiB memory cap on November 18–19, 2026.
Five of those samples appear in the 8- and 16-thread groups. The performance
averages cover completed pairs.

## Per-sample results

Each row represents one run of each version. Every row matched the stored
ExpansionHunter v5.0.0 digests for the VCF, JSON, BAM header, and decoded BAM
record multiset. The multiset check is independent of record order and
compression.

| Sample | Threads | v5 wall (min) | v1.1.0 wall (min) | Speedup | v5 peak RSS (GiB) | v1.1.0 peak RSS (GiB) | RSS reduction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ERR1955394 | 8 | 176.97 | 18.04 | 9.81× | 72.02 | 8.78 | 87.8% |
| ERR1955394 | 16 | 96.78 | 9.16 | 10.56× | 84.05 | 10.61 | 87.4% |
| ERR1955394 | 32 | 85.35 | 4.64 | 18.41× | 155.39 | 14.31 | 90.8% |
| ERR1955395 | 8 | 220.63 | 21.21 | 10.40× | 72.21 | 8.79 | 87.8% |
| ERR1955395 | 16 | 125.70 | 10.78 | 11.66× | 82.20 | 10.64 | 87.1% |
| ERR1955395 | 32 | 107.12 | 5.46 | 19.62× | 164.21 | 14.34 | 91.3% |
| ERR1955396 | 8 | 205.10 | 20.36 | 10.07× | 72.09 | 8.81 | 87.8% |
| ERR1955396 | 16 | 109.23 | 10.33 | 10.58× | 72.51 | 10.61 | 85.4% |
| ERR1955396 | 32 | 94.58 | 5.19 | 18.21× | 148.10 | 14.25 | 90.4% |
| ERR1955398 | 8 | 240.88 | 23.95 | 10.06× | 72.10 | 8.97 | 87.6% |
| ERR1955398 | 16 | 128.87 | 12.15 | 10.61× | 74.74 | 10.81 | 85.5% |
| ERR1955398 | 32 | 115.33 | 6.09 | 18.93× | 163.64 | 14.47 | 91.2% |
| SRR21753324 | 8 | 188.32 | 18.90 | 9.96× | 73.54 | 9.66 | 86.9% |
| SRR21753324 | 16 | 108.97 | 9.66 | 11.28× | 107.77 | 11.48 | 89.3% |
| SRR21753325 | 8 | 173.67 | 17.70 | 9.81× | 72.76 | 9.53 | 86.9% |
| SRR21753325 | 16 | 103.62 | 9.06 | 11.43× | 105.85 | 11.35 | 89.3% |
| SRR21753326 | 8 | 177.63 | 17.59 | 10.10× | 72.65 | 9.14 | 87.4% |
| SRR21753326 | 16 | 101.43 | 8.99 | 11.29× | 101.25 | 10.98 | 89.2% |
| SRR21753327 | 8 | 206.30 | 19.46 | 10.60× | 80.88 | 12.70 | 84.3% |
| SRR21753327 | 16 | 110.98 | 10.00 | 11.10× | 82.34 | 14.52 | 82.4% |
| SRR21753328 | 8 | 253.42 | 25.38 | 9.98× | 76.89 | 11.31 | 85.3% |
| SRR21753328 | 16 | 141.83 | 12.89 | 11.01× | 110.27 | 13.12 | 88.1% |

The archived GNU time records give these v5 measurements for the four
8-thread controls. The per-sample ratios and 8-thread summary use these
values before rounding.

| Sample | v5 wall (s) | v5 peak RSS (KiB) |
| --- | ---: | ---: |
| ERR1955394 | 10618 | 75516396 |
| ERR1955395 | 13238 | 75716152 |
| ERR1955396 | 12306 | 75592208 |
| ERR1955398 | 14453 | 75597428 |

## Method notes

- NuExpansionHunter v1.1.0 ran on November 23, 2026. Its executable SHA-256 was
  `536fef0e287122a554f71cd5a1af6be77e59d25ddb1dcc852c251af3788a4ed8`.
  The portable Release build used libdeflate, ISA-L, and jemalloc.
- The ExpansionHunter v5.0.0 wall times and peak RSS values come from the
  earlier matched runs reported on November 19, 2026.
- Both versions used streaming mode on eight AMD EPYC 9634 machines with the
  GangSTR hg38 v13 whole-genome catalog (832,152 loci) and the GRCh38 reference.
  Each v1.1.0 run used the same worker, input BAM, `--threads` value, and CPU
  mask as its v5 run. When setting up the comparisons, the v5 masks for 14 of
  the 22 pairs were reconstructed from the B1 runner’s fixed allocation rule.
  Completed B1 job tables confirm four of those masks, all for 8-thread
  controls. The other ten pairs still rely on the reconstructed masks. Each
  mask used both SMT siblings on 4, 8, or 16 physical cores for the 8-, 16-,
  or 32-thread groups, respectively. Both versions used `--aligner dag-aligner`
  and `--region-extension-length 1000`. The
  v1.1.0 runs were sequential on each worker under a 172 GiB memory cap.
- GNU time measured elapsed wall time for the ExpansionHunter invocation and
  maximum resident set size (peak RSS), converted from KiB to GiB. The runs
  used prepared BAMs and their indexes. The benchmark covers the streaming
  invocation, including catalog setup and output writing. Index construction
  has no separate timing in this comparison.
