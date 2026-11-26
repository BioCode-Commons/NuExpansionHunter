# NuExpansionHunter benchmarks

## NuExpansionHunter v1.1.0 versus ExpansionHunter v5.0.0

The benchmark covers 64 completed whole-genome streaming pairs across 23
genomes. At 4–16 threads, NuExpansionHunter v1.1.0 was **9.77–11.02× faster**
by geometric mean and used **86.8–88.3% less peak resident memory (RSS)** than
stock ExpansionHunter v5.0.0. Both programs used the GangSTR hg38 v13 catalog
(832,152 loci), the same `--threads` value, and the same logical CPUs within
each pair on AMD EPYC 9634 nodes with 188 GB RAM.

| Threads | Pairs | Geometric mean speedup (95% CI) | Speedup range | Median peak RSS, v5 → v1.1.0 (GiB) | Peak RSS reduction | CPU-time ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 8 | 9.77× (9.60–9.94) | 9.34–10.02× | 72.0 → 8.1 | 88.3% | 9.68× |
| 8 | 21 | 10.18× (9.99–10.37) | 9.53–11.17× | 72.2 → 9.0 | 87.2% | 9.99× |
| 16 | 23 | 11.02× (10.86–11.19) | 10.45–11.69× | 82.3 → 10.8 | 86.8% | 10.38× |
| 32 | 12 | 19.02× (18.74–19.30) | 18.21–19.62× | 156.7 → 14.3 | 90.8% | 13.75× |

At 32 threads, the stock build kept about 21–25 logical CPUs busy and hit the
172 GiB memory cap on 9 of 21 genomes. The 19.02× speedup covers the 12 that
completed. Rebuilding v5.0.0 with jemalloc and libdeflate nearly doubled its
speed on the two genomes that completed with both builds at 32 threads.
NuExpansionHunter was **9.7–10.5× faster at 32 threads** against the rebuild
on five genomes.

In a one-thread example (ERR1955399), wall time fell from 17.6 to 1.9 hours
(9.44×) and peak RSS from 73.5 to 7.1 GiB. All 64 headline pairs matched the
VCF, JSON, BAM header, and decoded BAM record multiset. Illumina’s official
v5.0.0 Linux binary also matched on both genomes tested.

See the [full benchmark report](docs/08_Benchmark.md) for per-genome results,
methods, HG002 agreement, TRExplorer, CRAM, seeking, and compute costs.

## NuExpansionHunter v1.1.0 versus v1.0.0

Both comparisons ran on November 23, 2026.

### Whole-genome results

The v1.1.0 binary used the portable Release build with libdeflate, ISA-L, and
jemalloc enabled. Each pair ran the GangSTR hg38 v13 whole-genome catalog
in streaming mode at `--threads 16`, on the same 16 logical CPUs (both SMT
siblings of 8 cores) of an AMD EPYC 9634 node under a 172 GiB memory cap.
The eight pairs ran on separate nodes, with v1.1.0 first in every pair.
The public SRA/ENA samples comprised four controls
(40.0–42.3× measured coverage) and four expansion-study samples (42.4–52.0×).
Every pair matched VCF and JSON bytes, BAM headers, and decoded BAM record
multisets. [Installation guide](docs/02_Installation.md)

Speedup is v1.0.0 wall time divided by v1.1.0 wall time. The geometric mean
was 1.367× across the eight samples (range 1.255–1.452×). Peak RSS remained
about 11 GiB.

| Sample | Cohort | Coverage | v1.0.0 wall (s) | v1.1.0 wall (s) | Speedup | v1.0.0 user + sys (s) | v1.1.0 user + sys (s) | v1.0.0 peak RSS (GiB) | v1.1.0 peak RSS (GiB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1955394 | Control | 41.8× | 723.3 | 548.8 | 1.318× | 11417.6 | 8602.8 | 10.54 | 10.60 |
| ERR1955398 | Control | 42.2× | 1012.7 | 726.6 | 1.394× | 16053.7 | 11448.4 | 10.74 | 10.80 |
| ERR1955395 | Control | 42.3× | 884.9 | 645.3 | 1.371× | 14010.3 | 10145.8 | 10.58 | 10.64 |
| ERR1955396 | Control | 40.0× | 836.4 | 619.6 | 1.350× | 13213.1 | 9721.4 | 10.54 | 10.61 |
| SRR7205175 | Expansion study | 42.4× | 1313.5 | 904.4 | 1.452× | 20881.7 | 14310.3 | 10.73 | 10.79 |
| SRR7205171 | Expansion study | 42.5× | 1409.4 | 973.9 | 1.447× | 22407.7 | 15414.3 | 10.97 | 11.04 |
| SRR7205174 | Expansion study | 52.0× | 712.3 | 524.8 | 1.357× | 11259.1 | 8238.3 | 11.62 | 11.69 |
| SRR21753326 | Expansion study | 43.1× | 676.6 | 539.3 | 1.255× | 10653.1 | 8438.2 | 10.89 | 10.98 |

### Chromosome subsets

The chromosome-subset comparison used portable Release builds of
NuExpansionHunter v1.1.0 and v1.0.0. All pairs used HG002 and ran on one AMD
EPYC 9634 machine. Both versions used the same inputs, `--threads` value, and
CPU set within each pair. The 16-thread pairs used both SMT siblings of 8
cores, and the 32-thread pairs used both siblings of 16 cores. The input BAM
was read into the page cache before each streaming run. The run order was
reversed for the second pair in each group. Speedup is the median of the
paired wall-time ratios (v1.0.0 divided by v1.1.0). All ten pairs matched VCF
and JSON bytes, BAM header digests, and decoded BAM record multiset digests.

| Catalog subset | Threads per arm | Pairs | v1.0.0 median wall (s) | v1.1.0 median wall (s) | Median paired speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| GangSTR hg38 v13, chr21–22 | 16 | 3 | 33.02 | 23.11 | 1.43× |
| GangSTR hg38 v13, chr21–22 | 32 | 3 | 17.15 | 12.19 | 1.40× |
| TRExplorer v2.1, chr22 | 16 | 2 | 257.52 | 137.97 | 1.87× |
| TRExplorer v2.1, chr22 | 32 | 2 | 129.43 | 70.72 | 1.83× |
