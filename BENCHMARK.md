# NuExpansionHunter benchmarks

## NuExpansionHunter v1.1.0 versus ExpansionHunter v5.0.0

Across 22 completed whole-genome comparisons, NuExpansionHunter v1.1.0
delivered mean speedups of **10.09–18.79×** and mean peak RSS reductions
of **86.9–90.9%** by thread count.
The runs used AMD EPYC 9634 CPUs, with both SMT siblings of 4, 8, or 16 cores
for 8, 16, or 32 threads, respectively.
All 22 matched the VCF, JSON, BAM header, and decoded BAM record multiset.

| Threads | Samples | v5 wall (min) | v1.1.0 wall (min) | Mean speedup | v5 peak RSS (GiB) | v1.1.0 peak RSS (GiB) | Mean RSS reduction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 9 | 204.8 | 20.3 | 10.09× | 73.9 | 9.7 | 86.9% |
| 16 | 9 | 114.2 | 10.3 | 11.06× | 91.2 | 11.6 | 87.1% |
| 32 | 4 | 100.6 | 5.3 | 18.79× | 157.8 | 14.3 | 90.9% |

Eight additional 32-thread v5 runs hit the 172 GiB memory cap and are excluded
from the averages. Per-sample results and methods are in the
[full benchmark report](docs/08_Benchmark.md).

## NuExpansionHunter v1.1.0 versus v1.0.0

Both comparisons ran on November 23, 2026.

### Whole-genome results

The v1.1.0 binary (SHA-256
`536fef0e287122a554f71cd5a1af6be77e59d25ddb1dcc852c251af3788a4ed8`)
used the portable Release build in the [installation guide](docs/02_Installation.md)
with libdeflate, ISA-L, and jemalloc enabled. Both versions ran the GangSTR
hg38 v13 whole-genome catalog in streaming mode at `--threads 16`, one arm at
a time on the same 16 logical CPUs (8 cores with SMT) of each of the eight AMD
EPYC 9634 machines used for the v5 comparison. The eight public SRA/ENA samples
comprised four 30× controls and four high-coverage patient samples. Every
pair matched VCF and JSON bytes, BAM headers, and decoded BAM record multisets.

Speedup is v1.0.0 wall time divided by v1.1.0 wall time. The geometric mean
was 1.367× across the eight samples (range 1.255–1.452×). Peak RSS remained
about 11 GiB.

| Sample | Cohort | v1.0.0 wall (s) | v1.1.0 wall (s) | Speedup | v1.0.0 user + sys (s) | v1.1.0 user + sys (s) | v1.0.0 peak RSS (GiB) | v1.1.0 peak RSS (GiB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1955394 | 30× control | 723.3 | 548.8 | 1.318× | 11417.6 | 8602.8 | 10.54 | 10.60 |
| ERR1955398 | 30× control | 1012.7 | 726.6 | 1.394× | 16053.7 | 11448.4 | 10.74 | 10.80 |
| ERR1955395 | 30× control | 884.9 | 645.3 | 1.371× | 14010.3 | 10145.8 | 10.58 | 10.64 |
| ERR1955396 | 30× control | 836.4 | 619.6 | 1.350× | 13213.1 | 9721.4 | 10.54 | 10.61 |
| SRR7205175 | High-coverage patient | 1313.5 | 904.4 | 1.452× | 20881.7 | 14310.3 | 10.73 | 10.79 |
| SRR7205171 | High-coverage patient | 1409.4 | 973.9 | 1.447× | 22407.7 | 15414.3 | 10.97 | 11.04 |
| SRR7205174 | High-coverage patient | 712.3 | 524.8 | 1.357× | 11259.1 | 8238.3 | 11.62 | 11.69 |
| SRR21753326 | High-coverage patient | 676.6 | 539.3 | 1.255× | 10653.1 | 8438.2 | 10.89 | 10.98 |

### Chromosome subsets

The v1.1.0 release gate compared the portable Release executable built from
this release's sources (SHA-256
`769a2df6423b56a26ce7127c9f5ce09da92822f8ee6fc4db650c08037da2b814`)
against the frozen portable NuExpansionHunter v1.0.0 binary. Both arms used
HG002, the same inputs and CPU set within each pair, and equal `--threads`
values. All pairs ran on one AMD EPYC 9634 host. The 16-thread pairs used 16
logical CPUs (8 cores, both SMT siblings), and the 32-thread pairs used 32
logical CPUs (16 cores, both SMT siblings). The gate warmed the input BAM
before each streaming run and reversed the arm order for the second pair in
each group. Speedup is the median of the paired v1.0.0 wall time divided by the
v1.1.0 wall time. All ten pairs matched VCF and JSON bytes, BAM header digests,
and decoded BAM record multiset digests.

| Catalog subset | Threads per arm | Pairs | v1.0.0 median wall (s) | v1.1.0 median wall (s) | Median paired speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| GangSTR hg38 v13, chr21–22 | 16 | 3 | 33.02 | 23.11 | 1.43× |
| GangSTR hg38 v13, chr21–22 | 32 | 3 | 17.15 | 12.19 | 1.40× |
| TRExplorer v2.1, chr22 | 16 | 2 | 257.52 | 137.97 | 1.87× |
| TRExplorer v2.1, chr22 | 32 | 2 | 129.43 | 70.72 | 1.83× |
