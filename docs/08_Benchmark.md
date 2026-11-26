# NuExpansionHunter v1.1.0 versus ExpansionHunter v5.0.0

The benchmark covers 64 completed GangSTR whole-genome streaming pairs on 23
genomes. At 4–16 threads, NuExpansionHunter v1.1.0 was 9.77–11.02× faster by
geometric mean and used 86.8–88.3% less peak resident memory (RSS) than stock
ExpansionHunter v5.0.0. All 64 pairs passed the four output identity checks.

## Summary

Speedup and CPU-time ratio are geometric means of the per-genome v5 divided
by v1.1.0 ratios. The 95% confidence intervals use Student’s *t* on the log
speedups, with the genome as the replicate. Peak RSS reduction is one minus
the reciprocal of the geometric mean of paired peak RSS ratios.

| Threads | Pairs | Geometric mean speedup (95% CI) | Speedup range | Median peak RSS, v5 → v1.1.0 (GiB) | Peak RSS reduction | CPU-time ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 8 | 9.77× (9.60–9.94) | 9.34–10.02× | 72.0 → 8.1 | 88.3% | 9.68× |
| 8 | 21 | 10.18× (9.99–10.37) | 9.53–11.17× | 72.2 → 9.0 | 87.2% | 9.99× |
| 16 | 23 | 11.02× (10.86–11.19) | 10.45–11.69× | 82.3 → 10.8 | 86.8% | 10.38× |
| 32 | 12 | 19.02× (18.74–19.30) | 18.21–19.62× | 156.7 → 14.3 | 90.8% | 13.75× |

At 32 threads, stock v5.0.0 finished 12 of 21 genomes (11 controls and one
expansion-study sample) and hit the 172 GiB memory cap on the other nine.
NuExpansionHunter completed all 21. The 19.02× speedup covers the 12 completed
pairs. On five genomes, the speedup against the [jemalloc and libdeflate build](#stronger-v500-baseline)
was 9.7–10.5× at 32 threads.

The sample mix varied by thread count. On the five genomes with completed
pairs at every headline thread count, geometric mean speedups against stock
v5.0.0 were 9.71×, 9.81×, 10.59×, and 18.89× at 4, 8, 16, and 32 threads,
respectively.

In the one-thread example (ERR1955399, 39.2× coverage), stock v5.0.0 took
63,354 s (17.6 h) and 73.5 GiB peak RSS. NuExpansionHunter took 6,712 s
(1.9 h) and 7.1 GiB, a 9.44× speedup and a 90.4% peak RSS reduction. The
two runs used one logical CPU each on separate cores in different cache
domains of the same node.

## Per-sample results

Each row is one completed pair. Coverage was measured from the aligned BAMs
and is shown to one decimal place.

| Sample | Coverage | Threads | v5 wall (s) | v1.1.0 wall (s) | Speedup | v5 peak RSS (GiB) | v1.1.0 peak RSS (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ERR1955394 | 41.8× | 4 | 20,978 | 2,150.86 | 9.75× | 71.77 | 7.86 |
| ERR1955394 | 41.8× | 8 | 10,618 | 1,082.64 | 9.81× | 72.02 | 8.78 |
| ERR1955394 | 41.8× | 16 | 5,807 | 549.78 | 10.56× | 84.05 | 10.61 |
| ERR1955394 | 41.8× | 32 | 5,121 | 278.21 | 18.41× | 155.39 | 14.31 |
| ERR1955395 | 42.3× | 8 | 13,238 | 1,272.89 | 10.40× | 72.21 | 8.79 |
| ERR1955395 | 42.3× | 16 | 7,542 | 647.08 | 11.66× | 82.20 | 10.64 |
| ERR1955395 | 42.3× | 32 | 6,427 | 327.54 | 19.62× | 164.21 | 14.34 |
| ERR1955396 | 40.0× | 8 | 12,306 | 1,221.69 | 10.07× | 72.09 | 8.81 |
| ERR1955396 | 40.0× | 16 | 6,554 | 619.66 | 10.58× | 72.51 | 10.61 |
| ERR1955396 | 40.0× | 32 | 5,675 | 311.64 | 18.21× | 148.10 | 14.25 |
| ERR1955398 | 42.2× | 8 | 14,453 | 1,436.90 | 10.06× | 72.10 | 8.97 |
| ERR1955398 | 42.2× | 16 | 7,732 | 728.93 | 10.61× | 74.74 | 10.81 |
| ERR1955398 | 42.2× | 32 | 6,920 | 365.54 | 18.93× | 163.64 | 14.47 |
| ERR1955399 | 39.2× | 4 | 20,468 | 2,191.62 | 9.34× | 71.46 | 7.93 |
| ERR1955399 | 39.2× | 8 | 10,426 | 1,093.70 | 9.53× | 71.54 | 8.84 |
| ERR1955399 | 39.2× | 16 | 5,888 | 553.04 | 10.65× | 81.78 | 10.64 |
| ERR1955399 | 39.2× | 32 | 5,314 | 275.74 | 19.27× | 142.45 | 14.33 |
| ERR1955400 | 39.4× | 4 | 21,120 | 2,176.41 | 9.70× | 71.42 | 7.85 |
| ERR1955400 | 39.4× | 8 | 10,704 | 1,092.02 | 9.80× | 71.51 | 8.75 |
| ERR1955400 | 39.4× | 16 | 5,820 | 556.82 | 10.45× | 83.04 | 10.57 |
| ERR1955400 | 39.4× | 32 | 5,201 | 280.17 | 18.56× | 151.91 | 14.14 |
| ERR1955401 | 41.1× | 4 | 23,135 | 2,331.79 | 9.92× | 71.81 | 8.06 |
| ERR1955401 | 41.1× | 8 | 11,717 | 1,173.52 | 9.98× | 71.89 | 8.97 |
| ERR1955401 | 41.1× | 16 | 6,448 | 594.27 | 10.85× | 76.85 | 10.77 |
| ERR1955401 | 41.1× | 32 | 5,648 | 298.09 | 18.95× | 152.30 | 14.42 |
| ERR1955402 | 40.8× | 8 | 11,493 | 1,149.73 | 10.00× | 72.01 | 8.89 |
| ERR1955402 | 40.8× | 16 | 6,465 | 585.31 | 11.05× | 82.26 | 10.70 |
| ERR1955402 | 40.8× | 32 | 5,609 | 294.12 | 19.07× | 154.96 | 14.36 |
| ERR1955403 | 42.5× | 16 | 7,368 | 675.26 | 10.91× | 79.66 | 10.80 |
| ERR1955405 | 45.8× | 8 | 14,150 | 1,383.04 | 10.23× | 73.07 | 9.29 |
| ERR1955405 | 45.8× | 16 | 7,716 | 704.21 | 10.96× | 82.41 | 11.11 |
| ERR1955405 | 45.8× | 32 | 6,913 | 356.72 | 19.38× | 170.74 | 14.78 |
| ERR1955406 | 38.8× | 8 | 11,626 | 1,107.32 | 10.50× | 71.33 | 8.45 |
| ERR1955406 | 38.8× | 16 | 6,313 | 565.39 | 11.17× | 84.41 | 10.26 |
| ERR1955406 | 38.8× | 32 | 5,603 | 286.84 | 19.53× | 158.74 | 13.92 |
| ERR1955409 | 40.8× | 8 | 11,069 | 1,130.80 | 9.79× | 71.89 | 8.91 |
| ERR1955409 | 40.8× | 16 | 6,142 | 574.94 | 10.68× | 88.79 | 10.71 |
| ERR1955409 | 40.8× | 32 | 5,532 | 289.50 | 19.11× | 157.92 | 14.36 |
| SRR21753323 | 46.4× | 4 | 23,126 | 2,352.27 | 9.83× | 73.90 | 9.27 |
| SRR21753323 | 46.4× | 8 | 11,672 | 1,177.06 | 9.92× | 74.59 | 10.17 |
| SRR21753323 | 46.4× | 16 | 6,273 | 600.03 | 10.45× | 75.58 | 11.91 |
| SRR21753323 | 46.4× | 32 | 5,764 | 299.38 | 19.25× | 163.45 | 15.53 |
| SRR21753324 | 50.1× | 8 | 11,299 | 1,133.96 | 9.96× | 73.54 | 9.66 |
| SRR21753324 | 50.1× | 16 | 6,538 | 579.86 | 11.28× | 107.77 | 11.48 |
| SRR21753325 | 44.9× | 4 | 20,591 | 2,121.90 | 9.70× | 72.49 | 8.64 |
| SRR21753325 | 44.9× | 8 | 10,420 | 1,061.99 | 9.81× | 72.76 | 9.53 |
| SRR21753325 | 44.9× | 16 | 6,217 | 543.73 | 11.43× | 105.85 | 11.35 |
| SRR21753326 | 43.1× | 4 | 21,054 | 2,100.26 | 10.02× | 72.25 | 8.22 |
| SRR21753326 | 43.1× | 8 | 10,658 | 1,055.67 | 10.10× | 72.65 | 9.14 |
| SRR21753326 | 43.1× | 16 | 6,086 | 539.20 | 11.29× | 101.25 | 10.98 |
| SRR21753327 | 76.7× | 8 | 12,378 | 1,167.58 | 10.60× | 80.88 | 12.70 |
| SRR21753327 | 76.7× | 16 | 6,659 | 599.98 | 11.10× | 82.34 | 14.52 |
| SRR21753328 | 63.1× | 4 | 30,111 | 3,045.72 | 9.89× | 76.47 | 10.39 |
| SRR21753328 | 63.1× | 8 | 15,205 | 1,522.84 | 9.98× | 76.89 | 11.31 |
| SRR21753328 | 63.1× | 16 | 8,510 | 773.16 | 11.01× | 110.27 | 13.12 |
| SRR7205168 | 40.9× | 8 | 15,148 | 1,394.09 | 10.87× | 72.57 | 8.94 |
| SRR7205168 | 40.9× | 16 | 8,266 | 707.19 | 11.69× | 84.75 | 10.76 |
| SRR7205171 | 42.5× | 8 | 21,013 | 1,935.75 | 10.86× | 72.63 | 9.20 |
| SRR7205171 | 42.5× | 16 | 10,962 | 973.96 | 11.26× | 76.04 | 11.04 |
| SRR7205174 | 52.0× | 8 | 10,757 | 1,021.64 | 10.53× | 74.11 | 9.89 |
| SRR7205174 | 52.0× | 16 | 5,828 | 523.33 | 11.14× | 86.97 | 11.70 |
| SRR7205175 | 42.4× | 8 | 19,912 | 1,782.92 | 11.17× | 72.09 | 8.97 |
| SRR7205175 | 42.4× | 16 | 10,526 | 903.46 | 11.65× | 78.23 | 10.78 |
| SRR7205176 | 41.4× | 16 | 10,016 | 886.92 | 11.29× | 73.59 | 10.71 |

## Method notes

- **Builds.** Stock ExpansionHunter v5.0.0 was built in Release mode with
  GCC 15.2 for the default x86-64 target, with a header include for compiler
  compatibility. It used bundled htslib 1.10.2, zlib, and the glibc allocator.
  NuExpansionHunter v1.1.0 used the portable x86-64-v2 Release configuration
  with GCC 15.2, run-time SIMD selection, libdeflate, ISA-L, and jemalloc.
  Builds from the `nueh-v1.1.0` release tag reproduced the benchmark outputs
  and ran within 0.10% of the benchmark binary's wall time in checks at 16 and
  32 threads. [Installation guide](02_Installation.md)
- **Inputs.** Both tools used the same prepared BAMs and indexes, GRCh38
  (no-alt), and the GangSTR hg38 v13 whole-genome catalog (832,152 loci).
  Measured coverage was 38.8–45.8× for the 12 Polaris HiSeq X controls
  (ENA PRJEB20654) and 40.9–76.7× for the 11 expansion-study samples
  (SRA PRJNA419676 and PRJNA885420).
- **Allocation.** Eight AMD EPYC 9634 nodes had 188 GB RAM and local NVMe.
  For each headline pair, both tools ran on the same node at equal `--threads`
  on the same recorded CPU mask, using both SMT siblings of 2, 4, 8, or 16
  cores at 4, 8, 16, or 32 threads. `--threads` sets the compute workers.
  NuExpansionHunter’s additional codec workers shared the same CPUs with
  default pool sizes. The nodes enforced a 172 GiB cgroup memory cap, which
  includes file cache, with swap disabled. [Thread options](03_Usage.md#optional-arguments)
- **Timing.** GNU time measured elapsed wall time, user plus system CPU time,
  and maximum resident set size (peak RSS, converted from KiB to GiB).
  Both tools ran with `--analysis-mode streaming --aligner dag-aligner
  --region-extension-length 1000` and the sample’s sex. Measurements cover
  the whole invocation, including catalog setup, input reading and output
  writing. Input alignment and BAM index construction were outside the
  timed invocation.
  Catalog indexes built within the invocation are included in wall time and
  have no separate timing. Each run loaded the catalog without a prepared
  catalog cache.
- **Run order and repeats.** The 4-thread runs alternated which tool ran first.
  The 8–32-thread pairs used archived stock measurements from November 18–20,
  2026, taken 3–6 days before the NuExpansionHunter runs and at lower scheduling
  priority. The CPU masks were recorded for all 56 archived baselines. Wall
  times differed by at most 1.05% over 12 stock repeats and 0.32% over nine
  NuExpansionHunter repeats. Each genome contributes one pair per thread count.
  Runs that hit the memory cap are excluded from speedup and RSS summaries.
  All eight stock runs attempted at 64 threads hit the cap.
- **Page cache.** Inputs were read before timing, though stock v5.0.0 still
  reread some input from NVMe at 16 and 32 threads. In four NuExpansionHunter
  runs starting with the BAM evicted from the page cache, wall time changed
  by less than 0.7% on this storage.

The measurements cover one processor model and local NVMe storage. Controls
and expansion-study samples ran on separate groups of nodes, so cohort and
node effects are confounded. The confidence intervals describe variation
among the measured genomes on this platform.

## Stronger v5.0.0 baseline

At 32 threads, stock v5.0.0 kept 21.3–24.6 of the 32 allocated logical CPUs
busy on completed runs. On the two genomes measured with both v5.0.0 builds
at 32 threads, rebuilding with jemalloc and libdeflate made v5.0.0
1.92–1.99× faster. This accounts for much of the 19.02× speedup against stock.
The rebuild also completed all three 32-thread retests selected from the nine
stock failures, using 66.6–74.8 GiB peak RSS.

Across the tests with the rebuilt baseline, NuExpansionHunter was 8.9–10.5×
faster. The ranges were 8.91–8.96× at 8 threads (two genomes), 9.14–9.22× at
16 threads (two genomes), and 9.7–10.5× at 32 threads (five genomes).

## Output identity and HG002 agreement

All 231 completed output sets in the benchmark identity series matched across 15
groups (199 from NuExpansionHunter runs and 32 from ExpansionHunter builds).
It includes repeated runs, library and SIMD variants, and BAM and CRAM input.
The four checks compare VCF and JSON bytes, the BAM header, and the multiset
of decoded BAM records regardless of order or compression. Illumina’s
unmodified v5.0.0 Linux release matched both the stock source build and
NuExpansionHunter on both genomes tested at 16 threads. It took 17–21% longer
than the stock source build on those genomes.

Against the GIAB HG002 GRCh38 tandem-repeat benchmark v1.0.1, 97.03% of
genotypes agreed exactly (both allele-length errors below 0.5 repeat unit)
and 99.03% agreed within one unit (both errors below 1.5 units) at 624,223
autosomal GangSTR loci fully contained in Tier 1 regions. Alleles were
matched by length, with the comparison window padded by one repeat unit on
each side. The primary analysis excludes 94 no-calls and 3,284 calls with
truth variants crossing the window. Including those 3,284 calls and counting
no-calls as disagreements gives 96.65% exact agreement across all 627,601
eligible loci.

Exact agreement dropped to 86.9% at the 69,100 scored loci where at least one
truth allele differed from the reference by half a repeat unit or more. The
two tools produced identical HG002 calls, so the agreement rates apply to
both. This test measures repeat length on one genome. Sequence, phasing, and
sensitivity for pathogenic expansions were outside its scope.

## Other workloads

| TRExplorer v2.1 catalog | Threads | NuExpansionHunter genomes | Wall time (h) | Peak RSS (GiB) |
| --- | --- | ---: | ---: | ---: |
| Full (5.66 million loci) | 16 and 32 | 4 | 1.8–2.8 at 16, 0.9–1.4 at 32 | 63.2–77.4 |
| 2-million-locus subset | 16 and 32 | 3 | 0.35–0.90 | 29.0–45.3 |

Stock v5.0.0 hit the 172 GiB cap while loading the full catalog on all 12
attempts across five genomes. The rebuild with jemalloc and libdeflate also
failed at load in its one attempt at 16 threads (944 s, 171.6 GiB peak RSS).
Stock v5.0.0 hit the same cap during streaming on all 10 attempts with the
2-million-locus subset. With the full catalog, NuExpansionHunter matched
its own output at 16 and 32 threads on all four genomes. These runs establish
that the catalog can be processed within the test budget. Genotype accuracy
across TRExplorer locus types was outside the scope of the comparison.

On GangSTR CRAM input, NuExpansionHunter was 7.44–7.45× faster than stock v5.0.0
on two genomes at 16 threads. NuExpansionHunter took 1.57–1.62× as long on
CRAM as on the corresponding BAMs. With the 31-locus default catalog in
seeking mode, NuExpansionHunter was 45–198× faster than v5.0.0 at 16 threads
on four genomes, with identical outputs. The reason for the variation among
genomes is unknown. In streaming mode with that catalog, both tools took
about 108 s on the one genome measured.

## Compute cost

At 16 threads the median CPU time was 26.7 CPU-hours per genome for stock
v5.0.0 and 2.6 for NuExpansionHunter. Projected compute cost was about
$1.78 and $0.10 per genome, respectively (17.8×), using Google Cloud C3D
on-demand prices in us-central1 on November 26, 2026. The model multiplies
the measured median wall time by the hourly rate for the least expensive C3D
instance with at least 16 vCPUs and RAM at least 10% above the largest
observed peak RSS. The instances were c3d-highmem-16 at $0.97976/hour for
stock v5.0.0 and c3d-highcpu-16 at $0.599696/hour for NuExpansionHunter. The
estimate covers compute only and assumes the measured wall times carry over
to the cloud instances.
[Google Cloud pricing](https://cloud.google.com/products/compute/pricing/general-purpose?hl=en)
