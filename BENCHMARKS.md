# Benchmark

[← NuExpansionHunter README](README.md#benchmarks)

## Happy Friday vs. dp-readprocessing on full HG002

We compared Happy Friday (`0888426`) with dp-readprocessing (`66200f6`) on the full HG002 BAM with the 174,293-locus Illumina catalog. Each core count used one pair of sequential runs on the same Oscar Genoa allocation with the same CPU affinity. Happy Friday ran second at 4 and 16 cores and first at 8 and 32 cores. Each binary used its own prepared catalog cache. Wall time includes cache loading and analyzer initialization and excludes cache preparation and output verification.

| Physical cores | Slurm job | dp-readprocessing wall (s) | Happy Friday wall (s) | Wall time saved | dp-readprocessing peak RSS (GiB) | Happy Friday peak RSS (GiB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 9052318 | 1354.772 | 1230.193 | 9.2% | 10.600 | 10.587 |
| 8 | 9052319 | 759.154 | 695.411 | 8.4% | 10.600 | 10.597 |
| 16 | 9052320 | 453.793 | 427.191 | 5.9% | 10.623 | 10.618 |
| 32 | 9052321 | 494.137 | 474.047 | 4.1% | 10.603 | 10.617 |

### Read processing, genotyping, and output

The sum covers these three logged phases. It excludes cache loading, analyzer initialization, and any time outside the phase timers.

| Physical cores | Version | Read processing (s) | Genotyping (s) | Output (s) | Sum (s) |
| ---: | :--- | ---: | ---: | ---: | ---: |
| 4 | dp-readprocessing | 1286.544 | 4.275 | 6.637 | 1297.456 |
| 4 | Happy Friday | 1156.150 | 4.259 | 6.672 | 1167.080 |
| 8 | dp-readprocessing | 698.027 | 3.545 | 6.780 | 708.351 |
| 8 | Happy Friday | 622.877 | 3.173 | 6.697 | 632.747 |
| 16 | dp-readprocessing | 403.006 | 0.942 | 6.764 | 410.712 |
| 16 | Happy Friday | 374.919 | 1.092 | 6.676 | 382.687 |
| 32 | dp-readprocessing | 442.890 | 0.500 | 7.570 | 450.961 |
| 32 | Happy Friday | 422.021 | 1.011 | 7.170 | 430.201 |

Happy Friday passed full output comparison against retained Illumina ExpansionHunter v5.0.0 at all four core counts. dp-readprocessing passed at 4, 8, and 16 cores. Its 32-core comparison was in progress at the last confirmed status. A passing comparison means byte-identical VCF and JSON and matching BAM headers and records after ignoring record order. Each result is a single paired observation.

## Multithread branch: complete HG002

The Multithread candidate at commit `91a359c49a9f50afced3bfea102604f2aef480d3` processed the complete HG002 BAM with the 174,293-locus Illumina catalog. Each core count ran once on Oscar; original EH timings are historical observations from separate runs. Per-sample wall time excludes one-time cache preparation and output comparison. Some candidate runs overlapped on a node, so the measurements are not isolated paired comparisons.

| Cores | Wall time, historical original → Multithread | Speedup | Allocated CPU-hours saved | Peak RSS, original → Multithread | Full output comparison |
| ---: | ---: | ---: | ---: | ---: | :--- |
| 4 | 3295.602 → 1706.411 s | **1.931×** | 48.22% | 15.338 → 10.596 GiB | **Output identical** |
| 8 | 1703.735 → 902.082 s | **1.889×** | 47.05% | 15.681 → 10.595 GiB | **Output identical** |
| 16 | 1248.945 → 576.285 s | **2.167×** | 53.86% | 29.584 → 10.616 GiB | **Output identical** |
| 32 | 1064.277 → 500.467 s | **2.127×** | 52.98% | 35.956 → 10.615 GiB | **Output identical** |

VCF and JSON were byte-identical to the retained original output; BAM headers and all 74,565,578 records, fields, tags and duplicates matched after ignoring record order. The correctness reference is distinct from the historical timing measurements. The candidate binary SHA-256 was `d730992506721e16d2a0eef1e451e0ea056670bbffd0af5b581564642f5a9f92`. The test build used explicit map initializers in two test files for GCC compatibility; product code was unchanged.

These measurements use default `--read-decode-threads 0`, so they do not measure the optional parallel read-preprocessing pool. Automatic decode and encode settings and the default BAM queue were used. Each row is a single observation; it does not establish run-to-run variation.

## Historical Jiapeng optimization results

The merged `Jiapeng_optimization` implementation was benchmarked on the complete HG002 BAM with the 174,293-locus Illumina catalog. Each core count ran once on Oscar. Original EH timings are historical measurements; the original was not rerun. The table excludes one-time catalog preparation and output-comparison time.

| Cores | Wall time, historical original → Jiapeng | Speedup | Allocated CPU-hours saved | Peak RSS, original → Jiapeng | Latest full output comparison |
| ---: | ---: | ---: | ---: | ---: | :--- |
| 4 | 3295.602 → 2168.563 s | **1.520×** | **34.20%** | 15.338 → 10.588 GiB | **Output identical** |
| 8 | 1703.735 → 1167.031 s | **1.460×** | **31.50%** | 15.681 → 10.593 GiB | **Output identical** |
| 16 | 1248.945 → 659.187 s | **1.895×** | **47.22%** | 29.584 → 10.613 GiB | **Output identical** |
| 32 | 1064.277 → 338.709 s | **3.142×** | **68.17%** | 35.956 → 10.614 GiB | **Output identical** |

The complete HG002 output comparison passed at every core count: VCF and JSON were byte-identical, and BAM headers plus all 74,565,578 records, fields, tags, and duplicate counts matched after ignoring record order only. The Jiapeng binary used a fresh compact-index cache; preparation took 443.993 seconds and produced a 9,191,803,561-byte archive. The benchmark used automatic decode/encode settings and the default unbounded BAM output queue. The main read, pair, query, and dispatch loop remains serial.

The comparison is against historical original EH and MMAP v2 timings, not a contemporaneous paired rerun. Relative to the accepted MMAP v2 observations, Jiapeng was 3.3%, 13.4%, 6.8%, and 44.5% faster at 4, 8, 16, and 32 cores respectively. These are single-run observations.

## MMAP v2 whole-catalog results

The latest binary, `catalog-mmap-v2` (SHA-256 `d6470f61b65421fd16387052d1c9df7864c2fd48ccf7c1b7f366bf10b42d5e07`), processed the complete HG002 BAM and the same 174,293-locus Illumina catalog described below. The timing table reports one run per core count on Oscar Genoa nodes. The original EH measurements were reused from the earlier benchmark; the original was not rerun. These are comparisons with a historical baseline, not new same-node pairs. The operating-system page cache was not cleared.

Per-sample wall time includes cache validation and mapping, analyzer/index initialization, read processing, genotyping, and output. One-time cache preparation and output comparison are excluded and reported separately. Allocated CPU-hours are physical cores × wall hours; process CPU-hours are user + system CPU seconds divided by 3600. RSS is reported in GiB (2³⁰ bytes).

| Cores | Wall time, historical original → MMAP | Speedup | Allocated CPU-hours, original → MMAP | CPU-hour saving | Process CPU-hours, original → MMAP | Peak RSS, original → MMAP | Latest full output comparison |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| 4 | 3295.602 → 2240.789 s | **1.471×** | 3.662 → 2.490 | **32.01%** | 3.390 → 2.300 | 15.338 → 15.196 GiB (−0.93%) | **Output identical*** |
| 8 | 1703.735 → 1348.256 s | **1.264×** | 3.786 → 2.996 | **20.86%** | 3.304 → 2.487 | 15.681 → 15.194 GiB (−3.10%) | **Output identical*** |
| 16 | 1248.945 → 703.918 s | **1.774×** | 5.551 → 3.129 | **43.64%** | 3.865 → 2.354 | 29.584 → 15.207 GiB (−48.60%) | **Output identical*** |
| 32 | 1064.277 → 489.479 s | **2.174×** | 9.460 → 4.351 | **54.01%** | 4.126 → 2.562 | 35.956 → 15.215 GiB (−57.68%) | **Output identical*** |

The candidate ran on `node2353`, `node2354`, `node2426`, and `node2406` for 4, 8, 16, and 32 cores respectively. Compared with the previous optimized version, the new 4-, 16-, and 32-core observations used less time, while the 8-core observation took **7.7% longer** (1348.256 vs 1251.675 seconds). This single-run comparison does not establish a uniform speed improvement from MMAP or isolate the cause of the 8-core difference.

### Per-sample phases

| Cores | Cache validation and mapping | Analyzer/index initialization | Read processing | Genotyping | Output |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 33.926 s | 29.469 s | 2158.812 s | 4.385 s | 10.630 s |
| 8 | 36.063 s | 57.441 s | 1235.432 s | 3.059 s | 11.871 s |
| 16 | 35.369 s | 48.325 s | 606.683 s | 1.227 s | 9.728 s |
| 32 | 35.047 s | 60.048 s | 381.358 s | 0.633 s | 9.691 s |

The first two columns correspond to the `catalog_graph_seconds` and `catalog_index_seconds` log fields. In cached runs, the first field measures loading and validating the prepared catalog rather than rebuilding graphs. These phase timers do not cover every setup and teardown operation, so their sum is slightly below total wall time.

### Cache preparation

Preparing this catalog took **490.755 seconds**, with **2.617 GiB peak RSS**, and produced a **13,402,416,541-byte archive** (12.482 GiB). Preparation used one application thread inside an eight-core allocation: **1.091 allocated CPU-hours** and **0.135 process CPU-hours**. This cost is paid once for a compatible catalog, reference, contig dictionary, index configuration, and builder source version; it is not included in the per-sample table. A first sample that also requires preparation pays both costs.

The archive maps graphs, strings, locus and variant metadata, 10-mer membership keys, 14-mer path indexes, and interval indexes. Sample evidence and alignment workspaces remain separate. MMAP does not impose a fixed RSS limit: file pages and mutable sample state both contribute to resident memory.

### Latest output validation

**Output identical on complete HG002 at 4, 8, 16, and 32 cores:** VCF/JSON were byte-identical, and BAM headers and all 74,565,578 records matched, including every field, tag, and duplicate count, ignoring record order only.

\* Complete outputs from all four core counts were compared locally with retained original-EH output from the earlier 32-core run `8574491`. That reference used the same complete input BAM, reference, catalog, sex, aligner, and region extension, and the original binary SHA-256 was `9605813238e33b4f359ce2498bd8e6524d6f704c77742e0a1abe420b70294da8`. The reference used for correctness is separate from the historical timing measurement used in the speedup table.

- VCF and JSON were byte-identical.
- BAM headers matched exactly.
- Both BAMs contained **74,565,578 records**. Every SAM field, tag value, and duplicate multiplicity matched after lossless sorting; only record order was ignored.

The 32-core comparison used the run reported in the timing table. Separate correctness runs at 4, 8, and 16 cores then passed the same full comparison against the retained original output. Their timings do not replace the measurements in the performance table. Each verified run contained 74,565,578 BAM records. Raw outputs were removed only after comparison passed; logs, measurements, comparison reports, and cleanup manifests remain on Oscar.

Before these runs, the latest source passed the EH unit suite and all 30 bundled GraphTools tests, including mapped-handle lifetime, archive round trips and corruption rejection, compact-index differential checks, and coordinate queries with remote mates and off-target regions. The example and a 10,000-read-pair HG002 subset with the complete catalog also passed exact VCF/JSON and BAM-content comparisons.

## Historical paired benchmark: stream-team3-bgzf-v1

The following results describe the previous optimized binary. Its full outputs were checked at every core count. They are retained here separately from the latest MMAP measurements.

We compare **Illumina ExpansionHunter v5.0.0** with **NuExpansionHunter** on the complete **HG002 human paired-end WGS** sample. The tested NuExpansionHunter binary is candidate `stream-team3-bgzf-v1` (SHA-256 `799bd667a8b8da55c40d41efdbdfafc65b03bb4a3db9e99b102df7dcb2e13da9`).

## Data and setup

| Setting | Value |
| :--- | :--- |
| Sample | GIAB HG002 / NA24385, PrecisionFDA Truth V2 |
| Sequencing | Illumina NovaSeq, PCR-free, nominal 35× WGS |
| Primary reads | **830,172,418**, comprising **415,086,209 pairs**; no downsampling |
| Input | Complete coordinate-sorted BAM, 46,033,519,192 bytes |
| Reference | `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` |
| Catalog | Illumina RepeatCatalogs hg38, revision `66e12774812319a58e19151e7ae054f5ecc88b93`, **174,293 loci** |
| Analysis | Streaming; male; DAG aligner; 1 kb region extension |
| Threads | 4, 8, 16, and 32 physical cores; no SMT |
| Platform | Brown Oscar `batch` partition; AMD EPYC 9684X (Genoa) nodes |
| Build | GCC 11.5.0, Release |
| Replication | One paired baseline/candidate run per core count, on the same node |

Streaming mode reads the complete BAM. The input also contains 5,088,387 supplementary alignment records, which are excluded from the primary-read count. No records are removed from the input.

The 4- and 8-core jobs were serialized to limit temporary-disk use. The operating-system page cache was not cleared. Each row is therefore one paired observation, not an estimate of run-to-run or node-to-node variance.

## Candidate

The candidate changes:

- batched mate retrieval in seeking mode;
- compact affine-alignment state with exact insertion traceback bits;
- rolling two-bit orientation and seed k-mers;
- encoded k-mer lookup with non-owning path views;
- ownership transfer for cached streaming mates;
- streaming locus queues allocated only while active;
- parallel BGZF compression with explicit writer completion and propagated errors; and
- initialized packed BAM sequence bytes.

The full-WGS benchmark uses streaming mode, so it does not measure the seeking-only batched mate path. The build pins HTSlib 1.10.2 and applies setup-rollback and asynchronous thread-pool/BGZF lifecycle repairs.

## Results

EH wall time covers ExpansionHunter execution only. Output comparison is separate. Allocated CPU-hours are physical cores × wall hours; process CPU-hours are GNU `time` user + system CPU seconds divided by 3600.

| Cores | EH wall time, original → candidate | Speedup | Allocated CPU-hours, original → candidate | CPU-hour saving | Process CPU-hours, original → candidate | Peak RSS, original → candidate | Output identical* |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| 4 | 3295.602 → 2313.767 s | **1.424×** | 3.662 → 2.571 | **29.79%** | 3.390 → 2.347 | 15.338 → 15.664 GiB (+2.12%) | Yes |
| 8 | 1703.735 → 1251.675 s | **1.361×** | 3.786 → 2.781 | **26.53%** | 3.304 → 2.366 | 15.681 → 15.674 GiB (−0.04%) | Yes |
| 16 | 1248.945 → 765.532 s | **1.631×** | 5.551 → 3.402 | **38.71%** | 3.865 → 2.356 | 29.584 → 15.683 GiB (−46.99%) | Yes |
| 32 | 1064.277 → 596.874 s | **1.783×** | 9.460 → 5.306 | **43.92%** | 4.126 → 2.499 | 35.956 → 15.701 GiB (−56.33%) | Yes |

\* VCF and JSON were byte-identical. BAM headers and records were identical after ignoring record order, as described below.

The speedup is not monotonic between 4 and 8 cores. In this single pair, the original scaled by 1.934× from 4 to 8 cores, while NuExpansionHunter scaled by 1.849×. The 16- and 32-core runs also show much lower peak memory with compact alignment state and demand-allocated locus queues. Repeated runs are needed before attributing the 4-to-8-core difference to a specific cause.

For this sample and catalog, all four runs used less wall time and fewer allocated CPU-hours. These single runs do not provide a population mean, error bars, or performance estimates for other machines, input technologies, catalogs, or samples.

## Output equivalence

At every core count:

- VCF and JSON matched byte for byte.
- BAM headers matched exactly.
- Both BAMs contained **74,565,578 records**.
- After lossless external sorting, every decompressed SAM field and tag value matched, and duplicate multiplicities were retained.
- Only BAM record order was ignored; raw BAM bytes were therefore not expected to match.

The comparator was tested with deliberately changed calls, tags, mapping quality, and duplicate counts. Comparison time is excluded from the EH timings. Raw VCF, JSON, and BAM outputs were removed only after the corresponding full comparison succeeded; measurement, comparison, status, and cleanup metadata were retained. The accepted 16- and 32-core metadata remain under the benchmark archive; the 4- and 8-core cleanup manifests record successful removal after comparison.

<details>
<summary><strong>Command</strong></summary>

Each binary uses the same inputs and options, with its core count substituted for `N` and a separate output prefix:

```sh
ExpansionHunter \
  --reads HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam \
  --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --variant-catalog genomewide.hg38.json \
  --output-prefix results/sample \
  --analysis-mode streaming \
  --threads N \
  --sex male \
  --aligner dag-aligner \
  --region-extension-length 1000 \
  --log-level info
```

</details>

## Validation before the full-WGS benchmark

The combined source passed the ExpansionHunter unit suite and all 30 bundled GraphTools tests. Targeted differential tests covered affine-alignment scores and traceback, rolling k-mer lookup, queue behavior, and cached-read ownership. Writer tests covered complete BAM content and failures during setup, record writing, and close. The official example and a paired human-read fixture matched at 2 and 8 threads before the full-WGS runs.

In this historical matrix, each core count passed its full-WGS output comparison. The latest MMAP validation scope is reported separately above.

## Build notes

Both versions use GCC 11.5.0, Release mode, and the same dependency versions. The following source-build fixes are included:

- download Boost 1.73.0 from `archives.boost.io`;
- keep the GoogleTest 1.10 `maybe-uninitialized` warning from stopping the GCC build;
- include `<limits>` explicitly in `GraphCoordinates.cpp`; and
- link the GCC runtime dynamically because Oscar does not provide its static development archive.

NuExpansionHunter additionally applies the two pinned HTSlib 1.10.2 lifecycle patches under [cmake/patches/](cmake/patches/) before compiling HTSlib. They do not change the HTSlib version.

## Sources

- [Google DeepVariant full-WGS inputs](https://github.com/google/deepvariant/blob/r1.6/scripts/inference_deeptrio.sh)
- [HG002 BAM](https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam?generation=1594940262604684) and [BAM index](https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam.bai?generation=1594940262931273)
- [GRCh38 no-alt reference](https://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz?generation=1594960641005742)
- [Illumina genome-wide hg38 STR catalog](https://github.com/Illumina/RepeatCatalogs/tree/66e12774812319a58e19151e7ae054f5ecc88b93/hg38)
