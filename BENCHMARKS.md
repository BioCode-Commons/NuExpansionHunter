# Benchmark

[← NuExpansionHunter README](README.md#benchmarks)

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

We accept a core count only after the full-WGS output comparison completes. A timing pair alone is not marked as passed.

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
