# Benchmark

[← NuExpansionHunter README](README.md#benchmarks)

We compare **Illumina ExpansionHunter v5.0.0** with **NuExpansionHunter** on the complete **HG002 human paired-end WGS** sample. NuExpansionHunter includes batched mate retrieval, a reusable DP row, and initialized BAM padding. This benchmark tests streaming mode; batched mate retrieval is used only in seeking mode.

## Data and setup

| Setting | Value |
| :--- | :--- |
| Sample | GIAB HG002 / NA24385, PrecisionFDA Truth V2 |
| Sequencing | Illumina NovaSeq, PCR-free, nominal 35× WGS |
| Primary reads | **830,172,418**, comprising **415,086,209 pairs**; no downsampling |
| Input | Complete coordinate-sorted BAM, 46,033,519,192 bytes |
| Reference | `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` |
| Catalog | Illumina RepeatCatalogs, 174,293 hg38 loci |
| Analysis | Streaming; 8, 16, 32, and 64 threads; male; DAG aligner; 1 kb region extension |
| Platform | Brown Oscar, `batch` partition; AMD Genoa nodes; one physical core per requested thread and 64 GiB allocated |
| Build | GCC 11.5.0, Release; matching compiler settings and libraries |
| Repeats | Three paired runs per thread count, each on one node: original/optimized, optimized/original, original/optimized |

Streaming mode reads the complete BAM. The input also contains 5,088,387 supplementary alignment records, which are excluded from the primary-read count. No records are removed from the input.

## Results

**Pending.** All four settings are submitted and await output-verification and storage dependencies. Each uses the complete HG002 input. The results will include the CPU model, node, and run settings.

| Threads | Median runtime: original → NuExpansionHunter | Speedup | Peak RAM: original → NuExpansionHunter | Output check |
| :--- | ---: | ---: | ---: | :--- |
| 8 | Pending | Pending | Pending | Pending |
| 16 | Pending | Pending | Pending | Pending |
| 32 | Pending | Pending | Pending | Pending |
| 64 | Pending | Pending | Pending | Pending |

We will report all three runtimes and the median and maximum peak RAM for each setting. Each speedup compares the two programs at the same thread count on the same node. All settings use the same CPU family; remaining hardware differences will be listed.

Small RAM increases are acceptable. Changes will be reported in MiB and percent. Runtime and peak RSS come from GNU `time`; output comparison runs separately and is excluded from runtime. The operating-system page cache is not cleared between runs.

<details>
<summary><strong>Command and output checks</strong></summary>

Both binaries use the same inputs and options, with separate output prefixes. This is the eight-thread command; the other rows set `--threads` to `16`, `32`, or `64`:

```sh
ExpansionHunter \
  --reads HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam \
  --reference GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --variant-catalog genomewide.hg38.json \
  --output-prefix results/sample \
  --analysis-mode streaming \
  --threads 8 \
  --sex male \
  --aligner dag-aligner \
  --region-extension-length 1000 \
  --log-level info
```

VCF and JSON must match byte for byte. BAM headers and every record must match in all fields, tag values, and duplicate counts; only record order is ignored. Raw BAM byte equality is reported separately.

The comparator is checked using deliberately changed calls, tags, mapping quality, and duplicate counts. All failed outputs are retained. After successful comparison, the first original/optimized outputs and all commands, measurements, and comparison records are kept.

</details>

## Build notes

Both versions use GCC 11.5.0, Release mode, and the same dependency versions. The following build fixes are applied to both:

- Download Boost 1.73.0 from `archives.boost.io`.
- Keep the GoogleTest 1.10 `maybe-uninitialized` warning from stopping the GCC build.
- Include `<limits>` explicitly in `GraphCoordinates.cpp`.
- Link the GCC runtime dynamically; Oscar does not provide its static development archive.

## Sources

- [Google DeepVariant full-WGS inputs](https://github.com/google/deepvariant/blob/r1.6/scripts/inference_deeptrio.sh)
- [HG002 BAM](https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam?generation=1594940262604684) and [BAM index](https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam.bai?generation=1594940262931273).
- [GRCh38 no-alt reference](https://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz?generation=1594960641005742).
- [Illumina genome-wide hg38 STR catalog](https://github.com/Illumina/RepeatCatalogs/tree/66e12774812319a58e19151e7ae054f5ecc88b93/hg38)
