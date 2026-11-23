# Usage

NuExpansionHunter requires a BAM or a CRAM file containing aligned reads
from a PCR-free WGS sample, a FASTA file with a reference genome assembly (which
must be the same as the one against which the reads were aligned), and a [variant catalog
file](04_VariantCatalogFiles.md).

NuExpansionHunter outputs a VCF file and a JSON file with variant genotypes and
other useful information along with a BAMlet containing alignments of reads that
overlap or located in close proximity to each variant. The VCF and JSON files
are largely equivalent, but the JSON file may be easier to parse
programmatically. Here is a template with the names of the required parameters.

```bash
ExpansionHunter --reads <BAM/CRAM file with aligned reads> \
                --reference <FASTA file with reference genome> \
                --variant-catalog <JSON file specifying variants to genotype> \
                --output-prefix <Prefix for the output files>
```

## Optional arguments

In addition to the required program options listed above, there are a number of
optional arguments.

* `--sex <arg>` Specifies sex of the sample; can be either `male` or `female`
  (default). This parameter only affects repeats on sex chromosomes.

* `--region-extension-length <int>` Specifies how far from on/off-target regions
   to search for informative reads. Set to 1000 by default.

* `--analysis-mode <mode>` Specify analysis mode, which can be either `seeking` or
  `streaming`. The default mode is `seeking`. See further description of analysis
   modes below.

* `--threads <int>` (`-n`) Specifies the number of compute workers. Set to 1 by
  default; must be at least 1. Codec pools and coordinator threads are additional,
  so this is not a total process thread limit.

* `--decode-threads auto|N` Specifies streaming input decompression workers.
  The default value is `auto`: `min(threads, max(4, threads / 4))` in streaming
  mode. `0` or `1` uses synchronous decompression. Seeking mode reads its input
  without decompression workers, and an explicit numeric value there is an error.

* `--encode-threads auto|N` Specifies BAM output compression workers.
  The default value is `auto`: `min(threads, max(2, threads / 4))` in streaming
  mode and `threads` in seeking mode. `0` uses synchronous compression. This
  setting is unused when BAMlet output is disabled.

* `--bamlet-queue-mb auto|N` Specifies the queued BAM-record budget in MiB.
  The default value is `auto`: 8 MiB, or 16 MiB with at least 16 compute
  workers, in both modes. Explicit `0` means unbounded.
  This is not a process RAM limit.

* `--catalog-cache PATH` Reads a prepared catalog from a read-only memory-mapped
  archive. No cache is selected by default. See catalog cache usage below.

* `--prepare-catalog-cache` Builds the archive named by `--catalog-cache` and
  exits before processing reads. Disabled by default.

* `--aligner <arg>` (`-a`) Specifies `dag-aligner` (default) or `path-aligner`.

* `--log-level <arg>` Specifies `trace`, `debug`, `info` (default), `warn`, or `error`.

* `--help` (`-h`) Prints the help message.

* `--version` (`-v`) Prints `Starting NuExpansionHunter v1.1.0 (based on ExpansionHunter v5.0.0)`
  to stderr and exits.

`auto` is a default value, not an option. In the formulas, `threads` is the
value of `--threads` and division is integer division. Explicit codec thread
counts and queue budgets must be non-negative integers. The short form of
`--analysis-mode` is `-m`.

The hidden internal option `--disable-bamlet-output` disables BAMlet output;
it is not listed by `--help` and is not a stable public interface.

Note that the full list of program options with brief explanations can be
obtained by running `ExpansionHunter --help`.

### Analysis modes

#### Seeking mode

Seeking mode retrieves reads by locus using the input alignment index. Use it
for a small number of loci. Streaming mode is the focus of the parallel input,
ordered worker pipeline and extension-trace cache.


#### Streaming mode

In streaming mode, the alignment file is read in a single pass and all variants are
analyzed during this reading operation. Streaming mode is recommended for the analysis
of large catalogs, but does require more memory as a function of catalog size. This mode
does not require that the BAM or CRAM file is sorted or indexed.

Indexed streaming requires the parallel reader and a suitable regular BAM with
a usable index. Streaming can fall back to sequential scanning. Realigned BAM
records in streaming mode are not in
ExpansionHunter v5.0.0's order: they have the same multiset of decoded records,
with scheduling-dependent order and compression level 1. See
[Performance and identity](07_Performance.md) for the complete limits.

### Catalog cache

Prepare a reusable catalog with the usual required inputs and an output prefix:

```bash
ExpansionHunter --reads "$READS" --reference "$REFERENCE" \
  --variant-catalog "$CATALOG" --output-prefix "$OUTPUT" \
  --analysis-mode streaming --prepare-catalog-cache --catalog-cache "$CACHE"
```

Preparation writes the cache and exits before processing reads. Its parent
directory must exist; preparation does not overwrite an existing file. To use
it, omit `--prepare-catalog-cache` and retain the usual input arguments:

```bash
ExpansionHunter --reads "$READS" --reference "$REFERENCE" \
  --variant-catalog "$CATALOG" --output-prefix "$OUTPUT" \
  --analysis-mode streaming --catalog-cache "$CACHE"
```

Loading validates the archive format and a signature over the program build
and inputs: the source-derived schema id, catalog and reference contents,
contig dictionary, region-extension length, and orientation and alignment k-mer
settings. The read data may differ across samples with a matching contig
dictionary. A cache belongs to the program's source: the schema id hashes every
source file, including comments, so a program built from different sources
rejects it. Prepare caches with the executable used for analysis, once per
catalog and program build, and reuse them for every sample. A cache mismatch is
an error; there is no fallback to parsing. Without `--catalog-cache`, the
catalog is built in memory from the input files.

## Running efficiently

* Use streaming mode (`--analysis-mode streaming`) for catalogs with more than a
  few loci. Seeking mode is meant for a handful of loci.
* Set `--threads` to the number of hardware threads. The automatic decode,
  encode and queue settings adapt to the worker count. The [benchmark](../BENCHMARK.md)
  uses these defaults with `--threads 8`, `--threads 16` and `--threads 32` on
  the same number of logical CPUs.
* For several samples with the same catalog, prepare a [catalog cache](#catalog-cache)
  once and pass it to every run. The run then skips catalog parsing and graph
  and index construction. The cache is memory-mapped read-only, so concurrent
  runs on one machine share it.
* Build with the optional libraries (libdeflate, ISA-L and jemalloc; see
  [Installation](02_Installation.md)). The benchmark builds include the optional libraries.
* Memory grows with `--threads` on large catalogs, because each worker keeps
  recently used locus indexes. Streaming also uses at most 64 MiB for the
  extension-trace table per pipeline. On a machine with little memory, lower
  `--threads` or split the catalog into parts and run them separately.
* Keep the BAM, its index and the reference on local disk.
