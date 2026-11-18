# NuExpansionHunter

There are a number of regions in the human genome consisting of repetitions of
short unit sequence (commonly a trimer). Such repeat regions can expand to a
size much larger than the read length and thereby cause a disease.
[Fragile X Syndrome](https://en.wikipedia.org/wiki/Fragile_X_syndrome),
[ALS](https://en.wikipedia.org/wiki/Amyotrophic_lateral_sclerosis), and
[Huntington's Disease](https://en.wikipedia.org/wiki/Huntington%27s_disease)
are well known examples.

NuExpansionHunter aims to estimate sizes of such repeats by performing a targeted
search through a BAM/CRAM file for reads that span, flank, and are fully
contained in each repeat.

NuExpansionHunter is a reimplementation version of Illumina ExpansionHunter v5.0.0. For every
input on which v5.0.0 has defined behaviour, it writes the same VCF and JSON bytes
and a realigned BAM with the same header and the same multiset of decoded records
(record order and compression differ), with the exceptions listed in
[Performance and identity](docs/07_Performance.md). See also the
[verification documentation](ehunter/verification/README.md).

In the current version of NuExpansionHunter, we focus only on the streaming mode of the 
original ExpansionHunter. For the seeking mode, we recommend users to run the original version by Illumina.

## Quick start

On Ubuntu, install the build packages
(`apt install build-essential git cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev`,
plus `libdeflate-dev libisal-dev libjemalloc-dev` for speed), then:

```bash
git clone https://github.com/BioCode-Commons/NuExpansionHunter.git && cd NuExpansionHunter
mkdir build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j"$(nproc)" && cd ..
example/run.sh  # compares a run on the example data with ExpansionHunter v5.0.0's output
build/install/bin/ExpansionHunter --reads sample.bam --reference genome.fa --variant-catalog catalog.json \
  --output-prefix sample --analysis-mode streaming --threads "$(nproc)"
```

For many samples, prepare the catalog cache once per catalog and program build,
then pass it on every run:

```bash
EH=build/install/bin/ExpansionHunter
"$EH" --reads sample1.bam --reference genome.fa --variant-catalog catalog.json \
  --output-prefix prepare --analysis-mode streaming --prepare-catalog-cache --catalog-cache catalog.cache
for bam in sample*.bam; do
  "$EH" --reads "$bam" --reference genome.fa --variant-catalog catalog.json --output-prefix "${bam%.bam}" \
    --analysis-mode streaming --threads "$(nproc)" --catalog-cache catalog.cache
done
```

The optimization work focuses on streaming mode (`--analysis-mode streaming`):
the parallel read input, worker pool, BAM output pipeline, compact catalog and
parallel genotyping exist only there. Seeking mode stays the default, as in
ExpansionHunter v5.0.0, and keeps v5.0.0's per-locus loop; it gains the faster
alignment and genotyping code, catalog loading and output writing.
See [Installation](docs/02_Installation.md), [Usage](docs/03_Usage.md) and
[Running efficiently](docs/03_Usage.md#running-efficiently) for details.

## License

NuExpansionHunter is provided under the terms and conditions of the
[Apache License Version 2.0](LICENSE.txt). See [NOTICE.txt](NOTICE.txt) for
attribution and [COPYRIGHT.txt](COPYRIGHT.txt) for third-party licenses.
NuExpansionHunter is not affiliated with or endorsed by Illumina.

## Documentation

Installation instructions, usage guide, and description of file formats are
contained in the docs folder:

* [Introduction](docs/01_Introduction.md)
* [Installation](docs/02_Installation.md)
* [Usage](docs/03_Usage.md)
* [Input variant catalogs](docs/04_VariantCatalogFiles.md)
* [Output JSON files](docs/05_OutputJsonFiles.md)
* [Output VCF files](docs/06_OutputVcfFiles.md)
* [Performance and identity](docs/07_Performance.md)

## Companion tools and resources

* [Bundled tools](tools/README.md)
* [Variant catalogs](variant_catalog/)
* [Example inputs and outputs](example/)

## Citation

The method is described in the following papers:

- Egor Dolzhenko, Joke van Vugt, Richard Shaw, Mitch Bekritsky, and others,
  [Detection of long repeat expansions from PCR-free whole-genome sequence data](http://genome.cshlp.org/content/27/11/1895),
  Genome Research 2017

- Egor Dolzhenko, Viraj Deshpande, Felix Schlesinger, Peter Krusche, Roman Petrovski, and others,
[ExpansionHunter: A sequence-graph based tool to analyze variation in short tandem repeat regions](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btz431/5499079),
Bioinformatics 2019
