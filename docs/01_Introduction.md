# A guide to NuExpansionHunter software

NuExpansionHunter is a tool for targeted genotyping of short tandem repeats and
flanking variants. It operates by performing a targeted search through a
BAM/CRAM file for reads that span, flank, and are fully contained in each
repeat. The program provides limited support for insertions,
deletions, and substitutions.

NuExpansionHunter is distributed as source (see [Installation](02_Installation.md));
in order to use it you need to (1) build the program from source, (2) obtain a set
of files specifying repeat regions of interest, and (3) have a BAM or a CRAM file
with alignments of reads from a PCR-free WGS sample.

NuExpansionHunter has two analysis modes (see [Usage](03_Usage.md)). The
optimization work focuses on streaming mode, which reads the alignment file in
one pass with parallel input, a worker pool, a BAM output pipeline, a compact
catalog and parallel genotyping. Seeking mode stays the default, as in
ExpansionHunter v5.0.0, and keeps v5.0.0's per-locus loop; it gains the faster
alignment and genotyping code, catalog loading and output writing. For large
catalogs, use `--analysis-mode streaming` with `--threads` set to the number of
hardware threads.

The following sections will help you to get started.

* [Installation](02_Installation.md)
* [Usage](03_Usage.md)
* [Input variant catalogs](04_VariantCatalogFiles.md)
* [Output JSON files](05_OutputJsonFiles.md)
* [Output VCF files](06_OutputVcfFiles.md)
