# Performance and identity

## What is optimized

The optimization work focuses on streaming mode. Items marked *(streaming)*
apply only to streaming mode; the others apply to both modes. Seeking mode
keeps ExpansionHunter v5.0.0's per-locus loop and index-based read retrieval.

* Vectorized affine dynamic programming (DP) with an early row stop.
* Exact budgeted sparse DP for low-cost extensions, with a measured search cap and skip history.
* Single-node and unique-ungapped fast paths proven by cost bounds.
* Prepared-graph cache.
* Frozen k-mer index layouts and SIMD probes.
* SIMD read decoding and BAM packing.
* Canonical-alignment fold that skips steps proven unchanged.
* *(streaming)* Parallel streaming pipeline with ordered per-locus reduction: pairing lanes, staging, a worker pool, a byte-bounded BAM sink and positional BGZF output.
* Flat evidence rows.
* *(streaming)* Compact frozen catalog. Its reusable cache (`--catalog-cache`) also serves seeking mode.
* *(streaming)* Parallel genotype grid.
* Parallel catalog loading and output writing.

## SIMD kernels

Every vector kernel has a scalar reference next to it and produces the same result bit for bit.
On x86-64 the program selects the widest supported tier at run time, so a portable build uses AVX2 or
AVX-512 where the CPU has them:

| Kernel | x86-64 tiers (run-time selection) | AArch64 |
| --- | --- | --- |
| DP row passes (deletion, match, consolidate, fused) | AVX2, SSE2 | NEON |
| Insertion scan | AVX-512BW, AVX2, scalar | NEON |
| Linear mismatch scan, read decoding, BAM packing | AVX2, scalar | NEON |
| Base-run matching, k-mer prefix probe | SSE2 | NEON |
| Orientation k-mer packing | SSSE3, scalar | NEON |

Three choices are deliberate. The DP row passes have no AVX-512 version: on Zen 4, 512-bit row passes were
slower than 256-bit ones for rows of a few hundred cells. The short per-read scans (mismatch scan, read decoding, BAM packing) have no
AVX-512 version: reads of 100-250 bases take only a few 32-byte AVX2 steps. They also have no SSE tier: CPUs
without AVX2 run them scalar, with identical output. NEON is part of AArch64 and is selected at compile time; it
has been tested under emulation, not on ARM hardware.

## Identity guarantee and its limits

For every input on which ExpansionHunter v5.0.0 has defined behaviour and which
stays inside the limits below, NuExpansionHunter writes the same VCF bytes,
the same JSON bytes, and a BAM with the same header and the same multiset of
decoded records as ExpansionHunter v5.0.0. See the
[verification documentation](../ehunter/verification/README.md) for the contracts.

1. **BAM file bytes and order.** Streaming output uses compression level 1 and
   fixed `0xff00`-byte BGZF blocks, with ISA-L when available. Record order
   depends on scheduling and differs from v5.0.0's per-locus input order.
   The header and decoded-record multiset are equal. Seeking mode builds
   records as v5.0.0 does and compresses them with htslib at its default
   level; with more than one thread, record order also depends on scheduling.
2. **Unused sequence nibble.** For an odd-length read, NuExpansionHunter writes
   zero to the unused low nibble of the last sequence byte; v5.0.0 leaves it
   uninitialized. Decoded records are equal.
3. **Size limits.** An evidence row's best-alignment index above 65,535 or a
   locus with more than 2^32 cells causes an error. v5.0.0 has neither limit.
4. **Graph consistency.** Freezing a catalog rejects a node whose reference
   projection length differs from its node length. Streaming always freezes
   the catalog.
5. **Error paths.** Seeking reports an uncatalogued read or inability to jump
   before decoding any mate. Failed BAM record writes throw, output files are
   checked on close, and an odd row count in fragment likelihoods throws.
6. **Negative contig indices.** Genome-mask queries return false for a negative
   contig index; v5.0.0 indexes out of bounds on that input.
7. **Locale.** Seed search uppercases using a locale captured once per aligner
   selector. Equality requires the global locale to remain constant during a run.
8. **Floating point.** Genotype-likelihood equality requires the same evaluation
   order, no fused multiply-add, and v5.0.0 built with the same compiler and math
   library. The CMake files add `-ffp-contract=off`; do not use `-ffast-math`.
9. **htslib coupling.** The in-place BAM reader depends on the BGZF and iterator
   layout of htslib 1.10.2, which the superbuild pins and patches as described in
   [Installation](02_Installation.md).
10. **Catalog caches.** The schema id ties each cache to the exact source tree.
    A cache from a different tree is rejected with an error; there is no fallback
    to parsing. Rebuild caches whenever the program is rebuilt.

## Resource controls

`--threads` sets compute workers; `--decode-threads` and `--encode-threads`
control additional codec pools. Their automatic formulas are in [Usage](03_Usage.md).
`--bamlet-queue-mb` bounds queued BAM bytes (automatically 8 MiB, or 16 MiB
with at least 16 compute workers, in both modes).
Explicit `0` is unbounded. The budget includes record structs, allocated record
buffers, batch pointer storage and compressed-block reservations, including the
consumer's batch. One oversized batch may enter an empty sink to ensure progress.
Worker-local packing buffers, codec scratch and bookkeeping are separate; this
is not a process RAM limit. Producers wait when the budget is full.
Catalog compaction releases construction storage; optional jemalloc returns
freed catalog memory. Unmatched mates and worker scratch are released before
genotyping, and analyzers are released as their findings are produced.

## Performance

Wall times on HG002 (GIAB, 35x PCR-free, aligned to GRCh38 no-alt), streaming mode, on an AMD EPYC 9634
(Zen 4) with the package's reference build: GCC 15.2, `-march=znver4 -ffp-contract=off`, libdeflate, ISA-L
and jemalloc. "Cores" are physical cores with both SMT threads. ExpansionHunter v5.0.0 ran with
`--threads` equal to the hardware threads; NuExpansionHunter ran with `--threads` 2.5 times the cores,
`--decode-threads` half and `--encode-threads` a quarter of the cores, and `--bamlet-queue-mb 256`.
Every run listed wrote the same VCF, JSON and BAM records as ExpansionHunter v5.0.0. NuExpansionHunter times
are medians of two or three runs (single runs for TRExplorer chr22 and the whole genome); PGO denotes the
profile-guided build of the same source, trained on a chr22 slice, the chr21-22 GangSTR set and one satellite
window. A dash marks a configuration that was not measured.

| Dataset (catalog on HG002) | Cores | ExpansionHunter v5.0.0 | NuExpansionHunter | NuExpansionHunter, PGO | Speed-up (no PGO) |
| --- | ---: | ---: | ---: | ---: | ---: |
| GangSTR hg38 v13, chr21-22 | 16 | 373.2 s | 17.0 s | 15.9 s | 22.0x |
| GangSTR hg38 v13, chr21-22 | 32 | — | 9.4 s | — | — |
| TRExplorer v2.1, chr22 | 16 | 1,121.0 s | 129.3 s | 118.9 s | 8.7x |
| TRExplorer v2.1, chr21-22 | 16 | 2,585.1 s | 293.8 s | 271.4 s | 8.8x |
| TRExplorer v2.1, chr21-22 | 32 | — | 148.6 s | — | — |
| GangSTR hg38 v13, whole genome | 16 | — | 383.9 s | 360.7 s | — |

Peak memory on the whole-genome GangSTR catalog at 16 cores is 16.4 GiB. A `-march` build runs only on CPUs
with the instructions it was compiled for; a portable build without `-march` is slower.
