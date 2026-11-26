# Performance and identity

## What is optimized

The optimization work focuses on streaming mode. Items marked *(streaming)*
apply only to streaming mode; the others apply to both modes. Seeking mode
keeps ExpansionHunter v5.0.0's per-locus loop and index-based read retrieval.

* Vectorized affine dynamic programming (DP), with certified linear-target bounds and an early row stop (K-38).
* Dense alignment by default. Adaptive budgeted sparse DP is available for differential verification.
* Single-node and unique-ungapped fast paths proven by cost bounds.
* Prepared-graph cache.
* *(streaming)* Bounded reuse of complete, ordered single-node extension traces across loci (K-37).
* Frozen k-mer index layouts and SIMD probes.
* SIMD read decoding and BAM packing.
* Direct selection from certified canonical seed products (K-39) and canonical-intersection step skipping.
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
is covered by emulated AArch64 tests. The extension cache, certified target bounds
and canonical product selection use the same eligibility rules in every SIMD tier.

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
11. **Extension traces.** Streaming reuse is local to one pipeline and requires
    the same read pair, physical mate, orientation, direction, complete query,
    and certified target prefix. Both loci independently satisfy the default-score,
    single-node and score-width guards. Accepted hits reproduce every ordered
    optimal traceback and remap its node IDs. A miss solves the receiving target.
12. **Alignment scope.** Certified linear bounds and canonical product selection
    preserve the optimal endpoints and canonical alignment consumed by genotyping.
    Internal second-best scores and nonoptimal matrix cells are outside this contract.

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
Catalog compaction releases construction storage. Optional jemalloc returns
freed catalog memory. Streaming releases unmatched mates and per-locus indexes
before genotyping. Some thread-local alignment buffers can stay allocated
until the worker exits. Analyzers are released as their findings are produced.

## Benchmarks

For whole-genome GangSTR streaming at 16 threads, NuExpansionHunter v1.1.0
took a median of 10 minutes with 10.8 GiB median peak RSS across 23 matched
pairs. Stock ExpansionHunter v5.0.0 took 1.8 hours with 82.3 GiB median peak
RSS. NuExpansionHunter’s largest peak at 16 threads was 14.5 GiB.

The full TRExplorer v2.1 catalog used 63.2–77.4 GiB of peak RSS at 16–32
threads. Some runs exceeded 64 GiB. Even the lowest measured peak left less
than 1 GiB on a 64 GiB machine for the operating system and other processes.
Set `--threads` based on the logical CPUs and RAM available to the job.
With limited RAM, lower `--threads` or split the catalog. [Usage](03_Usage.md#running-efficiently)
See the [benchmark report](08_Benchmark.md) for results by catalog and thread
count.
