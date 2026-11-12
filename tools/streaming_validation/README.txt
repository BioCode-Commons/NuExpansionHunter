OPT-IN EXPERIMENTAL FEATURE: benchmark on the intended host before enabling.

Bounded ordered sequence preprocessing for streaming

Use --analysis-mode streaming --read-decode-threads 1 to enable one additional
sequence preprocessing worker. The default is 0 (synchronous). This differs from
--decode-threads, which controls HTS decompression. --threads continues to mean
locus computation workers; it is NOT a total process thread budget. For a fixed
8-worker comparison use baseline --threads 8 and candidate --threads 7
--read-decode-threads 1, with identical explicit HTS decode/encode settings.
No automatic 32-core sizing or affinity policy is enabled.

The caller remains the single HTS reader and ordered dispatcher. It performs
primary/paired/coarse-mask filtering, copies only QNAME/packed bases/qualities
and routing metadata into two fixed slots per helper, and consumes decoded
records in input order. Each slot has at most 1024 records and 512 KiB raw data.
Workers perform pure base/quality decoding. Pairing, mate-coordinate semantics,
locus fanout, and per-locus evidence order remain on the original caller.
No shared concurrent mate map or reordered evidence is introduced.

Per-worker batch cache limit is 20,000,000 bytes (decimal MB); capacity accounting
includes both slots, raw/metadata/result vector capacities and owned strings.
This is not a process memory limit: allocator overhead, thread stacks, manager
controls, existing HTS buffers, mate/evidence queues and indexes are separate.
A single record exceeding a raw slot is decoded synchronously after preceding
batches drain, without a second large raw copy. The CLI logs conservative cache
high-water bytes. Storage grows with helper count, not catalog locus count.

Correctness: exact VCF and JSON bytes; exact BAM header and full BAM record byte
multiset, including duplicate counts. Only BAM record order may differ. The
validator externally sorts complete record bytes and compares actual bytes.
Stdout includes timestamps/telemetry and is not a result-file equality gate.

Run portable validation tools from this directory (Python3 standard library):
  python3 paired_streaming.py --baseline /path/base --candidate /path/new \
    --input /path/fixture --out /new/results --pairs 9 --compute 8 \
    --candidate-compute 7 --candidate-extra-json '["--read-decode-threads","1"]' \
    --candidate-preprocess-workers 1 --hts-decode 0 --encode 0 --cached
The output directory must not exist. Fixtures contain reference.fa, its .fai,
catalog.json, reads.bam, and optional manifest.json. --cached prepares one cache
per binary outside timed runs. --baseline-cache and --candidate-cache can reuse
matching caches with --cached; the application validates cache signatures.
Use --candidate-preprocess-workers 0 for the default-path control. All CPU-heavy
builds, fixture generation and benchmark processes must be serial, not parallel.
A discarded warmup plus 9 alternating AB/BA pairs records wall, user+system CPU,
and wait4 peak RSS. Reports preserve the original 1%/256KiB RSS noise allowance
and the separate authorized cache envelope. Retention is a human decision, not
an automatic conclusion from a fast pair or checksum.

Generate independent-coordinate synthetic loci with generate_catalog.py.
Compile sam_to_bam.c against the same HTSlib used by the application (or convert
its sorted SAM with samtools). Supply --bam-helper /path/sam_to_bam when using
the bundled helper. --loci 10000 --covered-loci 1000 --pairs-per-locus 50 gives
100,000 target records. These fixtures are reference-repeat-only; they are not
WGS/expanded-allele truth sets. --background-pairs-per-locus intentionally
separates frontend work by adding reads inside the coarse mask but outside
exact targets. Do not present that diagnostic as typical WGS throughput.

New tests cover quality thresholds/IUPAC, flags, coordinates, repeated names,
record/byte batch boundaries, filtering, empty inputs, oversized fallback,
error-prefix preservation and early destruction. Run the repository CTest suite.
Current host validation was macOS ARM64 with 10 reported logical CPUs. Linux
32-core throughput/affinity and AVX2 require separate hardware validation.

Observed local results (9 alternating pairs + discarded warmup, cached analysis;
percentages are medians of paired ratios, not ratios of marginal medians):
- 10k loci/2.004M frontend-stress records, compute8 vs compute7+helper1:
  wall -13.10%, CPU -2.04%, median peak RSS +0.28 MiB.
- Same input, compute8 vs compute8+helper1:
  wall -13.22%, CPU +2.42%, median peak RSS +0.52 MiB.
- 10k loci/100k target records, compute8 vs compute7+helper1:
  wall -9.41%, CPU -20.19%, median peak RSS +0.38 MiB.
- Same target input, compute7 vs compute7+helper1:
  wall +0.76%, CPU -0.26%, intervals cross parity: no demonstrated gain.
- Default0 frontend control: wall +0.91%, CPU -0.05%, intervals cross parity;
  RAM slightly lower. No proven speed gain or regression in this control.
All compared result files passed the exact equality policy. Compute-heavy CPU
savings in the fixed-budget case cannot be attributed solely to decoding: the
same-compute control is neutral. Preprocessing is an opt-in for a measured
frontend bottleneck, not an automatic recommendation for all datasets.
A 1M-target-record calibration had divergent warmup/formal timings; one pair
was not accepted as evidence of acceleration. No failed run was discarded and
relabelled as a successful timed pair.

The final unit run passes 139/140 EH tests and all30 GraphTools test programs.
The single PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime
failure is also present in the accepted baseline: graph edge-label JSON array
order [b,a] versus [a,b]. It is not claimed fixed by this change. Nine new decoder
tests pass under ASan/UBSan when those decoder TUs and tests are instrumented;
this does not establish that every linked dependency/the whole app is sanitized.

100,000-locus uncached check (4000 reads, 9 formal AB/BA pairs): all outputs
match, but the original process-RSS gates and the extra20MB envelope both FAIL.
Median peak RSS: 4,041,048,064 -> 4,108,648,448 bytes (+1.673%).
Maximum observed peak RSS: 4,107,698,176 -> 4,805,476,352 bytes (+697,778,176).
Paired wall +7.136% (interval crosses parity); CPU +5.306%. This sparse fixture
is dominated by existing catalog/index initialization, not read processing.
The two preliminary guards showed large RSS differences in opposite directions;
that does not invalidate or excuse the failed formal gate. The cause of the
large RSS variation remains unresolved. No stable large-catalog acceleration
or memory equivalence is claimed. This source is published for further target-
host benchmarking, with preprocessing disabled by default.

Focused ThreadSanitizer: all9 new decoder tests pass. Only decoder/helper/test
translation units were instrumented; this is not a whole-application race proof.
AVX2 and Linux dispatcher affinity are excluded from this feature.
