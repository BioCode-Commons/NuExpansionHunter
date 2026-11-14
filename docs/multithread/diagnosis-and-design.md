# Streaming parallel redesign: diagnosis and design

All reports and code comments for this change are in English.

## Scope and comparison point

The v6 implementation is frozen for benchmarking on `Jiapeng-multithread`,
created from `b00650d`; the `Jiapeng-RAMOPT` baseline branch is preserved.
The comparison binary is the previously validated RAMOPT implementation at
`b00650d297754f68d78778a441e2fcbf5096250b`. The work targets coordinate-sorted
BAM input, streaming analysis and prepared catalog caches. Measurements cover
cached sample analysis, including archive loading; catalog/index preparation
time and RAM are excluded. The primary workload has both many reads and many
repeat loci. The initial extra-memory reference is
20,000,000 bytes per requested compute worker. The user subsequently allowed
a modest excess, provided actual end-to-end peak RSS is reported honestly
and runtime does not regress. VCF and JSON bytes must match exactly; complete
BAM header and record bytes must match as a multiset, including duplicates.
Only BAM record order may change.

## What the published HG002 measurements establish

The public `BENCHMARKS.md` at revision `f3b1256` compares HappyFriday and
dp-readprocessing. It is an earlier comparison, distinct from this redesign's
RAMOPT baseline.
Both published versions used their own prepared catalog caches.

| HappyFriday phase | 16 workers (s) | 32 workers (s) | Difference (s) |
|---|---:|---:|---:|
| End to end | 427.191 | 474.047 | +46.856 |
| Read processing | 374.919 | 422.021 | +47.102 |
| Genotyping | 1.092 | 1.011 | -0.081 |
| Final output | 6.676 | 7.170 | +0.494 |
| Other elapsed time | 44.504 | 43.846 | -0.658 |

End-to-end time increases 10.97%; nearly all additional time is in read
processing. The dp-readprocessing baseline also slows down, by about 8.89%.
These measurements identify the affected phase but do not identify a unique
cause. The 16- and 32-worker runs were separate scheduler jobs; the available
report does not establish identical CPU/NUMA placement or otherwise identical
allocations across those two jobs. Version order also differs between them.

If automatic codec settings were used, the previous defaults configured
16+12+16=44 compute/decode/encode workers at T=16 and 32+12+32=76 at T=32,
plus input/output and codec coordination threads. Actual commands must be
checked before treating those counts as a property of the measured jobs.
The final-output timer excludes BAM production during read processing, so a
small final-output phase does not rule out an output bottleneck.

A source-level contention candidate was found in mapped graph traversal:
all locus graphs in a prepared catalog share an archive lifetime owner.
Owning adjacency ranges and their iterators repeatedly copy that shared owner.
That creates reference-count traffic on the same owner even for different
loci. New internal borrowed adjacency accessors avoid those copies while the
catalog/analyzer already guarantees graph lifetime. Public owning accessors
remain unchanged.

A completed controlled ownership ablation changed only the two borrowed
adjacency accessors to retain the archive owner. Both binaries otherwise used
the same v3 implementation, fixed codec pools of two decode and one
encode worker, and matching prepared caches. On the local 10,000-locus,
2,000,000-read fixture, three measured pairs per compute-worker setting gave:

| Compute workers | Median paired wall-time ratio, borrowed / retained | Median paired reduction |
|---|---:|---:|
| 4 | 0.793764 | 20.62% |
| 8 | 0.494501 | 50.55% |

Every measured pair matched exact JSON/VCF bytes and complete BAM header and
record-byte multisets. These figures isolate the ownership change within
v3; they are separate from the final v6 comparison against `b00650d` and do
not measure the later locale snapshot change. They support a concrete local
contention mechanism, without hardware-counter evidence or definitive
attribution to the original 32-worker HG002 run.
[Ownership-ablation measurements](evidence/owner-ablation-results.json) ·
[Ablation provenance](evidence/owner-ablation-manifest.json)

## Independent execution design

The new streaming execution layer uses its own fixed worker executor,
per-worker queues and work stealing. One stable workspace belongs to each
worker; an input publisher can execute queued work with a separate workspace
when admission would otherwise block. There is no scheduler queue or mutex
allocated for each catalog locus.

Ingress preserves primary-record admission, consecutive admitted-QNAME pairing
and the arriving record's routing metadata. It stores compact owned BAM fields
and defers sequence decoding until a routed pair is used. A shared pair is
decoded once, and each locus receives private mutable Read objects. Pending
raw payload is capped at 2,000,000 bytes; further unmatched records use exact
decoded-case sequences. Matching and missing-mate semantics are unchanged.
The pending table still grows with unmatched names; the raw-payload cap does
not bound all pending-read memory.

Alignment computations run independently in bounded chunks, including chunks
from the same locus. Reduction requires both computation completion and
completion of the preceding chunk for that locus. Workers do not wait for a
predecessor. Reduction replays statistics and variant evidence in encounter
order, preserving pair adjacency, tie handling, orientations and floating
accumulation order. Genotyping begins after EOF and complete reduction.

BAM records are created in reduction batches and transferred to an owned,
byte-bounded output sink. Its dedicated consumer can always drain independently
of compute-worker availability. Failed writes wake producers; all accepted
records retain RAII ownership through failure and shutdown. The same executor
is reused for initialization, alignment and genotyping. Pure sequence decoding
uses that executor rather than an additional preprocessing pool. Streaming
selectors also own an optional locale snapshot captured before alignment;
the scalar uppercase operation is preserved, and other callers continue using
the current global locale. In v6, each consumed item's temporary alignments
and oriented reads are released before the next item is consumed.

## Budget and diagnostics

`--threads` remains a compute-worker count, not a process thread limit.
Automatic streaming codec pools are capped at 4 decode and 2 encode workers;
explicit numeric settings remain supported. The default streaming BAM queue is
8 MiB; explicit zero remains unbounded. The legacy read-decode setting remains
accepted as a compatibility hint; decoding uses the shared executor.

The pipeline's byte charge is an admission estimate, not an allocator-exact
memory ceiling. Actual acceptance uses measured process RSS, with 20 MB per
compute worker retained as the original comparison reference and any excess
stated explicitly. Immutable catalogs, accumulated evidence, allocator
overhead, DP scratch and a single arbitrarily large input/alignment must be
distinguished from that estimate. Pending-read peaks, charged cache peaks, output resident
bytes and blocked output-submit time are logged separately. Summed worker
elapsed intervals include stalls and overlap; they are not CPU time.
Output-submit wait time is likewise summed across producers. Sink drain time
ends when its consumer joins and excludes the later HTS close/compression
flush; sink resident bytes are not a whole-process memory measurement.

## Local validation limits

The available machine is macOS ARM64 with 10 reported logical CPUs. Affinity
was not enforced, and exclusive host use or stable OS memory pressure was not
established. Local measurements cannot establish 16-to-32-core scaling on
the Genoa server. The runner in `tools/parallel_validation/` counterbalances
version and configuration order, records available topology/affinity data,
measures wall time, CPU time and peak RSS, and compares all outputs.
Its fixed-codec experiment and separate automatic-codec experiment answer
different questions. Real HG002 testing on the intended 32-core allocation
remains required before claiming that specific regression has been removed.

The final v6 primary comparison is complete: all three measured pairs were
faster and matched the full output policy. Median end-to-end time changed
from 184.0256 to 90.1127 seconds; median read-processing time changed from
169.5318 to 79.7250 seconds. Peak RSS increased in two of three pairs, and
the maximum observed peak rose from 3,839.406250 to 4,142.921875 MiB. This does
not establish a RAM reduction. All 13 stress cases and the independent timing recheck are complete. The
selected matrix contains 27 configurations and retains every observed slowdown
and RAM increase. The original clock-inconsistent configuration remains in
the evidence alongside its separately reported complete recheck. See the
[complete measurement matrix](measurement-overview.md) and
[regression audit](regression-audit.md).
[Final primary measurements](evidence/final-v6-primary.json)

Read-phase boundaries also differ: RAMOPT logs the phase before destroying
unmatched-read storage and the BAM reader, while v6 includes that cleanup.
Raw phase-time increases, especially for missing mates, therefore require
this boundary qualification alongside the end-to-end measurements. Phase
times are not adjusted or omitted to hide a regression. No further
optimization or memory profiling was performed during the final frozen
benchmark campaign.
