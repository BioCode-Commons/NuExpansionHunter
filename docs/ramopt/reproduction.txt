Streaming memory validation: exact index compression and per-locus storage

Scope and versions
------------------
The primary acceptance workload is coordinate-sorted BAM input, streaming
analysis, and a prepared catalog cache. Seeking mode is outside this work.
Owned/uncached runs and unsorted streaming inputs are supplementary diagnostics.
The target is approximately 50% less process memory without a runtime regression;
neither this README nor a smaller serialized index establishes that target for
every input. Consult the completed measurement records for each configuration.

The baseline is main at f3b1256. The frozen v1 candidate adds exact prefix-bucket
k-mer storage, compact eligible path payloads, fixed-size graph identity checks
in index archives, smaller per-locus queue metadata with pooled active queues,
and earlier release of streaming scratch/mates and completed analyzers.
Large or unsuitable indexes retain the previous full-width representation.
Existing baseline decoding, alignment and catalog optimizations are not new v1
changes. Read evidence is not capped, sampled, deduplicated, or discarded.

The local frozen snapshot is work/frozen-v1. Its manifest.json identifies the
base, changed product files and candidate executable SHA256; product.patch is
the product delta. The snapshot contains changed files, not a complete checkout.
Use complete source checkouts for the harness source-provenance arguments.

Candidate archives use format version 4. Earlier archives are rejected and must
be rebuilt with the corresponding executable and a new cache filename. Never
give the baseline and candidate the same cache. run_matrix.py prepares one
compatible cache per binary and fixture by default, outside measured runs.

A DP row-metadata compression experiment was reverted after local runtime
regression. It is not in v1. An isolated evidence-interning prototype was also
not included or validated as part of this product change.

Setup
-----
Run these commands from the candidate repository root. Choose a writable work
directory and fresh output paths; generators and benchmark runners generally
refuse to overwrite existing outputs. The local layout used sibling directories
work/NuExpansionHunter, work/baseline, work/build-baseline, work/build-candidate,
and work/frozen-v1. Substitute your own equivalent paths below.

  NU_MEM_WORK=/absolute/path/to/work
  NU_MEM_BASE="$NU_MEM_WORK/build-baseline/ExpansionHunter"
  NU_MEM_NEW="$NU_MEM_WORK/build-candidate/ExpansionHunter"

Build both complete application versions with identical compiler/build options
and dependencies, following the repository build instructions. The Python tools
use the standard library. Build tools/streaming_validation/sam_to_bam.c against
the same HTSlib as the application and place it at $NU_MEM_WORK/sam_to_bam.
It preserves records and builds a BAI, so its input must be coordinate sorted.
Build tools/memory_validation/sam_to_bam_streaming.c similarly and place it at
$NU_MEM_WORK/sam_to_bam_streaming. That converter preserves input order and does
not build an index; use it for the streaming stress generator, including its
unsorted cases. Do not substitute the indexed converter for unsorted input.

The candidate may be built by applying the delivered patch to f3b1256; no local
frozen-v1 snapshot is required. When an existing snapshot manifest is available,
optionally add --provenance-file /path/to/snapshot/manifest.json to run_matrix.py.
The harness records binary hashes and source information without that option.

Do all builds, fixture generation and independent validation outside timing
windows. Run benchmark processes serially, without other benchmarks, builds,
generation jobs or CPU-heavy background work. Do not run commands below in
parallel when collecting measurements.

Main catalog-size fixtures
--------------------------
These deterministic inputs have 150bp reads, ten reference CAG copies, mixed
qualities, one substitution per read and periodic missing quality fields.
The default seed is nueh-round5-v1. Extraction extension must remain 1000bp.

  python3 tools/streaming_validation/generate_catalog.py \
    --out "$NU_MEM_WORK/fixtures/covered10k" \
    --loci 10000 --covered-loci 10000 --pairs-per-locus 50 \
    --mismatch-every 1 --missing-quality-every 16 \
    --bam-helper "$NU_MEM_WORK/sam_to_bam"

  python3 tools/streaming_validation/generate_catalog.py \
    --out "$NU_MEM_WORK/fixtures/sparse100k" \
    --loci 100000 --covered-loci 100 --pairs-per-locus 20 \
    --mismatch-every 1 --missing-quality-every 16 \
    --bam-helper "$NU_MEM_WORK/sam_to_bam"

covered10k contains 1,000,000 records; sparse100k contains 4,000 records over a
100,000-locus catalog. They emphasize different sources of memory. In observed
covered10k outputs, estimated coverage is 8.1081, below the genotyping threshold.
Its LowDepth results do not establish that genotype likelihood computation was
exercised. Empty reported evidence counts also do not prove that no alignment
evidence accumulated before the low-coverage check.

Read and streaming stress fixtures
----------------------------------
The read generator writes four sorted cases: longread1000, higherror-indels150,
iupac-missingquality150 and expandedrepeat150. Each command below creates
200 loci and 25,600 records per case. Manifests describe exact changes and input
validation. The uncalibrated expansion case is not a genotype-positive assertion.

  python3 tools/memory_validation/generate_read_extremes.py \
    --out "$NU_MEM_WORK/fixtures/read-extremes" \
    --loci 200 --pairs-per-locus 64 \
    --sam-to-bam "$NU_MEM_WORK/sam_to_bam"

The following two commands reproduce the selected streaming stress scales.
shuffled and stripe-collisions are intentionally unsorted, supplementary cases.
far-apart, missing-mates, dense-fanout, hot-locus and empty-input are sorted.
Missing mates and empty input are intentional input conditions, not output
filtering. All generated input records are supplied to the application.

  python3 tools/memory_validation/generate_streaming_extremes.py \
    --out "$NU_MEM_WORK/fixtures/stream-extremes" \
    --sam-to-bam "$NU_MEM_WORK/sam_to_bam_streaming" \
    --only shuffled far-apart missing-mates --scale 8

  python3 tools/memory_validation/generate_streaming_extremes.py \
    --out "$NU_MEM_WORK/fixtures/stream-extremes" \
    --sam-to-bam "$NU_MEM_WORK/sam_to_bam_streaming" \
    --only dense-fanout hot-locus stripe-collisions empty-input --scale 1

For the longer sorted dense-fanout confirmation, generate scale 5 separately
and give the fixture a distinct basename. This creates 50,000 input records,
with each complete pair routed to 24 distinct locus IDs at the same coordinate.

  python3 tools/memory_validation/generate_streaming_extremes.py \
    --out "$NU_MEM_WORK/fixtures/stream-confirm" \
    --sam-to-bam "$NU_MEM_WORK/sam_to_bam_streaming" \
    --only dense-fanout --scale 5

  mv "$NU_MEM_WORK/fixtures/stream-confirm/dense-fanout" \
     "$NU_MEM_WORK/fixtures/stream-confirm/dense-fanout-large"

Genotype-positive supplemental fixtures
---------------------------------------
genotyped10k doubles target pairs to 100 per locus: 2,000,000 records. The
calibrated expansion case has 256 expanded-haplotype pairs plus eight exact
left-flank calibration pairs per locus: 105,600 records in total. Both reads of
each calibration pair lie wholly on the left flank, supplying a 300bp fragment
length observation. No biological truth claim is made for the called expansion.

  python3 tools/streaming_validation/generate_catalog.py \
    --out "$NU_MEM_WORK/fixtures/genotyped10k" \
    --loci 10000 --covered-loci 10000 --pairs-per-locus 100 \
    --mismatch-every 1 --missing-quality-every 16 \
    --bam-helper "$NU_MEM_WORK/sam_to_bam"

  python3 tools/memory_validation/generate_read_extremes.py \
    --out "$NU_MEM_WORK/fixtures/genotyped-read-extremes" \
    --only expandedrepeat150 --loci 200 --pairs-per-locus 256 \
    --calibration-pairs-per-locus 8 --sam-to-bam "$NU_MEM_WORK/sam_to_bam"

  mv "$NU_MEM_WORK/fixtures/genotyped-read-extremes/expandedrepeat150" \
     "$NU_MEM_WORK/fixtures/genotyped-read-extremes/expandedrepeat150-genotyped"

The positive expansion directory was renamed expandedrepeat150-genotyped to
distinguish it from the original expansion fixture and give run_matrix.py unique
input basenames. The generator's profile name inside its manifest remains
expandedrepeat150. Existing positive fixtures already use the renamed directory;
do not attempt to regenerate or rename them in place.

Non-timed candidate smoke checks observed calls at all 10,000 and all 200 loci,
respectively, with coverage 16.2162 and 22.0541 and fragment length 300. These
smoke observations are not paired performance results or genotype truth. Use the
validator below on every formal output to assert actual genotype execution.

Complete application measurements
---------------------------------
Example primary cached matrix, five measured pairs per configuration:

  python3 tools/frontend_validation/run_matrix.py \
    --baseline "$NU_MEM_BASE" --candidate "$NU_MEM_NEW" \
    --baseline-source "$NU_MEM_WORK/baseline" \
    --candidate-source "$NU_MEM_WORK/NuExpansionHunter" \
    --input "$NU_MEM_WORK/fixtures/covered10k" \
    --input "$NU_MEM_WORK/fixtures/sparse100k" \
    --input "$NU_MEM_WORK/fixtures/genotyped10k" \
    --out "$NU_MEM_WORK/reproduce-main-cached" \
    --threads 8 --preprocess 0 --decode 2 --encode 1 --pairs 5

Run the longer sorted dense-fanout confirmation with nine measured pairs and
both one and eight compute workers:

  python3 tools/frontend_validation/run_matrix.py \
    --baseline "$NU_MEM_BASE" --candidate "$NU_MEM_NEW" \
    --baseline-source "$NU_MEM_WORK/baseline" \
    --candidate-source "$NU_MEM_WORK/NuExpansionHunter" \
    --input "$NU_MEM_WORK/fixtures/stream-confirm/dense-fanout-large" \
    --out "$NU_MEM_WORK/reproduce-dense-confirm" \
    --threads 1 8 --preprocess 0 --decode 2 --encode 1 --pairs 9

For owned-index diagnostics, repeat with --uncached and a new output directory.
For additional explicit worker configurations use, for example, --threads 1 8.
The matrix also accepts --decode auto --encode auto; run that as a separate
configuration set and retain the resolved pool counts from the logs. Both sides
must receive the same requested thread policy. Worker counts are not a cap on
total process threads or a simulation of the corresponding physical-core count.

Example sorted stress matrix, three measured pairs per configuration:

  python3 tools/memory_validation/run_extremes.py \
    --baseline "$NU_MEM_BASE" --candidate "$NU_MEM_NEW" \
    --baseline-source "$NU_MEM_WORK/baseline" \
    --candidate-source "$NU_MEM_WORK/NuExpansionHunter" \
    --input-root "$NU_MEM_WORK/fixtures/read-extremes" \
    --input-root "$NU_MEM_WORK/fixtures/stream-extremes" \
    --input-root "$NU_MEM_WORK/fixtures/genotyped-read-extremes" \
    --select longread1000 higherror-indels150 iupac-missingquality150 \
      expandedrepeat150 far-apart missing-mates dense-fanout hot-locus \
      empty-input expandedrepeat150-genotyped \
    --out "$NU_MEM_WORK/reproduce-sorted-extremes" \
    --modes mapped --threads 1 8 --pairs 3 --timeout 900

Use a separate run_extremes invocation selecting shuffled stripe-collisions for
unsorted supplemental diagnostics. It automatically uses 2 and 8 compute workers
for stripe-collisions; this is intended to exercise a shared mutex stripe and
active-queue capacity, not to guarantee a particular scheduling interleaving.
run_extremes accepts --modes owned mapped and preserves failed or timed-out runs
in its outer results.json. A successful process exit alone is not an assertion
of the 50% memory target or the no-regression runtime requirement.

Measurement and equality policy
-------------------------------
Each configuration has one excluded warmup pair, then alternating AB/BA measured
pairs. The warmup is also checked for output equality. The runner uses
posix_spawn and wait4: wall time covers process start through exit, CPU time is
user plus system time, and memory is that process's maximum resident set size.
Darwin ru_maxrss is bytes; the helper converts Linux KiB to bytes. Reports retain
both the median of per-run peak RSS values and their maximum. These are process
high-water observations, not live-allocation totals or a universal memory bound.

Cache preparation, hashes and output comparison run outside timed samples;
cache loading, genotyping, output and final BAM drain remain in total wall time.
The read_processing phase includes reading, decode, pairing, dispatch and graph
alignment completion. It is neither a DP-only timer nor an end-to-end timer.
OS caches are not cleared, so these runs do not measure controlled cold-cache I/O.

Every baseline/candidate pair compares exact VCF and JSON bytes and the exact
BAM header plus full-record-byte multiset, preserving duplicate multiplicity.
Only BAM record order may differ. No BAM fields, quality values, tags, records or
outputs are removed to obtain a match. BAM compressed file hashes may differ
because order/compression differ; the complete decompressed record comparison
is the equality gate. A mismatch fails the run and retains the evidence.

Percent changes use the median of paired candidate/baseline ratios, which can
differ from the ratio of two marginal medians. Bootstrap intervals are descriptive
local evidence. An interval crossing 1 is inconclusive, not proof of no regression.
Short runs and three-pair samples provide limited timing precision. The generic
run_matrix reported RAM gate defaults to 1.05; it is not the separate 0.50 target
and a failed RAM-gate flag alone does not make that runner exit nonzero.
run_extremes records whether both median and maximum peak RSS meet that target.

Observed genotype validation and report generation
--------------------------------------------------
Validate each positive fixture's run subtree separately, including baseline,
candidate, warmup and measured outputs. The validator requires coverage >=10,
positive read and fragment lengths, actual repeat Genotype and confidence
interval fields, and VCF GT/REPCN agreement with JSON. --require-evidence also
requires reported repeat evidence. It exits nonzero on any failed assertion.

  python3 tools/memory_validation/validate_genotype_outputs.py \
    --expected-loci 10000 --require-evidence \
    --out "$NU_MEM_WORK/genotyped10k-validation.json" \
    "$NU_MEM_WORK/reproduce-main-cached/genotyped10k"

  python3 tools/memory_validation/validate_genotype_outputs.py \
    --expected-loci 200 --require-evidence \
    --out "$NU_MEM_WORK/expanded-positive-validation.json" \
    "$NU_MEM_WORK/reproduce-sorted-extremes/expandedrepeat150-genotyped-mapped"

Do not point this validator at a mixed root containing low-depth, empty, sparse
or differently sized fixtures. It deliberately does not filter out missing calls:
doing so would let an input that never genotyped pass vacuously. As a negative
control, covered10k must fail the same 10,000-locus assertion. Genotype validation
does not replace baseline/candidate full-output comparison or assert truth.

After the requested measurements finish, collect their actual records:

  python3 tools/memory_validation/report.py \
    --matrix "$NU_MEM_WORK/reproduce-main-cached/results.json" \
    --extremes "$NU_MEM_WORK/reproduce-sorted-extremes" \
    --out "$NU_MEM_WORK/reproduce-report"

--matrix takes a run_matrix results.json file; --extremes takes a run_extremes
output directory. Both options may be repeated. The report writes
benchmark-data.json and benchmark-comparison.html, retaining raw matrices,
preparation measurements and recorded stress failures. It is a collector, not
an independent completion gate: inspect complete/status fields and retained
failures before declaring acceptance. Do not silently omit inconvenient runs.

Interpretation limits and focused index checks
----------------------------------------------
The available measurements are synthetic diagnostics on one macOS ARM64 host.
They do not establish real-WGS, clinical, Linux, x86/SSE2, AVX2, multi-socket or
32-physical-core performance. Heap allocation, mapped-page residency, mate
distance, read length, catalog structure and retained evidence change the
fraction of process RSS that these optimizations can remove. The 50% target is
not promised for small catalogs or all stress inputs.

For independent exhaustive index/path checks, boundary/fallback cases, archive
lifetime and focused query timing, see tools/memory_validation/index_extremes.md.
Those probes measure exact index behavior and cache-hot queries; their timings
and serialized byte counts must not be substituted for application wall time or
process RSS. Retain executable/source/input hashes, command lines and the exact
output-comparison records when publishing any measured conclusion.
