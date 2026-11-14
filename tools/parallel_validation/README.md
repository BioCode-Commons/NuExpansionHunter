# Same-host streaming parallel validation

This runner compares the frozen RAMOPT baseline (`b00650d`) with the candidate
on the **same host and CPU allocation**, using streaming mode and prepared
catalog caches. It does not implement or change application scheduling.
Python 3.8 or later and Linux or macOS are required; no Python packages are
needed. The two executables must already be built.

Reserve one otherwise idle host with at least 32 physical cores for the server
matrix. Do not build, generate fixtures, hash unrelated large files, or run other
benchmarks concurrently. Keep CPU governor, SMT, NUMA policy, compiler settings,
storage, input data and clock/frequency policy unchanged between versions.
The runner does not set frequencies, clear the OS page cache or bind memory.
A 16-worker run and a 32-worker run in separate scheduler jobs are not a
same-allocation comparison merely because both jobs request the same node type.

## Server example

Inspect the allocated physical cores first, for example with
`lscpu -e=CPU,CORE,SOCKET,NODE,ONLINE`. The CPU list below is only an example;
replace it with 32 allowed logical CPU IDs representing 32 distinct physical
cores. All configurations use this **same complete mask**, including the
16-worker runs. Input, compression and coordination threads also share it.

```sh
python3 tools/parallel_validation/run.py \
  --baseline /builds/ramopt-b00650d/ExpansionHunter \
  --candidate /builds/candidate/ExpansionHunter \
  --baseline-source /sources/ramopt-b00650d \
  --candidate-source /sources/candidate \
  --reads /data/HG002.sorted.bam \
  --reference /data/GRCh38.fa \
  --variant-catalog /data/catalog.json \
  --threads 16 24 32 --pairs 5 \
  --cpu-list 0-31 \
  --decode 2 --encode 1 \
  --out /results/hg002-parallel-001
```

`--threads` retains its application meaning: compute workers, **not** a total
process thread limit or a physical-core count. Linux `--cpu-list` validates the
current allowed allocation and rejects SMT siblings. The runner records CPU,
socket/core, available NUMA-node and LLC-sharing information and checks that its
effective affinity does not change between runs. Without this option it records
the inherited Linux affinity but does not certify distinct physical cores.
macOS has no portable affinity enforcement here; local `--threads 1 4 8` runs
are useful for validation, not substitutes for the 32-core server result.

The defaults are five pairs at 16/24/32 workers, input decode 2 and output encode
1. These fixed codec settings avoid independently scaling the codec pools
while comparing the two versions.
`--decode auto --encode auto` is supported for a separate end-user-defaults
experiment; the resolved pools can then differ between binaries and must be
considered when interpreting results. The runner preserves each binary's BAM
queue default unless `--bamlet-queue-mb N` supplies a common explicit override.
`--read-decode-threads 0` is passed to both versions: the candidate decodes
sequences through its shared compute workers, without an extra preprocessing
pool. All resolved pool logs and new batch telemetry are retained.

## Caches and ordering

Separate baseline and candidate archives are prepared before measurements,
using the corresponding executable. Preparation is recorded but excluded from
summaries. Alternatively provide **both** `--baseline-cache` and
`--candidate-cache` with distinct existing paths. Each application validates its
own archive; the runner does not assume archives are compatible across builds.
Cache files and complete inputs are hashed outside timings. Archive loading
remains part of each application's measured wall/CPU time.

Each thread configuration receives one excluded warmup pair. Its complete
outputs are compared, and baseline warmups are also compared across thread
configurations. Measured rounds alternate baseline/candidate and
candidate/baseline. Thread-configuration order reverses each round and rotates
every two rounds, so one configuration is not always tested first. Five rounds
cannot perfectly balance every position; the exact schedule is saved. Runs are
strictly serial. No application process overlaps hashing or output comparison.

Application wall time, user+system CPU time and peak RSS come from the existing
`tools/streaming_validation/measure.py` `posix_spawn`/`wait4` implementation.
Linux RSS is converted from KiB; macOS reports bytes. A separate driver process
group makes interruption cleanup possible; driver startup is outside the
application timer. These measurements include the application's threads, not
the Python comparison process. Log-derived phase timings and batch/assistance
telemetry are additional diagnostics and can be absent in older binaries.
`alignment_worker_seconds` and `reduction_worker_seconds` sum elapsed intervals
across workers; they can overlap and include stalls. They are **not** measured
CPU time. Use `user_seconds + system_seconds` for process CPU consumption.

## Exact output checks and disk usage

The existing `tools/streaming_validation/compare_outputs.py` compares every VCF
and JSON byte, the complete BAM header, and the multiset of **complete BAM
record bytes**. External sorting preserves duplicate multiplicity. There is no
record sampling, tag filtering or hash-only acceptance. Compressed BAM bytes
may differ because record order is allowed to change; both compressed hashes
and canonical record hashes are retained as diagnostics.

By default all outputs are retained. Large real-data runs can use
`--discard-equal-output`: raw VCF/JSON/BAM files are deleted only after complete
equality and their hashes have been durably saved. The cross-configuration
warmup reference remains until every warmup comparison has finished. Logs,
measurement files, comparison results and caches remain. Unequal or incomplete
pairs are never deleted. Comparison scratch files are removed by the helper;
allow scratch space for external sorting of both BAM files. The default
comparison chunk size is 16 MiB per side, configurable with
`--compare-chunk-mib`; a single oversized BAM record can exceed that chunk size.

## Reading results

`results.json` is updated atomically after every preparation, application run,
comparison and cleanup. It includes excluded warmups, commands, binary/input
SHA-256, optional source Git status/diff and untracked-file hashes, affinity,
timestamps, stdout/stderr paths, resolved pools and phase/batch telemetry.
Source metadata is provenance, not proof that an executable was built from that
tree. Supply the frozen baseline executable; the runner does not fetch or build
`b00650d`.

For an isolated experiment with a different frozen baseline, pass a descriptive
`--baseline-label`, for example `--baseline-label "prototype5 memory control"`.
This changes metadata only. The default label is `RAMOPT b00650d`; a custom
label clears the expected-baseline-commit hint rather than attributing that
binary to `b00650d`. Executable SHA-256 and optional source provenance remain
authoritative; timing, comparisons and cleanup are unchanged.

`summary.json` contains paired wall/CPU/RSS ratios for each worker count, actual
RSS differences, and within-version wall scaling relative to the lowest tested
worker count using the same rounds. The configurable allowance defaults to
**20,000,000 bytes per compute worker** (`--extra-cache-mb-per-worker 20`), so the
candidate RSS limit for a paired T-worker run is:

```text
baseline observed peak RSS + T * 20,000,000 bytes
```

Median, maximum and every paired RSS difference are reported explicitly; this
is not a 1.05x allowance and is not a bound on other inputs. Review memory and
time together. Time fields report whether the median paired wall ratio and its
descriptive bootstrap upper limit are at most 1.0; no slowdown tolerance is
silently introduced. Five pairs and one host cannot prove a hardware-general
absence of regression, especially for short or noisy runs.

Exit 0 means all requested runs and complete output comparisons finished.
It **does not mean performance passed**: inspect the time and memory fields.
A command failure, output mismatch or interruption stops the run, preserves
partial evidence and sets `complete: false` / `all_requested_pairs_verified:
false`; no successful claim is emitted for a partial or unequal run. Use a new
`--out` directory for each invocation; automatic resumption is intentionally
unsupported.

## Optional Linux diagnostics

If the server still slows down above 16 workers, record topology and first
repeat on this one allocation. In a separate diagnostic experiment, wrap the
same application command saved in `results.json`, for example:

```sh
perf stat -e task-clock,cycles,instructions,context-switches,cpu-migrations,cache-misses \
  -- /builds/candidate/ExpansionHunter <the saved application arguments>
```

Use a fresh output prefix, the same CPU mask and the matching cache. Counter
availability/permissions vary by CPU and kernel; inspect `perf list` and record
unsupported events. Use `perf record`/annotation, `numastat -p` or suitable
memory-controller counters only as separate diagnostic runs, not mixed into
the paired timing matrix. Flat RSS alone says nothing about memory bandwidth,
NUMA traffic or cache-coherence contention. This tool provides no invented
counter values or diagnosis from worker count alone.

The runner itself has a tiny mock check, which creates no real genomic workload:

```sh
python3 tools/parallel_validation/test_runner.py
```

It checks alternating schedules, archive separation, complete reordered BAM
equality, duplicate-loss rejection, failure persistence and safe output cleanup.

## English Markdown and HTML reports

`report.py` reads only saved JSON metadata. It does not run applications, open
BAM/reference/catalog files, repeat comparisons or infer a successful outcome
from incomplete data. Each input label remains a separate comparison with its
own binary SHA-256; an owner-refcount ablation must not be pooled with a normal
baseline comparison.

```sh
python3 tools/parallel_validation/report.py \
  --result sparse20m=/results/sparse20m/results.json \
  --result owner-ablation=/results/owner-ablation/results.json \
  --label sparse20m="100k loci, 20M records" \
  --records sparse20m=20000000 --loci sparse20m=100000 \
  --details --out /reports/parallel-comparison-001
```

The new output directory contains `report.md`, standalone `index.html` and a
relative `data/` link for each compact case snapshot. Tables show requested and
completed pairs, observed slowdowns, wall confidence intervals when the saved
ratio set matches, median/maximum RSS, paired RSS differences and the recorded
decimal-MB-per-worker allowance. Warmups remain excluded from statistics but
must be verified; missing/duplicate rounds, failures and output mismatches are
explicit. Report generation itself is not a performance acceptance decision.
`--details` adds process CPU and logged read-phase tables; worker elapsed counters
are not treated as CPU time.

Counts can come from `--manifest LABEL=manifest.json`, an embedded manifest, or
`manifest.json` adjacent to the recorded reads path. `--records` and `--loci`
override those counts; `--label` changes only the display name. Missing counts
stay unknown rather than triggering a scan of large inputs. Default compact
snapshots omit absolute binary/source paths and logs. Use
`--include-raw-metadata` only when appropriate: it also copies full result JSON
and available manifest files, which can contain private paths. Existing report
directories are rejected so stale raw metadata cannot leak into a later report.
