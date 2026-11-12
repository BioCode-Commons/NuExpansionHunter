# Native cost-budget DP validation

This toolkit compares the production `PinnedDagAligner` wrapper on native ARM NEON or x86-64 AVX2. It checks complete ordered scores, paths and CIGAR outputs before timing. It does not require the full ExpansionHunter build or a genomic reference dataset.

The performance control is **published SIMD plus the same mapped `Graph::hasEdge` fix** as the candidate. This separates DP algorithm changes from removal of shared-owner contention. The script copies the baseline to the new output directory, verifies that `Graph.cpp` differs only by the exact reviewed `hasEdge` replacement, and applies that replacement. It also requires identical SIMD headers and verifies every GraphTools include file and the seven compiled source files, allowing changes only in Graph.cpp, PinnedDagAligner.hh and the new budget header. Both complete source manifests are recorded. It never changes either input checkout.

## Run

Prerequisites: Python 3.9+, Git, a C++11 GCC/Clang compiler, Boost headers, and a native ARM64 NEON or x86-64 AVX2 machine. Linux measurements use `wait4` RSS in KiB; macOS measurements normalize its byte units. Windows is not supported. Keep the machine otherwise idle during measurements.

Use the published SIMD commit `7fa386484fa7c6f783dbe98239ad77fd097d057b` for `BASELINE` and the budget-DP candidate checkout for `CANDIDATE`:

```sh
python3 CANDIDATE/tools/dp_budget_validation/run-native.py \
  --baseline BASELINE --candidate CANDIDATE --out RESULTS \
  --cxx c++
```

`RESULTS` must not exist and must be outside `BASELINE`. Boost is found through normal compiler include paths; add `--boost-include /path/to/include` if needed. `--cxx` accepts one executable path/name, not shell flags. The compiler and CPU details, Git revisions, source/probe/compiler/binary SHA256 hashes and exact build commands are saved in `results.json`. The expected published baseline commit is distinguished from the observed Git revision; archive snapshots remain usable and their source hashes are recorded. No global `-mavx2` flag is used: the production dispatch remains active.

On a native AVX2 host, add `--require-avx2`. A compiled runtime feature check must pass. The compiled helper architecture must match the host architecture. On macOS, x86-64 helpers reject Rosetta execution and unknown translation-query errors; the driver also records its translation check. An ARM64 helper matching an ARM64 host does not need the translation syscall: [Rosetta translates x86-64 processes, not ARM64 code](https://developer.apple.com/documentation/apple-silicon/about-the-rosetta-translation-environment). This lets native ARM validation run inside sandboxes that deny `sysctl`. The script never launches an emulator. On AVX2 it separately runs the independent insertion oracle with scalar, automatic-dispatch, and forced-AVX2 builds, compares every canonical output byte, and requires positive executed AVX2 block counters for both AVX2 builds. The force override is used only after the runtime feature check. ARM runs scalar and automatic NEON versions of that oracle.

`--prepare-only` validates and creates the fair control without compiling. `--skip-benchmarks` runs all correctness checks without timing. Build concurrency is controlled by `--jobs` (default 2); measurements are serialized. Reduced `--pairs` or `--target-seconds` values are for smoke tests and must not be presented as formal speed evidence.

## What is checked and measured

- 2,663 `Graph::hasEdge` checks comparing the published SIMD implementation directly against the candidate. Complete stdout, including exception messages, must match. Coverage includes owned/mapped membership, self-loops, unsorted and repeated successors, empty adjacency, invalid node IDs, malformed slice bounds/alignment, retained Graph and escaped-view ownership, and four-thread immutable reads. This separate correctness probe uses `-pthread`; the timed binaries retain their original flags.
- 9,380 complete pinned-wrapper cases, including owned/mapped graphs, prefix/suffix alignment, IUPAC, gaps, numeric guards, reuse, duplicate paths and more than ten tied paths.
- 28 explicit activation witnesses: 24 must succeed through sparse DP, including branch-spanning gaps and top-ten ties; four must use the specified fallback gates. Instrumented builds are excluded from timing.
- 24 low-error workloads × 64 calls: linear/join/repeat graphs, 150/250-base queries, and exactly 0/1/2/8 substitutions. Release control, release candidate and instrumented candidate must produce identical complete ordered output.
- Independent insertion recurrence/oracle checks, including narrowing, guard boundaries and packed traceback. Large binary oracle streams are byte-compared locally, then removed after their hashes are recorded.
- Nine ordinary workloads, including short queries, a long graph node, a short target, mixed direction, and difficult reads. These negative cases remain in the report even if they regress.
- Three selected one-error workloads at length 150: linear, join and repeat.

Each timed workload is calibrated toward 0.7 seconds for the faster binary, followed by one discarded warmup pair and nine alternating paired runs. Every timed stdout stream is byte-compared. Timed stdout contains result counts and a checksum over every returned operation; the separate preflight retains the full canonical alignment output. Reports include all pairs, wall time, CPU time, peak RSS and paired bootstrap intervals. A negative speedup is a regression. The toolkit does not fail merely because the candidate is slower, and does not discard those cases.

These are full-wrapper **single-process microbenchmarks**, including graph extraction/unrolling, DP, traceback and returned alignment parsing. They are not a multithreaded, full-application genotyping benchmark. Existing application-level validation remains necessary for production throughput claims. The pruned matrix is private to pinned alignment; its unused local second-best score is not promised to reproduce the generic dense-matrix second-best value.

## CI artifacts and optional workflow

`results.json` is the top-level index. Canonical text outputs, binaries and the generated control checkout remain local. The prepared CI workflow uploads JSON and text logs only; it excludes executable binaries, the copied checkout, and large canonical binary streams.

The workflow is delivered separately as `workflow-artifact/dp-native-avx2.yml`. Publishing tools does not install or execute it. Copy it to `.github/workflows/dp-native-avx2.yml` only with credentials authorized to update GitHub Actions workflows. It uses `contents: read`, Ubuntu 24.04, a 15-minute limit, the pinned SIMD baseline commit, and the current candidate commit. Automatic runs are limited to changes to that workflow on `Jiapeng-DP-SIMD`; manual `workflow_dispatch` is also enabled.

Official Node 24 action pins were verified against their release commit links and action metadata. GitHub retired the Node 20 runtime on [September 23, 2026](https://github.blog/changelog/2026-09-23-node-20-is-no-longer-available-in-github-actions/):

- [actions/checkout v7.0.1](https://github.com/actions/checkout/commit/3d3c42e5aac5ba805825da76410c181273ba90b1)
- [actions/upload-artifact v7.0.2](https://github.com/actions/upload-artifact/commit/cf430e030ddbb5b0abf93d22962f4752f3646cd9)
