# Native x86 streaming reproduction

This recipe measures the complete streaming process on a native Linux AVX2
host. It compares the published SIMD implementation **with the same mapped
`Graph::hasEdge` fix** against exact budget DP. Thus its percentage measures the
DP change after the graph ownership contention is removed from both variants.
The standalone correctness phase also checks actual AVX2 execution; it does
not provide the streaming timing result.

Use clean candidate and published SIMD checkouts, Python 3, and Ubuntu 24.04:

```sh
sudo apt-get update
sudo apt-get install --no-install-recommends -y \
  clang cmake ninja-build pkg-config libboost-dev \
  libboost-program-options-dev libboost-filesystem-dev libboost-system-dev \
  libgtest-dev libspdlog-dev libabsl-dev libhts-dev \
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev

python3 candidate/tools/dp_budget_validation/run-native.py \
  --baseline baseline --candidate candidate --out native-results \
  --cxx clang++ --jobs 2 --require-avx2 --skip-benchmarks

python3 candidate/tools/dp_budget_validation/native-streaming.py \
  --candidate candidate --correctness-results native-results/results.json \
  --out streaming-results --cc clang --cxx clang++ --jobs 2 --pairs 9
```

The baseline checkout must be commit
`7fa386484fa7c6f783dbe98239ad77fd097d057b`. Both output directories must be new
and outside the source checkouts. Allow 45 minutes for dependency installation,
standalone correctness, two complete builds, fixture generation and measurement
on a hosted runner. No benchmark processes overlap with compilation or with one
another. The script is a reproduction recipe; a successful run and its recorded
JSON are required before claiming a native x86 speedup.

## Build comparability

Both full applications use the same compiler, flags and system libraries.
Ubuntu's Boost 1.83 Math headers emit a preprocessor warning that future releases
will require C++14. Both builds use `-Wno-error=cpp` so the warning remains
visible without becoming an error under the repository's `-Werror`; the C++11
language standard and application source remain unchanged. Clang documents
[`-Wcpp` as a synonym for `-W#warnings`](https://releases.llvm.org/18.1.8/tools/clang/docs/DiagnosticsReference.html#wcpp).
Configure `ehunter/` directly instead of the root dependency superbuild.
Clang avoids the GraphTools GNU CMake path that adds global `-msse4.2` on x86;
the runner also checks generated compile commands for global SSE/AVX flags.
The AVX2 fallback retains function-specific dispatch, with no new SSE backend.
System shared `libhts.so` and `libcurl.so` are selected explicitly to avoid
accidentally selecting a static library with unresolved transitive dependencies.
The matching system CURL include directory is selected too, and CURL's config
package shortcut is disabled so a preinstalled `/usr/local` package cannot
override that choice. Ubuntu's [HTS development package](https://packages.ubuntu.com/en/noble/libhts-dev)
depends on the GnuTLS CURL development package used in the install command.
Package versions, compiler output, compile commands, source and binary hashes,
and dynamic dependencies are recorded. No application source is changed to
accommodate the system libraries.

The compiler flags suppress unknown-warning-option diagnostics so older Clang
versions can ignore the newer deprecated-literal-operator warning setting;
the repository's other warning-as-error checks remain active. The direct
application build declares CMake 3.13 and GraphTools declares 3.5, meeting
[CMake 4's minimum supported policy version](https://cmake.org/cmake/help/latest/command/cmake_minimum_required.html).
The older bundled GoogleTest download recipe is disabled by `BUILD_TESTS=OFF`;
system GTest is used. No policy-version override is needed for this recipe.

Ubuntu dependencies can differ from those in the ARM experiment. The Linux
comparison remains fair within this runner because both variants use identical
dependencies; do not compare Linux absolute seconds directly to ARM seconds.
The HTS decode and encode thread pools are disabled in this recipe.

## Workload and timing

The fixture is the same covered10k workload as the ARM run: 10,000 independent
synthetic loci, eight read pairs per locus, 160,000 reads of 150 bases, exactly
one substitution per read, mixed qualities and 10,000 missing-quality records.
The catalog, SAM, FASTA and FASTA-index SHA-256 hashes must match the ARM
fixture. Compressed BAM bytes need not match across HTS/zlib versions; a
separate parser validates all decoded input records against the reference.

Each application uses `min(7, max(1, logical_cpu_count - 1))` compute workers
and one preprocessing worker by default. The worker count is recorded and can
be explicitly selected with `--compute`. A small hosted runner may therefore
use fewer than the seven compute workers in the ARM experiment.

Each variant gets its own prepared catalog cache because the cache schema
incorporates source hashes. Cache preparation is recorded separately and is
outside the timing. After one discarded warmup pair, nine alternating AB/BA
pairs measure startup, input, alignment, genotyping and output through process
exit. Hashing, output comparisons and bootstrap statistics happen outside the
measured child processes.

Every warmup and measured pair must have exact VCF/JSON bytes and identical BAM
headers plus the complete BAM record-byte multiset, including duplicate
multiplicity. Only BAM record ordering may differ. The report includes paired
wall and CPU ratios, bootstrap intervals, and median/maximum process peak RSS.
The original RSS gates are retained: noise allowance is the larger of 256 KiB
and 1% of baseline RSS. There is no extra preprocessing-cache allowance for DP;
both variants have the same preprocessing configuration. A completed run is
not an automatic performance or memory acceptance decision: read each gate.

## CI artifact selection

Upload JSON and text evidence only; the complete generated fixture and BAM
outputs remain temporary runner files. Suggested artifact patterns:

```text
native-results/*.json
native-results/*.log
native-results/correctness/*.json
native-results/correctness/*.stderr
streaming-results/*.json
streaming-results/*.log
streaming-results/benchmark/results.json
streaming-results/benchmark/prepare-*-measure.json
streaming-results/benchmark/prepare-*-measure.stderr
streaming-results/benchmark/pair-*/measure.json
streaming-results/benchmark/pair-*/measure.stderr
```

Do not use `streaming-results/**/*.json`: that would also collect per-locus
genotype JSON outputs and add unnecessary artifact volume. The top-level
streaming report is `streaming-results/results.json`; full paired evidence is
`streaming-results/benchmark/results.json`. Source-helper provenance is in
`streaming-helper-provenance.json`.
