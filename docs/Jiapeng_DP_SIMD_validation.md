# Affine insertion DP: NEON and AVX2

Based on `Jiapeng-multithread@91a359c49a9f50afced3bfea102604f2aef480d3`.

The horizontal affine insertion recurrence is evaluated with a max-plus prefix
scan: eight cells at a time on ARM64 NEON, or sixteen cells on x86-64 AVX2.
The AVX2 helper has its own target attribute and runtime CPU detection; a global
`-mavx2` build flag is unnecessary. There is no separate SSE optimization path.

The fast paths require the default gap scores (`-8`, `-2`). Per-block lower bounds
of `INT16_MIN + 24` (NEON) and `INT16_MIN + 40` (AVX2) protect every original
scalar candidate from narrowing overflow, including dominated gap-opening terms.
Other scoring parameters, unsafe blocks, tails and unsupported CPUs retain the
original scalar recurrence. Scores, insertion traceback bits and tie order are
preserved. No additional heap buffer is introduced; the x86 matrix object remains
744 bytes in the tested build. Asymptotic time and space are unchanged.

## Measured ARM performance

Each result uses one discarded warmup pair followed by nine alternating AB/BA
pairs. Reductions use the median of paired candidate/baseline ratios; they are
not ratios of the two marginal medians. CPU time is user plus system time.

| Synthetic workload | Wall-time reduction | CPU-time reduction | Faster pairs | Median peak RSS increase |
| --- | ---: | ---: | ---: | ---: |
| 10,000 independent covered loci, 7 compute + 1 preprocessing worker, cached | 4.43% | 5.68% | 9/9 | 160 KiB |
| Repeated repository example, 1 compute worker, seeking, uncached | 20.88% | 19.74% | 9/9 | 96 KiB |
| Repeated repository example, 1 compute worker, streaming, uncached | 21.51% | 20.85% | 9/9 | 64 KiB |

The first fixture has 160,000 150-base reads, each with one substitution; it
contains the reference repeat allele only. The repeated example uses overlapping
copies of the same catalog coordinates. Neither represents real WGS throughput.
Local DP results were mixed: representative cases improved by about 9-37%, a
dense-branch case was inconclusive, a three-base query case regressed by 1.6%,
and a traceback-heavy case regressed by 3.0%. The latter cause is unconfirmed.

## Correctness and limitations

- Full DP differential suite: 1,578,916 cases, including complete observable
  scores, traceback predicates and selected complete ordered CIGAR paths.
- Direct insertion suite: 3,077,207 cases / 58,515,312 cells, compared against an
  independent scalar oracle and complete physical score/trace storage.
- ARM native and translated x86 auto, forced scalar and forced AVX2 builds
  produced identical canonical bytes. A separately instrumented x86 run executed
  2,444,621 AVX2 blocks; an uninstrumented forced-AVX2 release build also passed.
- Direct insertion ASan/UBSan checks passed, with container checking enabled and
  leak checking disabled. Focused boundary tests also passed. The general DP has
  an inherited query-padding access beyond vector logical size in
  `computeAlignPenalties`; this is not a full-application sanitizer-clean claim.
- GraphTools: 30/30 test programs passed after enabling `BUILD_TESTS`.
  Application: 139/140 tests passed; both baseline and candidate have the same
  `PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime` failure.
- End-to-end VCF and JSON matched complete bytes. Multithread BAM headers and
  full record-byte multisets matched, including duplicate counts; only record
  order was ignored. The one-worker examples also matched compressed BAM bytes.

All performance measurements were native ARM64 macOS. x86 correctness used
translation on that host; native Intel/AMD Linux performance and positive
runtime-dispatch behavior remain to be measured on the target server. Test-only
forced AVX2 bypasses CPU detection and must never be used on unsupported hardware.

## Reproduce the direct insertion checks

Prepare `BASE` at the baseline commit above and `NEW` at this branch. Run from
the repository root, using a new output directory:

```sh
python3 tools/dp_simd_validation/run.py \
  --baseline BASE --candidate NEW --out /tmp/eh-dp-simd-check
```

Add `--boost-include PATH` if Boost is not on the compiler search path, and
`--cxx g++` or `--cxx clang++` to choose the compiler. The runner builds native
release and ASan/UBSan probes and performs complete byte comparisons. On native
AVX2 hardware, add `--define GRAPHTOOLS_DP_TEST_INSTRUMENT_AVX2` and verify that
the candidate's `avx2_blocks` count is positive. Production builds should not
define any `GRAPHTOOLS_DP_TEST_*` macro.

The proposed exact cost-budget DP is a separate algorithm design. It is not
implemented or included in this branch, and none of the measurements above are
attributed to it.
