# Exact cost-budget pinned DP: streaming validation

This change reduces native ARM64 wall time by **56.77%** and CPU time by **59.93%** against the published SIMD version on the covered10k **streaming** workload. A control containing the same graph-membership fix isolates the DP change: wall time falls **52.41%**, CPU time **53.78%**. Both comparisons use nine serial pairs, all nine faster, with unchanged production outputs under the BAM record-order policy below. These are workload-specific measurements; the full probe table includes regressions.

The comparison baseline is [`Jiapeng-DP-SIMD@7fa3864`](https://github.com/BioCode-Commons/NuExpansionHunter/commit/7fa386484fa7c6f783dbe98239ad77fd097d057b). The portable native reproduction entry point is [`tools/dp_budget_validation/run-native.py`](../tools/dp_budget_validation/run-native.py); compact measured evidence is [`results-arm64.json`](../tools/dp_budget_validation/results-arm64.json).

## Algorithm and exactness

The input is a pinned query fragment of length `m`, an expanded target DAG with `B` base positions and `E` predecessor edges, and the existing base-matching/affine-gap scoring parameters. Seed finding and graph unrolling remain outside this matrix. The objective is the original maximum alignment score for the complete query, together with the original ordered optimal traceback results.

For the default match `+5`, mismatch `-4`, gap opening `-8`, extension `-2`, transform each score at consumed-query coordinate `i` to `C = 5*i - score`:

| Transition | Cost increment |
|---|---:|
| Match | 0 |
| Mismatch | 9 |
| Open insertion | 15 |
| Extend insertion | 7 |
| Open deletion | 10 |
| Extend deletion | 2 |

The legacy boundaries are preserved explicitly: virtual target row `H/E = 15*i`; query-zero column `H/F = 10*d_min`, with `E` unreachable there. An explicit infinity value represents unreachable states.

Budgets start at 9 and double, skipping budgets below a conservative endpoint length-cost bound. For each completed predecessor row, the finite `H` interval contains every finite `E/F` state. A target row therefore obtains all seeds from the predecessor intervals at the same query coordinate (deletion), shifted by one (match), and the legacy boundary. Extending the seed interval to the right by `floor(K/7)` conservatively covers all possible insertion descendants. This interval is intersected with the certified graph-distance corridor `max(0,d_min-floor(K/2)) ... min(m,d_max+floor(K/7))`; unreachable holes inside it remain represented.

All increments are nonnegative. Consequently every prefix of an alignment of cost at most `K` is retained. Once a **complete** budget pass contains a full-query endpoint, a discarded path of cost above `K` cannot beat it; all optimal ties are retained. The original `Aligner` traceback, tie priority, duplicate handling and top-10 behavior are reused. Exhausting the work allowance during a partial pass triggers dense fallback, even if an endpoint was already seen, so later tied endpoints cannot be lost.

The new matrix is restricted to `PinnedDagAligner`. The generic full-matrix implementation is unchanged. The private `secondBestScore` local can differ because higher-cost states need not exist in the sparse matrix; both pinned entry points ignore that value. Best scores, ordered raw optimal CIGARs and public results are the verified contract.

Non-default scoring, `m < 16`, `B < 8`, unsupported graph structure, or failure of `15*m + 10*B + 10 < 32768` use the original SIMD implementation. Failed sparse attempts also fall back. ARM NEON and runtime-dispatched AVX2 remain available, without a new SSE backend or a global `-mavx2` requirement.

## Complexity and memory

The dense recurrence costs `O(m*(B+E))` time and `O(m*B)` storage. A budget pass costs `O(B+E+A_K+R_K)`, where `A_K` is the total dynamic row width and `R_K` is the actual predecessor intersection work. Cumulative charged sparse recurrence work, including partial attempts, is capped at `m*(B+E)/2`. Including dense fallback, the worst-case asymptotic time remains `O(m*(B+E))`; this abstract operation cap is not a runtime guarantee. No general fixed-width band bound is asserted for arbitrary branching DAGs.

Sparse storage is proportional to the stored cells plus row metadata. Allocations are reused, and dense fallback buffers may coexist with sparse buffers. Worst-case retained storage remains `O(m*B)`. Reported RSS changes below are measured values, not a blanket memory-reduction claim.

## Why the graph-membership change is included

The first dynamic-support implementation sped up the DP stage but regressed the seven-worker application. Matching phase instrumentation localized the new cost to path construction. `Path::assertConnected` repeatedly calls mapped `Graph::hasEdge`, whose old implementation created owning `CatalogRange` iterators; these repeatedly copied the catalog's shared owner. Shortening DP exposed this shared-reference work much more strongly in parallel execution. A one-worker diagnostic already showed a 55.06% wall reduction before the graph fix, whereas the nine-pair seven-worker run regressed 17.93%.

The mapped membership check now searches its already-validated raw successor slice directly. The `Graph` retains the mapping throughout the call; no borrowed pointer escapes. Existing node-ID checks, bounds/alignment validation, and public owning range/iterator APIs remain intact. Unsorted and duplicate successor entries retain the same membership semantics. The fair SIMD control applies exactly this same change, so its paired comparison isolates the DP algorithm.

## Native ARM64 streaming measurements

Environment: ARM64, macOS 26.6.2, Apple Clang 21.0.0, release `-O3 -DNDEBUG`. Timed subprocesses run serially with one warmup pair and nine alternating AB/BA measured pairs. Microbench calibration targets 0.7 seconds for the faster binary. No instrumented binary is used for performance claims. CPU is user plus system time.

Reductions are `100*(1-median(candidate/baseline))` over pairs, not ratios of marginal medians. Intervals are descriptive 95% bootstrap intervals for the paired median ratio. Absolute times vary across separate batches; do not subtract their marginal medians to infer the standalone graph-fix effect.

The end-to-end dataset has 10,000 independently covered synthetic loci, 80,000 pairs / 160,000 reads of 150 bases, one substitution per read, 10,000 reads with missing quality, and only the reference 10-CAG repeat allele. Each streaming run uses seven compute workers plus one preprocessing worker and a compatible prepared catalog. Cache preparation is excluded from timing; loading that cache remains part of the timed application process. The end-to-end measurements cover the full ExpansionHunter subprocess: startup, read processing, alignment, genotyping and output writing. Fixture/cache preparation and post-run output comparison occur outside the timed process. It is a DP stress workload, not a representative clinical WGS benchmark.

| Baseline in this paired batch | Baseline median wall | V4 median wall | Paired wall reduction | Paired CPU reduction | Faster pairs |
|---|---:|---:|---:|---:|---:|
| Published SIMD | 5.112 s | 2.192 s | 56.77% | 59.93% | 9/9 |
| SIMD + same graph fix | 6.453 s | 3.147 s | 52.41% | 53.78% | 9/9 |

The percentage is the median paired reduction; it need not equal division of the two displayed marginal medians. The JSON preserves every measured pair and the excluded warmup, with wall, user CPU, system CPU, peak RSS and output-equality results.

| Case | Wall reduction | CPU reduction | Wall ratio 95% interval | Faster pairs | Median RSS change |
|---|---:|---:|---|---:|---:|
| `V4 vs published SIMD` | 56.77% | 59.93% | [0.424, 0.445] | 9/9 | -5840 KiB |
| `V4 vs SIMD + same graph fix` | 52.41% | 53.78% | [0.443, 0.533] | 9/9 | -6160 KiB |

VCF and JSON files in covered10k are byte-identical. BAM comparison covers the complete header and the multiset of complete record bytes, preserving duplicate counts; only record order is allowed to change. The covered10k RSS decrease needs no extra preprocessing-cache allowance.

The supplementary **dense32 streaming** check uses the repository example with 32 repeated identical-coordinate catalog entries, **one** compute worker, **no** prepared cache, and the published SIMD baseline. Median wall time is **3.555 → 2.612 s**, with a paired reduction of **26.54%** and CPU reduction of **26.20%**. All nine pairs are faster. VCF, JSON and the complete compressed BAM bytes match exactly, including order. Median peak RSS increases by **1,064,960 bytes (1,040 KiB)**, exceeding the local no-growth noise allowance; this RSS gate fails. This repeated-example workload is also not WGS.

### Why 26.54% and 56.77% differ

They describe two different streaming experiments: dense32 is a one-worker, uncached, repeated-coordinate example; covered10k has 10,000 independent loci, seven compute workers and a prepared catalog. Amdahl's law makes the overall improvement depend on the fraction of execution spent in accelerated work. Graph/query structure changes sparse-DP work, and caching and worker count change startup costs, shared-reference contention and the relative cost of other stages. The graph-membership fix also targets a specific mapped-catalog code path. The 52.41% algorithm-only result likewise belongs to the covered10k configuration.

Dense32 has not received the phase profiling used for covered10k, so these differences do not support assigning numerical causal shares to DP, I/O, initialization or parallel overhead. The 26.54% measurement does not represent a deterioration of the 56.77% experiment; it is the result for a different workload and configuration. The JSON retains the per-pair measurements, input/binary hashes, and previously collected nonfocus modes for review.

### Diagnostic pinned probes

The following probes measure the entire pinned call: graph extraction/unrolling, DP, ordered traceback and result construction. Both sides contain the graph fix. In the ordinary mixed probes, one in four queries is exact; the others have substitutions spaced 19 bases apart. Short-length and deliberately short-target controls are included.

| Case | Wall reduction | CPU reduction | Wall ratio 95% interval | Faster pairs | Median RSS change |
|---|---:|---:|---|---:|---:|
| `join150` | 59.96% | 60.04% | [0.394, 0.402] | 9/9 | -96 KiB |
| `linear150` | 69.50% | 69.51% | [0.304, 0.306] | 9/9 | -192 KiB |
| `linear250` | 64.37% | 64.32% | [0.349, 0.358] | 9/9 | -384 KiB |
| `long-node150` | 86.49% | 86.34% | [0.127, 0.165] | 9/9 | -1264 KiB |
| `repeat150` | 5.62% | 5.75% | [0.938, 0.951] | 9/9 | +96 KiB |
| `repeat250` | -21.82% | -21.95% | [1.166, 1.288] | 1/9 | +336 KiB |
| `short-target150` | -12.75% | -12.83% | [1.053, 1.185] | 0/9 | +176 KiB |
| `short16` | 22.34% | 22.18% | [0.733, 0.832] | 9/9 | +32 KiB |
| `short3` | -2.06% | -2.13% | [0.938, 1.094] | 4/9 | +32 KiB |

Low-error probes use valid graph walks with exactly the indicated one or two substitutions.

| Case | Wall reduction | CPU reduction | Wall ratio 95% interval | Faster pairs | Median RSS change |
|---|---:|---:|---|---:|---:|
| `lowerror-join150-e1` | 85.47% | 85.48% | [0.135, 0.155] | 9/9 | -256 KiB |
| `lowerror-join250-e2` | 90.95% | 90.75% | [0.086, 0.107] | 9/9 | -496 KiB |
| `lowerror-linear150-e1` | 89.14% | 89.16% | [0.096, 0.150] | 9/9 | -288 KiB |
| `lowerror-linear250-e2` | 90.79% | 90.81% | [0.082, 0.097] | 9/9 | -608 KiB |
| `lowerror-repeat150-e1` | 40.61% | 40.94% | [0.566, 0.618] | 9/9 | -144 KiB |
| `lowerror-repeat250-e2` | 44.26% | 44.12% | [0.538, 0.571] | 9/9 | -272 KiB |

**Regressions are retained:** mixed `repeat250` is 21.82% slower and increases median peak RSS by 336 KiB, exceeding the local RSS noise allowance. `short-target150` is 12.75% slower. `short3` has a 2.06% median slowdown with an interval crossing parity. These cases prevent any claim that every input gets faster or smaller.

## Validation scope

- Final V4: 9,380 native ARM pinned cases and 12,800 native Gapped cases have complete ordered outputs equal to the independent SIMD baseline.
- Final V4: 9,380 pinned cases each pass for translated x86_64 automatic dispatch and forced AVX2. These runs establish correctness, not native Intel/AMD speed.
- The unchanged DP header additionally passes 28 explicit activation/fallback witnesses per tested architecture (24 require sparse success), including indels, branch crossings and more than ten tied paths; 648 focused ASan/UBSan cases; 2,000 irregular-DAG cases; and 2,000 calls reusing one aligner across changing inputs, graphs, sparse success and fallback.
- The final fair low-error preflight compares 24 cases times 64 alignments in full.
- The graph fix passes 2,663 checks per build across three release and two ASan/UBSan builds: owned/mapped membership, duplicate/unsorted/empty adjacency, invalid IDs, malformed ranges, retained ownership, escaped views and immutable concurrent reads.
- An actual-input audit of the final DP header before the membership-only fix compares 160,000 calls and 655,466,727 canonical bytes per side: complete input context, best score, ordered pre-fix raw CIGAR operations and ordered public Path/Alignment outputs. Only the unused local `secondBestScore` is excluded and separately recorded. Both sides produce 160,000 CIGARs / 3,030,841 raw operations and 160,000 public results / 400,000 public operations.
- The final builds pass all 30 GraphTools test suites and 139 of 140 application cases. `PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime` is the same pre-existing failing application test as in the baseline. Focused sanitizer probes do not establish whole-application sanitizer cleanliness.

The earlier V2 topology-corridor implementation regressed covered10k wall time by 32.74%. V3 introduced dynamic finite support but, without the membership fix, still regressed seven-worker wall time by 17.93%. Their evidence remains recorded; the retained implementation is V4.

## Reproduction and provenance

Use the native harness README for compiler, Boost and architecture checks. It constructs a fair SIMD control, compares full outputs before timing, verifies actual AVX2 execution when requested, and runs serial paired probes. Native Intel/AMD performance is **not** included in this ARM report; translated timings must not be substituted for it.

The recorded source hashes and per-pair timings are in `results-arm64.json`. The unchanged SIMD header is `9f0b443bedda6fa778632002142f819b3d27f6b3b4de42c3c92ebca190add1bb`; the exact DP header is `59cac7722f2dbab69f081e5b9a2dcd6cb24416be40a0c2df68decc31eea723cd`; the mapped membership source is `b9152b5bd7217e3ea66475ba1f1b942b12009cdb63f53acc2abc6c58745f6f37`. Evidence paths in the JSON are relative to the budget-DP evidence bundle, with large canonical streams omitted from the compact repository summary.
