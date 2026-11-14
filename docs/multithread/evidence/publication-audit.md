# Publication preparation audit

Read-only audit on 2026-11-14, during the timed stress campaign. Branch:
`Jiapeng-multithread`. No product edits, builds, tests, binary hashes, or large
input scans were performed. The earlier intended-source manifest match was not
repeated. Only this note was written.

No new product correctness or missing-build-registration blocker was found.
Publication should use an explicit file list; the working tree contains
unexplained numbered source copies that must remain outside the commit.

## Local files to preserve and exclude

These paths are untracked. Their creation history is unknown; modification times
are metadata, not proof of origin. Do not stage them or delete them as part of
this review. The reference index is a generated example artifact.

| Path | Modification time (UTC) | Bytes |
| --- | --- | ---: |
| `ehunter/io/BamRecordSink 2.cpp` | 2026-11-14T18:14:28Z | 6109 |
| `ehunter/io/BamRecordSink 2.hh` | 2026-11-14T18:14:28Z | 2122 |
| `ehunter/io/BamletWriter 2.cpp` | 2026-11-14T18:14:28Z | 18701 |
| `ehunter/sample/HtsStreamingSampleAnalysis 2.cpp` | 2026-11-14T14:55:03Z | 19981 |
| `ehunter/sample/HtsStreamingSampleAnalysis 3.cpp` | 2026-11-14T18:28:56Z | 16895 |
| `ehunter/sample/StreamingPairSource 2.cpp` | 2026-11-14T18:12:18Z | 7272 |
| `ehunter/sample/StreamingPairSource 2.hh` | 2026-11-14T18:11:53Z | 4090 |
| `ehunter/tests/BamRecordSinkTest 2.cpp` | 2026-11-14T18:14:28Z | 13887 |
| `ehunter/tests/StreamingPairSourceTest 2.cpp` | 2026-11-14T18:15:49Z | 23238 |
| `example/input/reference.fa.fai` | 2026-11-14T19:10:45.517484Z | 22 |

The seven non-test numbered source/header copies are not named in executable
source lists. However, `ehunter/CMakeLists.txt:134` recursively includes `.cpp`
and header files when computing `PREPARED_SCHEMA_ID`. A future configure/build
in this local tree can therefore produce a different prepared-cache schema.
This is local workspace contamination: a clean clone that excludes these paths
is unaffected. It does not modify the already-built frozen benchmark binary.
No cleanup or rebuild was attempted during timing.

## Required tracked and new files

The archived v6 manifest names 23 intended product/build paths. All remain
present. Six are new files requiring explicit staging: `TaskExecutor.cpp/.hh`,
`BamRecordSink.cpp/.hh`, and `StreamingPairSource.cpp/.hh`. Their source/header
pairs are explicitly included in `ExpansionHunterLib`.

The four new test files (`BamRecordSinkTest.cpp`, `TaskExecutorTest.cpp`,
`StreamingPairSourceTest.cpp`, and `LocusParallelKernelTest.cpp`) are all included
in `UnitTests`. The modified `GraphTest.cpp` and `GappedAlignerTest.cpp` already
have graph test targets. Graph tests require `BUILD_TESTS=ON`, as before.

The Git index was empty at this review. For publication, stage intended tracked
changes and new paths explicitly using the archived manifest as the application
source list. Also include all four files in `tools/parallel_validation/`:
`run.py`, `test_runner.py`, `report.py`, and `README.md`. The report helper was
added after the frozen application manifest; it is documentation tooling, not a
new application-source delta. Avoid blanket `git add .` or `git add ehunter`.

## Reproduction and evidence

The validation README's example options correspond to the current runner and
report CLI definitions. Their required `measure.py` and `compare_outputs.py`
helpers are already tracked under `tools/streaming_validation/`. The inspected
runner/report documentation uses portable example paths and discloses prepared
cache generation outside timing, counterbalanced serial runs, complete output
comparisons including duplicate BAM records, wait4 peak RSS, decimal-MB memory
allowances, and local macOS versus same-allocation server limitations.

Report-stage metadata is explicitly partial while stress runs remain active.
Keep that status until all requested evidence is packaged; do not substitute a
successful primary run for unfinished stress coverage. The staged evidence also
correctly separates historical v4/v5 comparisons from final v6, discloses the
pre-existing PreparedCatalog edge-label-order test failure, and limits sanitizer
claims to instrumented components. Missing historical executor/graph exact
compile commands and contemporaneous source hashes remain documented gaps.

After the campaign, publish the completed English report and its relative data
and evidence links together. Final performance/RSS acceptance is outside this
lightweight preparation review.
