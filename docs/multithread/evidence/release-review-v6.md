# Frozen v6 release source review

Reviewed on 2026-11-14 while the primary benchmark was running. Branch:
`Jiapeng-multithread`; base commit:
`b00650d297754f68d78778a441e2fcbf5096250b`.

No correctness blocker was identified in this review. No product source was
edited, and no build, test, or benchmark was run. Benchmark acceptance and the
final performance/RSS report remain separate from this source review.

## Frozen-source verification

- All 35 current files listed in `prototype6-manifest.json` match their recorded
  SHA-256 hashes. This includes all 23 changed product/build files under
  `ehunter/`, excluding tests. The union of tracked changes and untracked product
  files contains no additional product file missing from that manifest.
- Comparing the v5 and v6 source manifests identifies exactly one product change:
  `ehunter/sample/HtsStreamingSampleAnalysis.cpp`. The three other changed
  manifest entries are the parallel-validation runner, its tests, and its README.
- The saved before/after source copies match the corresponding v5/v6 manifest
  hashes. Their diff contains one hunk: release each successfully consumed item's
  alignments, read, mate, and remaining shared input before consuming the next
  item. No additional product change is hidden in the v5-to-v6 delta.
- Reviewed pipeline SHA-256:
  `605a430d702c511090dc1d65d9edcf24d030f95a0b9a0bcfe47c0e93992c233f`.
  Isolated patch SHA-256:
  `6dbb7ffb7cb0c03e27765248816e44f902fd60c1a2e213d433889877d8d18154`.

## Early-release lifetime and ordering

The release at `HtsStreamingSampleAnalysis.cpp:293` follows a successful
`consumeMates` call. BAM records already own their serialized name, sequence,
quality, and graph-alignment tag. Retained repeat evidence contains `StrAlign`
values; locus and variant statistics retain values/counters. These consumers do
not retain references to the released `Read` or `GraphAlignment` objects.

The chunk's alignment computation has completed before reduction begins. The
change does not alter per-locus replay order, BAM batch finish, successor
publication, or admission-credit release. If consumption throws, the failing
item and subsequent items remain owned for unwinding; previously consumed items
have already produced owned evidence/output. If batch submission throws, the
existing failure and cleanup path still applies. No new dangling reference,
cross-thread mutation, or successful-output ordering change was found.

## Existing validation evidence inspected

- `unit8.log`: 239 tests ran; 238 passed. The single failure is
  `PreparedCatalog.WholeCatalogRoundTripAndIndependentHandleLifetime`; one other
  test remains disabled. The failure compares graph JSON edge-label arrays in
  different orders (`["b","a"]` versus `["a","b"]`). The same failure appears in
  `docs/ramopt/evidence/ctest-baseline.log`; it is pre-existing and unresolved.
- Passing coverage includes 19 executor tests, 26 pair-source tests, four kernel
  split/replay tests, 11 sink tests, and four BAM batch tests. These cover
  cancellation, capture ownership, publisher assistance, concurrent same-locus
  computation with ordered replay, output ownership, and failure propagation.
- `publisher-latency-focused-tests.log`: 45/45 focused tests passed, including
  both cancellation-polling tests added before v6.
- `graph-tests7.log`: 30/30 graph test programs passed before the v6 pipeline-only
  change; graph sources match the frozen manifests.
- Earlier executor sanitizer evidence reports 19 tests repeated 20 times under
  ASan/UBSan and separately 20 times under TSan, without sanitizer reports. This
  is focused component evidence, not a whole-pipeline sanitizer claim.

The active primary benchmark's final output comparisons, timing, and RSS results
were not inferred from these tests and are not claimed complete here.
