# Exact index boundary validation

`index_extremes_probe.cpp` builds deterministic in-process graph fixtures. It has no HTS dependency or external data download. Compile the same source once against the baseline GraphTools headers/library and once against the candidate headers/library, using identical compiler options, for example:

```sh
c++ -std=c++11 -O3 \
  -I /path/to/dependencies/include \
  -I /path/to/version/ehunter/thirdparty/graph-tools-master-f421f4c/include \
  tools/memory_validation/index_extremes_probe.cpp \
  /path/to/version-build/thirdparty/graph-tools-master-f421f4c/libgraphtools.a \
  -o /path/to/version-index-extremes
```

Run every exact audit first:

```sh
python3 tools/memory_validation/paired_index_extremes.py \
  --baseline /path/to/baseline-index-extremes \
  --candidate /path/to/candidate-index-extremes \
  --out /path/to/results/audits --audit-only
```

The audit compares every index key and its ordered paths, including duplicates, against independent path enumeration. It also compares complete diagnostic encodings between binaries, checks mapped/copy/move lifetime, and checks explicit bucket-boundary queries plus 100,000 deterministic pseudorandom queries per profile. The runner verifies that the candidate actually exercises the intended prefix representation or legacy fallback.

| Profiles | Boundary | Expected candidate representation |
| --- | --- | --- |
| `start2047`, `start2048` | Maximum single-node start position 2,047 vs. 2,048 | Prefix / fallback |
| `node3`, `node4` | Single-node path uses node ID 3 vs. 4 | Prefix / fallback |
| `dense14`, `dense10` | Homopolymer graph with many repeated path occurrences | Fallback |
| `skew14`, `skew10` | More than 32 keys share a prefix bucket | Fallback |
| `entries65535`, `entries65536` | Complete eight-base de Bruijn sequence, trimmed by one key or untrimmed | Fallback; 16-bit / 32-bit hash buckets |
| `complex14`, `complex10` | Four-node branching graph, two repeat loops, one ambiguous flank base | Prefix |
| `ambiguous14`, `ambiguous10` | Ambiguous repeat motif and flank bases; bounded expansion | Prefix |
| `tail14`, `tail10` | Last prefix bucket contains suffix 1; absent suffix 0 must not match padding | Prefix |
| `empty14` | No complete 14-mer exists | Fallback |
| `maximum16` | All-zero and all-one 32-bit keys, repeated path multiplicity | Fallback |

The dense graphs use 1 kb flanks. The largest sequence is 65,543 bases and has at most 65,536 distinct keys. Ambiguity is deliberately bounded to avoid exponential fixture construction.

A focused timing matrix can be run after correctness succeeds:

```sh
python3 tools/memory_validation/paired_index_extremes.py \
  --baseline /path/to/baseline-index-extremes \
  --candidate /path/to/candidate-index-extremes \
  --out /path/to/results/paired --pairs 5 --target-seconds 0.15 \
  --profiles start2047 start2048 node3 node4 dense14 skew14 entries65536 complex14
```

Each operation is calibrated toward 0.15 seconds for the faster binary, followed by one discarded warmup pair and five measured pairs with alternating execution order. Every timed stdout checksum must match. Reports preserve every pair, executable hashes, full-audit hashes, serialized record sizes, and representation selection. A ratio above one is a slowdown; the runner retains it rather than relaxing a gate.

Timing covers query loops only. Index construction, archive writing, query generation, and correctness checks are excluded. Hit queries use at most 4,096 deterministically shuffled keys. Miss queries are verified absent before timing. A complete 65,536-entry eight-base index has no valid absent key, so its miss operation is skipped. The 65,535-entry profile has exactly one valid missing key. The `paths` operation reads the path count and materializes the first path, rather than materializing every occurrence of a repeated key.

These are cache-hot index measurements, not full streaming throughput or process RSS. Record-byte reductions are serialized index storage measurements. Whole-application streaming benchmarks are required to judge end-to-end memory and runtime.

For sanitizer checks, compile both this probe and the candidate `KmerIndex.cpp` with `-fsanitize=address,undefined -fno-omit-frame-pointer`, then run each `audit PROFILE`. State explicitly whether supporting GraphTools objects were also instrumented. Unaligned SIMD loads, final-bucket padding, archive lifetime, and compact/full-path boundary transitions are exercised by these audits.
