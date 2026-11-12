Exact-output performance validation

These probes compare two separately built source trees with the same compiler,
Release flags and dependency versions. The reference used for the accompanying
report is commit 9b29510b63ddf39256f2649a1ca0b408b4fd6af9 on Jiapeng-OPT.
See validation-report.txt for measured results and known baseline test failures.

measure.py and end_to_end.py require macOS and Python 3. They use posix_spawn
and wait4: macOS ru_maxrss is in bytes. Do not use their RSS numbers on Linux
without converting Linux's KiB unit. Run benchmarks sequentially with no builds,
tests or other benchmark processes running concurrently. They do not drop caches.

Build both full programs using the repository's build instructions, in separate
build directories. Build each probe against that same variant's production
headers/library. Examples below assume the shell variables SOURCE, BUILD, DEPS,
PROBES and OUT have been set to appropriate absolute directories. SOURCE is the
repository root, BUILD contains libgraphtools.a below its thirdparty directory,
DEPS is the installation prefix containing Boost headers, and PROBES is this
directory. Repeat for baseline and candidate, changing OUT.

  mkdir -p "$OUT"
  c++ -std=c++11 -O3 -DNDEBUG -Wno-deprecated-builtins \
    -I"$SOURCE/ehunter/thirdparty/graph-tools-master-f421f4c/include" \
    -isystem "$DEPS/include" "$PROBES/dp-probe.cpp" -o "$OUT/dp-probe"

  c++ -std=c++11 -O3 -DNDEBUG \
    -Wno-error=deprecated-literal-operator -Wno-error=deprecated-declarations \
    -I"$SOURCE/ehunter" \
    -I"$SOURCE/ehunter/thirdparty/graph-tools-master-f421f4c/include" \
    -isystem "$DEPS/include" "$PROBES/string-probe.cpp" \
    "$SOURCE/ehunter/alignment/OrientationPredictor.cpp" \
    "$BUILD/thirdparty/graph-tools-master-f421f4c/libgraphtools.a" \
    -o "$OUT/string-probe"

For an additional dependency prefix, add another -isystem PREFIX/include.
The string oracle uses the frozen original suffix aligner with renamed classes;
the orientation oracle materializes substrings and both strands independently.

Correctness (replace executable paths):
  /baseline/string-probe check > baseline-string.json
  /candidate/string-probe check > candidate-string.json
  cmp baseline-string.json candidate-string.json
  /baseline/dp-probe verify --exhaustive --dump > baseline-dp.txt
  /candidate/dp-probe verify --exhaustive --dump > candidate-dp.txt
  cmp baseline-dp.txt candidate-dp.txt

Local timing examples (one warmup pair, then nine alternating AB/BA pairs):
  python3 "$PROBES/measure.py" --baseline /baseline/string-probe \
    --candidate /candidate/string-probe --out results/orientation \
    -- orientation mixed 6000
  python3 "$PROBES/measure.py" --baseline /baseline/dp-probe \
    --candidate /candidate/dp-probe --out results/dp \
    -- bench linear 150 160 60000
  python3 "$PROBES/measure.py" --baseline /baseline/dp-probe \
    --candidate /candidate/dp-probe --out results/memory -- memory 1024 25

The memory case stresses vector capacity history; its short duration is not
speed evidence. Increase repetitions when either timing median is below 0.5 s.
Use the other workloads in results.json to check representative regressions.

End-to-end (replace program/input paths):
  python3 "$PROBES/end_to_end.py" --baseline /baseline/ExpansionHunter \
    --candidate /candidate/ExpansionHunter --input /repo/example/input \
    --out results/e2e --copies 32 --pairs 9

The runner checks VCF, JSON and compressed BAM bytes, including record order.
It runs seeking and streaming, threads=1, encode-threads=0; streaming also uses
decode-threads=0. --copies 32 repeats the SAME genomic locus with distinct IDs;
it is not WGS, 32 independent regions, or 32 samples. --copies 1 --pairs 1 is
only a smoke check. --cached creates a separate compatible archive per binary
and excludes cache preparation from sample timing. Source fingerprints change
with the production patch, so old prepared catalog caches must be rebuilt.
Startup timestamps may differ in stdout; E2E correctness uses outputs_identical,
whereas the microbenchmark requires complete stdout equality.

Admission: at least nine pairs, both medians >=0.5 s, paired median wall ratio
<0.98, fixed-seed 10,000-resample bootstrap 95% upper bound <1. Inspect other
workloads for material regressions. Both median and maximum of peak RSS must
be <= their respective baseline value + max(256 KiB, 1%). This tolerance is an
empirical noise band, not a universal proof of zero additional allocation.

Default ASan container checks already fail in baseline computeAlignPenalties
because padded query reads exceed vector size while staying in capacity. The
limited sanitizer run in the report disables container-overflow and leak
detection; it must not be described as a clean default sanitizer pass.
