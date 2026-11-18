# Installation

NuExpansionHunter is distributed as source; no prebuilt binaries are provided.

Tested platforms and compilers:

* Linux x86-64 (Ubuntu) with GCC 11.5, 12.5, 13.4 and 15.2, using the commands
  below unchanged; the unit tests pass with each. The performance figures in
  [Performance](07_Performance.md) use GCC 15.2.
* Linux ARM64 (ARMv8-A): built with a cross compiler and tested under emulation.
  The NEON code has not been run on ARM hardware.

x86-64 builds need a CPU with SSE4.2 (x86-64-v2), which graph-tools requires.
macOS, Clang and other compilers are not tested.

## Building from source

Prerequisites:

* [GCC](https://gcc.gnu.org/) 11 or later (C++11).
* [CMake](https://cmake.org/) version 3.13.0 or above.
* Git: the superbuild applies the htslib patches with `git apply`.
* Development libraries for zlib, bzip2, xz and curl. On Ubuntu:
  `apt install build-essential git cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev`.

If the above prerequisites are satisfied, you are ready to
build the program. Note that during the build procedure, cmake will
attempt to download and install `abseil`, `boost`, `googletest`, `htslib`,
and `spdlog` so an active internet connection is required. Assuming
that the source code is contained in a directory `ExpansionHunter/`,
the build procedure can be initiated as follows:

```bash
cd ExpansionHunter
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

If all the above steps were successful, the `ExpansionHunter` executable can be found in:

    build/install/bin/ExpansionHunter

`make` also runs the unit tests. To check the build on the example data, run
`example/run.sh` from the source directory; it compares the result with the
stored ExpansionHunter v5.0.0 output (see `example/run.sh` for options).

If libcurl is installed outside the system paths, add
`-DCURL_INCLUDE_DIR=<directory containing curl/curl.h> -DCURL_LIBRARY=<path to libcurl.so>`
to the `cmake` command; the superbuild passes them to htslib and ExpansionHunter.

## CPU target

The `EH_TARGET` CMake option selects the instructions the compiler may use:

* `portable` (default): x86-64-v2 on x86-64, ARMv8-A on ARM64. The SIMD kernels
  still use AVX2 and AVX-512BW at run time on CPUs that have them.
* `native`: the CPU of the build machine.
* a CPU name accepted by `-march` on x86-64 or `-mcpu` on ARM64, for example
  `x86-64-v3`, `znver4`, `sapphirerapids` or `neoverse-v2`.

For example, `cmake .. -DCMAKE_BUILD_TYPE=Release -DEH_TARGET=native`. A binary
built for `native` or a named CPU runs only on CPUs with the same instructions.
The output does not depend on the target. When configuring `ehunter/`
directly, a `-march` or `-mcpu` flag in `CMAKE_CXX_FLAGS` takes the place of
`EH_TARGET`.

## Optional libraries and htslib

Three optional libraries make NuExpansionHunter faster. Install their
development packages before running `cmake` to use them. On Ubuntu:
`apt install libdeflate-dev libisal-dev libjemalloc-dev`. Without them, the
GangSTR chr21-22 run of [Performance](07_Performance.md) at 16 cores took 30%
longer and used 5% more memory.

| Library | Ubuntu package | Used for | Option |
| --- | --- | --- | --- |
| libdeflate | `libdeflate-dev` | htslib BGZF codec: input decompression, seeking-mode output | `EH_WITH_LIBDEFLATE` |
| ISA-L | `libisal-dev` | BAM compression in streaming mode | `EH_WITH_ISAL` |
| jemalloc | `libjemalloc-dev` | memory allocator; returns freed catalog memory after compaction | `EH_WITH_JEMALLOC` |

Each option takes `AUTO` (default), `ON` or `OFF`. `AUTO` uses a library when
its header and library are both found; a runtime-only package (for example
`libjemalloc2` without `libjemalloc-dev`) is not used. `ON` stops configuration
when the library is missing. `OFF` never uses it. The configure step prints one
line with the result, for example
`Optional libraries: libdeflate ON, ISA-L ON, jemalloc ON`.
Libraries outside the system paths can be given with `LIBDEFLATE_INCLUDE_DIR`
and `LIBDEFLATE_LIBRARY`, `EH_ISAL_INCLUDE_DIR` and `EH_ISAL_LIBRARY`, and
`EH_JEMALLOC_LIBRARY`. When configuring `ehunter/` directly, set
`EH_HTS_CODEC_LIBRARY` to the libdeflate library linked into htslib.
The superbuild reads `CPPFLAGS` from the environment when building htslib.

None of these libraries changes VCF, JSON or decoded BAM records.
htslib is pinned to 1.10.2: the in-place BAM record reader depends on its BGZF
and iterator layout. The superbuild applies three patches:

* `htslib-1.10.2-setup-rollback.patch`: rolls back partial thread-pool and BGZF setup on allocation or thread-creation failure.
* `htslib-1.10.2-runtime-lifecycle.patch`: drains and releases asynchronous jobs safely and propagates write, flush and close failures.
* `htslib-1.10.2-reuse-libdeflate.patch`: reuses per-thread libdeflate codec workspaces across BGZF blocks and releases them at thread exit.

## Optimized build

Run these commands from the source directory. Set `$BUILD` to an absolute build
path and `$CPU` to the target CPU name for `-march=<cpu>`. Build the dependencies
once with the superbuild, then configure `ehunter/` directly with
`-DCMAKE_PREFIX_PATH=<build>/install`: the superbuild does not forward
`CMAKE_CXX_FLAGS`.

```bash
cmake -S . -B "$BUILD" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD" --target Boost spdlog htslib googletest abseil --parallel
cmake -S ehunter -B "$BUILD/optimized" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$BUILD/install" \
  -DCMAKE_CXX_FLAGS="-march=$CPU -ffp-contract=off" -DEH_VERIFY=
cmake --build "$BUILD/optimized" --target ExpansionHunter --parallel
```

If htslib uses libdeflate, add `-DEH_HTS_CODEC_LIBRARY=<libdeflate library path>`
to the direct configuration command. Supply the ISA-L and jemalloc variables
there as needed. Keep the same compiler and optional libraries for every direct
configuration. The executable is `$BUILD/optimized/ExpansionHunter`.
The CMake files already add `-ffp-contract=off`; keep it when selecting CPU flags
and do not add `-ffast-math`. A `-march` binary runs only on CPUs with those
instructions.

## Profile-guided optimization (PGO)

Use the dependencies and CPU setting above. Set `$PGO` to a fresh absolute
profile directory. Both compiler passes must use the same build directory
(`$BUILD/optimized`) so the profile counts match the object paths.

1. Configure and build the instrumented executable:

   ```bash
   mkdir -p "$PGO/counts"
   cmake -S ehunter -B "$BUILD/optimized" \
     -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$BUILD/install" \
     -DCMAKE_CXX_FLAGS="-march=$CPU -ffp-contract=off -fprofile-generate=$PGO/counts -fprofile-update=atomic" \
     -DEH_VERIFY=
   cmake --build "$BUILD/optimized" --target ExpansionHunter --parallel
   ```

2. Train with three streaming runs, in order, accumulating counts in the same
   `$PGO/counts` directory. Use a repeat catalog on two chromosomes, a dense
   repeat catalog on one of those chromosomes, and a mixed catalog including
   satellite-pileup loci. For each run, set `$READS`, `$REFERENCE` and `$CATALOG`
   to the matching inputs and `$RUN` to `1`, `2` or `3`. Execute this command
   once for each input set without clearing the counts directory:

   ```bash
   "$BUILD/optimized/ExpansionHunter" \
     --reads "$READS" --reference "$REFERENCE" --variant-catalog "$CATALOG" \
     --output-prefix "$PGO/train-$RUN" --analysis-mode streaming \
     --sex male --aligner dag-aligner --region-extension-length 1000 \
     --log-level info --threads 40 --decode-threads 8 --encode-threads 4 \
     --bamlet-queue-mb 256
   ```

3. Reconfigure the same build directory to consume the accumulated counts and
   build the optimized executable:

   ```bash
   cmake -S ehunter -B "$BUILD/optimized" \
     -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$BUILD/install" \
     -DCMAKE_CXX_FLAGS="-march=$CPU -ffp-contract=off -fprofile-use=$PGO/counts -fprofile-correction -fprofile-partial-training" \
     -DEH_VERIFY=
   cmake --build "$BUILD/optimized" --target ExpansionHunter UnitTests --parallel
   "$BUILD/optimized/UnitTests"
   ```

## Verification builds

The `EH_VERIFY` CMake option selects verification implementations and checks.
Release builds leave it empty (`-DEH_VERIFY=`). The `install` target refuses a
verification build. See the [verification documentation](../ehunter/verification/README.md)
for supported settings and their contracts.

## Binary distribution

A binary that links htslib, libdeflate or ISA-L statically must ship their
notices; [COPYRIGHT.txt](../COPYRIGHT.txt) contains them.
