#
# Modified by NuExpansionHunter contributors, 2026.
# This file differs from ExpansionHunter v5.0.0: the build includes parallel processing and verification targets, hashes catalog schemas, selects the CPU target, and configures optional libraries and reproducible floating-point evaluation.
# Modifications are licensed under the Apache License, Version 2.0.
#
cmake_minimum_required(VERSION 3.13)
project(ExpansionHunter LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE} CACHE STRING
            "Choose the type of build (default: ${DEFAULT_CMAKE_BUILD_TYPE})" FORCE)
endif ()

message(STATUS "BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.73 REQUIRED COMPONENTS program_options filesystem system)


add_compile_options(-Wall -Werror -Wextra)
# GCC 14 and later raise these diagnostics in inherited code: unused-result in
# ExpansionHunter v5.0.0's graph-tools sources, cast-user-defined in sparsepp,
# array-bounds in v5.0.0's Cigar, stringop-overflow in Boost. Keep them as warnings.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
    add_compile_options(-Wno-error=unused-result -Wno-error=cast-user-defined -Wno-error=array-bounds
                        -Wno-error=stringop-overflow)
endif()
# GCC 13 and earlier report the boost::optional genotypes of ExpansionHunter v5.0.0's
# repeat and small-variant analyzers as maybe-uninitialized. Keep it a warning there.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
    set_source_files_properties(locus/RepeatAnalyzer.cpp locus/SmallVariantAnalyzer.cpp
                                PROPERTIES COMPILE_OPTIONS -Wno-error=maybe-uninitialized)
endif()
# Genotype likelihoods must be evaluated in source order, without fused
# multiply-add, whatever -march is chosen; this keeps them equal to v5.0.0's.
add_compile_options(-ffp-contract=off)

# Verification switches (see verification/README.md). Release builds leave this empty.
set(EH_VERIFY "" CACHE STRING "Verification switches, e.g. SIMD=SCALAR;BUDGET_SEARCH=NEVER; empty for a release build")
set(EH_VERIFY_VALUES_SIMD SCALAR SSE42 AVX2)
set(EH_VERIFY_VALUES_FAST_PATHS OFF)
set(EH_VERIFY_VALUES_BUDGET_SEARCH NEVER ALWAYS)
set(EH_VERIFY_VALUES_ROW_STOP OFF)
set(EH_VERIFY_VALUES_FUSED_ROWS OFF)
set(EH_VERIFY_VALUES_GRAPH_CACHE OFF)
set(EH_VERIFY_VALUES_FOLD_SKIP OFF)
set(EH_VERIFY_VALUES_LAYOUT REFERENCE)
set(EH_VERIFY_VALUES_SHADOW ON)
foreach(EH_VERIFY_SETTING IN LISTS EH_VERIFY)
    if(NOT EH_VERIFY_SETTING MATCHES "^([A-Z_]+)=([A-Z0-9]+)$")
        message(FATAL_ERROR "Unknown EH_VERIFY setting: ${EH_VERIFY_SETTING}")
    endif()
    set(EH_VERIFY_NAME ${CMAKE_MATCH_1})
    set(EH_VERIFY_VALUE ${CMAKE_MATCH_2})
    if(NOT DEFINED EH_VERIFY_VALUES_${EH_VERIFY_NAME})
        message(FATAL_ERROR "Unknown EH_VERIFY setting: ${EH_VERIFY_SETTING}")
    endif()
    if(NOT EH_VERIFY_VALUE IN_LIST EH_VERIFY_VALUES_${EH_VERIFY_NAME})
        message(FATAL_ERROR "EH_VERIFY: ${EH_VERIFY_NAME} must be one of ${EH_VERIFY_VALUES_${EH_VERIFY_NAME}}")
    endif()
    add_compile_definitions(EH_VERIFY_${EH_VERIFY_NAME}=EH_VERIFY_${EH_VERIFY_VALUE})
endforeach()
if(EH_VERIFY)
    message(WARNING "Verification build (EH_VERIFY=${EH_VERIFY}): not for production use")
endif()

# CPU target. The SIMD kernels choose AVX-512BW, AVX2 or SSSE3 at run time, so the
# portable build uses the wider instructions on CPUs that have them.
#   portable  x86-64-v2 on x86-64 (SSE4.2, which graph-tools also requires);
#             the compiler default (ARMv8-A) on ARM64. This is the default.
#   native    the build machine's CPU
#   <name>    a -march (x86-64) or -mcpu (ARM64) value, e.g. x86-64-v3, znver4, neoverse-v2
# A binary built for native or a named CPU runs only on that CPU family or newer.
# A -march or -mcpu flag in CMAKE_CXX_FLAGS takes the place of EH_TARGET.
set(EH_TARGET "portable" CACHE STRING "CPU target: portable, native, or a -march (x86-64) or -mcpu (ARM64) name")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
    set(EH_TARGET_FLAG "-mcpu=")
    set(EH_PORTABLE_TARGET "")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
    set(EH_TARGET_FLAG "-march=")
    set(EH_PORTABLE_TARGET "x86-64-v2")
else()
    set(EH_TARGET_FLAG "-march=")
    set(EH_PORTABLE_TARGET "")
endif()
if(EH_TARGET STREQUAL "portable")
    set(EH_TARGET_VALUE "${EH_PORTABLE_TARGET}")
else()
    set(EH_TARGET_VALUE "${EH_TARGET}")
endif()
if(CMAKE_CXX_FLAGS MATCHES "(^| )-m(arch|cpu)=")
    if(NOT EH_TARGET STREQUAL "portable")
        message(FATAL_ERROR "EH_TARGET=${EH_TARGET} and the CPU flag in CMAKE_CXX_FLAGS conflict; set only one")
    endif()
    message(STATUS "CPU target: from CMAKE_CXX_FLAGS")
elseif(EH_TARGET_VALUE)
    add_compile_options("${EH_TARGET_FLAG}${EH_TARGET_VALUE}")
    message(STATUS "CPU target: EH_TARGET=${EH_TARGET} (${EH_TARGET_FLAG}${EH_TARGET_VALUE})")
else()
    message(STATUS "CPU target: EH_TARGET=${EH_TARGET} (compiler default)")
endif()

enable_testing()

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)
find_package(spdlog REQUIRED)
find_package(absl REQUIRED)

find_library(htslib libhts.a)
find_library(htslib hts)

add_subdirectory(thirdparty/graph-tools-master-f421f4c)

add_library(ExpansionHunterLib
        alignment/AlignmentClassifier.hh alignment/AlignmentClassifier.cpp
        alignment/AlignmentFilters.hh alignment/AlignmentFilters.cpp
        alignment/ClassifierOfAlignmentsToVariant.hh alignment/ClassifierOfAlignmentsToVariant.cpp
        alignment/GraphVariantAlignmentStats.hh alignment/GraphVariantAlignmentStats.cpp
        alignment/GreedyAlignmentIntersector.hh alignment/GreedyAlignmentIntersector.cpp
        alignment/HighQualityBaseRunFinder.hh alignment/HighQualityBaseRunFinder.cpp
        alignment/OperationsOnAlignments.hh alignment/OperationsOnAlignments.cpp
        alignment/OrientationPredictor.hh alignment/OrientationPredictor.cpp
        alignment/SoftclippingAligner.hh alignment/SoftclippingAligner.cpp
        core/Common.hh core/Common.cpp
        core/ConcurrentQueue.hh
        core/CountTable.hh core/CountTable.cpp
        core/GenomicRegion.hh core/GenomicRegion.cpp
        core/HtsHelpers.hh core/HtsHelpers.cpp
        core/Parameters.hh core/Parameters.cpp
        core/Reference.hh core/Reference.cpp
        core/ReferenceContigInfo.hh core/ReferenceContigInfo.cpp
        core/LocusStats.hh core/LocusStats.cpp
        core/LogSum.hh
        core/OrderedBlocks.hh
        core/Read.hh core/Read.cpp
        core/ReadPairs.hh core/ReadPairs.cpp
        core/ReadSupportCalculator.hh core/ReadSupportCalculator.cpp
        core/ThreadPool.hh
        core/TaskExecutor.hh core/TaskExecutor.cpp
        core/VerifyPolicy.hh
        core/WeightedPurityCalculator.hh core/WeightedPurityCalculator.cpp
        genotyping/AlignMatrix.hh genotyping/AlignMatrix.cpp
        genotyping/AlignMatrixFiltering.hh genotyping/AlignMatrixFiltering.cpp
        genotyping/AlleleChecker.hh genotyping/AlleleChecker.cpp
        genotyping/FragLogliks.hh genotyping/FragLogliks.cpp
        genotyping/OneAlleleStrGenotyper.hh genotyping/OneAlleleStrGenotyper.cpp
        genotyping/RepeatGenotype.hh genotyping/RepeatGenotype.cpp
        genotyping/SmallVariantGenotype.hh genotyping/SmallVariantGenotype.cpp
        genotyping/SmallVariantGenotyper.hh genotyping/SmallVariantGenotyper.cpp
        genotyping/StrAlign.hh genotyping/StrAlign.cpp
        genotyping/StrGenotyper.hh genotyping/StrGenotyper.cpp
        genotyping/TwoAlleleStrGenotyper.hh genotyping/TwoAlleleStrGenotyper.cpp
        io/BamletWriter.hh io/BamletWriter.cpp
        io/BamRecordSink.hh io/BamRecordSink.cpp
        io/CatalogLoading.hh io/CatalogLoading.cpp
        io/PreparedCatalog.hh io/PreparedCatalog.cpp
        io/GraphBlueprint.hh io/GraphBlueprint.cpp
        io/JsonWriter.hh io/JsonWriter.cpp
        io/LocusSpecDecoding.hh io/LocusSpecDecoding.cpp
        io/ParameterLoading.hh io/ParameterLoading.cpp
        io/RegionGraph.hh io/RegionGraph.cpp
        io/SampleStats.hh io/SampleStats.cpp
        io/VcfHeader.hh io/VcfHeader.cpp
        io/VcfWriter.hh io/VcfWriter.cpp
        io/VcfWriterHelpers.hh io/VcfWriterHelpers.cpp
        sample/AnalyzerFinder.hh sample/AnalyzerFinder.cpp
        sample/GenomeMask.hh sample/GenomeMask.cpp
        sample/GenomeQueryCollection.hh sample/GenomeQueryCollection.cpp
        sample/HtsFileSeeker.hh sample/HtsFileSeeker.cpp
        sample/HtsFileStreamer.hh sample/HtsFileStreamer.cpp
        sample/HtsSeekingSampleAnalysis.hh sample/HtsSeekingSampleAnalysis.cpp
        sample/HtsStreamingReadPairQueue.hh sample/HtsStreamingReadPairQueue.cpp
        sample/HtsStreamingSampleAnalysis.hh sample/HtsStreamingSampleAnalysis.cpp
        sample/StreamingPairSource.hh sample/StreamingPairSource.cpp
        sample/ParallelPairSource.hh sample/ParallelPairSource.cpp
        sample/IndexBasedDepthEstimate.hh sample/IndexBasedDepthEstimate.cpp
        sample/MateExtractor.hh sample/MateExtractor.cpp
        )


set(CTPL_INCLUDE_DIRS thirdparty/ctpl/ctpl-0.0.2)
target_include_directories(ExpansionHunterLib PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${Boost_INCLUDE_DIRS}
        ${LIBLZMA_INCLUDE_DIRS}
        ${CURL_INCLUDE_DIRS}
        ${CTPL_INCLUDE_DIRS}
        )

# Set static linking of gcc standard libraries to simplify binary distribution
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(STATIC_FLAGS -static-libgcc -static-libstdc++)
endif()


target_link_libraries(ExpansionHunterLib PUBLIC
        ${STATIC_FLAGS}
        graphtools
        ${htslib}
        ${EH_HTS_CODEC_LIBRARY}
        ${Boost_LIBRARIES}
        ${LIBLZMA_LIBRARIES}
        ${CURL_LIBRARIES}
        ZLIB::ZLIB
        BZip2::BZip2
        Threads::Threads
        spdlog::spdlog
        absl::flat_hash_map
        absl::flat_hash_set)

# Prepared catalog semantics follow the source, not just the public EH version.
# Reject an archive produced by a different implementation even when its input
# paths and binary layout happen to match.
file(GLOB_RECURSE PREPARED_SCHEMA_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
list(FILTER PREPARED_SCHEMA_SOURCES EXCLUDE REGEX "/tests/|/verification/|Test\\.cpp$")
list(SORT PREPARED_SCHEMA_SOURCES)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${PREPARED_SCHEMA_SOURCES})
set(PREPARED_SCHEMA_PARTS "")
foreach(PREPARED_SOURCE IN LISTS PREPARED_SCHEMA_SOURCES)
    file(SHA256 "${PREPARED_SOURCE}" PREPARED_SOURCE_HASH)
    string(APPEND PREPARED_SCHEMA_PARTS "${PREPARED_SOURCE_HASH}")
endforeach()
string(SHA256 PREPARED_SCHEMA_ID "${PREPARED_SCHEMA_PARTS}")
set_property(SOURCE io/PreparedCatalog.cpp PROPERTY COMPILE_DEFINITIONS EH_PREPARED_SCHEMA_ID="${PREPARED_SCHEMA_ID}")

# Optional libraries: EH_WITH_ISAL and EH_WITH_JEMALLOC take AUTO, ON or OFF
# (see cmake/EhOptionalLibrary.cmake). Neither changes the output.
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/EhOptionalLibrary.cmake)
eh_optional_library(EH_WITH_ISAL isa-l/igzip_lib.h EH_ISAL_INCLUDE_DIR EH_ISAL_LIBRARY isal)
eh_optional_library(EH_WITH_JEMALLOC jemalloc/jemalloc.h EH_JEMALLOC_INCLUDE_DIR EH_JEMALLOC_LIBRARY jemalloc LINK_ONLY)
if(EH_HTS_CODEC_LIBRARY)
    set(EH_HTS_CODEC_SUMMARY ON)
else()
    set(EH_HTS_CODEC_SUMMARY OFF)
endif()
message(STATUS "Optional libraries: libdeflate (htslib) ${EH_HTS_CODEC_SUMMARY}, ISA-L ${EH_WITH_ISAL_FOUND}, jemalloc ${EH_WITH_JEMALLOC_FOUND}")
if(EH_WITH_ISAL_FOUND)
    target_compile_definitions(ExpansionHunterLib PRIVATE EH_HAVE_ISAL)
    target_include_directories(ExpansionHunterLib PRIVATE ${EH_ISAL_INCLUDE_DIR})
    target_link_libraries(ExpansionHunterLib PRIVATE ${EH_ISAL_LIBRARY})
endif()

add_executable(ExpansionHunter
        app/ExpansionHunter.cpp
        )
if(EH_WITH_JEMALLOC_FOUND)
    target_link_libraries(ExpansionHunter ExpansionHunterLib ${EH_JEMALLOC_LIBRARY})
else()
    target_link_libraries(ExpansionHunter ExpansionHunterLib)
endif()

add_executable(UnitTests
        tests/AlignMatrixTest.cpp
        tests/AlignmentClassifierTest.cpp
        tests/AlignmentSummaryTest.cpp
        tests/AlleleCheckerTest.cpp
        tests/BamletWriterTest.cpp
        tests/BamRecordSinkTest.cpp
        tests/CatalogIndexTest.cpp
        tests/ClassifierOfAlignmentsToVariantTest.cpp
        tests/ConcurrentQueueTest.cpp
        tests/CountTableTest.cpp
        tests/FragLogliksTest.cpp
        tests/GenomeMaskTest.cpp
        tests/GenomicRegionTest.cpp
        tests/GraphAlignmentOperationsTest.cpp
        tests/GraphBlueprintTest.cpp
        tests/GreedyAlignmentIntersectorTest.cpp
        tests/HighQualityBaseRunFinderTest.cpp
        tests/HtsHelpersTest.cpp
        tests/JsonWriterTest.cpp
        tests/LocusParallelKernelTest.cpp
        tests/LocusStatsTest.cpp
        tests/OrientationPredictorTest.cpp
        tests/ReadSupportCalculatorTest.cpp
        tests/ReadTest.cpp
        tests/RegionGraphTest.cpp
        tests/RepeatAnalyzerTest.cpp
        tests/RepeatGenotypeTest.cpp
        tests/SmallVariantGenotyperTest.cpp
        tests/SoftclippingAlignerTest.cpp
        tests/StrAlignTest.cpp
        tests/StreamingPairSourceTest.cpp
        tests/StrGenotyperTest.cpp
        tests/TaskExecutorTest.cpp
        tests/ThreadParametersTest.cpp
        tests/UnitTests.cpp
        tests/WeightedPurityCalculatorTest.cpp
        )
add_subdirectory(locus)

target_link_libraries(UnitTests ExpansionHunterLib GTest::GTest GTest::Main)

add_test(NAME UnitTests COMMAND UnitTests)

# Equivalence probes, re-hosted v5.0.0 references and verification drivers (see verification/README.md).
option(EH_BUILD_VERIFICATION "Build equivalence verification probes and drivers" OFF)
if(EH_BUILD_VERIFICATION)
    add_subdirectory(verification)
endif()

if(EH_VERIFY)
    install(CODE "message(FATAL_ERROR \"Refusing to install a verification build (EH_VERIFY=${EH_VERIFY})\")")
endif()
install(TARGETS ExpansionHunter UnitTests RUNTIME DESTINATION bin)
