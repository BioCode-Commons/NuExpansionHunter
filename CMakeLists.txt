#
# Modified by NuExpansionHunter contributors, 2026.
# This file differs from ExpansionHunter v5.0.0: the build uses the Boost archive mirror and GoogleTest compiler setting, applies htslib patches, selects optional libraries and the CPU target, pins every dependency by hash or commit, and passes libcurl to the htslib and ehunter builds.
# Modifications are licensed under the Apache License, Version 2.0.
#
cmake_minimum_required(VERSION 3.13)
project(ExpansionHunter)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check for required packages:
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_package(Git REQUIRED)

include(ExternalProject)
set(installDir ${CMAKE_CURRENT_BINARY_DIR}/install)
set(EH_TARGET "portable" CACHE STRING "CPU target: portable, native, or a -march (x86-64) or -mcpu (ARM64) name (see ehunter/CMakeLists.txt)")

# libcurl as found here (the system package, or CURL_INCLUDE_DIR and CURL_LIBRARY
# given on the command line) is used by both the htslib and the ehunter builds.
list(GET CURL_LIBRARIES 0 EH_CURL_LIBRARY)
get_filename_component(EH_CURL_LIBRARY_DIR "${EH_CURL_LIBRARY}" DIRECTORY)
set(CURL_CMAKE_ARGS "-DCURL_INCLUDE_DIR:PATH=${CURL_INCLUDE_DIRS}" "-DCURL_LIBRARY:FILEPATH=${EH_CURL_LIBRARY}")

# Optional libraries: EH_WITH_LIBDEFLATE, EH_WITH_ISAL and EH_WITH_JEMALLOC take
# AUTO (use the library when its development package is installed), ON or OFF.
include(${CMAKE_SOURCE_DIR}/cmake/EhOptionalLibrary.cmake)
eh_optional_library(EH_WITH_LIBDEFLATE libdeflate.h LIBDEFLATE_INCLUDE_DIR LIBDEFLATE_LIBRARY deflate)
eh_optional_library(EH_WITH_ISAL isa-l/igzip_lib.h EH_ISAL_INCLUDE_DIR EH_ISAL_LIBRARY isal)
eh_optional_library(EH_WITH_JEMALLOC jemalloc/jemalloc.h EH_JEMALLOC_INCLUDE_DIR EH_JEMALLOC_LIBRARY jemalloc LINK_ONLY)
message(STATUS "Optional libraries: libdeflate ${EH_WITH_LIBDEFLATE_FOUND}, ISA-L ${EH_WITH_ISAL_FOUND}, jemalloc ${EH_WITH_JEMALLOC_FOUND}")

# HTSlib's make-only build uses zlib unless its optional BGZF codec is enabled.
set(HTS_CPPFLAGS "$ENV{CPPFLAGS} -I${CURL_INCLUDE_DIRS}")
set(HTS_CC "CC=${CMAKE_C_COMPILER}")
if(EH_WITH_LIBDEFLATE_FOUND)
    set(HTS_CODEC_MAKE_ARGS ${HTS_CC}
        "CPPFLAGS=${HTS_CPPFLAGS} -DHAVE_LIBDEFLATE=1 -I${LIBDEFLATE_INCLUDE_DIR}"
        "LDFLAGS=$ENV{LDFLAGS} -L${EH_CURL_LIBRARY_DIR} -Wl,-rpath,${EH_CURL_LIBRARY_DIR}"
        "LIBS=-lz -lm -lbz2 -llzma -lcurl ${LIBDEFLATE_LIBRARY}")
    set(HTS_CODEC_CMAKE_ARG "-DEH_HTS_CODEC_LIBRARY:FILEPATH=${LIBDEFLATE_LIBRARY}")
else()
    set(HTS_CODEC_MAKE_ARGS ${HTS_CC}
        "CPPFLAGS=${HTS_CPPFLAGS}"
        "LDFLAGS=$ENV{LDFLAGS} -L${EH_CURL_LIBRARY_DIR} -Wl,-rpath,${EH_CURL_LIBRARY_DIR}")
    set(HTS_CODEC_CMAKE_ARG "-DEH_HTS_CODEC_LIBRARY:FILEPATH=")
endif()

# The ehunter build uses exactly the libraries selected here.
if(EH_WITH_ISAL_FOUND)
    set(ISA_L_CMAKE_ARGS -DEH_WITH_ISAL:STRING=ON "-DEH_ISAL_INCLUDE_DIR:PATH=${EH_ISAL_INCLUDE_DIR}"
                         "-DEH_ISAL_LIBRARY:FILEPATH=${EH_ISAL_LIBRARY}")
else()
    set(ISA_L_CMAKE_ARGS -DEH_WITH_ISAL:STRING=OFF)
endif()
if(EH_WITH_JEMALLOC_FOUND)
    set(JEMALLOC_CMAKE_ARGS -DEH_WITH_JEMALLOC:STRING=ON "-DEH_JEMALLOC_LIBRARY:FILEPATH=${EH_JEMALLOC_LIBRARY}")
    if(EH_JEMALLOC_INCLUDE_DIR)
        list(APPEND JEMALLOC_CMAKE_ARGS "-DEH_JEMALLOC_INCLUDE_DIR:PATH=${EH_JEMALLOC_INCLUDE_DIR}")
    endif()
else()
    set(JEMALLOC_CMAKE_ARGS -DEH_WITH_JEMALLOC:STRING=OFF)
endif()

ExternalProject_Add(htslib
	BUILD_IN_SOURCE YES
	GIT_REPOSITORY "https://github.com/samtools/htslib.git"
	GIT_TAG fd0f89554459b78c07303e2c8a42acacd6851b46 # tag 1.10.2
	UPDATE_COMMAND ""
	PATCH_COMMAND ${GIT_EXECUTABLE} apply ${CMAKE_SOURCE_DIR}/cmake/patches/htslib-1.10.2-setup-rollback.patch
		COMMAND ${GIT_EXECUTABLE} apply ${CMAKE_SOURCE_DIR}/cmake/patches/htslib-1.10.2-runtime-lifecycle.patch
		COMMAND ${GIT_EXECUTABLE} apply ${CMAKE_SOURCE_DIR}/cmake/patches/htslib-1.10.2-reuse-libdeflate.patch
	CONFIGURE_COMMAND ""
	BUILD_COMMAND $(MAKE) ${HTS_CODEC_MAKE_ARGS}
	INSTALL_COMMAND $(MAKE) ${HTS_CODEC_MAKE_ARGS} install prefix=${installDir}
	LOG_DOWNLOAD YES
)


# GCC 13 and later no longer include <cstdint> through other standard headers;
# spdlog 1.6.1, googletest 1.10.0 and abseil 20210324.2 rely on that.
set(EH_DEPENDENCY_CXX_FLAGS "")
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
    set(EH_DEPENDENCY_CXX_FLAGS "-include cstdint")
endif()

# Setup user config to force same c++ compiler in boost
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
    set(BOOST_UCONFIG "${CMAKE_BINARY_DIR}/user-config.jam")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        file(WRITE "${BOOST_UCONFIG}" "using gcc : : \"${CMAKE_CXX_COMPILER}\" ;\n")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        file(WRITE "${BOOST_UCONFIG}" "using clang : : \"${CMAKE_CXX_COMPILER}\" ;\n")
        set(B2_OPTIONS ${B2_OPTIONS} "toolset=clang")
    endif()
    set(BOOST_PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${BOOST_UCONFIG} tools/build/src/user-config.jam)
endif()

ExternalProject_Add(Boost
	BUILD_IN_SOURCE YES
	URL https://archives.boost.io/release/1.73.0/source/boost_1_73_0.tar.bz2
	URL_HASH SHA256=4eb3b8d442b426dc35346235c8733b5ae35ba431690e38c6a8263dce9fcbb402
	UPDATE_COMMAND ""
	PATCH_COMMAND ${BOOST_PATCH_COMMAND}
	CONFIGURE_COMMAND ./bootstrap.sh --prefix=${installDir}/lib
	BUILD_COMMAND ./b2 install -j8   --prefix=${installDir} --with-filesystem --with-system --with-program_options link=static ${B2_OPTIONS}
	INSTALL_COMMAND ""
)


ExternalProject_Add(spdlog
	GIT_REPOSITORY "https://github.com/gabime/spdlog.git"
	GIT_TAG 22a169bc319ac06948e7ee0be6b9b0ac81386604 # tag v1.6.1
	UPDATE_COMMAND ""
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_POLICY_VERSION_MINIMUM=3.5
		"-DCMAKE_CXX_FLAGS=${EH_DEPENDENCY_CXX_FLAGS}"
)


ExternalProject_Add(googletest
	GIT_REPOSITORY "https://github.com/google/googletest.git"
	GIT_TAG 703bd9caab50b139428cea1aaff9974ebee5742e # tag release-1.10.0
	UPDATE_COMMAND ""
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_POLICY_VERSION_MINIMUM=3.5
		"-DCMAKE_CXX_FLAGS=${EH_DEPENDENCY_CXX_FLAGS} -Wno-error=maybe-uninitialized"
)


# CMAKE_CXX_STANDARD is required for OS X 10.15
ExternalProject_Add(abseil
	GIT_REPOSITORY "https://github.com/abseil/abseil-cpp"
	GIT_TAG 278e0a071885a22dcd2fd1b5576cc44757299343 # tag 20210324.2
	UPDATE_COMMAND ""
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_CXX_STANDARD=11
		"-DCMAKE_CXX_FLAGS=${EH_DEPENDENCY_CXX_FLAGS}"
		-DBUILD_TESTING=OFF
)


ExternalProject_Add(ehunter
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/ehunter
	BUILD_ALWAYS YES
	TEST_COMMAND "ctest"
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_PREFIX_PATH:PATH=${installDir}
        ${HTS_CODEC_CMAKE_ARG} ${ISA_L_CMAKE_ARGS} ${JEMALLOC_CMAKE_ARGS} ${CURL_CMAKE_ARGS}
		-DEH_TARGET:STRING=${EH_TARGET}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${CMAKE_SOURCE_DIR}/ehunter
)


ExternalProject_Add_StepDependencies(ehunter configure
	Boost spdlog htslib googletest abseil)

