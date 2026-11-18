//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: analysis uses prepared catalogs, independent compute and codec pools, parallel text output, and explicit output-error propagation.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include <dlfcn.h>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// clang-format off
// Note that spdlog.h must be included before ostr.h
#include "spdlog/spdlog.h"
//#include "thirdparty/spdlog/include/spdlog/fmt/ostr.h"
// clang-format on

#include "app/Version.hh"
#include "core/Parameters.hh"
#include "core/OrderedBlocks.hh"
#include "core/VerifyPolicy.hh"
#include "io/BamletWriter.hh"
#include "io/CatalogLoading.hh"
#include "io/JsonWriter.hh"
#include "io/ParameterLoading.hh"
#include "io/PreparedCatalog.hh"
#include "io/SampleStats.hh"
#include "io/VcfWriter.hh"
#include "locus/VariantFindings.hh"
#include "sample/HtsSeekingSampleAnalysis.hh"
#include "sample/HtsStreamingSampleAnalysis.hh"

namespace spd = spdlog;

using namespace ehunter;

template <typename T> static void writeToFile(std::string fileName, T streamable)
{
    std::ofstream out;

    out.open(fileName.c_str());
    if (!out.is_open())
    {
        throw std::runtime_error("Failed to open " + fileName + " for writing (" + strerror(errno) + ")");
    }

    out << streamable;
    out.close();
    if (!out)
        throw std::runtime_error("Failed to write " + fileName);
}

void setLogLevel(LogLevel logLevel)
{
    switch (logLevel)
    {
    case LogLevel::kTrace:
        spdlog::set_level(spdlog::level::trace);
        break;
    case LogLevel::kDebug:
        spdlog::set_level(spdlog::level::debug);
        break;
    case LogLevel::kInfo:
        spdlog::set_level(spdlog::level::info);
        break;
    case LogLevel::kWarn:
        spdlog::set_level(spdlog::level::warn);
        break;
    case LogLevel::kError:
        spdlog::set_level(spdlog::level::err);
        break;
    }
}

// Return the pages of freed allocations to the system now. With jemalloc they
// would otherwise stay resident until its decay purges them.
static void returnFreedMemory()
{
    using Mallctl = int (*)(const char*, void*, size_t*, void*, size_t);
    const auto mallctl = reinterpret_cast<Mallctl>(dlsym(RTLD_DEFAULT, "mallctl"));
    if (mallctl)
        mallctl("arena.4096.purge", nullptr, nullptr, nullptr, 0);
}

int main(int argc, char** argv)
{
    spdlog::set_pattern("%Y-%m-%dT%H:%M:%S,[%v]");

    try
    {
        spdlog::info("Starting {}", kProgramVersion);
        if (verify::kAnySet)
            spdlog::warn("Verification build, not for production use: {}", verify::describeSettings());

        auto optionalProgramParameters = tryLoadingProgramParameters(argc, argv);
        if (!optionalProgramParameters)
        {
            return 0;
        }
        const ProgramParameters& params = *optionalProgramParameters;

        setLogLevel(params.logLevel());

        const SampleParameters& sampleParams = params.sample();

        spdlog::info("Analyzing sample {}", sampleParams.id());

        const InputPaths& inputPaths = params.inputPaths();

        spdlog::info("Initializing reference {}", inputPaths.reference());
        FastaReference reference(inputPaths.reference(), extractReferenceContigInfo(inputPaths.htsFile()));

        spdlog::info("Loading variant catalog from disk {}", inputPaths.catalog());
        HeuristicParameters heuristicParams = params.heuristics();
        RegionCatalog regionCatalog = (!params.catalogCache.empty() && !params.prepareCatalogCache)
            ? loadPreparedCatalog(params.catalogCache, inputPaths, heuristicParams, reference.contigInfo())
            : loadLocusCatalogFromDisk(inputPaths.catalog(), heuristicParams, reference, params.threadCount);

        if (params.prepareCatalogCache)
        {
            spdlog::info("Preparing reusable catalog k-mer indexes");
            prepareCatalogCache(
                params.catalogCache, regionCatalog, inputPaths, heuristicParams, reference.contigInfo());
            spdlog::info("Catalog k-mer preparation complete");
            return 0;
        }
        if (!params.catalogCache.empty())
        {
            heuristicParams.catalogIndexes = regionCatalog.archive();
        }

        const OutputPaths& outputPaths = params.outputPaths();

        locus::AlignWriterPtr bamletWriter;
        std::shared_ptr<BamletWriter> bamletOutput;
        if (params.disableBamletOutput)
        {
            bamletWriter.reset(new graphtools::BlankAlignmentWriter());
        }
        else
        {
            bamletOutput.reset(new BamletWriter(
                outputPaths.bamlet(), reference.contigInfo(), regionCatalog, params.encodeThreadCount,
                params.bamletQueueBytes, params.analysisMode() == AnalysisMode::kStreaming, params.threadCount));
            bamletWriter = bamletOutput;
        }

        SampleFindings sampleFindings;
        if (params.analysisMode() == AnalysisMode::kSeeking)
        {
            spdlog::info("Running sample analysis in seeking mode");
            sampleFindings = htsSeekingSampleAnalysis(
                inputPaths, sampleParams.sex(), heuristicParams, params.threadCount, regionCatalog, bamletWriter);
        }
        else
        {
            if (!regionCatalog.mappedBytes().data && !verify::kReferenceLayout)
            {
                // Hold the catalog in the compact frozen layout for the run. BAM
                // output has validated every locus mapping on the owned catalog.
                auto compact = compactCatalog(regionCatalog, static_cast<unsigned>(params.threadCount));
                regionCatalog.clear(static_cast<unsigned>(params.threadCount));
                regionCatalog = std::move(compact);
                returnFreedMemory();
            }
            spdlog::info("Running sample analysis in streaming mode");
            spdlog::info(
                "Thread pools: compute={}, decode={}, encode={} (I/O pools and coordinator threads are additional)",
                params.threadCount, params.decodeThreadCount, params.disableBamletOutput ? 0 : params.encodeThreadCount);
            sampleFindings = htsStreamingSampleAnalysis(
                inputPaths, sampleParams.sex(), heuristicParams, params.threadCount, regionCatalog, bamletWriter,
                params.decodeThreadCount);
        }

        spdlog::info("Writing output to disk");
        // VCF and JSON traverse the same findings in different orders. Sharing
        // workers lets either format use all available capacity as the other
        // enters serial sorting, header writing or its final file close.
        std::unique_ptr<ctpl::thread_pool> textWorkers;
        if (params.threadCount > 1 && sampleFindings.size() >= 1024)
            textWorkers.reset(new ctpl::thread_pool(params.threadCount));
        const unsigned jsonThreads = textWorkers ? params.threadCount : std::max(1, params.threadCount / 2);
        auto writeJson = [&] {
            JsonWriter writer(sampleParams, reference.contigInfo(), regionCatalog, sampleFindings, jsonThreads, textWorkers.get());
            writeToFile(outputPaths.json(), writer);
        };
        const auto policy = params.threadCount > 1 ? std::launch::async : std::launch::deferred;
        auto jsonOutput = std::async(policy, writeJson);
        VcfWriter vcfWriter(sampleParams.id(), reference, regionCatalog, sampleFindings,
                            textWorkers ? params.threadCount : std::max(1, params.threadCount - static_cast<int>(jsonThreads)),
                            textWorkers.get());
        writeToFile(outputPaths.vcf(), vcfWriter);
        jsonOutput.get();
        textWorkers.reset();
        if (bamletOutput)
            bamletOutput->finish();
        orderedBlocks(sampleFindings.size(), 1024, params.threadCount, [&](int, size_t begin, size_t end) {
            for (size_t index = begin; index < end; ++index)
            {
                // Destroy each locus's hash buckets on this worker as well as
                // its variants. clear() retains the buckets until process exit.
                auto& variants = sampleFindings[index].findingsForEachVariant;
                decltype(LocusFindings::findingsForEachVariant) empty;
                variants.swap(empty);
            }
            return true;
        }, [](bool) {});
        regionCatalog.clear(params.threadCount);
    }
    catch (const std::exception& e)
    {
        spdlog::error(e.what());
        return 1;
    }

    return 0;
}
