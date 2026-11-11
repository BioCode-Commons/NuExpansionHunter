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

#include <fstream>
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
#include <chrono>

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

int main(int argc, char** argv)
{
    spdlog::set_pattern("%Y-%m-%dT%H:%M:%S,[%v]");

    try
    {
        spdlog::info("Starting {}", kProgramVersion);

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
        const auto catalogStart = std::chrono::steady_clock::now();
        HeuristicParameters heuristicParams = params.heuristics();
        const RegionCatalog regionCatalog = (!params.catalogCache.empty() && !params.prepareCatalogCache)
            ? loadPreparedCatalog(params.catalogCache, inputPaths, heuristicParams, reference.contigInfo())
            : loadLocusCatalogFromDisk(inputPaths.catalog(), heuristicParams, reference);
        spdlog::info(
            "Phase catalog_graph_seconds={}",
            std::chrono::duration<double>(std::chrono::steady_clock::now() - catalogStart).count());

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
                params.bamletQueueBytes));
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
            spdlog::info("Running sample analysis in streaming mode");
            spdlog::info(
                "Thread pools: compute={}, decode={}, encode={} (I/O pools and coordinator threads are additional)",
                params.threadCount, params.decodeThreadCount, params.disableBamletOutput ? 0 : params.encodeThreadCount);
            sampleFindings = htsStreamingSampleAnalysis(
                inputPaths, sampleParams.sex(), heuristicParams, params.threadCount, regionCatalog, bamletWriter,
                params.decodeThreadCount, params.readDecodeThreadCount);
        }

        spdlog::info("Writing output to disk");
        const auto outputStart = std::chrono::steady_clock::now();
        VcfWriter vcfWriter(sampleParams.id(), reference, regionCatalog, sampleFindings);
        writeToFile(outputPaths.vcf(), vcfWriter);

        JsonWriter jsonWriter(sampleParams, reference.contigInfo(), regionCatalog, sampleFindings);
        writeToFile(outputPaths.json(), jsonWriter);
        if (bamletOutput)
            bamletOutput->finish();
        spdlog::info(
            "Phase output_seconds={}",
            std::chrono::duration<double>(std::chrono::steady_clock::now() - outputStart).count());
    }
    catch (const std::exception& e)
    {
        spdlog::error(e.what());
        return 1;
    }

    return 0;
}
