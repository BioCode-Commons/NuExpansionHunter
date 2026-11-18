//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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
// This file differs from ExpansionHunter v5.0.0: catalog loading supports parallel parsing of complete locus objects.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "io/CatalogLoading.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "core/Common.hh"
#include "core/OrderedBlocks.hh"
#include "core/Reference.hh"
#include "io/LocusSpecDecoding.hh"
#include "io/GraphBlueprint.hh"

using boost::optional;
using graphtools::NodeId;
using std::map;
using std::ostream;
using std::string;
using std::to_string;
using std::vector;

using Json = nlohmann::json;

namespace spd = spdlog;

namespace ehunter
{

static bool checkIfFieldExists(const Json& record, const string& fieldName)
{
    return record.find(fieldName) != record.end();
}

static void assertFieldExists(const Json& record, const string& fieldName)
{
    if (!checkIfFieldExists(record, fieldName))
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Field " + fieldName + " must be present in " + out.str());
    }
}

static void assertRecordIsArray(const Json& record)
{
    if (!record.is_array())
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Expected array but got this instead " + out.str());
    }
}

static void makeArray(Json& record)
{
    if (record.type() != Json::value_t::array)
    {
        record = Json::array({ record });
    }
}

static VariantTypeFromUser decodeVariantTypeFromUser(const string& encoding)
{
    if (encoding == "RareRepeat")
    {
        return VariantTypeFromUser::kRareRepeat;
    }
    if (encoding == "Repeat")
    {
        return VariantTypeFromUser::kCommonRepeat;
    }
    if (encoding == "SmallVariant")
    {
        return VariantTypeFromUser::kSmallVariant;
    }
    if (encoding == "SMN")
    {
        return VariantTypeFromUser::kSMN;
    }
    else
    {
        throw std::logic_error("Encountered invalid variant type: " + encoding);
    }
}

static vector<string> generateIds(const std::string& locusId, const Json& variantRegionEncodings)
{
    if (variantRegionEncodings.size() == 1)
    {
        return { locusId };
    }

    vector<string> variantIds;
    for (const auto& regionEncoding : variantRegionEncodings)
    {
        variantIds.push_back(locusId + "_" + regionEncoding.get<string>());
    }

    return variantIds;
}

static LocusDescriptionFromUser loadUserDescription(Json& locusJson, const ReferenceContigInfo& contigInfo)
{
    LocusDescriptionFromUser userDescription;

    assertFieldExists(locusJson, "LocusId");
    userDescription.locusId = locusJson["LocusId"].get<string>();

    assertFieldExists(locusJson, "ReferenceRegion");
    makeArray(locusJson["ReferenceRegion"]);
    for (const auto& encoding : locusJson["ReferenceRegion"])
    {
        GenomicRegion region = decode(contigInfo, encoding.get<string>());
        userDescription.referenceRegions.push_back(region);
    }

    assertFieldExists(locusJson, "LocusStructure");
    userDescription.locusStructure = locusJson["LocusStructure"].get<string>();

    assertFieldExists(locusJson, "VariantType");
    makeArray(locusJson["VariantType"]);
    for (const auto& encoding : locusJson["VariantType"])
    {
        userDescription.variantTypesFromUser.push_back(decodeVariantTypeFromUser(encoding.get<string>()));
    }

    if (checkIfFieldExists(locusJson, "TargetRegion"))
    {
        makeArray(locusJson["TargetRegion"]);
        for (const auto& locusEncoding : locusJson["TargetRegion"])
        {
            GenomicRegion region = decode(contigInfo, locusEncoding.get<string>());
            userDescription.targetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(locusJson, "VariantId"))
    {
        makeArray(locusJson["VariantId"]);
        for (const auto& variantId : locusJson["VariantId"])
        {
            userDescription.variantIds.push_back(variantId.get<string>());
        }
    }
    else
    {
        userDescription.variantIds = generateIds(userDescription.locusId, locusJson["ReferenceRegion"]);
    }

    if (checkIfFieldExists(locusJson, "OfftargetRegions"))
    {
        assertRecordIsArray(locusJson["OfftargetRegions"]);
        for (const auto& locusEncoding : locusJson["OfftargetRegions"])
        {
            GenomicRegion region = decode(contigInfo, locusEncoding.get<string>());
            userDescription.offtargetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(locusJson, "ErrorRate"))
    {
        userDescription.errorRate = locusJson["ErrorRate"].get<double>();
    }
    if (checkIfFieldExists(locusJson, "LikelihoodRatioThreshold"))
    {
        userDescription.likelihoodRatioThreshold = locusJson["LikelihoodRatioThreshold"].get<double>();
    }
    if (checkIfFieldExists(locusJson, "MinimalLocusCoverage"))
    {
        userDescription.minLocusCoverage = locusJson["MinimalLocusCoverage"].get<double>();
    }

    return userDescription;
}

// Locate whole objects in an array without constructing the whole catalog DOM.
// Workers parse each slice with the same JSON parser that handles the fallback.
static bool catalogObjects(const string& text, vector<std::pair<size_t, size_t>>& ranges)
{
    const auto next = [&](size_t pos) { return text.find_first_not_of(" \t\r\n", pos); };
    size_t pos = next(0);
    if (pos == string::npos || text[pos] != '[')
        return false;
    pos = next(pos + 1);
    if (pos < text.size() && text[pos] == ']')
        return next(pos + 1) == string::npos;
    while (pos < text.size() && text[pos] == '{')
    {
        const size_t begin = pos;
        size_t depth = 1;
        bool inString = false;
        while (++pos < text.size())
        {
            const char c = text[pos];
            if (inString)
            {
                if (c == '\\') ++pos;
                else if (c == '"') inString = false;
            }
            else if (c == '"') inString = true;
            else if (c == '{') ++depth;
            else if (c == '}' && --depth == 0) break;
        }
        if (pos >= text.size())
            return false;
        ranges.emplace_back(begin, pos + 1);
        pos = next(pos + 1);
        if (pos == string::npos)
            return false;
        if (text[pos] == ']')
            return next(pos + 1) == string::npos;
        if (text[pos] != ',')
            return false;
        pos = next(pos + 1);
    }
    return false;
}

// Contract [K-29]
// Ensures:  the sliced parallel load yields the same loci in the same order as the whole-file parse at one thread.
RegionCatalog loadLocusCatalogFromDisk(
    const string& catalogPath, const HeuristicParameters& heuristicParams, const Reference& reference, unsigned threads)
{
    std::ifstream inputStream(catalogPath.c_str());

    if (!inputStream.is_open())
    {
        throw std::runtime_error("Failed to open catalog file " + catalogPath);
    }

    Json catalogJson;
    string text;
    vector<std::pair<size_t, size_t>> ranges;
    bool useSlices = false;
    if (threads > 1)
    {
        char block[64 * 1024];
        while (inputStream.read(block, sizeof(block)) || inputStream.gcount())
            text.append(block, static_cast<size_t>(inputStream.gcount()));
        if (inputStream.bad())
            throw std::runtime_error("Failed to read catalog file " + catalogPath);
        useSlices = catalogObjects(text, ranges);
        if (!useSlices)
        {
            std::istringstream fallback(text);
            fallback >> catalogJson;
        }
    }
    else
        inputStream >> catalogJson;
    if (!useSlices)
        makeArray(catalogJson);

    const size_t count = useSlices ? ranges.size() : catalogJson.size();
    RegionCatalog catalog;
    catalog.reserve(count);
    const auto* fasta = dynamic_cast<const FastaReference*>(&reference);
    const size_t grain = 256;
    const size_t workers = fasta ? std::min<size_t>(std::max(1u, threads), (count + grain - 1) / grain) : 1;
    vector<std::unique_ptr<FastaReference>> handles(workers);
    auto decode = [&](int worker, size_t begin, size_t end) {
        if (workers > 1 && !handles[worker])
            handles[worker].reset(new FastaReference(fasta->path(), reference.contigInfo()));
        const Reference& local = handles[worker] ? *handles[worker] : reference;
        vector<LocusSpecification> loci;
        loci.reserve(end - begin);
        for (size_t i = begin; i < end; ++i)
        {
            Json record = useSlices ? Json::parse(text.cbegin() + ranges[i].first, text.cbegin() + ranges[i].second)
                                    : std::move(catalogJson[i]);
            auto description = loadUserDescription(record, local.contigInfo());
            loci.emplace_back(decodeLocusSpecification(description, local, heuristicParams));
        }
        return loci;
    };
    orderedBlocks(count, grain, workers, decode, [&](vector<LocusSpecification> loci) {
        for (auto& locus : loci)
            catalog.emplace_back(std::move(locus));
    });

    return catalog;
}

}
