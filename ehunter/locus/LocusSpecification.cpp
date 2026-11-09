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

#include "locus/LocusSpecification.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "core/Common.hh"
#include "core/Reference.hh"
#include "core/WeightedPurityCalculator.hh"

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
LocusSpecification::LocusSpecification(
    RegionId locusId, ChromType typeOfChromLocusLocatedOn, std::vector<GenomicRegion> targetReadExtractionRegions,
    graphtools::Graph regionGraph, NodeToRegionAssociation referenceRegions, GenotyperParameters genotyperParams)
    : owned_(std::make_shared<Owned>(Owned { std::move(locusId),
                                             typeOfChromLocusLocatedOn,
                                             std::move(targetReadExtractionRegions),
                                             {},
                                             std::move(regionGraph),
                                             {},
                                             std::move(referenceRegions),
                                             genotyperParams }))
{
}

LocusSpecification::LocusSpecification(graphtools::MappedBytes bytes, const MappedLocus* record)
    : mapped_(bytes)
    , record_(record)
    , mappedGraph_(graphtools::Graph::mapped(bytes, record->graph))
{
}

LocusSpecification::Owned& LocusSpecification::writable()
{
    if (record_)
        throw std::logic_error("Prepared locus is immutable");
    if (!owned_.unique())
        owned_ = std::make_shared<Owned>(*owned_);
    return *owned_;
}

graphtools::CatalogRange<VariantSpecification> LocusSpecification::variantSpecs() const
{
    if (!record_)
        return graphtools::CatalogRange<VariantSpecification>::fromVector(owned_->variants, owned_);
    return { mapped_.data,
             record_->variants.count,
             [](const void* p, size_t i, const void* aux, size_t extent)
             {
                 return VariantSpecification(
                     graphtools::MappedBytes { static_cast<const char*>(p), extent },
                     static_cast<const MappedVariant*>(aux) + i);
             },
             mapped_.array<MappedVariant>(record_->variants),
             mapped_.size,
             mapped_.owner };
}

graphtools::CatalogRange<std::pair<graphtools::NodeId, GenomicRegion>>
LocusSpecification::referenceProjectionOfNodes() const
{
    using Pair = std::pair<graphtools::NodeId, GenomicRegion>;
    if (!record_)
        return graphtools::CatalogRange<Pair>::fromContainer(owned_->projections, owned_);
    return { mapped_.array<MappedProjection>(record_->projections),
             record_->projections.count,
             [](const void* p, size_t i, const void*, size_t)
             {
                 const auto& item = static_cast<const MappedProjection*>(p)[i];
                 return Pair(item.node, item.region);
             },
             nullptr,
             0,
             mapped_.owner };
}

MappedLocus LocusSpecification::freeze(graphtools::BlobBuilder& out) const
{
    std::vector<GenomicRegion> targets = targetReadExtractionRegions(), offtargets = offtargetReadExtractionRegions();
    std::vector<MappedVariant> variants;
    for (const auto& variant : variantSpecs())
    {
        const std::vector<graphtools::NodeId> nodes = variant.nodes();
        const auto reference = variant.optionalRefNode();
        std::vector<graphtools::Slice> purity;
        if (variant.classification().subtype == VariantSubtype::kRareRepeat)
        {
            WeightedPurityCalculator calculator(std::string(regionGraph().nodeSeq(nodes.front())));
            for (const auto& permutation : calculator.permutations())
                purity.push_back(out.text(permutation));
        }
        variants.push_back({ out.text(variant.id()), out.append(nodes), variant.referenceLocus(),
                             variant.classification(), reference ? int64_t(*reference) : -1, out.append(purity) });
    }
    std::vector<MappedProjection> projections;
    for (const auto& item : referenceProjectionOfNodes())
    {
        if (regionGraph().nodeSeq(item.first).size() != uint64_t(item.second.length()))
            throw std::logic_error("Reference projection length does not match graph node");
        projections.push_back({ item.first, 0, item.second });
    }
    return { out.text(locusId()),  out.append(targets),     out.append(offtargets),      regionGraph().freeze(out),
             out.append(variants), out.append(projections), typeOfChromLocusLocatedOn(), 0,
             genotyperParameters() };
}

void LocusSpecification::addVariantSpecification(
    std::string id, VariantClassification classification, GenomicRegion referenceLocus, vector<NodeId> nodes,
    optional<NodeId> refNode)
{
    writable().variants.emplace_back(
        std::move(id), classification, std::move(referenceLocus), std::move(nodes), refNode);
}

VariantSpecification LocusSpecification::getVariantSpecById(const std::string& variantSpecId) const
{
    for (const auto& variantSpec : variantSpecs())
    {
        if (variantSpec.id() == variantSpecId)
        {
            return variantSpec;
        }
    }

    throw std::logic_error("There is no variant " + variantSpecId + " in locus " + locusId());
}

bool LocusSpecification::requiresGenomeWideDepth() const
{
    for (const auto& variantSpec : variantSpecs())
    {
        if (variantSpec.classification().subtype == VariantSubtype::kSMN)
        {
            return true;
        }
    }

    return false;
}

}
