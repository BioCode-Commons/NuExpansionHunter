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

#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"
#include "thirdparty/json/json.hpp"

#include "core/Common.hh"
#include "core/GenomicRegion.hh"
#include "core/Reference.hh"
#include "locus/VariantSpecification.hh"

namespace ehunter
{

using RegionId = std::string;
using NodeToRegionAssociation = std::unordered_map<graphtools::NodeId, GenomicRegion>;

struct MappedProjection
{
    graphtools::NodeId node;
    uint32_t reserved;
    GenomicRegion region;
};
struct MappedLocus
{
    graphtools::Slice id, targets, offtargets, graph, variants, projections;
    ChromType chromosomeType;
    uint32_t reserved;
    GenotyperParameters genotyper;
};

class LocusSpecification
{
public:
    LocusSpecification(graphtools::MappedBytes bytes, const MappedLocus* record);
    LocusSpecification(
        RegionId locusId, ChromType typeOfChromLocusLocatedOn, std::vector<GenomicRegion> targetReadExtractionRegions,
        graphtools::Graph regionGraph, NodeToRegionAssociation referenceRegions, GenotyperParameters genotyperParams);

    graphtools::TextView locusId() const
    {
        return record_ ? mapped_.text(record_->id) : graphtools::TextView(owned_->id);
    }
    ChromType typeOfChromLocusLocatedOn() const { return record_ ? record_->chromosomeType : owned_->chromosomeType; }
    /*
     * List of all regions in the reference this graph describes
     * i.e. where we expect relevant reads to align
     */
    graphtools::CatalogRange<GenomicRegion> targetReadExtractionRegions() const
    {
        if (record_)
            return graphtools::CatalogRange<GenomicRegion>::fromArray(
                mapped_.array<GenomicRegion>(record_->targets), record_->targets.count, mapped_.owner);
        return graphtools::CatalogRange<GenomicRegion>::fromVector(owned_->targets, owned_);
    }
    /*
     * List of regions that additional relevant reads might be found
     * Require filtering or special considerations
     */
    graphtools::CatalogRange<GenomicRegion> offtargetReadExtractionRegions() const
    {
        if (record_)
            return graphtools::CatalogRange<GenomicRegion>::fromArray(
                mapped_.array<GenomicRegion>(record_->offtargets), record_->offtargets.count, mapped_.owner);
        return graphtools::CatalogRange<GenomicRegion>::fromVector(owned_->offtargets, owned_);
    }
    void setOfftargetReadExtractionRegions(const std::vector<GenomicRegion>& offtargetReadExtractionRegions)
    {
        writable().offtargets = offtargetReadExtractionRegions;
    }

    const graphtools::Graph& regionGraph() const& { return record_ ? mappedGraph_ : owned_->graph; }
    graphtools::Graph regionGraph() const&& { return record_ ? mappedGraph_ : owned_->graph; }
    const GenotyperParameters& genotyperParameters() const& { return record_ ? record_->genotyper : owned_->genotyper; }
    GenotyperParameters genotyperParameters() const&& { return record_ ? record_->genotyper : owned_->genotyper; }
    bool isPrepared() const { return record_ != nullptr; }
    graphtools::CatalogRange<VariantSpecification> variantSpecs() const;

    void addVariantSpecification(
        std::string id, VariantClassification classification, GenomicRegion referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode);

    VariantSpecification getVariantSpecById(const std::string& variantSpecId) const;

    graphtools::CatalogRange<std::pair<graphtools::NodeId, GenomicRegion>> referenceProjectionOfNodes() const;
    MappedLocus freeze(graphtools::BlobBuilder& out) const;

    bool requiresGenomeWideDepth() const;

private:
    struct Owned
    {
        std::string id;
        ChromType chromosomeType;
        std::vector<GenomicRegion> targets, offtargets;
        graphtools::Graph graph;
        std::vector<VariantSpecification> variants;
        NodeToRegionAssociation projections;
        GenotyperParameters genotyper;
    };
    Owned& writable();
    std::shared_ptr<Owned> owned_;
    graphtools::MappedBytes mapped_;
    const MappedLocus* record_ = nullptr;
    graphtools::Graph mappedGraph_;
};

// The mapped mode keeps just a descriptor: no vector of locus/graph objects is
// reconstructed. A small value handle is produced when a locus is visited.
class RegionCatalog
{
public:
    RegionCatalog() = default;
    RegionCatalog(std::initializer_list<LocusSpecification> loci)
    {
        for (const auto& locus : loci)
            push_back(locus);
    }
    RegionCatalog(
        graphtools::MappedBytes bytes, graphtools::Slice loci,
        std::shared_ptr<const graphtools::KmerIndexArchive> owner, graphtools::Slice ids, graphtools::Slice intervals)
        : mapped_(bytes)
        , loci_(loci)
        , ids_(ids)
        , intervals_(intervals)
        , owner_(std::move(owner))
    {
        mapped_.owner = owner_;
        mapped_.array<MappedLocus>(loci_);
    }
    size_t size() const { return mapped_.data ? loci_.count : owned_.size(); }
    bool empty() const { return size() == 0; }
    LocusSpecification operator[](size_t index) const
    {
        if (index >= size())
            throw std::out_of_range("Catalog locus index");
        return mapped_.data ? LocusSpecification(mapped_, mapped_.array<MappedLocus>(loci_) + index) : owned_[index];
    }
    auto begin() const -> graphtools::CatalogRange<LocusSpecification>::iterator { return range().begin(); }
    auto end() const -> graphtools::CatalogRange<LocusSpecification>::iterator { return range().end(); }
    void reserve(size_t size)
    {
        requireOwned();
        owned_.reserve(size);
    }
    void push_back(LocusSpecification locus)
    {
        requireOwned();
        ownedIds_.emplace(std::string(locus.locusId()), owned_.size());
        owned_.push_back(std::move(locus));
    }
    void emplace_back(LocusSpecification locus) { push_back(std::move(locus)); }
    graphtools::MappedBytes mappedBytes() const { return mapped_; }
    graphtools::Slice preparedIntervals() const { return intervals_; }
    size_t findLocusIndex(graphtools::TextView id) const
    {
        if (!mapped_.data)
            return ownedIds_.at(std::string(id));
        if (!ids_.count || (ids_.count & (ids_.count - 1)))
            throw std::runtime_error("Invalid mapped locus table");
        const auto* table = mapped_.array<uint64_t>(ids_);
        uint64_t bucket = graphtools::hashText(id) & (ids_.count - 1);
        for (uint64_t probe = 0; probe < ids_.count; ++probe)
        {
            const auto slot = table[bucket];
            if (!slot)
                throw std::out_of_range("Unknown catalog locus");
            if (slot > size())
                throw std::runtime_error("Invalid mapped locus ID");
            if (mapped_.text(mapped_.array<MappedLocus>(loci_)[slot - 1].id) == id)
                return slot - 1;
            bucket = (bucket + 1) & (ids_.count - 1);
        }
        throw std::runtime_error("Mapped locus table has no empty slot");
    }
    const std::shared_ptr<const graphtools::KmerIndexArchive>& archive() const { return owner_; }
    size_t nodeCount(size_t locus) const
    {
        if (!mapped_.data)
            return owned_.at(locus).regionGraph().numNodes();
        const auto& record = mapped_.array<MappedLocus>(loci_)[locus];
        return mapped_.array<graphtools::MappedGraph>(record.graph)->nodes.count;
    }
    boost::optional<GenomicRegion> referenceRegion(size_t locus, graphtools::NodeId node) const
    {
        if (mapped_.data)
        {
            const auto slice = mapped_.array<MappedLocus>(loci_)[locus].projections;
            const auto* values = mapped_.array<MappedProjection>(slice);
            for (uint64_t i = 0; i < slice.count; ++i)
                if (values[i].node == node)
                    return values[i].region;
        }
        else
            for (const auto& item : owned_.at(locus).referenceProjectionOfNodes())
                if (item.first == node)
                    return item.second;
        return boost::none;
    }

private:
    graphtools::CatalogRange<LocusSpecification> range() const
    {
        return { this, size(), [](const void* p, size_t i, const void*, size_t)
                 { return (*static_cast<const RegionCatalog*>(p))[i]; } };
    }
    void requireOwned() const
    {
        if (mapped_.data)
            throw std::logic_error("Prepared catalog is immutable");
    }
    std::vector<LocusSpecification> owned_;
    std::unordered_map<std::string, size_t> ownedIds_;
    graphtools::MappedBytes mapped_;
    graphtools::Slice loci_ { 0, 0 };
    graphtools::Slice ids_ { 0, 0 }, intervals_ { 0, 0 };
    std::shared_ptr<const graphtools::KmerIndexArchive> owner_;
};

}
