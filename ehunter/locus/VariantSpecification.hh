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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"

#include "core/Common.hh"
#include "core/GenomicRegion.hh"
#include "core/Parameters.hh"

namespace ehunter
{

enum class VariantType
{
    kRepeat,
    kSmallVariant
};

enum class VariantSubtype
{
    kCommonRepeat,
    kRareRepeat,
    kInsertion,
    kDeletion,
    kSwap,
    kSMN
};

struct VariantClassification
{
    VariantClassification(VariantType type, VariantSubtype subtype)
        : type(type)
        , subtype(subtype)
    {
    }

    bool operator==(const VariantClassification& other) const
    {
        return (other.type == type && other.subtype == subtype);
    }

    VariantType type;
    VariantSubtype subtype;
};

struct MappedVariant
{
    graphtools::Slice id, nodes;
    GenomicRegion referenceLocus;
    VariantClassification classification;
    int64_t refNode;
    graphtools::Slice purity;
};

class VariantSpecification
{
public:
    VariantSpecification(graphtools::MappedBytes bytes, const MappedVariant* record)
        : mapped_(std::move(bytes))
        , record_(record)
    {
    }
    VariantSpecification(
        std::string id, VariantClassification classification, GenomicRegion referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode)
        : owned_(std::make_shared<Owned>(
              Owned { std::move(id), classification, referenceLocus, std::move(nodes), optionalRefNode }))
    {
        assertConsistency();
    }

    graphtools::TextView id() const { return record_ ? mapped_.text(record_->id) : graphtools::TextView(owned_->id); }
    VariantClassification classification() const { return record_ ? record_->classification : owned_->classification; }
    const GenomicRegion& referenceLocus() const& { return record_ ? record_->referenceLocus : owned_->referenceLocus; }
    GenomicRegion referenceLocus() const&& { return record_ ? record_->referenceLocus : owned_->referenceLocus; }
    void retainOwner(std::shared_ptr<const void> owner)
    {
        if (record_)
            mapped_.owner = std::move(owner);
    }
    void requirePreparedAddress(const void* address, size_t length) const
    {
        const auto start = reinterpret_cast<uintptr_t>(mapped_.data);
        const auto value = reinterpret_cast<uintptr_t>(address);
        if (!record_ || !mapped_.owner || value < start || value - start > mapped_.size
            || length > mapped_.size - (value - start))
            throw std::logic_error("Prepared analyzer parameters must belong to its retained catalog mapping");
    }
    graphtools::CatalogRange<graphtools::TextView> purityPermutations() const
    {
        if (!record_)
            return {};
        return { mapped_.data,
                 record_->purity.count,
                 [](const void* base, size_t i, const void* aux, size_t)
                 {
                     const auto slice = static_cast<const graphtools::Slice*>(aux)[i];
                     return graphtools::TextView(static_cast<const char*>(base) + slice.offset, slice.count);
                 },
                 mapped_.array<graphtools::Slice>(record_->purity),
                 mapped_.size,
                 mapped_.owner };
    }
    graphtools::CatalogRange<graphtools::NodeId> nodes() const
    {
        if (record_)
            return graphtools::CatalogRange<graphtools::NodeId>::fromArray(
                mapped_.array<graphtools::NodeId>(record_->nodes), record_->nodes.count, mapped_.owner);
        return graphtools::CatalogRange<graphtools::NodeId>::fromVector(owned_->nodes, owned_);
    }
    boost::optional<graphtools::NodeId> optionalRefNode() const
    {
        if (!record_)
            return owned_->optionalRefNode;
        return record_->refNode < 0 ? boost::optional<graphtools::NodeId>()
                                    : boost::optional<graphtools::NodeId>(record_->refNode);
    }

    bool operator==(const VariantSpecification& other) const
    {
        return id() == other.id() && classification() == other.classification() && nodes() == other.nodes();
    }

    void assertConsistency() const;

private:
    struct Owned
    {
        std::string id;
        VariantClassification classification;
        GenomicRegion referenceLocus;
        std::vector<graphtools::NodeId> nodes;
        boost::optional<graphtools::NodeId> optionalRefNode;
    };
    std::shared_ptr<const Owned> owned_;
    graphtools::MappedBytes mapped_;
    const MappedVariant* record_ = nullptr;
};

std::ostream& operator<<(std::ostream& out, VariantType type);
std::ostream& operator<<(std::ostream& out, VariantSubtype subtype);
std::ostream& operator<<(std::ostream& out, VariantClassification classification);
std::ostream& operator<<(std::ostream& out, const VariantSpecification& variantSpec);

}
