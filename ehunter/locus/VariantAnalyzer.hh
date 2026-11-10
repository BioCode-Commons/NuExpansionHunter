//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <memory>
#include <string>
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "core/Common.hh"
#include "core/LocusStats.hh"
#include "core/Parameters.hh"
#include "core/Read.hh"
#include "locus/VariantFindings.hh"
#include "locus/VariantSpecification.hh"

namespace ehunter
{

class VariantAnalyzer
{
public:
    VariantAnalyzer(
        std::string variantId, const graphtools::Graph& graph, std::vector<graphtools::NodeId> nodeIds,
        GenotyperParameters genotyperParams)
        : variantId_(std::move(variantId))
        , graph_(graph)
        , nodeIds_(std::move(nodeIds))
        , genotyperParams_(genotyperParams)
    {
    }
    virtual ~VariantAnalyzer() = default;
    VariantAnalyzer(VariantSpecification variant, const graphtools::Graph& graph, const GenotyperParameters& params)
        : graph_(graph)
        , preparedVariant_(std::move(variant))
        , preparedGenotyper_(&params)
    {
        preparedVariant_->requirePreparedAddress(&params, sizeof(params));
    }

    virtual void processMates(
        const Read& read, const graphtools::GraphAlignment& readAlignment, const Read& mate,
        const graphtools::GraphAlignment& mateAlignment)
        = 0;

    virtual std::unique_ptr<VariantFindings> analyze(const LocusStats& stats) = 0;

    graphtools::TextView variantId() const
    {
        return preparedVariant_ ? preparedVariant_->id() : graphtools::TextView(variantId_);
    }
    const graphtools::Graph& graph() const { return graph_; }
    graphtools::CatalogRange<graphtools::NodeId> nodeIds() const
    {
        return preparedVariant_ ? preparedVariant_->nodes()
                                : graphtools::CatalogRange<graphtools::NodeId>::fromVector(nodeIds_);
    }
    const GenotyperParameters& genotyperParameters() const
    {
        return preparedGenotyper_ ? *preparedGenotyper_ : genotyperParams_;
    }

protected:
    std::string variantId_;
    const graphtools::Graph& graph_;
    std::vector<graphtools::NodeId> nodeIds_;
    GenotyperParameters genotyperParams_;
    boost::optional<VariantSpecification> preparedVariant_;
    const GenotyperParameters* preparedGenotyper_ = nullptr;
};

}
