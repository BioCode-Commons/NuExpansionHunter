//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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

#include <boost/optional.hpp>

#include "alignment/OrientationPredictor.hh"
#include "core/Parameters.hh"
#include "core/Read.hh"

#include "graphalign/GappedAligner.hh"
#include "graphio/AlignmentWriter.hh"

namespace ehunter
{
namespace locus
{

class LocusAligner
{
public:
    using GraphPtr = const graphtools::Graph*;
    using Align = graphtools::GraphAlignment;
    using OptionalAlign = boost::optional<Align>;
    using AlignedPair = std::pair<OptionalAlign, OptionalAlign>;
    using AlignWriterPtr = std::shared_ptr<graphtools::AlignmentWriter>;

    LocusAligner(
        graphtools::TextView locusId, GraphPtr graph, const HeuristicParameters& params, AlignWriterPtr writer,
        uint64_t archiveSlot = 0);

    /// \param[in,out] alignerSelector A per-thread alignment workspace which mutates during alignment
    ///
    AlignedPair align(Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector);

    // Computes the same pair result without writing output. Concurrent calls
    // require distinct mutable reads and independently constructed selectors.
    // Returned alignments borrow the graph supplied to this aligner.
    AlignedPair alignMates(Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector) const;

    // Emit a result with the reads in their post-alignment orientation. The
    // caller must honor the writer's lifetime and serialization contract.
    void writeAlignedPair(const Read& read, const Read* mate, const AlignedPair& alignedPair) const;

private:
    OptionalAlign align(Read& read, graphtools::AlignerSelector& alignerSelector) const;

    graphtools::TextView locusId_;
    graphtools::GappedGraphAligner aligner_;
    OrientationPredictor orientationPredictor_;
    AlignWriterPtr writer_;
};

}
}
