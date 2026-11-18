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
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: locus evidence is separated from alignment computation and supports prepared catalogs and parallel genotyping.
// Modifications are licensed under the Apache License, Version 2.0.
//

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphalign/GappedAligner.hh"
#include "graphio/AlignmentWriter.hh"

#include "core/GenomicRegion.hh"
#include "core/LocusStats.hh"
#include "core/Read.hh"
#include "locus/IrrPairFinder.hh"
#include "locus/LocusAligner.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"
#include "locus/VariantAnalyzer.hh"

namespace ehunter
{
namespace locus
{

// Regions of the reference genome that can contain reads that originated in a given locus are partitioned into target
// and offtarget regions. Target regions typically consist of the reference region of the locus and possibly other
// highly-similar regions where reads typically misalign. Offtarget regions are regions where certain kinds of
// relevant reads might occasionally misalign and that require special handling (usually for efficiency reasons).
enum class RegionType
{
    kTarget,
    kOfftarget
};

using AlignWriterPtr = std::shared_ptr<graphtools::AlignmentWriter>;

class LocusEvidence
{
public:
    using Align = graphtools::GraphAlignment;
    using AlignedPair = LocusAligner::AlignedPair;
    using Node = graphtools::NodeId;

    explicit LocusEvidence(LocusSpecification spec);

    graphtools::TextView locusId() const { return locusSpec_.locusId(); }
    const LocusSpecification& locusSpec() const { return locusSpec_; }

    void consumeEvidence(const Read& read, const Read* mate, RegionType regionType, const AlignedPair& alignedPair);
    LocusFindings analyze(Sex sampleSex, boost::optional<double> genomeWideDepth, ctpl::thread_pool* gridWorkers = nullptr);

    const boost::optional<IrrPairFinder>& irrPairFinder() const { return irrPairFinder_; }
    // Evidence rows of all variants; orders the genotyping work.
    uint64_t evidenceRows() const
    {
        uint64_t rows = 0;
        for (const auto& analyzer : variantAnalyzers_)
            rows += analyzer->evidenceRows();
        return rows;
    }
    void addIrrPairFinder(std::string motif);

    void addRepeatAnalyzer(std::string variantId, Node nodeId);
    void addSmallVariantAnalyzer(
        std::string variantId, VariantSubtype subtype, std::vector<Node> nodes, boost::optional<Node> refNode);

private:
    void consumeOntargetMates(const Read& read, const Read* mate, const AlignedPair& alignedPair);
    void processOfftargetMates(const Read& read, const Read& mate);
    void runVariantAnalysis(const Read& read, const Align& readAlign, const Read& mate, const Align& mateAlign);

    LocusSpecification locusSpec_;
    LocusStatsCalculator statsCalc_;
    boost::optional<IrrPairFinder> irrPairFinder_;
    std::vector<std::unique_ptr<VariantAnalyzer>> variantAnalyzers_;
};

class LocusAnalyzer : public LocusEvidence
{
public:
    LocusAnalyzer(LocusSpecification locusSpec, const HeuristicParameters& params, AlignWriterPtr writer,
                  uint64_t archiveSlot = 0);

    void processMates(Read& read, Read* mate, RegionType regionType, graphtools::AlignerSelector& alignerSelector);

    // Concurrent computation is safe only with task-private reads and a
    // separately constructed selector per worker. No evidence or output is
    // changed here; offtarget computation leaves the reads unchanged.
    AlignedPair computeMates(
        Read& read, Read* mate, RegionType regionType, graphtools::AlignerSelector& alignerSelector) const;

    // Consume in original input order for this locus, using the same oriented
    // reads and region type passed to computeMates. The analyzer and its graph
    // must outlive all outstanding computations and returned alignments.
    void consumeMates(Read& read, Read* mate, RegionType regionType, const AlignedPair& alignedPair);
    void writeAlignedMates(const Read& read, const Read* mate, RegionType regionType, const AlignedPair& alignedPair) const;
private:
    LocusAligner aligner_;
};

}
}
