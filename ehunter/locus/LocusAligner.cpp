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
// This file differs from ExpansionHunter v5.0.0: alignment computation and output are separate and shared reads are copied only for reorientation.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "locus/LocusAligner.hh"

#include "alignment/AlignmentFilters.hh"
#include "alignment/OperationsOnAlignments.hh"

namespace ehunter
{
namespace locus
{

LocusAligner::LocusAligner(
    graphtools::TextView locusId, GraphPtr graph, const HeuristicParameters& params, AlignWriterPtr writer,
    uint64_t archiveSlot)
    : locusId_(std::move(locusId))
    , aligner_(
          graph, params.kmerLenForAlignment(), params.paddingLength(), params.seedAffixTrimLength(),
          params.catalogIndexes, archiveSlot)
    , orientationPredictor_(
          graph, params.orientationPredictorKmerLen(), params.orientationPredictorMinKmerCount(), params.catalogIndexes,
          archiveSlot + 1)
    , writer_(std::move(writer))
{
}

LocusAligner::AlignedPair LocusAligner::align(Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector)
{
    auto alignedPair = alignMates(read, mate, alignerSelector);
    writeAlignedPair(read, mate, alignedPair);
    return alignedPair;
}

static LocusAligner::AlignedPair filterPair(size_t readLength, LocusAligner::OptionalAlign readAlign,
    LocusAligner::OptionalAlign mateAlign)
{
    int numMatchingBases = static_cast<int>(static_cast<double>(readLength) / 7.5);
    numMatchingBases = std::max(numMatchingBases, 10);
    LinearAlignmentParameters parameters;
    const int kMinNonRepeatAlignmentScore = numMatchingBases * parameters.matchScore;

    if (!checkIfLocallyPlacedReadPair(readAlign, mateAlign, kMinNonRepeatAlignmentScore))
    {
        return { boost::none, boost::none };
    }

    return { std::move(readAlign), std::move(mateAlign) };
}

void LocusAligner::writeAlignedPair(const Read& read, const Read* mate, const AlignedPair& alignedPair) const
{
    if (alignedPair.first && alignedPair.second)
    {
        writer_->write(
            locusId_, read.fragmentId(), read.sequence(), read.isFirstMate(), read.isReversed(), read.isReversed(),
            *alignedPair.first);
        writer_->write(
            locusId_, mate->fragmentId(), mate->sequence(), mate->isFirstMate(), mate->isReversed(), mate->isReversed(),
            *alignedPair.second);
    }
}

template <typename Writable>
LocusAligner::OptionalAlign LocusAligner::alignRead(
    const Read& read, Writable writable, graphtools::AlignerSelector& alignerSelector) const
{
    OrientationPrediction predictedOrientation = orientationPredictor_.predict(read.sequence());
    const Read* oriented = &read;
    if (predictedOrientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
    {
        Read& reversed = writable();
        reversed.reverseComplement();
        oriented = &reversed;
    }
    else if (predictedOrientation == OrientationPrediction::kDoesNotAlign)
    {
        return {};
    }

    auto readAligns = aligner_.align(oriented->sequence(), alignerSelector);
    if (readAligns.empty())
    {
        return {};
    }

    if (readAligns.size() == 1)
    {
        return std::move(readAligns.front());
    }

    return computeCanonicalAlignment(readAligns);
}

LocusAligner::AlignedPair LocusAligner::alignMates(
    Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector) const
{
    auto first = alignRead(read, [&]() -> Read& { return read; }, alignerSelector);
    auto second = mate ? alignRead(*mate, [&]() -> Read& { return *mate; }, alignerSelector) : boost::none;
    return filterPair(read.sequence().size(), std::move(first), std::move(second));
}

LocusAligner::AlignedPair LocusAligner::alignSharedMates(const Read& read, const Read* mate,
    boost::optional<Read>& reversedRead, boost::optional<Read>& reversedMate,
    graphtools::AlignerSelector& alignerSelector) const
{
    auto first = alignRead(read, [&]() -> Read& { reversedRead.emplace(read); return *reversedRead; }, alignerSelector);
    auto second = mate ? alignRead(*mate, [&]() -> Read& {
        reversedMate.emplace(*mate);
        return *reversedMate;
    }, alignerSelector) : boost::none;
    return filterPair(read.sequence().size(), std::move(first), std::move(second));
}

}
}
