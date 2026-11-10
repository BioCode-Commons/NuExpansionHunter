//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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
#include <unordered_map>
#include <vector>

#include "thirdparty/intervaltree/IntervalTree.h"

#include "core/Read.hh"
#include "locus/LocusAnalyzer.hh"
#include "sample/GenomeMask.hh"
#include "sample/PreparedIntervals.hh"

namespace ehunter
{

// Specifies which mates should be processed with a given locus analyzer
enum class AnalyzerInputType
{
    kReadOnly,
    kMateOnly,
    kBothReads
};

// Stores information needed to properly pass reads to the analyzer
struct AnalyzerBundle
{
    AnalyzerBundle(locus::RegionType regionType, const size_t initLocusAnalyzerIndex)
        : regionType(regionType)
        , inputType(AnalyzerInputType::kBothReads)
        , locusIndex(initLocusAnalyzerIndex)
    {
    }

    locus::RegionType regionType;
    AnalyzerInputType inputType;
    size_t locusIndex;
};

void processAnalyzerBundleReadPair(
    locus::LocusAnalyzer& locusAnalyzer, locus::RegionType regionType, AnalyzerInputType inputType, Read& read,
    Read& mate, graphtools::AlignerSelector& alignerSelector);

// Enables retrieval of appropriate locus analyzers by genomic coordinates of read alignments
class AnalyzerFinder
{
public:
    AnalyzerFinder(
        std::vector<std::unique_ptr<locus::LocusAnalyzer>>& locusAnalyzers, const RegionCatalog* catalog = nullptr);

    // Retrieves analyzers appropriate for the given read pair
    std::vector<AnalyzerBundle> query(
        int32_t readContigId, int64_t readStart, int64_t readEnd, int32_t mateContigId, int64_t mateStart,
        int64_t mateEnd) const;

    // Retrieves analyzers appropriate for the given read
    std::vector<AnalyzerBundle> query(int32_t readContigId, int64_t readStart, int64_t readEnd) const;

    // Current alignments advance in coordinate order; mate coordinates need not.
    std::vector<AnalyzerBundle> queryStreaming(
        int32_t readContigId, int64_t readStart, int64_t readEnd, int32_t mateContigId, int64_t mateStart,
        int64_t mateEnd);

private:
    using AnalyzerIntervalTree = IntervalTree<std::size_t, AnalyzerBundle>;
    using AnalyzerIntervalTrees = std::unordered_map<int32_t, AnalyzerIntervalTree>;

    AnalyzerIntervalTrees intervalTrees_;
    graphtools::MappedBytes prepared_;
    graphtools::Slice preparedRoots_ { 0, 0 };
    const PreparedContigIntervals* preparedRoot(int32_t contig) const;
    void queryPrepared(graphtools::Slice tree, uint64_t start, uint64_t end, std::vector<AnalyzerBundle>& result) const;
    struct SweepInterval
    {
        std::size_t start, stop, rank;
        AnalyzerBundle bundle;
    };
    struct Sweep
    {
        std::vector<SweepInterval> intervals;
        graphtools::CatalogRange<PreparedInterval> prepared;
        std::vector<std::size_t> active;
        std::size_t next = 0;
        int64_t position = -1;
        size_t size() const { return prepared.empty() ? intervals.size() : prepared.size(); }
        SweepInterval at(size_t i) const
        {
            if (prepared.empty())
                return intervals[i];
            const auto interval = prepared[i];
            return { size_t(interval.start), size_t(interval.stop), size_t(interval.rank),
                     AnalyzerBundle(static_cast<locus::RegionType>(interval.regionType), interval.locus) };
        }
    };
    std::unordered_map<int32_t, Sweep> sweeps_;
};

}
