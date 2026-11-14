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

#include "sample/AnalyzerFinder.hh"

using boost::optional;
using ehunter::locus::LocusAnalyzer;
using ehunter::locus::RegionType;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

namespace ehunter
{

namespace
{
bool areMatesNearby(int32_t readContigId, int64_t readPosition, int32_t mateContigId, int64_t matePosition)
{
    const int kMaxMateDistance = 1000;
    return ((readContigId == mateContigId) && (std::abs(readPosition - matePosition) < kMaxMateDistance));
}

inline RegionType coalesceRegionTypes(RegionType readRegionType, RegionType mateRegionType)
{
    if (readRegionType == RegionType::kTarget || mateRegionType == RegionType::kTarget)
    {
        return RegionType::kTarget;
    }

    return RegionType::kOfftarget;
}

/// Remove items from \p bundles which refer to a LocusAnalyzer already found in \p commonBundles
///
void filterOutCommonBundles(const vector<AnalyzerBundle>& commonBundles, vector<AnalyzerBundle>& bundles)
{
    if (commonBundles.size() == bundles.size())
    {
        bundles.clear();
    }
    else
    {
        auto nextRemaining = bundles.begin();
        for (auto& bundle : bundles)
        {
            bool bundleAssigned(false);
            for (const auto& commonBundle : commonBundles)
            {
                if (bundle.locusIndex == commonBundle.locusIndex)
                {
                    bundleAssigned = true;
                    break;
                }
            }
            if (not bundleAssigned)
            {
                *nextRemaining++ = bundle;
            }
        }
        bundles.erase(nextRemaining, bundles.end());
    }
}

/// \param[in,out] readBundles Used to find common bundles, updated to contain any remaining bundles on return
///
/// \param[in,out] mateBundles Used to find common bundles, updated to contain any remaining bundles on return
///
vector<AnalyzerBundle> coalesceCommonBundles(vector<AnalyzerBundle>& readBundles, vector<AnalyzerBundle>& mateBundles)
{
    vector<AnalyzerBundle> commonBundles;

    for (const auto& readBundle : readBundles)
    {
        for (const auto& mateBundle : mateBundles)
        {
            if (readBundle.locusIndex == mateBundle.locusIndex)
            {
                commonBundles.emplace_back(readBundle);
                commonBundles.back().regionType = coalesceRegionTypes(readBundle.regionType, mateBundle.regionType);
                break;
            }
        }
    }

    if (not commonBundles.empty())
    {
        filterOutCommonBundles(commonBundles, readBundles);
        filterOutCommonBundles(commonBundles, mateBundles);
    }

    return commonBundles;
}

/// We ignore nearby pairs where one mate is inside and one mate is outside of the offtarget region
///
/// \param[in,out] bundles Coalesced bundles are appended to this structure
///
void coalesceBundlesForNearbyMates(
    const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles,
    vector<AnalyzerBundle>& bundles)
{
    for (const auto& bundle : readBundles)
    {
        if (bundle.regionType == RegionType::kTarget)
        {
            bundles.push_back(bundle);
            bundles.back().inputType = AnalyzerInputType::kReadOnly;
        }
    }

    for (const auto& bundle : mateBundles)
    {
        if (bundle.regionType == RegionType::kTarget)
        {
            bundles.push_back(bundle);
            bundles.back().inputType = AnalyzerInputType::kMateOnly;
        }
    }
}

/// \param[in,out] bundles Coalesced bundles are appended to this structure
///
void coalesceBundlesForFarawayMates(
    const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles,
    vector<AnalyzerBundle>& bundles)
{
    for (const auto& bundle : readBundles)
    {
        bundles.push_back(bundle);
        bundles.back().inputType = AnalyzerInputType::kBothReads;
    }

    for (const auto& bundle : mateBundles)
    {
        bundles.push_back(bundle);
        bundles.back().inputType = AnalyzerInputType::kBothReads;
    }
}
}

void processAnalyzerBundleReadPair(
    locus::LocusAnalyzer& locusAnalyzer, locus::RegionType regionType, AnalyzerInputType inputType, Read& read,
    Read& mate, graphtools::AlignerSelector& alignerSelector)
{

    switch (inputType)
    {
    case AnalyzerInputType::kBothReads:
        locusAnalyzer.processMates(read, &mate, regionType, alignerSelector);
        break;
    case AnalyzerInputType::kReadOnly:
        locusAnalyzer.processMates(read, nullptr, regionType, alignerSelector);
        break;
    case AnalyzerInputType::kMateOnly:
        locusAnalyzer.processMates(mate, nullptr, regionType, alignerSelector);
        break;
    }
}

AnalyzerFinder::AnalyzerFinder(vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers, const RegionCatalog* catalog)
{
    if (catalog && catalog->mappedBytes().data)
    {
        prepared_ = catalog->mappedBytes();
        preparedRoots_ = catalog->preparedIntervals();
        prepared_.array<PreparedContigIntervals>(preparedRoots_);
        return;
    }
    using IntervalWithLocusTypeAndAnalyzer = Interval<std::size_t, AnalyzerBundle>;

    unordered_map<int32_t, vector<IntervalWithLocusTypeAndAnalyzer>> contigToIntervals;

    const unsigned locusAnalzerCount(locusAnalyzers.size());
    for (unsigned locusAnalyzerIndex(0); locusAnalyzerIndex < locusAnalzerCount; ++locusAnalyzerIndex)
    {
        const auto& locusAnalyzer(locusAnalyzers[locusAnalyzerIndex]);
        const LocusSpecification& locusSpec = locusAnalyzer->locusSpec();
        for (const auto& region : locusSpec.targetReadExtractionRegions())
        {
            AnalyzerBundle bundle(RegionType::kTarget, locusAnalyzerIndex);
            contigToIntervals[region.contigIndex()].emplace_back(region.start(), region.end(), bundle);
        }

        for (const auto& region : locusSpec.offtargetReadExtractionRegions())
        {
            AnalyzerBundle bundle(RegionType::kOfftarget, locusAnalyzerIndex);
            contigToIntervals[region.contigIndex()].emplace_back(region.start(), region.end(), bundle);
        }
    }

    for (auto& contigAndIntervals : contigToIntervals)
    {
        int32_t contigIndex = contigAndIntervals.first;
        auto intervals = contigAndIntervals.second;

        intervalTrees_.emplace(std::make_pair(contigIndex, AnalyzerIntervalTree(std::move(intervals))));
    }
}

vector<AnalyzerBundle> AnalyzerFinder::query(int32_t contigIndex, int64_t start, int64_t end) const
{
    if (prepared_.data)
    {
        vector<AnalyzerBundle> result;
        const auto root = preparedRoot(contigIndex);
        if (root && start >= 0 && end >= 0)
            queryPrepared(root->tree, start, end, result);
        return result;
    }
    const auto contigTreeIterator = intervalTrees_.find(contigIndex);
    if (contigTreeIterator == intervalTrees_.end())
    {
        return vector<AnalyzerBundle>();
    }

    vector<AnalyzerBundle> analyzerBundles;
    contigTreeIterator->second.visit_overlapping(
        start, end,
        [&](const AnalyzerIntervalTree::interval& intervalWithBundle)
        {
            const bool isReadContainedInInterval = static_cast<int64_t>(intervalWithBundle.start) <= start
                && end <= static_cast<int64_t>(intervalWithBundle.stop);
            if (isReadContainedInInterval)
            {
                analyzerBundles.push_back(intervalWithBundle.value);
            }
        });

    return analyzerBundles;
}

const PreparedContigIntervals* AnalyzerFinder::preparedRoot(int32_t contig) const
{
    if (contig < 0 || uint64_t(contig) >= preparedRoots_.count)
        return nullptr;
    return prepared_.array<PreparedContigIntervals>(preparedRoots_) + contig;
}

void AnalyzerFinder::queryPrepared(
    graphtools::Slice tree, uint64_t start, uint64_t end, std::vector<AnalyzerBundle>& result) const
{
    if (!tree.count)
        return;
    if (tree.count != 1)
        throw std::runtime_error("Invalid prepared interval tree");
    const auto& node = *prepared_.array<PreparedIntervalNode>(tree);
    const auto* intervals = prepared_.array<PreparedInterval>(node.intervals);
    if (node.intervals.count && end >= intervals[0].start)
        for (uint64_t i = 0; i < node.intervals.count; ++i)
        {
            // Prepared node buckets retain IntervalTree's increasing-start order.
            if (intervals[i].start > start)
                break;
            if (end <= intervals[i].stop)
                result.emplace_back(static_cast<locus::RegionType>(intervals[i].regionType), intervals[i].locus);
        }
    if (start <= node.center)
        queryPrepared(node.left, start, end, result);
    if (end >= node.center)
        queryPrepared(node.right, start, end, result);
}

static vector<AnalyzerBundle> coalesceQueryResults(
    vector<AnalyzerBundle> readAnalyzerBundles, vector<AnalyzerBundle> mateAnalyzerBundles, int32_t readContigId,
    int64_t readStart, int32_t mateContigId, int64_t mateStart)
{
    vector<AnalyzerBundle> bundles = coalesceCommonBundles(readAnalyzerBundles, mateAnalyzerBundles);

    if ((not readAnalyzerBundles.empty()) or (not mateAnalyzerBundles.empty()))
    {
        if (areMatesNearby(readContigId, readStart, mateContigId, mateStart))
        {
            coalesceBundlesForNearbyMates(readAnalyzerBundles, mateAnalyzerBundles, bundles);
        }
        else
        {
            coalesceBundlesForFarawayMates(readAnalyzerBundles, mateAnalyzerBundles, bundles);
        }
    }
    return bundles;
}

static vector<AnalyzerBundle> coalesceStreamingQueryResults(
    vector<AnalyzerBundle> readBundles, vector<AnalyzerBundle> mateBundles, int32_t readContigId,
    int64_t readStart, int32_t mateContigId, int64_t mateStart)
{
    if (readBundles.empty() && mateBundles.empty())
        return readBundles;
    // Reuse the query result allocation when only one side has candidates.
    if (readBundles.empty() || mateBundles.empty())
    {
        const bool readIsEmpty = readBundles.empty();
        auto& bundles = readIsEmpty ? mateBundles : readBundles;
        const bool nearby = areMatesNearby(readContigId, readStart, mateContigId, mateStart);
        auto next = bundles.begin();
        for (auto bundle : bundles)
        {
            if (nearby && bundle.regionType != RegionType::kTarget)
                continue;
            bundle.inputType = nearby
                ? (readIsEmpty ? AnalyzerInputType::kMateOnly : AnalyzerInputType::kReadOnly)
                : AnalyzerInputType::kBothReads;
            *next++ = bundle;
        }
        bundles.erase(next, bundles.end());
        return std::move(bundles);
    }
    if (readBundles.size() == 1 && mateBundles.size() == 1
        && readBundles.front().locusIndex == mateBundles.front().locusIndex)
    {
        readBundles.front().regionType
            = coalesceRegionTypes(readBundles.front().regionType, mateBundles.front().regionType);
        return readBundles;
    }
    // The general path retains its first-match, duplicate and ordering rules.
    return coalesceQueryResults(
        std::move(readBundles), std::move(mateBundles), readContigId, readStart, mateContigId, mateStart);
}

vector<AnalyzerBundle> AnalyzerFinder::query(
    int32_t readContigId, int64_t readStart, int64_t readEnd, int32_t mateContigId, int64_t mateStart,
    int64_t mateEnd) const
{
    return coalesceQueryResults(
        query(readContigId, readStart, readEnd), query(mateContigId, mateStart, mateEnd), readContigId, readStart,
        mateContigId, mateStart);
}

vector<AnalyzerBundle> AnalyzerFinder::queryStreaming(
    int32_t readContigId, int64_t readStart, int64_t readEnd, int32_t mateContigId, int64_t mateStart, int64_t mateEnd)
{
    vector<AnalyzerBundle> readBundles;
    const auto tree = intervalTrees_.find(readContigId);
    const auto root = prepared_.data ? preparedRoot(readContigId) : nullptr;
    if ((root || tree != intervalTrees_.end()) && readStart >= 0)
    {
        auto sweepIt = sweeps_.find(readContigId);
        if (sweepIt == sweeps_.end())
        {
            Sweep sweep;
            if (root)
            {
                sweep.prepared = graphtools::CatalogRange<PreparedInterval>::fromArray(
                    prepared_.array<PreparedInterval>(root->sorted), root->sorted.count);
            }
            else
            {
                // Assign the original tree traversal rank before sorting. This preserves
                // bundle order, including duplicate intervals, for every query.
                tree->second.visit_near(
                    0, std::numeric_limits<std::size_t>::max(),
                    [&](const AnalyzerIntervalTree::interval& interval)
                    {
                        sweep.intervals.push_back(
                            { interval.start, interval.stop, sweep.intervals.size(), interval.value });
                    });
                std::sort(
                    sweep.intervals.begin(), sweep.intervals.end(),
                    [](const SweepInterval& a, const SweepInterval& b) { return a.start < b.start; });
            }
            sweepIt = sweeps_.emplace(readContigId, std::move(sweep)).first;
        }
        auto& sweep = sweepIt->second;
        if (readStart < sweep.position)
            readBundles = query(readContigId, readStart, readEnd);
        else
        {
            sweep.position = readStart;
            while (sweep.next < sweep.size() && sweep.at(sweep.next).start <= std::size_t(readStart))
                sweep.active.push_back(sweep.next++);
            auto& matching = sweep.matching;
            matching.clear();
            auto nextActive = sweep.active.begin();
            for (auto index : sweep.active)
            {
                const auto stop = sweep.at(index).stop;
                if (stop < std::size_t(readStart))
                    continue;
                *nextActive++ = index;
                if (readEnd >= 0 && stop >= std::size_t(readEnd))
                    matching.push_back(index);
            }
            sweep.active.erase(nextActive, sweep.active.end());
            if (matching.size() > 1)
                std::sort(
                    matching.begin(), matching.end(),
                    [&](size_t a, size_t b) { return sweep.at(a).rank < sweep.at(b).rank; });
            readBundles.reserve(matching.size());
            for (auto index : matching)
                readBundles.push_back(sweep.at(index).bundle);
        }
    }
    else
        readBundles = query(readContigId, readStart, readEnd);
    // Never expire mate or offtarget information using the current read position.
    return coalesceStreamingQueryResults(
        std::move(readBundles), query(mateContigId, mateStart, mateEnd), readContigId, readStart, mateContigId,
        mateStart);
}

}
