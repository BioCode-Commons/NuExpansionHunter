// Modified for NuExpansionHunter (2026): batch mate retrieval while preserving query results.
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

#include "sample/MateExtractor.hh"

#include <map>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "core/HtsHelpers.hh"

namespace ehunter
{

using boost::optional;
using std::string;

namespace htshelpers
{
MateExtractor::MateExtractor(const string& htsFilePath, const std::string& htsReferencePath)
    : htsFilePath_(htsFilePath)
    , htsReferencePath_(htsReferencePath)
    , contigInfo_({})
{
    openFile();
    loadHeader();
    loadIndex();
    htsAlignmentPtr_ = bam_init1();
}

MateExtractor::~MateExtractor()
{
    bam_destroy1(htsAlignmentPtr_);
    htsAlignmentPtr_ = nullptr;

    hts_idx_destroy(htsIndexPtr_);
    htsIndexPtr_ = nullptr;

    bam_hdr_destroy(htsHeaderPtr_);
    htsHeaderPtr_ = nullptr;

    sam_close(htsFilePtr_);
    htsFilePtr_ = nullptr;
}

void MateExtractor::openFile()
{
    htsFilePtr_ = sam_open(htsFilePath_.c_str(), "r");

    if (!htsFilePtr_)
    {
        throw std::runtime_error("Failed to read BAM file " + htsFilePath_);
    }

    // Required step for parsing of some CRAMs
    if (hts_set_fai_filename(htsFilePtr_, htsReferencePath_.c_str()) != 0)
    {
        throw std::runtime_error("Failed to set index of: " + htsReferencePath_);
    }
}

void MateExtractor::loadHeader()
{
    htsHeaderPtr_ = sam_hdr_read(htsFilePtr_);

    if (!htsHeaderPtr_)
    {
        throw std::runtime_error("Failed to read header of " + htsFilePath_);
    }

    contigInfo_ = htshelpers::decodeContigInfo(htsHeaderPtr_);
}

void MateExtractor::loadIndex()
{
    htsIndexPtr_ = sam_index_load(htsFilePtr_, htsFilePath_.c_str());

    if (!htsIndexPtr_)
    {
        throw std::runtime_error("Failed to read index of " + htsFilePath_);
    }
}

optional<Read> MateExtractor::extractMate(
    const Read& read, const LinearAlignmentStats& alignmentStats, LinearAlignmentStats& mateStats)
{
    const int32_t searchRegionContigIndex
        = alignmentStats.isMateMapped ? alignmentStats.mateChromId : alignmentStats.chromId;

    const int32_t searchRegionStart = alignmentStats.isMateMapped ? alignmentStats.matePos : alignmentStats.pos;
    const int32_t searchRegionEnd = searchRegionStart + 1;

    hts_itr_t* htsRegionPtr_
        = sam_itr_queryi(htsIndexPtr_, searchRegionContigIndex, searchRegionStart, searchRegionEnd);

    if (!htsRegionPtr_)
    {
        const string& contigName = contigInfo_.getContigName(searchRegionContigIndex);
        const string regionEncoding
            = contigName + ":" + std::to_string(searchRegionStart) + "-" + std::to_string(searchRegionEnd);

        throw std::logic_error("Unable to jump to " + regionEncoding + " to recover a mate");
    }

    while (sam_itr_next(htsFilePtr_, htsRegionPtr_, htsAlignmentPtr_) >= 0)
    {
        const bool isSecondaryAlignment = htsAlignmentPtr_->core.flag & BAM_FSECONDARY;
        const bool isSupplementaryAlignment = htsAlignmentPtr_->core.flag & BAM_FSUPPLEMENTARY;
        const bool isPrimaryAlignment = !(isSecondaryAlignment || isSupplementaryAlignment);
        if (!isPrimaryAlignment)
        {
            continue;
        }

        Read putativeMate = htshelpers::decodeRead(htsAlignmentPtr_);

        const bool belongToSameFragment = read.fragmentId() == putativeMate.fragmentId();
        const bool formProperPair = read.mateNumber() != putativeMate.mateNumber();

        if (belongToSameFragment && formProperPair)
        {
            mateStats = decodeAlignmentStats(htsAlignmentPtr_);
            hts_itr_destroy(htsRegionPtr_);
            return putativeMate;
        }
    }
    hts_itr_destroy(htsRegionPtr_);

    return optional<Read>();
}

void MateExtractor::extractMates(std::vector<MateQuery>& queries)
{
    using Coordinate = std::pair<int32_t, hts_pos_t>;
    using PendingQueries = std::multimap<Coordinate, MateQuery*>;
    PendingQueries pending;
    std::unordered_multimap<string, PendingQueries::iterator> byFragment;
    std::vector<string> regions;

    for (auto& query : queries)
    {
        const auto& stats = *query.alignmentStats;
        const int32_t tid = stats.isMateMapped ? stats.mateChromId : stats.chromId;
        const int32_t start = stats.isMateMapped ? stats.matePos : stats.pos;

        // Negative values denote HTSlib's non-coordinate queries. Keep their
        // single-query semantics; they cannot be represented as genomic points.
        if (tid < 0 || start < 0)
        {
            query.mate = extractMate(*query.read, stats, query.mateStats);
            continue;
        }

        // Validate each query with the same API as extractMate. In particular,
        // the multi-region API must not silently skip an unindexed reference.
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> check(
            sam_itr_queryi(htsIndexPtr_, tid, start, static_cast<hts_pos_t>(start) + 1), &hts_itr_destroy);
        if (!check)
        {
            const string region = contigInfo_.getContigName(tid) + ":" + std::to_string(start) + "-"
                + std::to_string(static_cast<hts_pos_t>(start) + 1);
            throw std::logic_error("Unable to jump to " + region + " to recover a mate");
        }

        const auto entry = pending.emplace(Coordinate(tid, start), &query);
        byFragment.emplace(query.read->fragmentId(), entry);
    }

    if (pending.empty())
    {
        return;
    }

    // Query the exact union of the original one-base intervals. Quoting the
    // reference name also supports names containing colons.
    for (auto it = pending.begin(); it != pending.end(); it = pending.upper_bound(it->first))
    {
        const auto& coordinate = it->first;
        const string position = std::to_string(coordinate.second + 1);
        regions.push_back("{" + contigInfo_.getContigName(coordinate.first) + "}:" + position + "-" + position);
    }
    std::vector<char*> regionPointers;
    for (auto& region : regions)
    {
        regionPointers.push_back(&region[0]);
    }
    // Single-region CRAM queries set a decoder-level range, even when the
    // iterator is only used for validation. Clear it before the multi-region
    // iterator applies its own interval union. This uses the public HTSlib API.
    {
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> reset(
            sam_itr_queryi(htsIndexPtr_, HTS_IDX_START, 0, 0), &hts_itr_destroy);
        if (!reset)
        {
            throw std::logic_error("Unable to reset input before recovering mates");
        }
    }
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_regarray(htsIndexPtr_, htsHeaderPtr_, regionPointers.data(), regionPointers.size()), &hts_itr_destroy);
    if (!iterator)
    {
        throw std::logic_error("Unable to query regions to recover mates");
    }

    while (!pending.empty() && sam_itr_next(htsFilePtr_, iterator.get(), htsAlignmentPtr_) >= 0)
    {
        if (!isPrimaryAlignment(htsAlignmentPtr_))
        {
            continue;
        }
        const int32_t tid = htsAlignmentPtr_->core.tid;
        const hts_pos_t begin = htsAlignmentPtr_->core.pos;
        const hts_pos_t end = bam_endpos(htsAlignmentPtr_);
        const auto active = pending.lower_bound(Coordinate(tid, begin));
        if (active == pending.end() || active->first.first != tid || active->first.second >= end)
        {
            continue;
        }

        // Decode every primary record that an unfinished original query would
        // inspect, preserving sequence/quality decoding and input validation.
        Read mate = decodeRead(htsAlignmentPtr_);
        const auto matches = byFragment.equal_range(mate.fragmentId());
        auto it = matches.first;
        while (it != matches.second)
        {
            const auto position = it->second->first;
            auto& query = *it->second->second;
            if (position.first == tid && position.second >= begin && position.second < end
                && query.read->mateNumber() != mate.mateNumber())
            {
                // HTSlib returns records in file order, so this is the same
                // first matching primary mate as the individual query.
                query.mate = mate;
                query.mateStats = decodeAlignmentStats(htsAlignmentPtr_);
                pending.erase(it->second);
                it = byFragment.erase(it);
            }
            else
            {
                ++it;
            }
        }
    }
}

}

}
