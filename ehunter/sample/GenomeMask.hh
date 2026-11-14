//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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

#include "graphutils/MappedData.hh"
#include "sample/PreparedIntervals.hh"
#include <string>
#include <vector>

namespace ehunter
{

/// Binned Bit-mask covering the genome
/// Used to mark target regions for fast read screening
class GenomeMask
{
public:
    // Accept an unsigned coordinate so computing a bin never shifts a negative
    // signed value. Query callers reject negative coordinates first.
    static constexpr uint64_t binIndex(uint64_t position) { return position >> 10; }

    void addRegion(int32_t contigId, int64_t start, int64_t stop);
    bool query(int32_t contigId, int64_t pos) const
    {
        if (contigId < 0 || pos < 0)
            return false;
        if (prepared_.data)
        {
            if (uint64_t(contigId) >= roots_.count)
                return false;
            const auto mask = prepared_.array<PreparedContigIntervals>(roots_)[contigId].mask;
            const auto bin = binIndex(static_cast<uint64_t>(pos));
            return bin / 64 < mask.count
                && (prepared_.array<uint64_t>(mask)[bin / 64] & (uint64_t(1) << (bin % 64)));
        }
        if (contigId >= static_cast<int>(mask_.size()))
            return false;
        const contigMask& cmask = mask_[contigId];
        const size_t bin = binIndex(static_cast<uint64_t>(pos));
        return bin < cmask.size() && cmask[bin];
    }
    void usePrepared(graphtools::MappedBytes bytes, graphtools::Slice roots)
    {
        prepared_ = bytes;
        roots_ = roots;
    }

private:
    using contigMask = std::vector<bool>;
    std::vector<contigMask> mask_;
    graphtools::MappedBytes prepared_;
    graphtools::Slice roots_ { 0, 0 };
};

}
