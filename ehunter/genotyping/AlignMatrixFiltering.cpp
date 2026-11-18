//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
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
// This file differs from ExpansionHunter v5.0.0: filtering compacts surviving read pairs in the flat alignment matrix.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "genotyping/AlignMatrixFiltering.hh"

#include <algorithm>
#include <utility>

namespace ehunter
{
namespace strgt
{

void filter(AlignMatrix& aligns)
{
    assert(aligns.numReads() % 2 == 0);
    int numStrReads = 0;
    int numStrReadsWithIndels = 0;
    int longestStrSize = 0;
    int longestStrSizeWithoutIndels = 0;
    for (int readIndex = 0; readIndex != aligns.numReads(); ++readIndex)
    {
        const auto& topAlign = aligns.getBestAlign(readIndex);

        longestStrSize = std::max(longestStrSize, topAlign.numMotifs());

        if (topAlign.type() != StrAlign::Type::kOutside)
        {
            ++numStrReads;
            if (topAlign.numIndels() > 0)
            {
                ++numStrReadsWithIndels;
            }
            else
            {
                longestStrSizeWithoutIndels = std::max(longestStrSizeWithoutIndels, topAlign.numMotifs());
            }
        }
    }

    const int maxOutliers = std::max(1, static_cast<int>(numStrReads * 0.20));
    if (numStrReadsWithIndels == 0 or maxOutliers < numStrReadsWithIndels)
    {
        return;
    }

    const double increase = (static_cast<double>(longestStrSize) - longestStrSizeWithoutIndels) / longestStrSize;

    if (increase < 0.1)
    {
        return;
    }

    // Read pairs produce an even row count, but public remove() can leave an
    // odd count in a release build. For odd counts, use v5.0.0's reverse
    // traversal, including its handling of a surviving unpaired tail.
    if (aligns.numReads() % 2 != 0)
    {
        for (int readIndex = aligns.numReads() - 1; readIndex != -1; --readIndex)
        {
            if (aligns.getBestAlign(readIndex).numIndels() == 0)
                continue;
            if (readIndex % 2 == 0)
            {
                aligns.remove(readIndex + 1);
                aligns.remove(readIndex);
            }
            else
            {
                aligns.remove(readIndex);
                aligns.remove(readIndex - 1);
                --readIndex;
            }
        }
        return;
    }

    // Retain surviving pairs in their original order, compacting rows forward
    // in place: a destination never passes its source.
    int retained = 0;
    size_t written = 0;
    const int numReads = aligns.numReads();
    for (int readIndex = 0; readIndex != numReads; readIndex += 2)
    {
        if (aligns.getBestAlign(readIndex).numIndels() != 0 || aligns.getBestAlign(readIndex + 1).numIndels() != 0)
        {
            continue;
        }

        const uint32_t begin = aligns.rowBegin_[readIndex];
        const uint32_t end = aligns.rowBegin_[readIndex + 2];
        const uint32_t mateBegin = aligns.rowBegin_[readIndex + 1];
        if (written != begin)
            std::copy(aligns.cells_.begin() + begin, aligns.cells_.begin() + end, aligns.cells_.begin() + written);
        aligns.rowBegin_[retained] = static_cast<uint32_t>(written);
        aligns.rowBegin_[retained + 1] = static_cast<uint32_t>(written + (mateBegin - begin));
        aligns.bestIndex_[retained] = aligns.bestIndex_[readIndex];
        aligns.bestIndex_[retained + 1] = aligns.bestIndex_[readIndex + 1];
        written += end - begin;
        retained += 2;
    }
    aligns.cells_.erase(aligns.cells_.begin() + written, aligns.cells_.end());
    aligns.bestIndex_.erase(aligns.bestIndex_.begin() + retained, aligns.bestIndex_.end());
    if (retained == 0)
        aligns.rowBegin_.clear();
    else
    {
        aligns.rowBegin_[retained] = static_cast<uint32_t>(written);
        aligns.rowBegin_.erase(aligns.rowBegin_.begin() + retained + 1, aligns.rowBegin_.end());
    }
}

CountTable countAligns(StrAlign::Type alignType, const AlignMatrix& aligns)
{
    CountTable countTable;
    for (int readIndex = 0; readIndex != aligns.numReads(); ++readIndex)
    {
        const auto& align = aligns.getBestAlign(readIndex);
        if (align.type() == alignType)
        {
            countTable.incrementCountOf(align.numMotifs());
        }
    }

    return countTable;
}

}
}
