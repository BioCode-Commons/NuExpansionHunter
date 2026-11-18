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
// This file differs from ExpansionHunter v5.0.0: alignment rows share a flat cell array with row offsets and best-alignment indexes.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "genotyping/AlignMatrix.hh"

#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <numeric>

using graphtools::GraphAlignment;
using std::vector;

namespace ehunter
{
namespace strgt
{

AlignMatrix::AlignMatrix(int strNode)
    : strNode_(strNode)
    , alignmentCalculator_(strNode)
{
}

StrAlign AlignMatrix::getAlign(int readIndex, int alleleSize) const
{
    if (readIndex >= numReads())
    {
        throw std::runtime_error("Encountered invalid alignment matrix index " + std::to_string(readIndex));
    }

    const auto begin = rowBegin_[readIndex];
    const auto end = rowBegin_[readIndex + 1];
    if (alleleSize < static_cast<int>(end - begin))
    {
        return cells_[begin + alleleSize];
    }
    else
    {
        assert(end != begin);
        return cells_[end - 1];
    }
}

void AlignMatrix::closeRow(size_t first, size_t bestIndex)
{
    if (cells_.size() > std::numeric_limits<uint32_t>::max())
        throw std::length_error("Alignment matrix exceeds its cell index range");
    if (bestIndex > std::numeric_limits<uint16_t>::max())
        throw std::length_error("Alignment matrix row exceeds its best-alignment index range");
    // assign() rather than push_back() on the empty array: GCC 11 reports a
    // false stringop-overflow in Boost's reallocation for the latter.
    if (rowBegin_.empty())
        rowBegin_.assign(1, static_cast<uint32_t>(first));
    rowBegin_.push_back(static_cast<uint32_t>(cells_.size()));
    bestIndex_.push_back(static_cast<uint16_t>(bestIndex));
}

// Contract [K-25]
// Requires: a row's best index fits 16 bits and the locus holds fewer than 2^32 cells (closeRow throws otherwise).
// Ensures:  under the abstraction that maps (cells, row offsets, best index) to v5.0.0's vector of rows, every
//           operation gives the same abstract state and return value as v5's AlignMatrix.
// Reference: verification/reference/genotyping/AlignMatrix (v5.0.0 text).
void AlignMatrix::add(const GraphAlignment& read, const GraphAlignment& mate)
{
    const int numMotifsInRead = std::count(read.path().begin(), read.path().end(), strNode_);
    const int numMorifsInMate = std::count(mate.path().begin(), mate.path().end(), strNode_);

    if (numMotifsInRead != 0 || numMorifsInMate != 0)
    {
        add(read);
        add(mate);
    }
}

void AlignMatrix::add(const graphtools::GraphAlignment& graphAlign)
{
    const int numMotifsInAlign = std::count(graphAlign.path().begin(), graphAlign.path().end(), strNode_);
    StrAlign alignToMostConsistentAllele = alignmentCalculator_.findConsistentAlignment(numMotifsInAlign, graphAlign);

    const size_t first = cells_.size();
    try
    {
        for (int numMotifs = numMotifsInAlign - 1; numMotifs != -1; --numMotifs)
        {
            StrAlign align = alignmentCalculator_.findConsistentAlignment(numMotifs, graphAlign);
            cells_.emplace_back(align);
        }
        std::reverse(cells_.begin() + first, cells_.end());
        cells_.emplace_back(alignToMostConsistentAllele);
        StrAlign previousAlign = cells_.back();
        for (int numMotifs = numMotifsInAlign + 1;; ++numMotifs)
        {
            StrAlign align = alignmentCalculator_.findConsistentAlignment(numMotifs, graphAlign);
            if (align.type() == previousAlign.type() && align.score() == previousAlign.score())
            {
                break;
            }
            cells_.emplace_back(align);
            previousAlign = align;
        }
        closeRow(first, static_cast<size_t>(numMotifsInAlign));
    }
    catch (...)
    {
        cells_.erase(cells_.begin() + first, cells_.end());
        throw;
    }
}

int AlignMatrix::getMaxMotifCount() const
{
    int maxMotifCount = 0;
    for (int readIndex = 0; readIndex != numReads(); ++readIndex)
    {
        const int maxMotifsInAlign = rowSize(readIndex) - 1;
        if (maxMotifsInAlign > maxMotifCount)
        {
            maxMotifCount = maxMotifsInAlign;
        }
    }

    return maxMotifCount;
}

std::ostream& operator<<(std::ostream& out, const AlignMatrix& matrix)
{
    const unsigned alignCount(matrix.numReads());
    std::vector<unsigned> alignIndexes(alignCount);
    std::iota(alignIndexes.begin(), alignIndexes.end(), 0);

    // Sort alignIndexes to canonicalize the alignmatrix output
    std::sort(
        alignIndexes.begin(), alignIndexes.end(),
        [&matrix](const unsigned ai, const unsigned bi) -> bool
        {
            const StrAlign* a = matrix.rowData(ai);
            const StrAlign* b = matrix.rowData(bi);
            const unsigned aSize = matrix.rowSize(ai);
            const unsigned bSize = matrix.rowSize(bi);
            if (aSize < bSize)
            {
                return true;
            }
            if (bSize < aSize)
            {
                return false;
            }
            for (unsigned i(0); i < aSize; ++i)
            {
                if (a[i] < b[i])
                {
                    return true;
                }
                if (b[i] < a[i])
                {
                    return false;
                }
            }
            return false;
        });

    auto dumpStrAlign = [&out](const StrAlign& align)
    {
        out << "(";
        switch (align.type())
        {
        case StrAlign::Type::kOutside:
            out << "O";
            break;
        case StrAlign::Type::kInRepeat:
            out << "I";
            break;
        case StrAlign::Type::kSpanning:
            out << "S";
            break;
        case StrAlign::Type::kFlanking:
            out << "F";
            break;
        }
        out << "," << align.numMotifs() << "," << align.score() << "), ";
    };

    for (unsigned origAlignIndex(0); origAlignIndex < alignCount; ++origAlignIndex)
    {
        const unsigned alignIndex(alignIndexes[origAlignIndex]);
        const StrAlign* aligns = matrix.rowData(alignIndex);
        for (int i = 0; i != matrix.rowSize(alignIndex); ++i)
        {
            dumpStrAlign(aligns[i]);
        }
        out << "\n";
    }

    return out;
}

void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs)
{
    // Find highest-scoring in-repeat read
    const int longIrrLowerBound = static_cast<int>(0.90 * maxMotifsInRead);

    int topIrrIndex = -1;
    int topIrrScore = -1;
    for (int readIndex = 0; readIndex != alignMatrix.numReads(); ++readIndex)
    {
        const StrAlign align = alignMatrix.getBestAlign(readIndex);
        const bool isLongIrr = align.type() == StrAlign::Type::kInRepeat && align.numMotifs() >= longIrrLowerBound;

        if (isLongIrr && align.score() > topIrrScore)
        {
            topIrrScore = align.score();
            topIrrIndex = readIndex;
        }
    }

    if (topIrrIndex == -1)
    {
        return;
    }

    const vector<StrAlign> irrAligns(
        alignMatrix.rowData(topIrrIndex), alignMatrix.rowData(topIrrIndex) + alignMatrix.rowSize(topIrrIndex));
    const size_t irrBest = alignMatrix.bestIndex_[topIrrIndex];

    for (int newIrrIndex = 0; newIrrIndex != 2 * numIrrPairs; ++newIrrIndex)
    {
        const size_t first = alignMatrix.cells_.size();
        alignMatrix.cells_.insert(alignMatrix.cells_.end(), irrAligns.begin(), irrAligns.end());
        alignMatrix.closeRow(first, irrBest);
    }
}

void AlignMatrix::remove(int readIndex)
{
    assert(readIndex < numReads());
    const uint32_t begin = rowBegin_[readIndex];
    const uint32_t length = rowBegin_[readIndex + 1] - begin;
    cells_.erase(cells_.begin() + begin, cells_.begin() + begin + length);
    rowBegin_.erase(rowBegin_.begin() + readIndex + 1);
    for (size_t row = readIndex + 1; row < rowBegin_.size(); ++row)
        rowBegin_[row] -= length;
    bestIndex_.erase(bestIndex_.begin() + readIndex);
    if (bestIndex_.empty())
        rowBegin_.clear();
}

}
}
