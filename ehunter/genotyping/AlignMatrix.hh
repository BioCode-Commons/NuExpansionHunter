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

#pragma once

#include <cstdint>
#include <iosfwd>
#include <limits>
#include <vector>

#include <boost/container/vector.hpp>

#include "graphalign/GraphAlignment.hh"

#include "genotyping/StrAlign.hh"

namespace ehunter
{
namespace strgt
{

// One row of alignment summaries per read, for allele sizes 0, 1, ... The last
// element of a row stands for every larger allele size. Rows are stored back
// to back in one array; a read's best alignment is the row element at the
// read's own motif count.
class AlignMatrix
{
public:
    explicit AlignMatrix(int strNode);
    int numReads() const { return static_cast<int>(bestIndex_.size()); }
    void add(const graphtools::GraphAlignment& read, const graphtools::GraphAlignment& mate);
    void remove(int readIndex);
    StrAlign getAlign(int readIndex, int alleleSize) const;
    StrAlign getBestAlign(int readIndex) const { return cells_[rowBegin_[readIndex] + bestIndex_[readIndex]]; }
    // Row access: rowSize(i) summaries starting at rowData(i).
    int rowSize(int readIndex) const { return static_cast<int>(rowBegin_[readIndex + 1] - rowBegin_[readIndex]); }
    const StrAlign* rowData(int readIndex) const { return cells_.data() + rowBegin_[readIndex]; }
    int getMaxMotifCount() const;

    friend void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs);
    friend void filter(AlignMatrix& aligns);

private:
    // Grow by half instead of doubling: millions of loci keep their rows until
    // genotyping, so unused capacity is retained memory.
    template <typename T>
    using Array = boost::container::vector<
        T, void, typename boost::container::vector_options<
                     boost::container::growth_factor<boost::container::growth_factor_50>>::type>;

    void add(const graphtools::GraphAlignment& graphAlign);
    // Appends a row of cells [first, cells_.size()) with its best index.
    void closeRow(size_t first, size_t bestIndex);

    int strNode_;
    ConsistentAlignmentCalculator alignmentCalculator_;
    Array<StrAlign> cells_;
    // Row i is cells_[rowBegin_[i], rowBegin_[i + 1]); empty when there are no rows.
    Array<uint32_t> rowBegin_;
    Array<uint16_t> bestIndex_;
};

std::ostream& operator<<(std::ostream& out, const AlignMatrix& matrix);

void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs);

}
}
