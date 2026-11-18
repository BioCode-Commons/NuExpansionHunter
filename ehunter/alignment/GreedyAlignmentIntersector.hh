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
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: alignment intersection takes a path as its second operand and counts repeat visits without temporary lists.
// Modifications are licensed under the Apache License, Version 2.0.
//

#pragma once

#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"

namespace ehunter
{

// Intersects an alignment with a path. The second operand is a path, not an alignment: a step of the
// canonical fold depends only on the running alignment and the path of the next alignment (K-21).
// For a GraphAlignment a, intersecting with a.path() is v5.0.0's intersection with a, because every
// alignment has one node alignment per path node.
class GreedyAlignmentIntersector
{
public:
    GreedyAlignmentIntersector(const graphtools::GraphAlignment& firstAlignment, const graphtools::Path& secondPath)
        : firstAlignment_(firstAlignment)
        , firstPath_(firstAlignment.path())
        , secondPath_(secondPath)
    {
        initialize();
    }

    boost::optional<graphtools::GraphAlignment> intersect();

private:
    void initialize();
    bool tryAdvancingIndexesToCommonNode();
    bool checkIfAlignmentEndReached(int firstPathIndex, int secondPathIndex);
    bool checkIfCommonNodeIsLoop() const;
    void advanceIndexesToMatchRemainingIterations();
    void advanceIndexesToLastCommonNode();
    void computeIntersectionEnds();
    bool checkIfIntersectionIsConsistent() const;
    boost::optional<graphtools::GraphAlignment> softclipFirstAlignmentToIntersection() const;

    const graphtools::GraphAlignment& firstAlignment_;
    const graphtools::Path& firstPath_;
    const graphtools::Path& secondPath_;

    int nodeIndexOfIntersectionStartOnFirstPath_;
    int nodeIndexOfIntersectionStartOnSecondPath_;
    int nodeIndexOfIntersectionEndOnFirstPath_ = -1;
    int nodeIndexOfIntersectionEndOnSecondPath_ = -1;
    int intersectionStart_ = -1;
    int intersectionEnd_ = -1;
};

}
