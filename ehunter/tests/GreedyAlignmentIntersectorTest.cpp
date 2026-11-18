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
// This file differs from ExpansionHunter v5.0.0: tests cover path-based intersection and canonical-fold equivalence.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "alignment/GreedyAlignmentIntersector.hh"
#include "alignment/OperationsOnAlignments.hh"

#include <list>
#include <vector>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignment.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

#include "io/RegionGraph.hh"
#include "locus/LocusSpecification.hh"

using graphtools::Graph;
using graphtools::Path;

using namespace ehunter;

TEST(CanonicalAlignmentTransfer, PreservesInputsAndShrunkIntersection)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));
    const auto first = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]2[2M]", &graph);
    const auto second = decodeGraphAlignment(1, "1[2M]1[3M]2[1M]", &graph);
    const std::list<graphtools::GraphAlignment> input = { first, second };
    const auto expected = decodeGraphAlignment(1, "1[3S2M]1[3M]2[1M1S]", &graph);
    EXPECT_EQ(expected, computeCanonicalAlignment(input));
    EXPECT_EQ(first, input.front());
    EXPECT_EQ(second, input.back());
    EXPECT_EQ(first, computeCanonicalAlignment({ first }));
    const auto disjoint = decodeGraphAlignment(0, "4[4M]", &graph);
    EXPECT_EQ(first, computeCanonicalAlignment({ first, disjoint }));
}

TEST(IntersectingPaths, IntersectionStartsAtLoopNode_PathsIntersected)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));

    {
        const auto firstAlignment = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]2[2M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(1, "1[2M]1[3M]2[1M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment.path());
        const auto intersection = alignmentIntersector.intersect();

        const auto expectedAlignment = decodeGraphAlignment(1, "1[3S2M]1[3M]2[1M1S]", &graph);
        EXPECT_EQ(expectedAlignment, *intersection);
    }

    {
        const auto firstAlignment = decodeGraphAlignment(1, "1[2M]1[3M]2[1M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]2[3M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment.path());
        const auto intersection = alignmentIntersector.intersect();

        const auto expectedAlignment = decodeGraphAlignment(1, "1[2M]1[3M]2[1M]", &graph);
        EXPECT_EQ(expectedAlignment, *intersection);
    }

    {
        const auto firstAlignment = decodeGraphAlignment(2, "0[2M]1[1M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(2, "1[1M]1[3M]2[2M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment.path());
        const auto intersection = alignmentIntersector.intersect();

        const auto expectedAlignment = decodeGraphAlignment(0, "1[2S1M]", &graph);
        EXPECT_EQ(expectedAlignment, *intersection);
    }
}

TEST(IntersectingPaths, IntersectionStartsAtRegularNode_PathsIntersected)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));

    const auto firstAlignment = decodeGraphAlignment(1, "0[3M]1[3M]1[3M]1[3M]1[3M]", &graph);
    const auto secondAlignment = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]1[2M]", &graph);

    GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment.path());
    const auto intersection = alignmentIntersector.intersect();

    const auto expectedAlignment = decodeGraphAlignment(2, "0[1S2M]1[3M]1[3M]1[2M4S]", &graph);
    EXPECT_EQ(expectedAlignment, *intersection);
}

TEST(IntersectingPaths, NonintersectingPaths_HandleledProperly)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));

    {
        const auto firstAlignment = decodeGraphAlignment(1, "0[3M]1[3M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(0, "3[3M]4[2M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment.path());

        EXPECT_FALSE(alignmentIntersector.intersect());
    }

    {
        const auto firstAlignment = decodeGraphAlignment(1, "0[3M]1[3M]2[2M]", &graph);
        const auto secondAlignment = decodeGraphAlignment(2, "2[4M]3[3M]4[2M]", &graph);

        GreedyAlignmentIntersector alignmentIntersector(firstAlignment, secondAlignment.path());

        EXPECT_FALSE(alignmentIntersector.intersect());
    }
}

TEST(CanonicalAlignmentFold, RepeatedPathsGiveTheReferenceResult)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CAG)*CAACAG(CCG)*CCTT"));
    // Same path, different per-node alignments; unchanged and changing steps;
    // an element with a disjoint path (empty intersection) at the end.
    const auto a = decodeGraphAlignment(2, "0[2M]1[3M]1[3M]2[2M]", &graph);
    const auto aMismatch = decodeGraphAlignment(2, "0[2M]1[1M1X1M]1[3M]2[2M]", &graph);
    const auto b = decodeGraphAlignment(1, "1[2M]1[3M]2[1M]", &graph);
    const auto bMismatch = decodeGraphAlignment(1, "1[2M]1[1X2M]2[1M]", &graph);
    const auto disjoint = decodeGraphAlignment(0, "4[4M]", &graph);
    const std::vector<std::list<graphtools::GraphAlignment>> lists = {
        { a, a, aMismatch, aMismatch },
        { a, aMismatch, b, bMismatch, b },
        { b, bMismatch, a, aMismatch },
        { a, aMismatch, aMismatch, disjoint },
        { a, aMismatch, b, b, b, disjoint },
    };
    for (const auto& list : lists)
        EXPECT_EQ(computeCanonicalAlignmentReference(list), computeCanonicalAlignment(list));
    // Ten times ten ties on two paths.
    std::list<graphtools::GraphAlignment> ties;
    for (int i = 0; i < 50; ++i)
        ties.push_back(i % 2 ? aMismatch : a);
    for (int i = 0; i < 50; ++i)
        ties.push_back(i % 2 ? bMismatch : b);
    EXPECT_EQ(computeCanonicalAlignmentReference(ties), computeCanonicalAlignment(ties));
}

// The first sorted element of a prefix x suffix product is not built from the
// smallest suffix: merging the seed's match into a suffix that starts with a
// match changes the order of two suffixes with the same path.
TEST(CanonicalAlignmentFactoring, SeamReordersSuffixCandidates)
{
    using graphtools::Alignment;
    const Alignment seed(0, "10M");
    const Alignment deletionFirst(10, "1D5M");
    const Alignment matchFirst(10, "2M1D3M");
    EXPECT_TRUE(matchFirst < deletionFirst);
    EXPECT_TRUE(graphtools::mergeAlignments(seed, deletionFirst) < graphtools::mergeAlignments(seed, matchFirst));
}

// The intersection does not factor into a prefix part and a suffix part: one
// extra repeat copy on the right of the second alignment clips the left end of
// the canonical alignment, although the two alignments agree on the left.
TEST(CanonicalAlignmentFactoring, RightCopyNumberChangesLeftEnd)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    const auto moreCopies = decodeGraphAlignment(0, "1[1M]1[1M]1[1M]2[3M]", &graph);
    const auto fewerCopies = decodeGraphAlignment(0, "1[1M]1[1M]2[1I3M]", &graph);
    ASSERT_TRUE(moreCopies < fewerCopies);
    const auto canonical = computeCanonicalAlignment({ moreCopies, fewerCopies });
    EXPECT_EQ(computeCanonicalAlignmentReference({ moreCopies, fewerCopies }), canonical);
    EXPECT_FALSE(moreCopies == canonical);
    // The left end lost its first base to a soft clip, or the first copy altogether.
    EXPECT_TRUE(canonical[0].numClipped() > 0 || canonical.path().numNodes() < moreCopies.path().numNodes());
}
