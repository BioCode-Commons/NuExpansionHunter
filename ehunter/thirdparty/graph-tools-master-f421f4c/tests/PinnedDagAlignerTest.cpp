//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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

#include <stdio.h>
#include <string.h>

#include "graphalign/DagAlignerAffine.hh"
#include "graphalign/GappedAligner.hh"
#include "graphalign/PinnedDagAligner.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"

#include "gtest/gtest.h"

using std::string;
using namespace testing;
using namespace graphtools;
using namespace graphalign;
using namespace graphalign::dagAligner;

template <typename T> std::string toString(const T& obj)
{
    std::stringstream ss;
    ss << obj;
    return ss.str();
}

TEST(SimpleGraph, AlignSuffix)
{
    Graph graph = makeSwapGraph("AAAA", "C", "T", "GGGG");
    graph.addEdge(1, 1);
    Path seed(
        &graph, 1,
        {
            0,
        },
        3);
    PinnedDagAligner dag_pinned_aligner(1, -1, 0, -2);
    int32_t top_dag_score = INT32_MIN;
    std::list<std::pair<Path, Alignment>> res = dag_pinned_aligner.prefixAlign(seed, "ACGG", 8, top_dag_score);
    EXPECT_EQ(4, top_dag_score);
    EXPECT_EQ("(0@1)-(1)-(3@2)", toString(res.front().first));
    EXPECT_EQ("4M", res.front().second.generateCigar());
}

TEST(SimpleGraph, AlignPrefix)
{
    Graph graph = makeSwapGraph("AAAA", "C", "T", "GGGG");
    graph.addEdge(1, 1);
    Path seed(&graph, 0, { 2, 3 }, 2);
    PinnedDagAligner dag_pinned_aligner(1, -1, 0, -2);
    int32_t top_dag_score = INT32_MIN;
    std::list<std::pair<Path, Alignment>> res = dag_pinned_aligner.suffixAlign(seed, "AA", 8, top_dag_score);
    EXPECT_EQ(2, top_dag_score);
    EXPECT_EQ("(0@2)-(2)-(3@2)", toString(res.front().first));
    EXPECT_EQ("2M", res.front().second.generateCigar());
}

TEST(SingleNodeExtension, KeepsEntireRemainingNodeInBothDirections)
{
    Graph graph(1);
    graph.setNodeSeq(0, "ACGTACGT");
    graph.addEdge(0, 0);
    PinnedDagAligner aligner(5, -4, -8, -2);
    int score = INT32_MIN;
    // The extraction distance limits graph discovery, not bases within the
    // final node. A query may consume more than this distance.
    auto forward = aligner.prefixAlign(Path(&graph, 0, { 0 }, 2), "GTACGT", 2, score);
    ASSERT_EQ(1u, forward.size());
    EXPECT_EQ(30, score);
    EXPECT_EQ("(0@0)-(0@8)", toString(forward.front().first));
    EXPECT_EQ("6M", forward.front().second.generateCigar());
    auto reverse = aligner.suffixAlign(Path(&graph, 6, { 0 }, 8), "ACGTAC", 2, score);
    ASSERT_EQ(1u, reverse.size());
    EXPECT_EQ(30, score);
    EXPECT_EQ("(0@0)-(0@8)", toString(reverse.front().first));
    EXPECT_EQ("6M", reverse.front().second.generateCigar());
}

TEST(SingleNodeExtension, PreservesRepeatBoundaryCrossing)
{
    Graph graph(1);
    graph.setNodeSeq(0, "ACGT");
    graph.addEdge(0, 0);
    PinnedDagAligner aligner(5, -4, -8, -2);
    int score = INT32_MIN;
    auto forward = aligner.prefixAlign(Path(&graph, 0, { 0 }, 4), "ACGTAC", 6, score);
    ASSERT_EQ(1u, forward.size());
    EXPECT_EQ(30, score);
    EXPECT_EQ("ACGTACGTAC", forward.front().first.seq());
    EXPECT_EQ("6M", forward.front().second.generateCigar());
    auto reverse = aligner.suffixAlign(Path(&graph, 0, { 0 }, 4), "GTACGT", 6, score);
    ASSERT_EQ(1u, reverse.size());
    EXPECT_EQ(30, score);
    EXPECT_EQ("GTACGTACGT", reverse.front().first.seq());
    EXPECT_EQ("6M", reverse.front().second.generateCigar());
}

TEST(EdgeMapLinearMetadata, OptionalRootEdgeHasIdenticalConnectivity)
{
    for (const int length : { 1, 2, 15, 16, 17, 31, 32, 33, 150, 1000 })
    {
        const EdgeMap compact({ { length, length } }, { 7 });
        const EdgeMap general({ { -1, 0 }, { length, length } }, { 7 });
        EXPECT_EQ(toString(general), toString(compact));
        EXPECT_TRUE(compact.hasSingleLinearNode(length));
        EXPECT_FALSE(compact.hasSingleLinearNode(length + 1));
        EXPECT_FALSE(general.hasSingleLinearNode(length));
        for (int offset = 0; offset < length; ++offset)
        {
            EXPECT_EQ(7u, compact.getNodeId(offset));
            ASSERT_EQ(1, std::distance(compact.prevNodesBegin(offset), compact.prevNodesEnd(offset)));
            EXPECT_EQ(offset - 1, *compact.prevNodesBegin(offset));
        }
    }
}
