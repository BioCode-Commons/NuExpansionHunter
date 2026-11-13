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

TEST(SingleNodeExtension, TwoInternalMismatchesKeepUniqueAlignmentAndMergeAdjacentOperations)
{
    Graph graph(1);
    graph.setNodeSeq(0, string(155, 'A'));
    PinnedDagAligner aligner(5, -4, -8, -2);
    for (const bool adjacent : { false, true })
    {
        string query(150, 'A');
        query[1] = query[adjacent ? 2 : 148] = 'C';
        const string expectedCigar = adjacent ? "1M2X147M" : "1M1X146M1X1M";
        int score = INT32_MIN;
        const auto forward = aligner.prefixAlign(Path(&graph, 0, { 0 }, 0), query, 150, score);
        ASSERT_EQ(1u, forward.size());
        EXPECT_EQ(732, score);
        EXPECT_EQ(Path(&graph, 0, { 0 }, 150), forward.front().first);
        EXPECT_EQ(expectedCigar, forward.front().second.generateCigar());
        const auto reverse = aligner.suffixAlign(Path(&graph, 155, { 0 }, 155), query, 150, score);
        ASSERT_EQ(1u, reverse.size());
        EXPECT_EQ(732, score);
        EXPECT_EQ(Path(&graph, 5, { 0 }, 155), reverse.front().first);
        EXPECT_EQ(expectedCigar, reverse.front().second.generateCigar());
    }
}

TEST(SingleNodeExtension, TwoMismatchesKeepAllCost18DeletionTies)
{
    Graph graph(1);
    string target(21, 'A');
    target[8] = 'C';
    graph.setNodeSeq(0, target);
    string query = target;
    query.erase(3, 5);
    PinnedDagAligner aligner(5, -4, -8, -2);
    int score = INT32_MIN;
    const auto results = aligner.prefixAlign(Path(&graph, 0, { 0 }, 0), query, 16, score);
    ASSERT_EQ(10u, results.size());
    EXPECT_EQ(62, score);
    EXPECT_EQ("1D1M1D1M3D14M", results.front().second.generateCigar());
    EXPECT_EQ(Path(&graph, 0, { 0 }, 21), results.front().first);
    EXPECT_EQ("3M1X4M1X7M", results.back().second.generateCigar());
    EXPECT_EQ(Path(&graph, 0, { 0 }, 16), results.back().first);
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

TEST(SingleNodeExtension, PreservesOperationsAndReferenceStartInBothDirections)
{
    Graph graph(3);
    graph.setNodeSeq(0, "AAAAAAAA");
    const string nodeSequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    graph.setNodeSeq(1, nodeSequence);
    graph.setNodeSeq(2, "CCCCCCCC");
    graph.addEdge(0, 1);
    graph.addEdge(1, 1);
    graph.addEdge(1, 2);
    PinnedDagAligner aligner(5, -4, -8, -2);
    const Path forwardSeed(&graph, 2, { 0, 1 }, 4);
    const Path reverseSeed(&graph, 60, { 1, 2 }, 5);

    for (const int mismatch : { -1, 0, 15, 31 })
    {
        string forwardQuery = nodeSequence.substr(4, 32);
        string reverseQuery = nodeSequence.substr(28, 32);
        string cigar = "32M";
        if (mismatch >= 0)
        {
            // An unrecognized byte is a mismatch, including on signed-char
            // platforms. Endpoint mismatches must not produce zero-length M.
            forwardQuery[mismatch] = reverseQuery[mismatch] = static_cast<char>(0xff);
            cigar = mismatch ? std::to_string(mismatch) + "M" : "";
            cigar += "1X";
            if (mismatch != 31)
                cigar += std::to_string(31 - mismatch) + "M";
        }
        int score = INT32_MIN;
        // The query can consume more than the extraction distance: the entire
        // last node remains part of the discovered graph in both directions.
        const auto forward = aligner.prefixAlign(forwardSeed, forwardQuery, 1, score);
        ASSERT_EQ(1u, forward.size());
        EXPECT_EQ(mismatch < 0 ? 160 : 151, score);
        EXPECT_EQ(Path(&graph, 2, { 0, 1 }, mismatch == 31 ? 35 : 36), forward.front().first);
        EXPECT_EQ(10u, forward.front().second.referenceStart());
        EXPECT_EQ(mismatch == 31 ? "31M1S" : cigar, forward.front().second.generateCigar());

        const auto reverse = aligner.suffixAlign(reverseSeed, reverseQuery, 1, score);
        ASSERT_EQ(1u, reverse.size());
        EXPECT_EQ(mismatch < 0 ? 160 : 151, score);
        EXPECT_EQ(Path(&graph, mismatch == 0 ? 29 : 28, { 1, 2 }, 5), reverse.front().first);
        EXPECT_EQ(0u, reverse.front().second.referenceStart());
        EXPECT_EQ(mismatch == 0 ? "1S31M" : cigar, reverse.front().second.generateCigar());
    }
}

TEST(SingleNodeExtension, NonDefaultGapCostsKeepOptimalGappedAlignment)
{
    Graph graph(1);
    string target(64, 'A');
    target[16] = 'C';
    graph.setNodeSeq(0, target);
    PinnedDagAligner aligner(5, -4, -1, -1);
    int score = INT32_MIN;
    const auto result = aligner.prefixAlign(Path(&graph, 0, { 0 }, 0), string(32, 'A'), 32, score);
    ASSERT_FALSE(result.empty());
    EXPECT_EQ(158, score);
    EXPECT_EQ("16M1D16M", result.front().second.generateCigar());
}

TEST(SingleNodeExtension, TerminalNodesAllowDiscoveryPastTheirEnds)
{
    Graph graph(2);
    graph.setNodeSeq(0, string(64, 'A'));
    graph.setNodeSeq(1, string(64, 'C'));
    graph.addEdge(0, 1);
    PinnedDagAligner aligner(5, -4, -8, -2);
    int score = INT32_MIN;
    string forwardQuery(32, 'C');
    forwardQuery[0] = 'A';
    const auto forward = aligner.prefixAlign(Path(&graph, 5, { 0, 1 }, 32), forwardQuery, 200, score);
    ASSERT_EQ(1u, forward.size());
    EXPECT_EQ(151, score);
    EXPECT_EQ(Path(&graph, 5, { 0, 1 }, 64), forward.front().first);
    EXPECT_EQ(91u, forward.front().second.referenceStart());
    EXPECT_EQ("1X31M", forward.front().second.generateCigar());

    string reverseQuery(32, 'A');
    reverseQuery[31] = 'C';
    const auto reverse = aligner.suffixAlign(Path(&graph, 32, { 0, 1 }, 5), reverseQuery, 200, score);
    ASSERT_EQ(1u, reverse.size());
    EXPECT_EQ(151, score);
    EXPECT_EQ(Path(&graph, 0, { 0, 1 }, 5), reverse.front().first);
    EXPECT_EQ(0u, reverse.front().second.referenceStart());
    EXPECT_EQ("31M1X", reverse.front().second.generateCigar());
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

namespace
{
string repeatShortcutReference(size_t length, uint32_t state)
{
    string result(length, 'A');
    for (char& base : result)
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        base = "ACGT"[state % 4];
    }
    return result;
}
}

// These 150 bp Path/score/CIGAR expectations were independently checked against
// unmodified Illumina v5.0.0, including its empty first-unrolled-node handling.
TEST(UniqueUngappedExtension, PreservesRepeatExitAndEmptySeedBoundaries)
{
    Graph graph(3);
    const string left = repeatShortcutReference(300, 13619);
    const string right = repeatShortcutReference(300, 7619);
    graph.setNodeSeq(0, left);
    graph.setNodeSeq(1, "CA");
    graph.setNodeSeq(2, right);
    graph.addEdge(0, 1);
    graph.addEdge(1, 1);
    graph.addEdge(1, 2);
    string repeat;
    for (int i = 0; i < 50; ++i)
        repeat += "CA";
    PinnedDagAligner aligner(5, -4, -8, -2);
    for (int boundary : { 0, 1, 2 })
        for (bool reverse : { false, true })
            for (int mismatch : { -1, 75 })
            {
                SCOPED_TRACE(toString(boundary) + "/" + toString(reverse) + "/" + toString(mismatch));
                string query = reverse ? left + repeat + string("CA").substr(0, boundary)
                                       : string("CA").substr(boundary) + repeat + right;
                query = reverse ? query.substr(query.size() - 150) : query.substr(0, 150);
                if (mismatch >= 0)
                    query[mismatch] = 'N';
                const Path seed(&graph, boundary, { 1 }, boundary);
                int score = INT32_MIN;
                const auto results = reverse ? aligner.suffixAlign(seed, query, 160, score)
                                             : aligner.prefixAlign(seed, query, 160, score);
                ASSERT_EQ(1u, results.size());
                std::vector<NodeId> nodes(51, 1);
                if (reverse)
                    nodes.insert(nodes.begin(), 0);
                else
                    nodes.push_back(2);
                EXPECT_EQ(Path(&graph, reverse ? 250 + boundary : boundary, nodes,
                               reverse ? boundary : 48 + boundary), results.front().first);
                EXPECT_EQ(mismatch < 0 ? 750 : 741, score);
                EXPECT_EQ(mismatch < 0 ? "150M" : "75M1X74M", results.front().second.generateCigar());
                EXPECT_EQ(0u, results.front().second.referenceStart());
            }
}

TEST(UniqueUngappedExtension, KeepsOrderedTiesAndReplacesWorseTiesWithAnExactPath)
{
    const string query = repeatShortcutReference(150, 459);
    for (bool exactBranch : { false, true })
    {
        Graph graph(exactBranch ? 5 : 4);
        graph.setNodeSeq(0, "A");
        const NodeId terminal = exactBranch ? 4 : 3;
        graph.setNodeSeq(terminal, "G");
        for (NodeId node = 1; node < terminal; ++node)
        {
            string sequence = query;
            if (node < 3)
            {
                const size_t mismatch = node == 1 ? 75 : 80;
                sequence[mismatch] = sequence[mismatch] == 'A' ? 'C' : 'A';
            }
            graph.setNodeSeq(node, sequence);
            graph.addEdge(0, node);
            graph.addEdge(node, terminal);
        }
        PinnedDagAligner aligner(5, -4, -8, -2);
        int score = INT32_MIN;
        const auto results = aligner.prefixAlign(Path(&graph, 1, { 0 }, 1), query, 160, score);
        ASSERT_EQ(exactBranch ? 1u : 2u, results.size());
        EXPECT_EQ(exactBranch ? 750 : 741, score);
        NodeId expectedNode = exactBranch ? 3 : 1;
        for (const auto& result : results)
        {
            EXPECT_EQ(Path(&graph, 1, { 0, expectedNode }, 150), result.first);
            EXPECT_EQ(exactBranch ? "150M" : expectedNode == 1 ? "75M1X74M" : "80M1X69M",
                      result.second.generateCigar());
            EXPECT_EQ(0u, result.second.referenceStart());
            ++expectedNode;
        }
    }
}

TEST(EdgeMapNodeSpans, PreservesDuplicateIncomingEdgesAndTheirOrder)
{
    const std::vector<std::pair<int, int>> edges = {
        { 1, 2 }, { 0, 2 }, { -1, 2 }, { 1, 2 }, { 3, 4 }, { 2, 4 }, { 6, 6 }
    };
    // Arrays produced by the original v5 EdgeMap constructor.
    const EdgeMap expected({ 7, 7, 4, 4, 9, 9 }, { 0, 1, 2, 6, 7, 9, 10 },
                           { -1, 0, 1, 0, -1, 1, 2, 3, 2, 4 });
    EXPECT_EQ(toString(expected), toString(EdgeMap(edges, { 7, 4, 9 })));
    auto explicitRoot = edges;
    explicitRoot.insert(explicitRoot.begin(), { -1, 0 });
    EXPECT_EQ(toString(expected), toString(EdgeMap(explicitRoot, { 7, 4, 9 })));
}

TEST(EdgeMapNodeSpans, RetainsGeneralConstructionForUnusualOffsets)
{
    const EdgeMap negativePredecessor({ { -2, 2 }, { 4, 4 } }, { 5, 9 });
    const EdgeMap expectedNegative({ 5, 5, 9, 9 }, { 0, 1, 2, 3, 4 }, { -1, 0, -2, 2 });
    EXPECT_EQ(toString(expectedNegative), toString(negativePredecessor));
    const EdgeMap duplicateRoot({ { -1, 0 }, { -1, 0 }, { 4, 4 } }, { 7, 9 });
    const EdgeMap expectedRoot({ 7, 9, 9, 9 }, { 0, 1, 2, 3, 4 }, { -1, -1, 1, 2 });
    EXPECT_EQ(toString(expectedRoot), toString(duplicateRoot));
}

TEST(PreparedGraphReuse, SequenceMutationAssignmentAndReusedAddressesKeepFreshResults)
{
    PinnedDagAligner aligner(5, -4, -8, -2);
    string query(150, 'A');
    query[45] = query[90] = 'C';
    auto verify = [&](Graph& graph, bool changed)
    {
        int score = INT32_MIN;
        const auto results = aligner.prefixAlign(Path(&graph, 0, { 0 }, 0), query, 160, score);
        ASSERT_EQ(1u, results.size());
        EXPECT_EQ(changed ? 723 : 732, score);
        EXPECT_EQ(changed ? "45M1X14M1X29M1X59M" : "45M1X44M1X59M",
                  results.front().second.generateCigar());
        EXPECT_EQ(Path(&graph, 0, { 0 }, 150), results.front().first);
    };
    Graph graph(1);
    graph.setNodeSeq(0, string(180, 'A'));
    verify(graph, false);
    verify(graph, false);
    Graph saved = graph;
    string changed(180, 'A');
    changed[60] = 'C';
    graph.setNodeSeq(0, changed);
    verify(graph, true);
    verify(saved, false);
    graph = saved;
    verify(graph, false);
    alignas(Graph) unsigned char storage[sizeof(Graph)];
    for (bool modified : { false, true })
    {
        Graph* placed = new (storage) Graph(1);
        placed->setNodeSeq(0, modified ? changed : string(180, 'A'));
        verify(*placed, modified);
        placed->~Graph();
    }
}

TEST(PreparedGraphReuse, MutatingACopiedGraphDoesNotChangeTheOriginalPreparation)
{
    Graph graph(3);
    graph.setNodeSeq(0, string(75, 'A'));
    graph.setNodeSeq(1, string(100, 'A'));
    string alternate(100, 'A');
    alternate[15] = alternate[45] = 'C';
    graph.setNodeSeq(2, alternate);
    graph.addEdge(0, 1);
    string query(150, 'A');
    for (int mismatch : { 20, 45, 90, 120 })
        query[mismatch] = 'C';
    PinnedDagAligner aligner(5, -4, -8, -2);
    auto verify = [&](Graph& current, bool connected)
    {
        int score = INT32_MIN;
        const auto results = aligner.prefixAlign(Path(&current, 0, { 0 }, 0), query, 160, score);
        ASSERT_EQ(1u, results.size());
        EXPECT_EQ(connected ? 732 : 714, score);
        EXPECT_EQ(connected ? "20M1X24M1X104M" : "20M1X24M1X44M1X29M1X29M",
                  results.front().second.generateCigar());
        EXPECT_EQ(Path(&current, 0, { 0, NodeId(connected ? 2 : 1) }, 75), results.front().first);
    };
    verify(graph, false);
    Graph copied = graph;
    verify(copied, false);
    copied.addEdge(0, 2);
    verify(copied, true);
    verify(graph, false);
    graph = copied;
    verify(graph, true);
}

// Repeated high-cost calls may choose a wider initial DP search. A subsequent
// cost-18 call must keep the original optimum and ordered result, too.
TEST(PinnedBudgetGuidance, KeepsExactResultsWhenRecentAlignmentCostsChange)
{
    Graph graph(2);
    const string sequence = repeatShortcutReference(220, 81723);
    graph.setNodeSeq(0, sequence.substr(0, 75));
    graph.setNodeSeq(1, sequence.substr(75));
    graph.addEdge(0, 1);
    for (bool reverse : { false, true })
    {
        PinnedDagAligner aligner(5, -4, -8, -2);
        const Path seed = reverse ? Path(&graph, 105, { 1 }, 105) : Path(&graph, 0, { 0 }, 0);
        const string exact = sequence.substr(reverse ? 30 : 0, 150);
        for (int mismatches : { 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 2, 2 })
        {
            SCOPED_TRACE(toString(reverse) + "/" + toString(mismatches));
            string query = exact;
            for (int k = 0; k < mismatches; ++k)
                query[20 + 30 * k] = 'N';
            int score = INT32_MIN;
            const auto results = reverse ? aligner.suffixAlign(seed, query, 180, score)
                                         : aligner.prefixAlign(seed, query, 180, score);
            ASSERT_EQ(1u, results.size());
            EXPECT_EQ(750 - 9 * mismatches, score);
            EXPECT_EQ(Path(&graph, reverse ? 30 : 0, { 0, 1 }, reverse ? 105 : 75), results.front().first);
            EXPECT_EQ(mismatches == 4 ? "20M1X29M1X29M1X29M1X39M" : "20M1X29M1X99M",
                      results.front().second.generateCigar());
            EXPECT_EQ(0u, results.front().second.referenceStart());
        }
    }
}
