//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
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

#include "graphcore/Graph.hh"

#include <limits>
#include <memory>

#include "graphalign/GappedAligner.hh"
#include "graphcore/GraphBuilders.hh"
#include "gtest/gtest.h"

using std::set;
using std::string;
using std::vector;

using namespace graphtools;

TEST(GraphConstruction, TypicalNodeCount_GraphConstructed)
{
    Graph graph(3);
    ASSERT_EQ(3ul, graph.numNodes());
}

TEST(NodeNameManipulation, TypicalNode_NodeNameSet)
{
    Graph graph(3);
    graph.setNodeName(1, "LF");
    ASSERT_EQ("LF", graph.nodeName(1));
}

TEST(NodeNameManipulation, NonexistingNode_ExceptionRaised)
{
    Graph graph(1);
    ASSERT_ANY_THROW(graph.setNodeName(1, "LF"));
    ASSERT_ANY_THROW(graph.nodeName(1));
}

TEST(NodeSequenceManipulation, TypicalSequence_SequenceSet)
{
    Graph graph(3);
    graph.setNodeSeq(1, "ATT");
    EXPECT_EQ("ATT", graph.nodeSeq(1));
}

TEST(NodeSequenceManipulation, DegenerateSequence_SequenceExpansionObtained)
{
    Graph graph(3);
    graph.setNodeSeq(1, "WC");
    vector<string> expected_expansion = { "AC", "TC" };
    EXPECT_EQ(expected_expansion, graph.nodeSeqExpansion(1));
}

TEST(NodeSequenceManipulation, NonexistingNode_ExceptionRaised)
{
    Graph graph(3);
    EXPECT_ANY_THROW(graph.setNodeSeq(4, "ATT"));
    EXPECT_ANY_THROW(graph.nodeSeq(4));
    EXPECT_ANY_THROW(graph.nodeSeqExpansion(4));
}

TEST(NodeSequenceManipulation, EmptySequence_ExceptionRaised)
{
    Graph graph(3);
    EXPECT_ANY_THROW(graph.setNodeSeq(1, ""));
}

TEST(AddingEdges, TypicalEdge_EdgeAdded)
{
    Graph graph(3);
    graph.addEdge(0, 1);
    graph.addEdge(0, 0);
    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(0, 0));
}

TEST(AddingEdges, EdgeBreakingTopologicalOrder_ExceptionRaised)
{
    Graph graph(3);
    EXPECT_ANY_THROW(graph.addEdge(2, 1));
}

TEST(AddingEdges, EdgesBetweenNonexistingNodes_ExceptionRaised)
{
    Graph graph(4);
    EXPECT_ANY_THROW(graph.addEdge(-1, 2));
    EXPECT_ANY_THROW(graph.addEdge(1, 4));
    EXPECT_ANY_THROW(graph.addEdge(4, 5));
}

TEST(AddingEdges, EdgesThatAlreadyExist_ExceptionRaised)
{
    Graph graph(4);
    graph.addEdge(1, 2);
    EXPECT_ANY_THROW(graph.addEdge(1, 2));
}

TEST(CheckingIfEdgesExist, EdgesBetweenNonexistingNodes_ExceptionRaised)
{
    Graph graph(4);
    EXPECT_ANY_THROW(graph.hasEdge(-1, 2));
    EXPECT_ANY_THROW(graph.hasEdge(1, 4));
    EXPECT_ANY_THROW(graph.hasEdge(4, 5));
}

TEST(EdgeLabelManipulation, TypicalEdges_EdgesLabeled)
{
    Graph graph(4);
    graph.addEdge(0, 2);
    graph.addLabelToEdge(0, 2, "ref");
    graph.addLabelToEdge(0, 2, "alt");
    Labels expected_labels = { "ref", "alt" };
    ASSERT_EQ(expected_labels, graph.edgeLabels(0, 2));
}

TEST(EdgeLabelManipulation, NonexistingEdges_ExceptionRaised)
{
    Graph graph(4);
    EXPECT_ANY_THROW(graph.addLabelToEdge(0, 1, "ref"));
    EXPECT_ANY_THROW(graph.addLabelToEdge(0, 4, "ref"));
    EXPECT_ANY_THROW(graph.edgeLabels(0, 1));
}

TEST(GettingNodeNeighbors, TypicalNode_SuccessorsFound)
{
    Graph graph(4);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(0, 3);
    graph.addEdge(2, 3);

    const set<NodeId> expected_successors = { 1, 2, 3 };
    ASSERT_EQ(expected_successors, graph.successors(0));
    ASSERT_TRUE(graph.successors(1).empty());
}

TEST(GettingNodeNeighbors, LoopAtNode_SuccessorsFound)
{
    Graph graph(4);
    graph.addEdge(0, 0);
    graph.addEdge(0, 1);

    const set<NodeId> expected_predecessors = { 0, 1 };
    ASSERT_EQ(expected_predecessors, graph.successors(0));
}

TEST(GettingNodeNeighbors, TypicalNode_PredecessorsFound)
{
    Graph graph(4);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(0, 3);
    graph.addEdge(2, 3);

    const set<NodeId> expected_predecessors = { 0, 2 };
    ASSERT_EQ(expected_predecessors, graph.predecessors(3));
}

TEST(GettingNodeNeighbors, NeighborsOfNonexistingNode_ExceptionRaised)
{
    Graph graph(4);
    EXPECT_ANY_THROW(graph.successors(4));
    EXPECT_ANY_THROW(graph.predecessors(-1));
}

namespace
{
Graph makeTraversalGraph()
{
    Graph graph(6);
    for (NodeId node = 0; node < graph.numNodes(); ++node)
        graph.setNodeSeq(node, "ACG");
    // Deliberately insert out of order: traversal order must remain sorted.
    for (const auto& edge : vector<NodeIdPair> { { 0, 4 }, { 0, 0 }, { 0, 2 }, { 0, 1 },
                                              { 2, 4 }, { 1, 4 }, { 4, 4 } })
        graph.addEdge(edge.first, edge.second);
    return graph;
}

Graph mapGraph(const Graph& graph, std::shared_ptr<string>& storage)
{
    BlobBuilder builder;
    const auto record = graph.freeze(builder);
    storage = std::make_shared<string>(std::move(builder.bytes));
    return Graph::mapped(MappedBytes(storage->data(), storage->size(), storage), record);
}
}

TEST(BorrowedGraphTraversal, OwnedAndMappedNeighborsPreserveOrderInBothDirections)
{
    const auto owned = makeTraversalGraph();
    std::shared_ptr<string> storage;
    const auto mapped = mapGraph(owned, storage);
    const vector<vector<NodeId>> expectedSuccessors { { 0, 1, 2, 4 }, { 4 }, { 4 }, {}, { 4 }, {} };
    const vector<vector<NodeId>> expectedPredecessors { { 0 }, { 0 }, { 0 }, {}, { 0, 1, 2, 4 }, {} };

    for (const auto* graph : { &owned, &mapped })
    {
        SCOPED_TRACE(graph->isMapped());
        const ReverseGraph reverse(*graph);
        for (NodeId node = 0; node < graph->numNodes(); ++node)
        {
            SCOPED_TRACE(node);
            EXPECT_EQ(expectedSuccessors[node], vector<NodeId>(graph->borrowedSuccessors(node)));
            EXPECT_EQ(expectedPredecessors[node], vector<NodeId>(graph->borrowedPredecessors(node)));
            EXPECT_EQ(expectedSuccessors[node], vector<NodeId>(graph->successors(node)));
            EXPECT_EQ(expectedPredecessors[node], vector<NodeId>(graph->predecessors(node)));
            EXPECT_EQ(expectedPredecessors[node], vector<NodeId>(reverse.borrowedSuccessors(node)));
            EXPECT_EQ(expectedSuccessors[node], vector<NodeId>(reverse.borrowedPredecessors(node)));
        }
        for (NodeId invalid : { NodeId(graph->numNodes()), std::numeric_limits<NodeId>::max() })
        {
            EXPECT_THROW(graph->borrowedSuccessors(invalid), std::logic_error);
            EXPECT_THROW(graph->borrowedPredecessors(invalid), std::logic_error);
            EXPECT_THROW(reverse.borrowedSuccessors(invalid), std::logic_error);
            EXPECT_THROW(reverse.borrowedPredecessors(invalid), std::logic_error);
        }
    }
}

TEST(BorrowedGraphTraversal, MappedBorrowedRangesAndIteratorsDoNotRetainStorage)
{
    std::shared_ptr<string> storage;
    const auto graph = mapGraph(makeTraversalGraph(), storage);
    const auto ownersBeforeTraversal = storage.use_count();
    const auto successors = graph.borrowedSuccessors(0);
    const auto predecessors = graph.borrowedPredecessors(4);
    auto next = successors.begin();
    const auto end = successors.end();
    auto previous = predecessors.rbegin();
    EXPECT_EQ(ownersBeforeTraversal, storage.use_count());
    EXPECT_EQ(0u, *next);
    EXPECT_NE(end, ++next);
    EXPECT_EQ(1u, *next);
    EXPECT_EQ(4u, *previous);
    EXPECT_EQ(2u, *++previous);
    EXPECT_EQ(ownersBeforeTraversal, storage.use_count());
}

TEST(BorrowedGraphTraversal, CustomAndDerivedGraphTypesKeepTheirOwnAdjacencySemantics)
{
    struct CustomGraph
    {
        vector<NodeId> next { 3, 1 }, previous { 2, 0 };
        const vector<NodeId>& successors(NodeId) const { return next; }
        const vector<NodeId>& predecessors(NodeId) const { return previous; }
    };
    struct DerivedGraph : Graph
    {
        DerivedGraph() : Graph(5) { addEdge(0, 1); }
        vector<NodeId> successors(NodeId) const { return { 4, 2 }; }
        vector<NodeId> predecessors(NodeId) const { return { 3, 0 }; }
    };
    const CustomGraph custom;
    const DerivedGraph derived;
    // The fallback also preserves reference-returning APIs; no borrowed API is required.
    EXPECT_EQ(&custom.next, &detail::borrowedSuccessors(custom, 0));
    EXPECT_EQ(&custom.previous, &detail::borrowedPredecessors(custom, 0));
    EXPECT_EQ((vector<NodeId> { 4, 2 }), detail::borrowedSuccessors(derived, 0));
    EXPECT_EQ((vector<NodeId> { 3, 0 }), detail::borrowedPredecessors(derived, 0));
}

TEST(MappedGraphEdges, ExistingMissingLoopAndInvalidEdgesMatchOwnedGraph)
{
    const auto owned = makeTraversalGraph();
    std::shared_ptr<string> storage;
    const auto mapped = mapGraph(owned, storage);
    const set<NodeIdPair> expected { { 0, 0 }, { 0, 1 }, { 0, 2 }, { 0, 4 }, { 1, 4 }, { 2, 4 }, { 4, 4 } };
    for (const auto* graph : { &owned, &mapped })
    {
        SCOPED_TRACE(graph->isMapped());
        for (NodeId from = 0; from < graph->numNodes(); ++from)
            for (NodeId to = 0; to < graph->numNodes(); ++to)
                EXPECT_EQ(expected.count({ from, to }) != 0, graph->hasEdge(from, to));
        for (NodeId invalid : { NodeId(graph->numNodes()), std::numeric_limits<NodeId>::max() })
        {
            EXPECT_THROW(graph->hasEdge(invalid, 0), std::logic_error);
            EXPECT_THROW(graph->hasEdge(0, invalid), std::logic_error);
            EXPECT_THROW(graph->hasEdge(invalid, invalid), std::logic_error);
        }
    }
}

TEST(MappedGraphEdges, MalformedAdjacencySlicesAreRejectedBeforeTraversal)
{
    BlobBuilder builder;
    const auto record = makeTraversalGraph().freeze(builder);
    const MappedBytes bytes(builder.bytes.data(), builder.bytes.size());
    const auto root = *bytes.array<MappedGraph>(record);
    const auto originalNode = *bytes.array<MappedNode>(root.nodes);
    const vector<Slice> invalidSlices { { uint64_t(builder.bytes.size()) + 8, 1 },
                                       { root.nodes.offset, std::numeric_limits<uint64_t>::max() },
                                       { 1, 1 } };
    for (const auto& invalid : invalidSlices)
    {
        SCOPED_TRACE(invalid.offset);
        auto node = originalNode;
        node.successors = invalid;
        node.predecessors = invalid;
        builder.replace(root.nodes.offset, node);
        const auto graph = Graph::mapped(MappedBytes(builder.bytes.data(), builder.bytes.size()), record);
        EXPECT_THROW(graph.successors(0), std::runtime_error);
        EXPECT_THROW(graph.predecessors(0), std::runtime_error);
        EXPECT_THROW(graph.borrowedSuccessors(0), std::runtime_error);
        EXPECT_THROW(graph.borrowedPredecessors(0), std::runtime_error);
        EXPECT_THROW(graph.hasEdge(0, 0), std::runtime_error);
    }
}

TEST(MappedGraphLifetime, OwningNeighborsAndTemporaryRangeIteratorsOutliveSourceGraph)
{
    CatalogRange<NodeId> successors;
    auto previous = CatalogRange<NodeId>().begin();
    auto previousEnd = CatalogRange<NodeId>().end();
    std::weak_ptr<string> storageLifetime;
    {
        std::shared_ptr<string> storage;
        const auto graph = mapGraph(makeTraversalGraph(), storage);
        storageLifetime = storage;
        successors = graph.successors(0);
        previous = graph.predecessors(4).begin();
        previousEnd = graph.predecessors(4).end();
    }
    ASSERT_FALSE(storageLifetime.expired());
    EXPECT_EQ((vector<NodeId> { 0, 1, 2, 4 }), vector<NodeId>(successors));
    successors = CatalogRange<NodeId>();
    ASSERT_FALSE(storageLifetime.expired()); // The iterators retain their temporary range's storage.
    EXPECT_EQ((vector<NodeId> { 0, 1, 2, 4 }), vector<NodeId>(previous, previousEnd));
    previous = CatalogRange<NodeId>().begin();
    previousEnd = CatalogRange<NodeId>().end();
    EXPECT_TRUE(storageLifetime.expired());
}

TEST(MappedGraphAlignment, RepeatAndDeletionAlignmentsMatchOwnedGraphIncludingOrder)
{
    const vector<Graph> graphs { makeStrGraph("ATATTA", "CG", "TATTT"),
                                 makeDeletionGraph("ACGTAC", "GG", "TACGTA"),
                                 makeStrGraph("AAG", "CGG", "CTT") };
    const vector<vector<string>> queries { { "TTCGCTAT", "ATTCGCTAT", "TACGCGTAT", "NNNNNN" },
                                          { "ACGTACGGTACGTA", "ACGTACTACGTA", "ACGTACTTCGTA" },
                                          { "GCGGC", "AAGCGGCGGCTT" } };
    for (size_t graphIndex = 0; graphIndex < graphs.size(); ++graphIndex)
    {
        SCOPED_TRACE(graphIndex);
        const auto& owned = graphs[graphIndex];
        std::shared_ptr<string> storage;
        const auto mapped = mapGraph(owned, storage);
        const GappedGraphAligner ownedAligner(&owned, 3, 2, 0), mappedAligner(&mapped, 3, 2, 0);
        for (auto type : { AlignerType::PATH_ALIGNER, AlignerType::DAG_ALIGNER })
        {
            SCOPED_TRACE(static_cast<int>(type));
            AlignerSelector ownedSelector(type), mappedSelector(type);
            for (const auto& query : queries[graphIndex])
            {
                SCOPED_TRACE(query);
                const auto expected = ownedAligner.align(query, ownedSelector);
                const auto observed = mappedAligner.align(query, mappedSelector);
                // Rejection of the deliberately mutated/ambiguous queries is also an output to preserve.
                if (query != "NNNNNN" && query != "ACGTACTTCGTA")
                    ASSERT_FALSE(expected.empty());
                ASSERT_EQ(expected.size(), observed.size());
                auto actual = observed.begin();
                for (const auto& alignment : expected)
                {
                    // Path equality includes the graph address, so compare every value explicitly.
                    EXPECT_EQ(alignment.path().encode(), actual->path().encode());
                    EXPECT_EQ(alignment.alignments(), actual->alignments());
                    EXPECT_EQ(alignment.generateCigar(), actual->generateCigar());
                    ++actual;
                }
            }
        }
    }
}
