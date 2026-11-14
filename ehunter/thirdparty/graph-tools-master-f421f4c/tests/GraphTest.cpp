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

TEST(NodeSequenceManipulation, LengthOnlyAccessPreservesOwnedMappedAndReverseSemantics)
{
    Graph graph(2);
    graph.setNodeSeq(0, "ACGTN");
    graph.setNodeSeq(1, string(1024, 'A'));
    graph.addEdge(0, 1);
    BlobBuilder blob;
    const auto record = graph.freeze(blob);
    MappedBytes bytes(blob.bytes.data(), blob.bytes.size());
    Graph mapped = Graph::mapped(bytes, record);
    ReverseGraph reverse(graph), mappedReverse(mapped);
    for (NodeId node = 0; node < graph.numNodes(); ++node)
    {
        EXPECT_EQ(graph.nodeSeq(node).size(), graph.nodeSeqLength(node));
        EXPECT_EQ(reverse.nodeSeq(node).size(), reverse.nodeSeqLength(node));
        EXPECT_EQ(mapped.nodeSeq(node).size(), mapped.nodeSeqLength(node));
        EXPECT_EQ(mappedReverse.nodeSeq(node).size(), mappedReverse.nodeSeqLength(node));
    }
    EXPECT_THROW(graph.nodeSeqLength(2), std::logic_error);
    EXPECT_THROW(mappedReverse.nodeSeqLength(2), std::logic_error);
    graph.setNodeSeq(0, "C");
    EXPECT_EQ(reverse.nodeSeqLength(0), 1u);

    const auto nodes = bytes.array<MappedGraph>(record)->nodes;
    auto corrupt = bytes.array<MappedNode>(nodes)[0];
    corrupt.sequence = { bytes.size, 1 };
    blob.replace(nodes.offset, corrupt);
    EXPECT_THROW(mapped.nodeSeq(0), std::runtime_error);
    EXPECT_THROW(mapped.nodeSeqLength(0), std::runtime_error);
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

TEST(GettingNodeNeighbors, BorrowedRangesPreserveOwnedMappedAndReverseTraversal)
{
    Graph graph(4);
    for (NodeId node = 0; node < graph.numNodes(); ++node)
        graph.setNodeSeq(node, "ACG");
    graph.addEdge(0, 0);
    graph.addEdge(0, 1);
    graph.addEdge(0, 3);
    graph.addEdge(1, 3);
    graph.addEdge(3, 3);
    BlobBuilder blob;
    const auto record = graph.freeze(blob);
    Graph mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);

    for (const Graph* current : { &graph, &mapped })
    {
        ReverseGraph reverse(*current);
        for (NodeId node = 0; node < current->numNodes(); ++node)
        {
            EXPECT_EQ(current->successors(node), current->successorsBorrowed(node));
            EXPECT_EQ(current->predecessors(node), current->predecessorsBorrowed(node));
            EXPECT_EQ(reverse.successors(node), reverse.successorsBorrowed(node));
            EXPECT_EQ(reverse.predecessors(node), reverse.predecessorsBorrowed(node));
            EXPECT_EQ(current->predecessorsBorrowed(node), reverse.successorsBorrowed(node));
            EXPECT_EQ(current->successorsBorrowed(node), reverse.predecessorsBorrowed(node));
            const auto borrowed = current->successorsBorrowed(node);
            vector<NodeId> reversed(borrowed.rbegin(), borrowed.rend());
            vector<NodeId> expected = current->successors(node);
            std::reverse(expected.begin(), expected.end());
            EXPECT_EQ(expected, reversed);
        }
        EXPECT_THROW(current->successorsBorrowed(4), std::logic_error);
        EXPECT_THROW(current->predecessorsBorrowed(-1), std::logic_error);
        EXPECT_THROW(reverse.successorsBorrowed(4), std::logic_error);
        EXPECT_THROW(reverse.predecessorsBorrowed(-1), std::logic_error);
    }

    const auto nodes = MappedBytes(blob.bytes.data(), blob.bytes.size()).array<MappedGraph>(record)->nodes;
    auto corrupt = MappedBytes(blob.bytes.data(), blob.bytes.size()).array<MappedNode>(nodes)[0];
    corrupt.successors = { blob.bytes.size(), 1 };
    corrupt.predecessors = { blob.bytes.size(), 1 };
    blob.replace(nodes.offset, corrupt);
    EXPECT_THROW(mapped.successors(0), std::runtime_error);
    EXPECT_THROW(mapped.successorsBorrowed(0), std::runtime_error);
    EXPECT_THROW(mapped.predecessors(0), std::runtime_error);
    EXPECT_THROW(mapped.predecessorsBorrowed(0), std::runtime_error);
}

TEST(GettingNodeNeighbors, BorrowedRangesAndIteratorsDoNotRetainMappedOwner)
{
    Graph graph(3);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(1, 2);
    auto owner = std::make_shared<BlobBuilder>();
    const auto record = graph.freeze(*owner);
    Graph mapped = Graph::mapped(MappedBytes(owner->bytes.data(), owner->bytes.size(), owner), record);
    const auto retained = owner.use_count();

    const auto successors = mapped.successorsBorrowed(0);
    const auto predecessors = mapped.predecessorsBorrowed(2);
    const auto copied = successors;
    auto first = copied.begin();
    const auto end = copied.end();
    const auto beforeIncrement = first++;
    const auto found = predecessors.find(1);
    const auto reversed = ReverseGraph(mapped).successorsBorrowed(2);
    const auto reverseFirst = reversed.begin();
    EXPECT_EQ(retained, owner.use_count());
    EXPECT_EQ(1u, *beforeIncrement);
    EXPECT_EQ(2u, *first);
    EXPECT_NE(first, end);
    EXPECT_EQ(1u, *found);
    EXPECT_EQ(0u, *reverseFirst);

    const auto owning = mapped.successors(0);
    EXPECT_EQ(retained + 1, owner.use_count());
    const auto owningFirst = owning.begin();
    const auto owningEnd = owning.end();
    EXPECT_EQ(retained + 3, owner.use_count());
    EXPECT_NE(owningFirst, owningEnd);
}

TEST(GettingNodeNeighbors, OwningMappedRangesAndIteratorsOutliveGraph)
{
    CatalogRange<NodeId> successors, predecessors;
    std::weak_ptr<BlobBuilder> lifetime;
    {
        Graph graph(3);
        graph.addEdge(0, 1);
        graph.addEdge(0, 2);
        graph.addEdge(1, 2);
        auto owner = std::make_shared<BlobBuilder>();
        lifetime = owner;
        const auto record = graph.freeze(*owner);
        Graph mapped = Graph::mapped(MappedBytes(owner->bytes.data(), owner->bytes.size(), owner), record);
        successors = mapped.successors(0);
        predecessors = ReverseGraph(mapped).successors(2);
    }
    ASSERT_FALSE(lifetime.expired());
    EXPECT_EQ((vector<NodeId> { 1, 2 }), vector<NodeId>(successors));
    EXPECT_EQ((vector<NodeId> { 0, 1 }), vector<NodeId>(predecessors));
    {
        const auto firstSuccessor = successors.begin(), endSuccessor = successors.end();
        const auto firstPredecessor = predecessors.begin(), endPredecessor = predecessors.end();
        successors = CatalogRange<NodeId>();
        predecessors = CatalogRange<NodeId>();
        ASSERT_FALSE(lifetime.expired());
        EXPECT_EQ((vector<NodeId> { 1, 2 }), vector<NodeId>(firstSuccessor, endSuccessor));
        EXPECT_EQ((vector<NodeId> { 0, 1 }), vector<NodeId>(firstPredecessor, endPredecessor));
    }
    EXPECT_TRUE(lifetime.expired());
}
