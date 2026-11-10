//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>

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

#pragma once

#include <algorithm>
#include <cstdint>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graphutils/MappedData.hh"
#include "graphutils/PairHashing.hh"

namespace graphtools
{

using NodeId = uint32_t;
using NodeIdPair = std::pair<NodeId, NodeId>;
using SortedLabels = std::set<std::string>;
using Labels = std::unordered_set<std::string>;
using AdjacencyList = std::vector<std::set<NodeId>>;

struct Node
{
    std::string name;
    std::string sequence;
    std::vector<std::string> sequence_expansion;
};

struct MappedNode
{
    Slice name, sequence, expansions, successors, predecessors;
};
struct MappedEdge
{
    NodeId from, to;
    Slice labels;
};
struct MappedGraph
{
    Slice nodes, edges, id;
};

/**
 * Sequence graph that can hold degenerate nucleotide sequences
 */
class Graph
{
public:
    explicit Graph(size_t num_nodes = 0, std::string const& id = "")
        : graphId_(id)
    {
        init(num_nodes);
    }

    virtual ~Graph() = default;
    bool isMapped() const { return mapped_.data != nullptr; }

    size_t numNodes() const { return mapped_.data ? mappedGraph_->nodes.count : nodes_.size(); }
    size_t numEdges() const { return mapped_.data ? mappedGraph_->edges.count : edge_labels_.size(); }
    TextView nodeName(NodeId node_id) const;
    void setNodeName(NodeId node_id, const std::string& node_name);
    TextView nodeSeq(NodeId node_id) const;
    void setNodeSeq(NodeId node_id, const std::string& sequence);
    CatalogRange<TextView> nodeSeqExpansion(NodeId node_id) const;
    void addEdge(NodeId source_node_id, NodeId sink_node_id);
    bool hasEdge(NodeId source_node_id, NodeId sink_node_id) const;

    SortedLabels allLabels() const;
    CatalogRange<TextView> edgeLabels(NodeId source_node_id, NodeId sink_node_id) const;
    // All edges in the graph with this label
    std::set<NodeIdPair> edgesWithLabel(const std::string& label) const;
    /**
     * All label to an existing edge
     * @throws if source -> sink is not an edge in the graph
     */
    void addLabelToEdge(NodeId source_node_id, NodeId sink_node_id, const std::string& label);
    // Remove this label from all edges
    void eraseLabel(const std::string& label);

    // this cannot be const if graphs are to be assigned. Currently this happens in unit tests.
    TextView graphId() const { return mapped_.data ? mapped_.text(mappedGraph_->id) : TextView(graphId_); }
    void setGraphId(const std::string& id)
    {
        assertMutable();
        graphId_ = id;
    }
    CatalogRange<NodeId> successors(NodeId node_id) const;
    CatalogRange<NodeId> predecessors(NodeId node_id) const;
    Slice freeze(BlobBuilder& out) const;
    static Graph mapped(MappedBytes bytes, Slice record);

private:
    void init(size_t num_nodes);
    void assertNodeExists(NodeId node_id) const;
    void assertEdgeExists(NodeIdPair edge) const;
    void assertMutable() const;
    const MappedNode& mappedNode(NodeId node) const { return mapped_.array<MappedNode>(mappedGraph_->nodes)[node]; }
    MappedBytes mapped_;
    const MappedGraph* mappedGraph_ = nullptr;
    std::string graphId_;
    std::vector<Node> nodes_;
    std::unordered_map<NodeIdPair, Labels> edge_labels_;
    AdjacencyList adjacency_list_;
    AdjacencyList reverse_adjacency_list_;
};

class ReverseGraph
{
    const Graph& graph_;

public:
    explicit ReverseGraph(const Graph& graph)
        : graph_(graph)
    {
    }

    size_t numNodes() const { return graph_.numNodes(); }
    size_t numEdges() const { return graph_.numEdges(); }
    TextView nodeName(NodeId nodeId) const { return graph_.nodeName(nodeId); }
    const std::string nodeSeq(NodeId nodeId) const
    {
        std::string ret = graph_.nodeSeq(nodeId);
        std::reverse(ret.begin(), ret.end());
        return ret;
    }
    bool hasEdge(NodeId sourceNodeId, NodeId sinkNodeId) const { return graph_.hasEdge(sourceNodeId, sinkNodeId); }

    SortedLabels allLabels() const { return graph_.allLabels(); }
    CatalogRange<TextView> edgeLabels(NodeId sourceNodeId, NodeId sinkNodeId) const
    {
        return graph_.edgeLabels(sourceNodeId, sinkNodeId);
    }

    // All edges in the graph with this label
    std::set<NodeIdPair> edgesWithLabel(const std::string& label) const { return graph_.edgesWithLabel(label); }

    CatalogRange<NodeId> successors(NodeId nodeId) const { return graph_.predecessors(nodeId); }
    CatalogRange<NodeId> predecessors(NodeId nodeId) const { return graph_.successors(nodeId); }
};
}
