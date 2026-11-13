//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include <stdexcept>

#include "graphutils/SequenceOperations.hh"

using std::set;
using std::string;
using std::to_string;
using std::vector;

namespace graphtools
{

void Graph::init(size_t numnodes_)
{
    nodes_.resize(numnodes_);
    adjacency_list_.resize(numnodes_);
    reverse_adjacency_list_.resize(numnodes_);
}

void Graph::assertNodeExists(NodeId node_id) const
{
    if (node_id >= numNodes())
    {
        throw std::logic_error("Node with id " + to_string(node_id) + " does not exist");
    }
}

void Graph::assertEdgeExists(NodeIdPair node_id_pair) const
{
    if (!hasEdge(node_id_pair.first, node_id_pair.second))
    {
        throw std::logic_error(
            "There is no edge between " + to_string(node_id_pair.first) + " and " + to_string(node_id_pair.second));
    }
}

void assertValidSequence(string const& seq)
{
    if (seq.empty())
    {
        throw std::logic_error("Invalid node sequence " + seq);
    }
}

TextView Graph::nodeName(NodeId node_id) const
{
    assertNodeExists(node_id);
    return mapped_.data ? mapped_.text(mappedNode(node_id).name) : TextView(nodes_[node_id].name);
}
void Graph::setNodeName(NodeId node_id, const std::string& node_name)
{
    assertMutable();
    assertNodeExists(node_id);
    nodes_[node_id].name = node_name;
}

TextView Graph::nodeSeq(const NodeId node_id) const
{
    assertNodeExists(node_id);
    return mapped_.data ? mapped_.text(mappedNode(node_id).sequence) : TextView(nodes_[node_id].sequence);
}

size_t Graph::nodeSeqLength(const NodeId node_id) const
{
    assertNodeExists(node_id);
    if (!mapped_.data)
        return nodes_[node_id].sequence.size();
    const auto sequence = mappedNode(node_id).sequence;
    // Preserve the same mapped bounds checks as nodeSeq().
    mapped_.array<char>(sequence);
    return sequence.count;
}

CatalogRange<TextView> Graph::nodeSeqExpansion(NodeId node_id) const
{
    assertNodeExists(node_id);
    if (mapped_.data)
    {
        const auto strings = mappedNode(node_id).expansions;
        return CatalogRange<TextView>(
            mapped_.data, strings.count,
            [](const void* base, size_t i, const void* aux, size_t)
            {
                const auto slice = static_cast<const Slice*>(aux)[i];
                return TextView(static_cast<const char*>(base) + slice.offset, slice.count);
            },
            mapped_.array<Slice>(strings), mapped_.size, mapped_.owner);
    }
    const auto& strings = nodes_[node_id].sequence_expansion;
    return CatalogRange<TextView>(
        &strings, strings.size(), [](const void* base, size_t i, const void*, size_t)
        { return TextView((*static_cast<const vector<string>*>(base))[i]); });
}

void Graph::setNodeSeq(NodeId node_id, const string& sequence)
{
    assertMutable();
    assertNodeExists(node_id);
    assertValidSequence(sequence);
    nodes_[node_id].sequence = sequence;
    expandReferenceSequence(sequence, nodes_[node_id].sequence_expansion);
}

void Graph::addEdge(NodeId source_id, NodeId sink_id)
{
    assertMutable();
    assertNodeExists(source_id);
    assertNodeExists(sink_id);

    const string edge_encoding = "(" + to_string(source_id) + " ," + to_string(sink_id) + ")";
    if (hasEdge(source_id, sink_id))
    {
        throw std::logic_error("Graph already contains edge " + edge_encoding);
    }

    if (source_id > sink_id)
    {
        throw std::logic_error("Edge " + edge_encoding + " breaks topological order");
    }

    NodeIdPair node_id_pair = std::make_pair(source_id, sink_id);
    edge_labels_[node_id_pair];

    adjacency_list_[source_id].emplace(sink_id);
    reverse_adjacency_list_[sink_id].emplace(source_id);
}

bool Graph::hasEdge(NodeId source_id, NodeId sink_id) const
{
    assertNodeExists(source_id);
    assertNodeExists(sink_id);
    if (mapped_.data)
    {
        const auto slice = mappedNode(source_id).successors;
        const auto* ids = mapped_.array<NodeId>(slice);
        return std::find(ids, ids + slice.count, sink_id) != ids + slice.count;
    }
    return edge_labels_.find(NodeIdPair(source_id, sink_id)) != edge_labels_.end();
}

void Graph::addLabelToEdge(NodeId source_id, NodeId sink_id, const std::string& label)
{
    assertMutable();
    NodeIdPair node_pair = std::make_pair(source_id, sink_id);
    assertEdgeExists(node_pair);
    edge_labels_[node_pair].insert(label);
}

SortedLabels Graph::allLabels() const
{
    SortedLabels labels;
    if (mapped_.data)
    {
        for (const auto& edge :
             CatalogRange<MappedEdge>::fromArray(mapped_.array<MappedEdge>(mappedGraph_->edges), numEdges()))
            for (const auto& label : edgeLabels(edge.from, edge.to))
                labels.insert(label);
        return labels;
    }
    for (const auto& single_edge_labels_ : edge_labels_)
    {
        labels.insert(single_edge_labels_.second.begin(), single_edge_labels_.second.end());
    }
    return labels;
}

CatalogRange<TextView> Graph::edgeLabels(NodeId source_id, NodeId sink_id) const
{
    if (mapped_.data)
    {
        for (size_t i = 0; i < numEdges(); ++i)
        {
            const auto& edge = mapped_.array<MappedEdge>(mappedGraph_->edges)[i];
            if (edge.from != source_id || edge.to != sink_id)
                continue;
            return { mapped_.data,
                     edge.labels.count,
                     [](const void* base, size_t i, const void* aux, size_t)
                     {
                         const auto slice = static_cast<const Slice*>(aux)[i];
                         return TextView(static_cast<const char*>(base) + slice.offset, slice.count);
                     },
                     mapped_.array<Slice>(edge.labels),
                     mapped_.size,
                     mapped_.owner };
        }
        throw std::out_of_range("Missing mapped graph edge");
    }
    NodeIdPair node_id_pair = std::make_pair(source_id, sink_id);
    const auto& labels = edge_labels_.at(node_id_pair);
    return { &labels, labels.size(), [](const void* base, size_t i, const void*, size_t)
             {
                 auto it = static_cast<const Labels*>(base)->begin();
                 std::advance(it, i);
                 return TextView(*it);
             } };
}

std::set<NodeIdPair> Graph::edgesWithLabel(const std::string& label) const
{
    std::set<NodeIdPair> result;
    if (mapped_.data)
    {
        for (size_t i = 0; i < numEdges(); ++i)
        {
            const auto& edge = mapped_.array<MappedEdge>(mappedGraph_->edges)[i];
            if (edgeLabels(edge.from, edge.to).count(label))
                result.emplace(edge.from, edge.to);
        }
        return result;
    }
    for (const auto& pair : edge_labels_)
    {
        if (pair.second.count(label))
        {
            result.insert(pair.first);
        }
    }
    return result;
}

void Graph::eraseLabel(const std::string& label)
{
    assertMutable();
    for (auto& pair : edge_labels_)
    {
        pair.second.erase(label);
    }
}

CatalogRange<NodeId> Graph::successors(NodeId node_id) const
{
    assertNodeExists(node_id);
    if (mapped_.data)
    {
        const auto slice = mappedNode(node_id).successors;
        return CatalogRange<NodeId>::fromArray(mapped_.array<NodeId>(slice), slice.count, mapped_.owner);
    }
    return CatalogRange<NodeId>::fromContainer(adjacency_list_[node_id]);
}

CatalogRange<NodeId> Graph::predecessors(NodeId node_id) const
{
    assertNodeExists(node_id);
    if (mapped_.data)
    {
        const auto slice = mappedNode(node_id).predecessors;
        return CatalogRange<NodeId>::fromArray(mapped_.array<NodeId>(slice), slice.count, mapped_.owner);
    }
    return CatalogRange<NodeId>::fromContainer(reverse_adjacency_list_[node_id]);
}

void Graph::assertMutable() const
{
    if (mapped_.data)
        throw std::logic_error("A prepared graph is immutable");
}

Slice Graph::freeze(BlobBuilder& out) const
{
    std::vector<MappedNode> nodes;
    std::vector<MappedEdge> edges;
    for (NodeId node = 0; node < numNodes(); ++node)
    {
        std::vector<Slice> expansions;
        for (auto text : nodeSeqExpansion(node))
            expansions.push_back(out.text(text));
        std::vector<NodeId> next = successors(node), previous = predecessors(node);
        nodes.push_back({ out.text(nodeName(node)), out.text(nodeSeq(node)), out.append(expansions), out.append(next),
                          out.append(previous) });
        for (auto successor : next)
        {
            std::vector<Slice> labels;
            for (const auto& label : edgeLabels(node, successor))
                labels.push_back(out.text(label));
            edges.push_back({ node, successor, out.append(labels) });
        }
    }
    return out.append(MappedGraph { out.append(nodes), out.append(edges), out.text(graphId()) });
}

Graph Graph::mapped(MappedBytes bytes, Slice record)
{
    if (record.count != 1)
        throw std::runtime_error("Invalid mapped graph record");
    Graph graph;
    graph.mapped_ = bytes;
    graph.mappedGraph_ = bytes.array<MappedGraph>(record);
    bytes.array<MappedNode>(graph.mappedGraph_->nodes);
    bytes.array<MappedEdge>(graph.mappedGraph_->edges);
    return graph;
}
}
