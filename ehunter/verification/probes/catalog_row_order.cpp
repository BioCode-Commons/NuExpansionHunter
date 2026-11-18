//
// NuExpansionHunter
// Copyright 2026 NuExpansionHunter contributors
//
// This file is part of NuExpansionHunter, a modified version of
// Illumina ExpansionHunter v5.0.0 (Copyright 2016-2019 Illumina, Inc.).
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

// Precondition check for K-04/K-05: every graph that ExpansionHunter aligns gives the dense matrix
// ordered DAG rows. For each locus of the given catalogs, as loaded by the program, this runs the
// edge discovery that PinnedDagAligner::prepareGraph uses, forward (prefix alignment) and on the
// reverse graph (suffix alignment), from every node, at several start offsets and extension lengths.
// Each edge list must go from a lower to a higher row and be sorted by row (the precondition). The
// summary also counts edge maps with hasOrderedNodeSpans, which selects the reset that skips clearing
// interior cells. Edge lists that discovery rejects with an exception never reach the matrix.
//
// usage: catalog_row_order <reads.bam> <reference.fa> <catalog.json>...
// Prints one summary line per catalog; exits nonzero on any violation.

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/Parameters.hh"
#include "core/Reference.hh"
#include "graphalign/PinnedDagAligner.hh"
#include "graphalign/dagAligner/Details.hh"
#include "graphcore/Graph.hh"
#include "io/CatalogLoading.hh"
#include "io/SampleStats.hh"

using graphtools::Graph;
using graphtools::NodeId;
using graphtools::ReverseGraph;

template <typename Tag, typename Tag::Member member> struct Access
{
    friend typename Tag::Member access(Tag) { return member; }
};
template <typename GraphT>
using Discover = void (*)(const GraphT&, NodeId, std::size_t, std::size_t, std::vector<int>&,
    std::vector<std::pair<int, int>>&, std::string&, std::vector<NodeId>&);
struct ForwardTag { using Member = Discover<Graph>; friend Member access(ForwardTag); };
struct ReverseTag { using Member = Discover<ReverseGraph>; friend Member access(ReverseTag); };
template struct Access<ForwardTag, &graphtools::PinnedDagAligner::bfsDiscoverEdges<Graph>>;
template struct Access<ReverseTag, &graphtools::PinnedDagAligner::bfsDiscoverEdges<ReverseGraph>>;

struct Counts
{
    uint64_t loci = 0, preparations = 0, ordered = 0, violations = 0, spans = 0, rejected = 0;
};

template <typename GraphT>
void check(Discover<GraphT> discover, const GraphT& graph, NodeId node, std::size_t offset, std::size_t extension,
    Counts& counts, const std::string& locus)
{
    std::vector<int> nodeIds;
    std::vector<std::pair<int, int>> edges;
    std::string target;
    std::vector<NodeId> originalIds;
    try
    {
        discover(graph, node, offset, extension, nodeIds, edges, target, originalIds);
    }
    catch (const std::exception&)
    {
        ++counts.rejected;
        return;
    }
    // prepareGraph appends the length marker and builds the EdgeMap only for a nonempty target.
    edges.emplace_back(int(target.length()), int(target.length()));
    if (target.empty())
        return;
    ++counts.preparations;
    bool ok = true;
    int previousRow = 0;
    for (std::size_t i = 0; i + 1 < edges.size(); ++i)
    {
        const auto& edge = edges[i];
        ok = ok && edge.first >= -1 && edge.first < edge.second && edge.second >= previousRow
            && edge.second <= int(target.length());
        previousRow = edge.second;
    }
    const graphalign::dagAligner::EdgeMap edgeMap(edges, nodeIds);
    if (edgeMap.hasOrderedNodeSpans(target.length()))
        ++counts.spans;
    if (ok)
        ++counts.ordered;
    else if (counts.violations++ < 10)
        std::cout << "violation: locus " << locus << " node " << node << " offset " << offset << " extension "
                  << extension << '\n';
}

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cerr << "usage: catalog_row_order <reads.bam> <reference.fa> <catalog.json>...\n";
        return 2;
    }
    using namespace ehunter;
    const FastaReference reference(argv[2], extractReferenceContigInfo(argv[1]));
    const HeuristicParameters heuristics(1000, 20, false, graphtools::AlignerType::DAG_ALIGNER);
    const std::size_t extensions[] = { 1, 2, 16, 64, 150, 151, 250, 300, 1000, 3000 };
    bool allOrdered = true;
    for (int c = 3; c < argc; ++c)
    {
        const RegionCatalog catalog = loadLocusCatalogFromDisk(argv[c], heuristics, reference);
        Counts counts;
        for (std::size_t i = 0; i < catalog.size(); ++i)
        {
            const LocusSpecification locus = catalog[i];
            const Graph& graph = locus.regionGraph();
            const ReverseGraph reverse(graph);
            const std::string id(locus.locusId().data(), locus.locusId().size());
            ++counts.loci;
            for (NodeId node = 0; node < graph.numNodes(); ++node)
            {
                const std::size_t length = graph.nodeSeqLength(node);
                std::vector<std::size_t> offsets = { 0, 1, length / 2, length ? length - 1 : 0, length };
                std::sort(offsets.begin(), offsets.end());
                offsets.erase(std::unique(offsets.begin(), offsets.end()), offsets.end());
                for (std::size_t offset : offsets)
                {
                    if (offset > length)
                        continue;
                    for (std::size_t extension : extensions)
                    {
                        check(access(ForwardTag()), graph, node, offset, extension, counts, id);
                        check(access(ReverseTag()), reverse, node, offset, extension, counts, id);
                    }
                }
            }
        }
        std::string name(argv[c]);
        name = name.substr(name.find_last_of('/') + 1);
        std::cout << "catalog " << name << ": loci " << counts.loci << ", preparations " << counts.preparations
                  << ", ordered rows " << counts.ordered << ", violations " << counts.violations
                  << ", ordered node spans " << counts.spans
                  << ", rejected by discovery " << counts.rejected << '\n';
        allOrdered = allOrdered && counts.violations == 0 && counts.ordered == counts.preparations;
    }
    return allOrdered ? 0 : 1;
}
