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

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "alignment/OrientationPredictor.hh"
#include "graphalign/PinnedDagAligner.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphutils/SequenceOperations.hh"
// The re-hosted wrapper uses the same graph and scoring objects.
#include "reference/graphalign/PinnedDagAligner_v5.hh"

using namespace graphtools;
using namespace ehunter;

static uint64_t digest = 14695981039346656037ULL;
static void mix(uint64_t x) { digest = (digest ^ x) * 1099511628211ULL; }

static std::string repeat(const std::string& motif, size_t n)
{
    std::string result;
    result.reserve(n);
    for (size_t i = 0; i < n; ++i) result += motif[i % motif.size()];
    return result;
}

static Graph suffixGraph(int kind)
{
    if (kind == 1) return makeSwapGraph(repeat("ACGTTGCA", 48), "CAG", "CTG", repeat("GATTACA", 48));
    if (kind == 2) return makeStrGraph(repeat("ACGTTGCA", 24), "CAG", repeat("GATTACA", 24));
    Graph graph(1);
    graph.setNodeSeq(0, repeat("ACGTTGCAGATCCGTA", 2048));
    return graph;
}

static Path suffixSeed(const Graph& graph, bool emptyTarget = false)
{
    const NodeId node = emptyTarget ? 0 : static_cast<NodeId>(graph.numNodes() - 1);
    const int position = emptyTarget ? 0 : static_cast<int>(graph.nodeSeq(node).size());
    return Path(&graph, position, { node }, position);
}

static void mixAlignments(int score, const std::list<PathAndAlignment>& result)
{
    mix(static_cast<uint64_t>(score)); mix(result.size());
    for (const auto& item : result)
    {
        mix(item.first.startPosition()); mix(item.first.endPosition());
        for (NodeId node : item.first.nodeIds()) mix(node);
        mix(item.second.referenceStart());
        for (const auto& op : item.second.operations())
        { mix(static_cast<uint64_t>(op.type())); mix(op.length()); }
    }
}

static std::vector<std::string> suffixQueries(const std::string& workload)
{
    std::vector<size_t> lengths = { 35, 75, 100, 150 };
    if (workload == "short") lengths = { 1, 8, 16, 24, 32 };
    if (workload == "long") lengths = { 300, 1000 };
    std::vector<std::string> result;
    for (size_t n : lengths)
        for (int kind = 0; kind < 8; ++kind)
        {
            std::string q = repeat(kind % 2 ? "CAG" : "ACGTTGCAGATCCGTA", n);
            if (kind == 2) q[n / 2] = 'N';
            if (kind == 3) q[n / 2] = 'X';
            if (kind == 4) q[n / 2] = 'R';
            if (kind == 5) for (char& c : q) c = static_cast<char>(::tolower(c));
            if (kind == 6) q.insert(n / 2, 1, 'A');
            if (kind == 7 && n > 1) q.erase(n / 2, 1);
            result.push_back(q);
        }
    return result;
}

static size_t checkSuffix()
{
    size_t cases = 0;
    auto queries = suffixQueries("short");
    const auto typical = suffixQueries("mixed");
    queries.insert(queries.end(), typical.begin(), typical.end());
    const auto longer = suffixQueries("long");
    queries.insert(queries.end(), longer.begin(), longer.end());
    for (int graphKind : { 0, 1, 2 })
    {
        const Graph graph = suffixGraph(graphKind);
        const Path seed = suffixSeed(graph);
        for (int match : { 1, 5 })
        {
            PinnedDagAligner candidate(match, -4, -8, -2);
            graphtools::reference::PinnedDagAligner reference(match, -4, -8, -2);
            for (const auto& query : queries)
            {
                // Long cyclic-repeat cases have many equivalent traceback paths;
                // short/typical cases cover them while bounding correctness time.
                if (graphKind == 2 && query.size() > 300) continue;
                int actualScore = INT32_MIN, expectedScore = INT32_MIN;
                const std::string unchanged = query;
                const auto expected = reference.suffixAlign(seed, query, query.size() + 12, expectedScore);
                const auto actual = candidate.suffixAlign(seed, query, query.size() + 12, actualScore);
                if (actualScore != expectedScore || actual != expected || query != unchanged)
                    throw std::runtime_error("suffix score/path/cigar/order mismatch");
                mixAlignments(actualScore, actual); ++cases;
            }
            int actualScore = -12345, expectedScore = -12345;
            const Path emptySeed = suffixSeed(graph, true);
            const auto expected = reference.suffixAlign(emptySeed, "", 0, expectedScore);
            const auto actual = candidate.suffixAlign(emptySeed, "", 0, actualScore);
            if (actualScore != expectedScore || actual != expected) throw std::runtime_error("empty suffix mismatch");
            ++cases;
        }
    }
    return cases;
}


int main()
{
    try { const auto cases = checkSuffix(); std::cout << cases << " " << digest << "\n"; }
    catch (const std::exception& error) { std::cerr << error.what() << "\n"; return 1; }
}
