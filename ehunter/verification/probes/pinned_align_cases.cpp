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

// Verification build: -DEH_VERIFY="BUDGET_SEARCH=NEVER;FAST_PATHS=OFF".
// Exercise actual prefixAlign/suffixAlign calls, including graph expansion,
// DP, traceback, graph-CIGAR parsing, and returned results. Owned and mapped
// graphs use the same production path. Never interpret this as a getter bench.
// Build this identical source with each independent production tree.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "PinnedAlignmentOracle.hh"

using namespace graphtools;

struct Digest
{
    uint64_t hash = UINT64_C(1469598103934665603), words = 0;
    bool dump = false;
    void add(uint64_t value)
    {
        ++words;
        if (dump) std::cout << value << '\n';
        for (int i = 0; i < 8; ++i) { hash ^= (value >> (8 * i)) & 255; hash *= UINT64_C(1099511628211); }
    }
    void text(const std::string& value)
    {
        add(value.size());
        for (unsigned char c : value) add(c);
    }
};

uint32_t randomWord(uint32_t& state)
{
    state ^= state << 13; state ^= state >> 17; state ^= state << 5;
    return state;
}

std::string dna(size_t size, uint32_t state)
{
    std::string sequence(size, 'A');
    for (char& base : sequence) base = "ACGT"[randomWord(state) % 4];
    return sequence;
}

Graph makeGraph(const std::string& shape, size_t queryLength, size_t extension, size_t nodeLength = 0)
{
    if (shape == "linear")
    {
        Graph graph(1);
        graph.setNodeSeq(0, dna(nodeLength ? nodeLength : queryLength + extension + 32, 813));
        return graph;
    }
    if (shape == "join")
    {
        Graph graph(4);
        graph.setNodeSeq(0, "ACGTACTG");
        const size_t middleLength = nodeLength ? nodeLength : std::max<size_t>(8, queryLength > 16 ? queryLength - 16 : 8);
        graph.setNodeSeq(1, dna(middleLength, 421));
        graph.setNodeSeq(2, dna(middleLength, 911));
        graph.setNodeSeq(3, "TTGCACGA");
        graph.addEdge(0, 1); graph.addEdge(0, 2); graph.addEdge(1, 3); graph.addEdge(2, 3);
        return graph;
    }
    if (shape == "repeat")
    {
        Graph graph(3);
        graph.setNodeSeq(0, "ACGTACTG");
        graph.setNodeSeq(1, nodeLength ? dna(nodeLength, 941) : "AG");
        graph.setNodeSeq(2, "TTGCACGA");
        graph.addEdge(0, 1); graph.addEdge(1, 1); graph.addEdge(1, 2); graph.addEdge(0, 2);
        return graph;
    }
    throw std::invalid_argument("unknown graph shape");
}

std::string referenceWalk(const Graph& graph, const std::string& shape, size_t queryLength, unsigned alternate)
{
    if (shape == "linear") return std::string(graph.nodeSeq(0));
    std::string sequence(graph.nodeSeq(0));
    if (shape == "join") sequence += std::string(graph.nodeSeq(1 + alternate % 2));
    else
        while (sequence.size() < queryLength + 8) sequence += graph.nodeSeq(1);
    sequence += std::string(graph.nodeSeq(graph.numNodes() - 1));
    return sequence;
}

std::string queryFor(const Graph& graph, const std::string& shape, size_t queryLength, bool suffix, unsigned sample)
{
    auto walk = referenceWalk(graph, shape, queryLength, sample);
    // An intentionally short non-repeat graph exercises clipped query tails.
    if (walk.size() < queryLength) walk += dna(queryLength - walk.size(), 631);
    std::string query = suffix ? walk.substr(walk.size() - queryLength) : walk.substr(0, queryLength);
    // Every fourth query is exact; others have sparse substitutions. Construction
    // is outside benchmark timing, and cases exercise multi-operation CIGARs.
    if (sample % 4)
        for (size_t i = 5 + sample % 11; i < query.size(); i += 19)
            query[i] = query[i] == 'A' ? 'C' : 'A';
    return query;
}

Path makeSeed(const Graph& graph, bool suffix, unsigned seedMode)
{
    const NodeId node = suffix ? graph.numNodes() - 1 : 0;
    const size_t length = graph.nodeSeq(node).size();
    const int seedLength = seedMode == 2 ? length : seedMode ? std::min<size_t>(2, length) : 0;
    const int end = suffix ? length : seedLength;
    return Path(&graph, end - seedLength, {node}, end);
}

void inspectAlignment(const Alignment& alignment, Digest& out)
{
    out.add(alignment.referenceStart()); out.add(alignment.referenceLength()); out.add(alignment.queryLength());
    out.add(alignment.numOperations()); out.add(alignment.numMatched()); out.add(alignment.numMismatched());
    out.add(alignment.numClipped()); out.add(alignment.numInserted()); out.add(alignment.numDeleted());
    out.text(alignment.generateCigar());
    for (const auto& operation : alignment)
    {
        out.add(static_cast<int>(operation.type())); out.add(operation.length());
        out.add(operation.referenceLength()); out.add(operation.queryLength()); out.text(operation.generateCigar());
    }
}

void verifyCase(CheckedPinnedAligner& aligner, const Graph& graph, const std::string& shape,
                int scoring, int length, int extension, bool suffix, unsigned seedMode,
                unsigned sample, const std::string& storage, size_t nodeLength, Digest& out)
{
    if (out.dump) std::cout << "CASE pinned " << storage << ' ' << shape << ' ' << scoring << ' '
        << length << ' ' << extension << ' ' << suffix << ' ' << seedMode << ' ' << sample << ' ' << nodeLength << '\n';
    const auto query = queryFor(graph, shape, length, suffix, sample);
    const auto seed = makeSeed(graph, suffix, seedMode);
    out.text(storage); out.text(shape); out.add(scoring); out.add(length); out.add(extension);
    out.add(suffix); out.add(seedMode); out.add(sample); out.add(nodeLength); out.text(query);
    int score = -773;
    try
    {
        const auto results = suffix ? aligner.suffixAlign(seed, query, extension, score)
                                    : aligner.prefixAlign(seed, query, extension, score);
        out.add(1); out.add(score); out.add(results.size());
        for (const auto& result : results)
        {
            out.text(result.first.encode()); out.text(result.first.seq());
            inspectAlignment(result.second, out);
        }
    }
    catch (const std::exception& error) { out.add(0); out.add(score); out.text(error.what()); }
}

uint64_t verifyPinned(Digest& out)
{
    uint64_t cases = 0;
    const int lengths[] = {1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 150, 250, 8};
    for (const std::string shape : {"linear", "join", "repeat"})
        for (int scoring : {0, 1})
        {
            CheckedPinnedAligner aligner(scoring ? 3 : 5, scoring ? -2 : -4, scoring ? -5 : -8, -2);
            for (int length : lengths)
            {
                Graph owned = makeGraph(shape, length, length + 8);
                BlobBuilder blob;
                const auto record = owned.freeze(blob);
                const Graph mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);
                for (const std::string storage : {"owned", "mapped"})
                {
                    const Graph& graph = storage == "owned" ? owned : mapped;
                    for (int extension : {0, 1, length + 8})
                        for (bool suffix : {false, true})
                            for (unsigned seedMode : {0u, 1u})
                                for (unsigned sample : {0u, 1u, 3u})
                                {
                                    verifyCase(aligner, graph, shape, scoring, length, extension, suffix,
                                               seedMode, sample, storage, 0, out);
                                    ++cases;
                                }
                }
            }
        }
    // Override actual node length independently of read length. This includes
    // single-base nodes, SSO-sized sequences, heap strings and long branches.
    for (const std::string shape : {"linear", "join", "repeat"})
        for (size_t nodeLength : {size_t(1), size_t(22), size_t(256), size_t(2048)})
            for (int length : {1, 16, 150})
            {
                Graph owned = makeGraph(shape, length, length + 8, nodeLength);
                BlobBuilder blob;
                const auto record = owned.freeze(blob);
                const Graph mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);
                CheckedPinnedAligner aligner(5, -4, -8, -2);
                for (const std::string storage : {"owned", "mapped"})
                {
                    const Graph& graph = storage == "owned" ? owned : mapped;
                    for (bool suffix : {false, true})
                        for (unsigned seedMode : {0u, 1u, 2u})
                            for (unsigned sample : {0u, 3u})
                            {
                                verifyCase(aligner, graph, shape, 0, length, length + 8, suffix,
                                           seedMode, sample, storage, nodeLength, out);
                                ++cases;
                            }
                }
            }
    // Equal-score graph alternatives exercise output order and tie expansion.
    for (const std::string shape : {"linear", "join", "repeat"})
        for (int length : {1, 7, 16})
        {
            Graph owned = makeGraph(shape, length, length + 2, 3);
            BlobBuilder blob;
            const auto record = owned.freeze(blob);
            const Graph mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);
            CheckedPinnedAligner aligner(0, 0, 0, 0);
            for (const std::string storage : {"owned", "mapped"})
            {
                const Graph& graph = storage == "owned" ? owned : mapped;
                for (bool suffix : {false, true})
                    for (unsigned seedMode : {0u, 1u, 2u})
                    {
                        verifyCase(aligner, graph, shape, 2, length, length + 2, suffix,
                                   seedMode, 0, storage, 3, out);
                        ++cases;
                    }
            }
        }
    return cases;
}

void inspectExplicit(CheckedPinnedAligner& aligner, const Graph& graph, const std::string& query,
                     int extension, bool suffix, unsigned seedMode, const std::string& tag, Digest& out)
{
    if (out.dump) std::cout << "CASE explicit " << tag << ' ' << query.size() << ' ' << extension
        << ' ' << suffix << ' ' << seedMode << '\n';
    out.text(tag); out.text(query); out.add(extension); out.add(suffix); out.add(seedMode);
    const auto seed = makeSeed(graph, suffix, seedMode);
    int score = -773;
    try
    {
        const auto results = suffix ? aligner.suffixAlign(seed, query, extension, score)
                                    : aligner.prefixAlign(seed, query, extension, score);
        out.add(1); out.add(score); out.add(results.size());
        for (const auto& result : results)
        {
            out.text(result.first.encode()); out.text(result.first.seq());
            out.add(result.first.startPosition()); out.add(result.first.endPosition());
            out.add(result.first.numNodes());
            for (NodeId node : result.first.nodeIds()) out.add(node);
            inspectAlignment(result.second, out);
        }
    }
    catch (const std::exception& error) { out.add(0); out.add(score); out.text(error.what()); }
}

uint64_t verifyAlgorithmCases(Digest& out)
{
    uint64_t cases = 0;
    const auto inspect = [&](const Graph& owned, const std::vector<std::string>& queries,
                             const std::string& tag, int padding, bool alternateScoring)
    {
        BlobBuilder blob;
        const auto record = owned.freeze(blob);
        const Graph mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);
        for (const std::string storage : {"owned", "mapped"})
            for (int scoring = 0; scoring < (alternateScoring ? 3 : 1); ++scoring)
            {
                const Graph& graph = storage == "owned" ? owned : mapped;
                CheckedPinnedAligner aligner(scoring == 1 ? 4 : 5, scoring == 1 ? -3 : -4,
                                        scoring == 1 ? -7 : -8, scoring == 2 ? -1 : -2);
                for (const auto& query : queries)
                    for (bool suffix : {false, true})
                        for (unsigned seedMode : {0u, 1u, 2u})
                        {
                            inspectExplicit(aligner, graph, query, int(query.size()) + padding,
                                suffix, seedMode, tag + "/" + storage + "/" + std::to_string(scoring), out);
                            ++cases;
                        }
            }
    };

    // More equal-score graph paths than the production top-ten traceback cap.
    // These cases use default scores; zero-score cases cover fallback separately.
    for (int alternatives : {2, 10, 11, 20})
    {
        Graph graph(alternatives + 2);
        graph.setNodeSeq(0, "ACGTG");
        graph.setNodeSeq(alternatives + 1, "TTGCA");
        for (int node = 1; node <= alternatives; ++node)
        {
            graph.setNodeSeq(node, "CAGCAGCAG");
            graph.addEdge(0, node); graph.addEdge(node, alternatives + 1);
        }
        inspect(graph, {"ACGTGCAGCAGCAGTTGCA", "ACGTGCAGCAGTTGCA", "ACGTGCAGCAGCAGCAGTTGCA",
                        "ACGTGCAGCNGCAGTTGCA", "CAGCAGCAG", ""},
                "many-ties-" + std::to_string(alternatives), 8, false);
    }

    // Repeat-length alternatives can merge into the same downstream position.
    // This tests graph-distance bounds and the ordering of repeated copy IDs.
    {
        Graph graph(5);
        graph.setNodeSeq(0, "ACGTG"); graph.setNodeSeq(1, "CAG"); graph.setNodeSeq(2, "CAACAG");
        graph.setNodeSeq(3, "CCG"); graph.setNodeSeq(4, "TTGCA");
        for (const auto& edge : {std::make_pair(0,1), {1,1}, {1,2}, {0,2}, {2,3}, {3,3}, {3,4}, {2,4}})
            graph.addEdge(edge.first, edge.second);
        std::vector<std::string> queries;
        for (int c : {0, 1, 2, 5, 15})
            for (int d : {0, 1, 4})
            {
                std::string q = "ACGTG";
                for (int i = 0; i < c; ++i) q += "CAG";
                q += "CAACAG";
                for (int i = 0; i < d; ++i) q += "CCG";
                q += "TTGCA";
                queries.push_back(q);
                if (q.size() > 8) { q[7] = 'N'; queries.push_back(q); }
            }
        inspect(graph, queries, "two-repeats", 10, true);
    }

    // Asymmetric long insertions/deletions, leading/trailing gaps and IUPAC.
    {
        Graph graph(1); const std::string reference = dna(180, 9041); graph.setNodeSeq(0, reference);
        std::vector<std::string> queries = {reference, reference.substr(0, 150), "", "A", "N"};
        for (int gap : {1, 3, 9, 25, 60})
        {
            std::string q = reference; q.insert(80, dna(gap, 533)); queries.push_back(q);
            q = reference; q.erase(60, gap); queries.push_back(q);
            queries.push_back(dna(gap, 431) + reference);
            queries.push_back(reference + dna(gap, 223));
        }
        for (const char base : std::string("NRYKMSWBDHV"))
        {
            std::string q = reference; q[75] = base; queries.push_back(q);
        }
        std::string lower = reference;
        for (char& base : lower) base = char(base - 'A' + 'a');
        queries.push_back(lower); queries.push_back(std::string(180, 'N'));
        inspect(graph, queries, "gaps-iupac", 12, true);
        Graph ambiguous(1); ambiguous.setNodeSeq(0, "ACGTRYSWKMBDHVNACGT");
        inspect(ambiguous, {"ACGTACGTACGTACGTAC", "ACGTRYSWKMBDHVNACGT", "NNNNNNNNNNNNNNNNNN"},
                "ambiguous-target", 5, false);
    }

    // Explicit overflow-guard/fallback dimensions; preserve reference exceptions.
    for (int length : {1200, 1300, 1320, 2200, 3300})
    {
        Graph graph(1); const std::string reference = dna(length + 8, 381); graph.setNodeSeq(0, reference);
        CheckedPinnedAligner aligner(5, -4, -8, -2);
        std::string q = reference.substr(0, length);
        for (bool suffix : {false, true})
        {
            inspectExplicit(aligner, graph, q, length + 8, suffix, 0, "long-default", out); ++cases;
            q[length / 2] = q[length / 2] == 'A' ? 'C' : 'A';
            inspectExplicit(aligner, graph, q, length + 8, suffix, 0, "long-mutated", out); ++cases;
        }
    }
    return cases;
}

void benchmark(const std::string& storage, const std::string& direction, const std::string& shape,
               int queryLength, int extension, int nodeLength, int repetitions)
{
    if (queryLength <= 0 || extension <= 0 || nodeLength < 0 || repetitions <= 0)
        throw std::invalid_argument("positive workload sizes required");
    if (direction != "prefix" && direction != "suffix" && direction != "mixed")
        throw std::invalid_argument("unknown direction");
    if (storage != "owned" && storage != "mapped") throw std::invalid_argument("unknown storage");
    const Graph owned = makeGraph(shape, queryLength, extension, nodeLength);
    BlobBuilder blob;
    Graph mapped;
    if (storage == "mapped")
    {
        const auto record = owned.freeze(blob);
        mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);
    }
    const Graph& graph = storage == "owned" ? owned : mapped;
    const auto prefixSeed = makeSeed(graph, false, false), suffixSeed = makeSeed(graph, true, false);
    std::vector<std::string> prefixQueries, suffixQueries;
    for (unsigned sample = 0; sample < 64; ++sample)
    {
        prefixQueries.push_back(queryFor(graph, shape, queryLength, false, sample));
        suffixQueries.push_back(queryFor(graph, shape, queryLength, true, sample));
    }
    CheckedPinnedAligner aligner(5, -4, -8, -2);
    uint64_t checksum = 0, prefixResults = 0, suffixResults = 0, operations = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int repetition = 0; repetition < repetitions; ++repetition)
    {
        const bool suffix = direction == "suffix" || (direction == "mixed" && repetition % 2);
        const unsigned sample = repetition % 64;
        int score = -773;
        const auto results = suffix ? aligner.suffixAlign(suffixSeed, suffixQueries[sample], extension, score)
                                    : aligner.prefixAlign(prefixSeed, prefixQueries[sample], extension, score);
        // Require a real completed alignment and consume every path/operation.
        // This measures the full production alignment call, not only a getter.
        if (results.empty()) throw std::runtime_error("benchmark returned no alignment");
        (suffix ? suffixResults : prefixResults) += results.size();
        checksum = checksum * 33 + static_cast<uint64_t>(score);
        for (const auto& result : results)
        {
            const auto& path = result.first;
            const auto& alignment = result.second;
            checksum = checksum * 33 + path.startPosition(); checksum = checksum * 33 + path.endPosition();
            checksum = checksum * 33 + path.numNodes();
            for (NodeId node : path.nodeIds()) checksum = checksum * 33 + node;
            checksum = checksum * 33 + alignment.referenceStart();
            for (const auto& operation : alignment)
            {
                checksum = checksum * 33 + static_cast<int>(operation.type());
                checksum = checksum * 33 + operation.length(); ++operations;
            }
        }
    }
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (!operations) throw std::runtime_error("benchmark returned no operations");
    std::cout << "{\"mode\":\"pinned\",\"storage\":\"" << storage << "\",\"node_length\":" << nodeLength << ",\"direction\":\"" << direction << "\",\"shape\":\"" << shape
              << "\",\"query_length\":" << queryLength << ",\"extension\":" << extension
              << ",\"repetitions\":" << repetitions << ",\"prefix_results\":" << prefixResults
              << ",\"suffix_results\":" << suffixResults << ",\"operations\":" << operations
              << ",\"checksum\":" << checksum << "}\n";
    std::cerr << "{\"kernel_seconds\":" << elapsed << "}\n";
}


std::string lowErrorQuery(const Graph& graph, const std::string& shape, size_t length,
                          bool suffix, unsigned sample, unsigned errors)
{
    const auto walk = referenceWalk(graph, shape, length, sample);
    if (walk.size() < length || errors > length) throw std::invalid_argument("invalid low-error fixture");
    const auto exact = suffix ? walk.substr(walk.size() - length) : walk.substr(0, length);
    auto query = exact;
    for (unsigned index = 0; index < errors; ++index)
    {
        const size_t position = (sample * 7 + index * (length / errors)) % length;
        query[position] = query[position] == 'A' ? 'C' : 'A';
    }
    unsigned differences = 0;
    for (size_t i = 0; i < length; ++i) differences += query[i] != exact[i];
    if (differences != errors) throw std::runtime_error("fixture error count mismatch");
    return query;
}

void lowErrorBenchmark(const std::string& storage, const std::string& direction, const std::string& shape,
               int queryLength, int extension, int nodeLength, int errors, int repetitions, bool dump)
{
    if (queryLength <= 0 || extension <= 0 || nodeLength < 0 || repetitions <= 0)
        throw std::invalid_argument("positive workload sizes required");
    if (direction != "prefix" && direction != "suffix" && direction != "mixed")
        throw std::invalid_argument("unknown direction");
    if (storage != "owned" && storage != "mapped") throw std::invalid_argument("unknown storage");
    const Graph owned = makeGraph(shape, queryLength, extension, nodeLength);
    BlobBuilder blob;
    Graph mapped;
    if (storage == "mapped")
    {
        const auto record = owned.freeze(blob);
        mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);
    }
    const Graph& graph = storage == "owned" ? owned : mapped;
    const auto prefixSeed = makeSeed(graph, false, false), suffixSeed = makeSeed(graph, true, false);
    std::vector<std::string> prefixQueries, suffixQueries;
    for (unsigned sample = 0; sample < 64; ++sample)
    {
        prefixQueries.push_back(lowErrorQuery(graph, shape, queryLength, false, sample, errors));
        suffixQueries.push_back(lowErrorQuery(graph, shape, queryLength, true, sample, errors));
    }
    CheckedPinnedAligner aligner(5, -4, -8, -2);
    Digest output; output.dump = dump;
    uint64_t checksum = 0, prefixResults = 0, suffixResults = 0, operations = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int repetition = 0; repetition < repetitions; ++repetition)
    {
        const bool suffix = direction == "suffix" || (direction == "mixed" && repetition % 2);
        const unsigned sample = repetition % 64;
        int score = -773;
        const auto results = suffix ? aligner.suffixAlign(suffixSeed, suffixQueries[sample], extension, score)
                                    : aligner.prefixAlign(prefixSeed, prefixQueries[sample], extension, score);
        if (dump)
        {
            std::cout << "CASE lowerror " << shape << ' ' << queryLength << ' ' << errors << ' ' << suffix << ' ' << sample << '\n';
            output.add(score); output.add(results.size());
            for (const auto& result : results)
            {
                output.text(result.first.encode()); output.text(result.first.seq());
                output.add(result.first.startPosition()); output.add(result.first.endPosition());
                for (NodeId node : result.first.nodeIds()) output.add(node);
                inspectAlignment(result.second, output);
            }
        }
        // Require a real completed alignment and consume every path/operation.
        // This measures the full production alignment call, not only a getter.
        if (results.empty()) throw std::runtime_error("benchmark returned no alignment");
        (suffix ? suffixResults : prefixResults) += results.size();
        checksum = checksum * 33 + static_cast<uint64_t>(score);
        for (const auto& result : results)
        {
            const auto& path = result.first;
            const auto& alignment = result.second;
            checksum = checksum * 33 + path.startPosition(); checksum = checksum * 33 + path.endPosition();
            checksum = checksum * 33 + path.numNodes();
            for (NodeId node : path.nodeIds()) checksum = checksum * 33 + node;
            checksum = checksum * 33 + alignment.referenceStart();
            for (const auto& operation : alignment)
            {
                checksum = checksum * 33 + static_cast<int>(operation.type());
                checksum = checksum * 33 + operation.length(); ++operations;
            }
        }
    }
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (!operations) throw std::runtime_error("benchmark returned no operations");
    std::cout << "{\"mode\":\"pinned\",\"storage\":\"" << storage << "\",\"node_length\":" << nodeLength << ",\"direction\":\"" << direction << "\",\"shape\":\"" << shape
              << "\",\"query_length\":" << queryLength << ",\"extension\":" << extension
              << ",\"repetitions\":" << repetitions << ",\"prefix_results\":" << prefixResults
              << ",\"suffix_results\":" << suffixResults << ",\"operations\":" << operations
              << ",\"checksum\":" << checksum << "}\n";
    std::cerr << "{\"kernel_seconds\":" << elapsed << "}\n";
}

struct ProofCase
{
    std::string name;
    Graph graph;
    std::string query;
    bool suffix = false;
    int match = 5, mismatch = -4, gapOpen = -8, gapExt = -2;
    bool requireSparse = true, requireAttempt = true;
    int inserted = 0, deleted = 0, expectedResults = 0;
};

Graph proofGraph(int alternatives, bool suffix)
{
    Graph graph(alternatives + 2);
    graph.setNodeSeq(0, dna(suffix ? 512 : 48, 91823));
    graph.setNodeSeq(alternatives + 1, dna(suffix ? 48 : 512, 72149));
    for (int branch = 1; branch <= alternatives; ++branch)
    {
        graph.setNodeSeq(branch, dna(48, 11813));
        graph.addEdge(0, branch); graph.addEdge(branch, alternatives + 1);
    }
    return graph;
}

ProofCase graphCase(int alternatives, bool suffix, int insertion, int deletion)
{
    ProofCase c;
    c.name = "branches-" + std::to_string(alternatives) + (suffix ? "-suffix" : "-prefix")
        + "-i" + std::to_string(insertion) + "-d" + std::to_string(deletion);
    c.graph = proofGraph(alternatives, suffix); c.suffix = suffix;
    const std::string walk = std::string(c.graph.nodeSeq(0)) + std::string(c.graph.nodeSeq(1))
        + std::string(c.graph.nodeSeq(alternatives + 1));
    const size_t length = 128 - insertion + deletion;
    c.query = suffix ? walk.substr(walk.size() - length) : walk.substr(0, length);
    // Delete across the branch/merge boundary; insertion lies on that edge.
    // In suffix orientation the corresponding boundary is length-96.
    const size_t boundary = suffix ? length - 96 : 96;
    if (deletion) c.query.erase(boundary - deletion / 2, deletion);
    if (insertion) c.query.insert(deletion ? (suffix ? boundary + 16 : boundary - 16) : boundary, dna(insertion, 76233));
    if (c.query.size() != 128) throw std::runtime_error("invalid proof query length");
    c.inserted = insertion ? 1 : 0; c.deleted = deletion ? 1 : 0;
    c.expectedResults = alternatives > 10 ? 10 : (insertion || deletion ? 0 : alternatives);
    return c;
}

void checkCaseOutcome(const ProofCase& c)
{
    using namespace graphalign::dagAligner;
    std::vector<int> begin(c.graph.numNodes()), end(c.graph.numNodes()), ids, predecessors;
    std::vector<size_t> index;
    std::string target;
    std::vector<NodeId> order;
    for (NodeId node = 0; node < c.graph.numNodes(); ++node) order.push_back(node);
    if (c.suffix) std::reverse(order.begin(), order.end());
    for (NodeId node : order)
    {
        std::string sequence(c.graph.nodeSeq(node));
        if (c.suffix) std::reverse(sequence.begin(), sequence.end());
        begin[node] = target.size(); target += sequence; end[node] = target.size() - 1;
    }
    for (NodeId node : order)
        for (int row = begin[node]; row <= end[node]; ++row)
        {
            ids.push_back(node); index.push_back(predecessors.size());
            if (row != begin[node]) predecessors.push_back(row - 1);
            else
            {
                const auto incoming = c.suffix ? c.graph.successors(node) : c.graph.predecessors(node);
                if (incoming.empty()) predecessors.push_back(-1);
                for (NodeId parent : incoming) predecessors.push_back(end[parent]);
            }
        }
    index.push_back(predecessors.size());
    const EdgeMap edges(ids, index, predecessors);
    std::string query = c.query;
    if (c.suffix) std::reverse(query.begin(), query.end());
    PinnedBudgetAlignMatrix matrix(BaseMatchingPenaltyMatrix(c.match, c.mismatch), c.gapOpen, c.gapExt);
    BudgetSearch search;
    matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges, 0, 0, &search);
    if (!c.requireAttempt && search.outcome != BudgetSearch::kNotAttempted)
        throw std::runtime_error(c.name + ": ineligible search was attempted");
    if (c.requireAttempt && search.outcome != BudgetSearch::kFound && search.outcome != BudgetSearch::kFailed)
        throw std::runtime_error(c.name + ": search did not expose an attempted outcome");
    if (c.requireSparse && !graphtools::verify::kBudgetSearchNever && search.outcome != BudgetSearch::kFound)
        throw std::runtime_error(c.name + ": sparse witness fell back");
    if (c.requireAttempt && graphtools::verify::kBudgetSearchNever && search.outcome != BudgetSearch::kFailed)
        throw std::runtime_error(c.name + ": dense policy did not force fallback");
}

void runProof(const ProofCase& c, Digest& out)
{
    checkCaseOutcome(c);
    CheckedPinnedAligner aligner(c.match, c.mismatch, c.gapOpen, c.gapExt);
    const auto seed = makeSeed(c.graph, c.suffix, 0);
    int score = -773;
    const auto results = c.suffix ? aligner.suffixAlign(seed, c.query, c.query.size() + 24, score)
                                  : aligner.prefixAlign(seed, c.query, c.query.size() + 24, score);
    if (results.empty()) throw std::runtime_error(c.name + ": no result");
    if (c.expectedResults && int(results.size()) != c.expectedResults)
        throw std::runtime_error(c.name + ": wrong expected number of tied results: " + std::to_string(results.size()));
    std::cout << "CASE activation " << c.name << '\n';
    out.text(c.name); out.text(c.query); out.add(score); out.add(results.size());
    unsigned inserted = 0, deleted = 0;
    for (const auto& result : results)
    {
        const auto& path = result.first;
        out.text(path.encode()); out.text(path.seq()); out.add(path.startPosition()); out.add(path.endPosition());
        out.add(path.numNodes()); for (NodeId node : path.nodeIds()) out.add(node);
        inspectAlignment(result.second, out);
        inserted += result.second.numInserted(); deleted += result.second.numDeleted();
        if (c.graph.numNodes() > 1 && path.numNodes() < 3) throw std::runtime_error(c.name + ": result does not traverse branch");
    }
    if (inserted < unsigned(c.inserted) || deleted < unsigned(c.deleted))
        throw std::runtime_error(c.name + ": intended insertion/deletion not present in optimum");
}

int activationCases()
{
    try
    {
        Digest out; out.dump = true;
        unsigned cases = 0;
        for (bool suffix : {false, true})
            for (int alternatives : {2, 20})
                for (const auto& gap : {std::make_pair(0, 0), {1, 0}, {3, 0}, {0, 3}, {0, 8}, {1, 3}})
                { runProof(graphCase(alternatives, suffix, gap.first, gap.second), out); ++cases; }
        for (int gate = 0; gate < 4; ++gate)
        {
            ProofCase c; c.name = std::string("fallback-") + (gate == 0 ? "short-query" : gate == 1 ? "score" : gate == 2 ? "numeric-guard" : "exhausted-budget");
            const int length = gate == 0 ? 15 : gate == 2 ? 1500 : 128;
            c.graph = makeGraph("linear", length, 24);
            c.query = std::string(c.graph.nodeSeq(0)).substr(0, length);
            c.requireSparse = false; c.requireAttempt = gate == 3;
            if (gate == 1) c.match = 4;
            if (gate == 3) c.query.assign(length, 'N');
            runProof(c, out); ++cases;
        }
        std::cout << "{\"cases\":" << cases << ",\"words\":" << out.words << ",\"checksum\":" << out.hash << "}\n";
        return 0;
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}


void checkSearchOutcomes()
{
    using namespace graphalign::dagAligner;
    // Raw row-order edges prevent the compact single-node representation.
    std::vector<int> ids(128, 0), predecessors;
    std::vector<size_t> index;
    for (int row = 0; row < 128; ++row) { index.push_back(predecessors.size()); predecessors.push_back(row - 1); }
    index.push_back(predecessors.size());
    const EdgeMap edges(ids, index, predecessors);
    const std::string target = dna(128, 71231), query = target.substr(0, 64);
    PinnedBudgetAlignMatrix matrix(BaseMatchingPenaltyMatrix(5, -4), -8, -2);
    BudgetSearch search;
    matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges, 0, 0, &search);
    if (search.outcome != (graphtools::verify::kBudgetSearchNever ? BudgetSearch::kFailed : BudgetSearch::kFound)) throw std::runtime_error("Sparse activation witness did not find an endpoint");
    search.skip = true;
    matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges, 0, 0, &search);
    const auto expected = graphtools::verify::kBudgetSearchAlways ? BudgetSearch::kFound : BudgetSearch::kSkipped;
    if (search.outcome != expected) throw std::runtime_error("Skip outcome differs");
    const std::string shortQuery = query.substr(0, 15);
    matrix.init(shortQuery.begin(), shortQuery.end(), target.begin(), target.end(), edges, 0, 0, &search);
    if (search.outcome != BudgetSearch::kNotAttempted) throw std::runtime_error("Eligibility outcome differs");
}
int main()
{
    try
    {
        checkSearchOutcomes();
        Digest out; out.dump = true;
        const auto count = verifyPinned(out) + verifyAlgorithmCases(out);
        std::cout << "pinned-cases " << count << '\n';
        for (const auto& shape : {"linear", "join", "repeat"})
            for (const auto& storage : {"owned", "mapped"})
                for (int errors : {0, 1, 2, 8})
                    lowErrorBenchmark(storage, "mixed", shape, 128, 152, 0, errors, 64, true);
        return activationCases();
    }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
