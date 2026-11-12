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

#include "graphalign/PinnedDagAligner.hh"

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

void verifyCase(PinnedDagAligner& aligner, const Graph& graph, const std::string& shape,
                int scoring, int length, int extension, bool suffix, unsigned seedMode,
                unsigned sample, const std::string& storage, size_t nodeLength, Digest& out)
{
    if (out.dump) std::cout << "CASE inherited " << storage << ' ' << shape << ' ' << scoring << ' '
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

uint64_t verify(Digest& out)
{
    uint64_t cases = 0;
    const int lengths[] = {1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 150, 250, 8};
    for (const std::string shape : {"linear", "join", "repeat"})
        for (int scoring : {0, 1})
        {
            PinnedDagAligner aligner(scoring ? 3 : 5, scoring ? -2 : -4, scoring ? -5 : -8, -2);
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
                PinnedDagAligner aligner(5, -4, -8, -2);
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
            PinnedDagAligner aligner(0, 0, 0, 0);
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

void inspectExplicit(PinnedDagAligner& aligner, const Graph& graph, const std::string& query,
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
                PinnedDagAligner aligner(scoring == 1 ? 4 : 5, scoring == 1 ? -3 : -4,
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
    // Unlike the inherited zero-score fallback cases, these use default scores.
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

    // Explicit overflow-guard/fallback dimensions; preserve legacy exceptions.
    for (int length : {1200, 1300, 1320, 2200, 3300})
    {
        Graph graph(1); const std::string reference = dna(length + 8, 381); graph.setNodeSeq(0, reference);
        PinnedDagAligner aligner(5, -4, -8, -2);
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
    PinnedDagAligner aligner(5, -4, -8, -2);
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


void reportBudgetStats()
{
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
    const auto& s = graphalign::dagAligner::PinnedBudgetAlignMatrix::stats();
    std::cerr << "{\"budget_dp\":{\"calls\":" << s.calls << ",\"attempts\":" << s.attempts
              << ",\"succeeded\":" << s.succeeded << ",\"fallbacks\":" << s.fallbacks
              << ",\"dense_cells\":" << s.denseCells << ",\"band_cells\":" << s.bandCells << "}}\n";
#endif
}

int main(int argc, char** argv)
{
    try
    {
        if (argc >= 2 && std::string(argv[1]) == "verify")
        {
            Digest out;
            if (argc == 3 && std::string(argv[2]) == "--dump") out.dump = true;
            else if (argc != 2) throw std::invalid_argument("unknown verify flag");
            const uint64_t cases = verify(out) + verifyAlgorithmCases(out);
            std::cout << "{\"cases\":" << cases << ",\"words\":" << out.words << ",\"checksum\":" << out.hash << "}\n";
            reportBudgetStats();
            return 0;
        }
        if (argc == 9 && std::string(argv[1]) == "bench")
        { benchmark(argv[2], argv[3], argv[4], std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stoi(argv[8])); return 0; }
        throw std::invalid_argument("usage: graph-build-probe verify [--dump] | bench owned|mapped prefix|suffix|mixed linear|join|repeat QLEN EXTENSION NODE_LENGTH REP");
    }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
