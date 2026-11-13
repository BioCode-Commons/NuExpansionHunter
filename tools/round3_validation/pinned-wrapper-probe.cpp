// Compile this identical source against baseline and candidate GraphTools.
// Compare complete stdout, including result order, exceptions and score state.
// This is a differential correctness probe, not a microbenchmark.
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "graphalign/PinnedDagAligner.hh"

using namespace graphtools;

namespace
{
uint64_t cases = 0;

std::string hex(const std::string& text)
{
    const char* digits = "0123456789abcdef";
    std::string output;
    output.reserve(2 * text.size());
    for (unsigned char c : text)
    {
        output += digits[c >> 4];
        output += digits[c & 15];
    }
    return output;
}

std::string dna(size_t size, uint32_t state)
{
    std::string sequence(size, 'A');
    for (char& base : sequence)
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        base = "ACGT"[state % 4];
    }
    return sequence;
}

void inspect(PinnedDagAligner& aligner, const Path& seed, const std::string& query, size_t extension,
             bool suffix, const std::string& tag)
{
    std::cout << "CASE " << ++cases << ' ' << tag << ' ' << suffix << ' ' << extension << ' '
              << seed.encode() << ' ' << hex(query) << '\n';
    int score = -773;
    try
    {
        const auto results = suffix ? aligner.suffixAlign(seed, query, extension, score)
                                    : aligner.prefixAlign(seed, query, extension, score);
        std::cout << "OK " << score << ' ' << results.size() << '\n';
        for (const auto& result : results)
        {
            const Path& path = result.first;
            const Alignment& alignment = result.second;
            std::cout << "PATH " << path.encode() << ' ' << hex(path.seq()) << ' '
                      << path.startPosition() << ' ' << path.endPosition() << ' ' << path.numNodes();
            for (NodeId node : path.nodeIds())
                std::cout << ' ' << node;
            std::cout << "\nALIGN " << alignment.referenceStart() << ' ' << alignment.referenceLength() << ' '
                      << alignment.queryLength() << ' ' << alignment.numOperations() << ' '
                      << alignment.numMatched() << ' ' << alignment.numMismatched() << ' '
                      << alignment.numClipped() << ' ' << alignment.numInserted() << ' '
                      << alignment.numDeleted() << ' ' << alignment.generateCigar() << '\n';
            for (const auto& operation : alignment)
                std::cout << "OP " << int(operation.type()) << ' ' << operation.length() << ' '
                          << operation.referenceLength() << ' ' << operation.queryLength() << ' '
                          << operation.generateCigar() << '\n';
        }
    }
    catch (const std::exception& error)
    {
        std::cout << "EXCEPTION " << score << ' ' << hex(error.what()) << '\n';
    }
}

std::vector<std::string> mutations(const std::string& exact)
{
    std::vector<std::string> queries { exact };
    for (size_t position : { size_t(0), exact.size() / 2, exact.size() - 1 })
    {
        std::string query = exact;
        query[position] = static_cast<char>(0xff);
        queries.push_back(query);
    }
    std::string query = exact;
    query[0] = query[query.size() - 1] = 'N';
    queries.push_back(query);
    query = exact;
    query.insert(query.size() / 2, "G");
    queries.push_back(query);
    query = exact;
    query.erase(query.size() / 2, 1);
    queries.push_back(query);
    return queries;
}

template <typename Function> void bothStorage(const Graph& owned, Function function)
{
    BlobBuilder blob;
    const auto record = owned.freeze(blob);
    const Graph mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);
    function(owned, "owned");
    function(mapped, "mapped");
}

void terminalCases()
{
    Graph owned(2);
    owned.setNodeSeq(0, dna(384, 913));
    owned.setNodeSeq(1, dna(384, 531));
    owned.addEdge(0, 1);
    bothStorage(owned, [](const Graph& graph, const std::string& storage) {
        PinnedDagAligner aligner(5, -4, -8, -2);
        for (bool suffix : { false, true })
        {
            const Path seed = suffix ? Path(&graph, 352, { 0, 1 }, 7) : Path(&graph, 7, { 0, 1 }, 32);
            const std::string node = graph.nodeSeq(suffix ? 0 : 1);
            for (size_t length : { 16u, 17u, 31u, 32u, 33u, 63u, 64u, 150u })
            {
                const std::string exact = node.substr(suffix ? 352 - length : 32, length);
                for (size_t extension : { size_t(1), length, size_t(352), size_t(353), size_t(552) })
                    for (const auto& query : mutations(exact))
                        inspect(aligner, seed, query, extension, suffix, "terminal/" + storage);
            }
            // An empty query and an empty terminal target preserve score and
            // exception behavior even after a successful compact-wrapper call.
            inspect(aligner, seed, "", 400, suffix, "empty-query/" + storage);
            const Path boundary = suffix ? Path(&graph, 0, { 0, 1 }, 7) : Path(&graph, 7, { 0, 1 }, 384);
            inspect(aligner, boundary, "ACGTACGTACGTACGT", 400, suffix, "empty-target/" + storage);
        }
    });
}

void rangeCases()
{
    // Two equal repeat branches feed another repeat, a join, and an optional
    // downstream repeat. Reverse discovery must retain its ordered traversal
    // even though it follows the original NodeIds in the opposite direction.
    Graph owned(7);
    owned.setNodeSeq(0, "ACGTGATCGTACGATTACGATCGT");
    owned.setNodeSeq(1, "CA");
    owned.setNodeSeq(2, "CA");
    owned.setNodeSeq(3, "CAG");
    owned.setNodeSeq(4, "TACG");
    owned.setNodeSeq(5, "CCG");
    owned.setNodeSeq(6, "TAGCATCGATACGTAGCTAGTCGA");
    for (const auto edge : { std::make_pair(0, 1), { 0, 2 }, { 0, 3 }, { 0, 4 }, { 1, 1 }, { 2, 2 },
                             { 1, 3 }, { 2, 3 }, { 3, 3 }, { 3, 4 }, { 1, 4 }, { 2, 4 },
                             { 4, 5 }, { 5, 5 }, { 5, 6 }, { 4, 6 } })
        owned.addEdge(edge.first, edge.second);
    bothStorage(owned, [](const Graph& graph, const std::string& storage) {
        for (int scoring : { 0, 1 })
        {
            PinnedDagAligner aligner(5, -4, scoring ? -1 : -8, scoring ? -1 : -2);
            for (int a : { 0, 1, 7 })
                for (int b : { 0, 2 })
                    for (int c : { 0, 3 })
                    {
                        std::string exact = graph.nodeSeq(0);
                        for (int i = 0; i < a; ++i) exact += "CA";
                        for (int i = 0; i < b; ++i) exact += "CAG";
                        exact += "TACG";
                        for (int i = 0; i < c; ++i) exact += "CCG";
                        exact += graph.nodeSeq(6);
                        for (bool suffix : { false, true })
                            for (int seedMode : { 0, 1, 2 })
                            {
                                const int tailLength = graph.nodeSeqLength(6);
                                const Path seed = suffix
                                    ? (seedMode == 2 ? Path(&graph, 3, { 5, 6 }, 4)
                                                    : Path(&graph, tailLength - (seedMode ? 2 : 0), { 6 }, tailLength))
                                    : (seedMode == 2 ? Path(&graph, 4, { 0, 1 }, 0)
                                                    : Path(&graph, 0, { 0 }, seedMode ? 2 : 0));
                                std::string substitution = exact;
                                substitution[substitution.size() / 2] = 'N';
                                std::string deletion = exact;
                                deletion.erase(deletion.size() / 2, 1);
                                for (const auto& query : { exact, substitution, deletion })
                                    inspect(aligner, seed, query, query.size() + 12, suffix,
                                        "ranges/" + storage + "/" + std::to_string(scoring));
                            }
                    }
        }
    });
}

void boundaryAndScoringCases()
{
    for (int targetLength : { 16, 17, 3251, 3252 })
    {
        Graph owned(1);
        owned.setNodeSeq(0, dna(targetLength, 81239));
        bothStorage(owned, [&](const Graph& graph, const std::string& storage) {
            // 16*15 + targetLength*10 + 10 is just below/above 32768.
            for (const auto scores : { std::array<int, 4> { 5, -4, -8, -2 },
                                       std::array<int, 4> { 5, -4, -1, -1 },
                                       std::array<int, 4> { 65541, 65532, 65528, 65534 } })
            {
                PinnedDagAligner aligner(scores[0], scores[1], scores[2], scores[3]);
                for (bool suffix : { false, true })
                {
                    const Path seed(&graph, suffix ? targetLength : 0, { 0 }, suffix ? targetLength : 0);
                    std::string query = graph.nodeSeq(0).substr(suffix ? targetLength - 16 : 0, 16);
                    query[15] = 'N';
                    for (size_t extension : { size_t(0), size_t(1), size_t(targetLength + 12) })
                        inspect(aligner, seed, query, extension, suffix,
                            "score-bound/" + storage + "/" + std::to_string(scores[0]) + "/"
                                + std::to_string(scores[2]));
                }
            }
        });
    }
    Graph oracle(1);
    std::string target(64, 'A');
    target[16] = 'C';
    oracle.setNodeSeq(0, target);
    bothStorage(oracle, [](const Graph& graph, const std::string& storage) {
        PinnedDagAligner aligner(5, -4, -1, -1);
        inspect(aligner, Path(&graph, 0, { 0 }, 0), std::string(32, 'A'), 32, false,
            "nondefault-unit-test-oracle/" + storage);
    });
}

void smallGraphBoundaryCases()
{
    for (int nodeCount : { 63, 64, 65 })
    {
        Graph owned(nodeCount);
        for (int node = 0; node != nodeCount; ++node)
            owned.setNodeSeq(node, "A");
        const NodeId repeat = nodeCount - 2, tail = nodeCount - 1;
        owned.setNodeSeq(0, dna(16, 7241));
        owned.setNodeSeq(1, dna(40, 523));
        owned.setNodeSeq(2, "TT");
        owned.setNodeSeq(repeat, "CAG");
        owned.setNodeSeq(tail, dna(48, 842));
        owned.addEdge(0, 1);
        owned.addEdge(0, 2);
        owned.addEdge(1, repeat);
        owned.addEdge(2, repeat);
        owned.addEdge(repeat, repeat);
        owned.addEdge(repeat, tail);
        // The first branch enqueues the join at distance 40; the second
        // improves it to 2 while the first queue entry is still pending.
        // This exercises duplicate FIFO entries and offset updates. IDs at
        // 63/64 also exercise the stack-array boundary and large-graph fallback.
        bothStorage(owned, [&](const Graph& graph, const std::string& storage) {
            for (int scoring : { 0, 1 })
            {
                PinnedDagAligner aligner(5, -4, scoring ? -1 : -8, scoring ? -1 : -2);
                const std::string walk = std::string(graph.nodeSeq(0)) + "TTCAGCAGCAGCAG"
                    + std::string(graph.nodeSeq(tail));
                for (bool suffix : { false, true })
                    for (int seedMode : { 0, 1 })
                    {
                        const Path seed = suffix ? Path(&graph, seedMode ? 0 : 48, { tail }, 48)
                                                 : Path(&graph, 0, { 0 }, seedMode ? 16 : 0);
                        const std::string exact = !seedMode ? walk : suffix ? walk.substr(0, walk.size() - 48)
                                                                                    : walk.substr(16);
                        std::string mismatch = exact;
                        mismatch[mismatch.size() / 2] = 'N';
                        std::string deletion = exact;
                        deletion.erase(deletion.size() / 2, 1);
                        for (const auto& query : { exact, mismatch, deletion })
                            for (size_t extension : { size_t(1), query.size() + 12 })
                                inspect(aligner, seed, query, extension, suffix,
                                    "small-graph-boundary/" + storage + "/" + std::to_string(nodeCount)
                                        + "/" + std::to_string(scoring));
                    }
            }
        });
    }
}
}

int main()
{
    try
    {
        terminalCases();
        rangeCases();
        boundaryAndScoringCases();
        smallGraphBoundaryCases();
        std::cout << "TOTAL " << cases << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
