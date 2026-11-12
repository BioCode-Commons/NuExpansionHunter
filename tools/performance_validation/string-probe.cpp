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
// Frozen original implementation, with only its class names changed. It uses
// the same graph/scoring objects as the selected production implementation.
#include "reference-pinned-dag.hh"

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

static int referenceCount(const std::string& query, const KmerIndex& index)
{
    int count = 0;
    const size_t k = index.kmerLength();
    for (size_t pos = 0; pos + k <= query.size();)
    {
        std::string word = query.substr(pos, k);
        std::transform(word.begin(), word.end(), word.begin(), ::toupper);
        if (index.contains(word)) { ++count; pos += k; }
        else ++pos;
    }
    return count;
}

static OrientationPrediction referencePredict(const std::string& q, const KmerIndex& index, int threshold)
{
    const int f = referenceCount(q, index);
    const int r = referenceCount(reverseComplement(q), index);
    if (std::max(f, r) < threshold) return OrientationPrediction::kDoesNotAlign;
    return f >= r ? OrientationPrediction::kAlignsInOriginalOrientation
                  : OrientationPrediction::kAlignsInReverseComplementOrientation;
}

static std::vector<std::string> reads(const std::string& workload)
{
    std::vector<size_t> lengths = { 35, 75, 100, 150, 250, 300 };
    if (workload == "short") lengths = { 1, 9, 10, 15, 16, 17, 25 };
    if (workload == "long") lengths = { 500, 1000, 5000 };
    std::mt19937 rng(78191);
    std::vector<std::string> result;
    const std::string alphabet = "ACGTacgtRYSWKMBDHVNryswkmbdhvnXx?";
    for (size_t n : lengths)
        for (int i = 0; i < 32; ++i)
        {
            const int kind = workload == "forward" ? 0 : workload == "reverse" ? 1
                : workload == "invalid" ? 3 : i % 4;
            std::string q = repeat(kind == 1 ? "CTG" : "CAG", n);
            if (kind == 2)
            {
                for (size_t j = 0; j < n; j += 13) q[j] = "ACGTN"[rng() % 5];
                for (size_t j = 1; j < n; j += 5) q[j] = static_cast<char>(::tolower(q[j]));
            }
            if (kind == 3) for (char& c : q) c = alphabet[rng() % alphabet.size()];
            result.push_back(q);
        }
    return result;
}

static size_t checkOrientation()
{
    size_t cases = 0;
    const Graph graph = makeStrGraph("ACGTTGCA", "CAG", "GATCCGTA");
    auto queries = reads("short");
    const auto regular = reads("mixed");
    queries.insert(queries.end(), regular.begin(), regular.end());
    queries.push_back("");
    for (int byte = 0; byte < 128; ++byte)
    {
        std::string q = repeat("ctgCTG", 51);
        q[17] = static_cast<char>(byte);
        queries.push_back(q);
    }
    for (int k : { 1, 2, 3, 10, 16 })
    {
        const KmerIndex index(graph, k, true);
        for (int threshold : { 0, 1, 3, 50 })
        {
            OrientationPredictor predictor(&graph, k, threshold);
            for (const auto& q : queries)
            {
                const auto expected = referencePredict(q, index, threshold);
                const auto actual = predictor.predict(q);
                if (expected != actual) throw std::runtime_error("orientation mismatch k=" + std::to_string(k));
                mix(static_cast<uint64_t>(actual));
                ++cases;
            }
        }
    }
    return cases;
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
            ReferencePinnedDagAligner reference(match, -4, -8, -2);
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

int main(int argc, char** argv)
{
    try
    {
        const std::string mode = argc > 1 ? argv[1] : "check";
        if (mode == "check")
        {
            const size_t orientation = checkOrientation();
            const size_t suffix = checkSuffix();
            std::cout << "{\"orientation_cases\":" << orientation << ",\"suffix_cases\":" << suffix
                      << ",\"checksum\":" << digest << "}\n";
            return 0;
        }
        const std::string workload = argc > 2 ? argv[2] : "mixed";
        const size_t repetitions = argc > 3 ? std::stoull(argv[3]) : 100;
        size_t calls = 0;
        double elapsed = 0;
        if (mode == "orientation")
        {
            const auto graph = makeStrGraph("ACGTTGCA", "CAG", "GATCCGTA");
            const OrientationPredictor predictor(&graph, 10, 3);
            const auto queries = reads(workload);
            const auto start = std::chrono::steady_clock::now();
            for (size_t rep = 0; rep < repetitions; ++rep)
                for (const auto& q : queries) { mix(static_cast<uint64_t>(predictor.predict(q))); ++calls; }
            elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        }
        else if (mode == "suffix")
        {
            const auto graph = suffixGraph(workload == "branch" ? 2 : 0);
            const auto seed = suffixSeed(graph);
            PinnedDagAligner aligner(5, -4, -8, -2);
            const auto queries = suffixQueries(workload);
            const auto start = std::chrono::steady_clock::now();
            for (size_t rep = 0; rep < repetitions; ++rep)
                for (const auto& q : queries)
                {
                    int score = INT32_MIN;
                    const auto result = aligner.suffixAlign(seed, q, q.size() + 12, score);
                    mixAlignments(score, result); ++calls;
                }
            elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        }
        else throw std::runtime_error("mode must be check, orientation, or suffix");
        std::cout << "{\"mode\":\"" << mode << "\",\"workload\":\"" << workload << "\",\"calls\":" << calls
                  << ",\"checksum\":" << digest << "}\n";
        std::cerr << "loop_seconds=" << elapsed << '\n';
    }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
