// Build this unchanged source separately against each candidate's production
// include directory. `verify --dump` emits exact canonical observable results;
// compare the files byte-for-byte. `bench` is for sequential paired runs only.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "graphalign/DagAlignerAffine.hh"
#include "graphalign/dagAligner/BaseMatchingPenaltyMatrix.hh"

namespace da = graphalign::dagAligner;
using Score = da::Score;
using Penalty = da::BaseMatchingPenaltyMatrix;
template <bool penalize> using Matrix = da::AffineAlignMatrixVectorized<Penalty, penalize>;
template <bool penalize> using Aligner = da::Aligner<Matrix<penalize>, false>;

struct Digest
{
    uint64_t hash = UINT64_C(1469598103934665603);
    uint64_t words = 0;
    bool dump = false;
    void add(uint64_t x)
    {
        ++words;
        if (dump) std::cout << x << '\n';
        for (int i = 0; i < 8; ++i)
        {
            hash ^= (x >> (8 * i)) & 255;
            hash *= UINT64_C(1099511628211);
        }
    }
    void text(const std::string& s)
    {
        add(s.size());
        for (unsigned char c : s) add(c);
    }
};

uint32_t nextRandom(uint32_t& state)
{
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

std::string sequence(int len, uint32_t seed, const std::string& alphabet)
{
    std::string ret;
    ret.reserve(len);
    for (int i = 0; i < len; ++i) ret.push_back(alphabet[nextRandom(seed) % alphabet.size()]);
    return ret;
}

da::EdgeMap makeEdges(int len, int shape)
{
    std::vector<int> nodeIds;
    std::vector<std::size_t> index(1, 0);
    std::vector<int> predecessors;
    for (int t = 0; t < len; ++t)
    {
        nodeIds.push_back(t / 4);
        if (shape == 5)
        {
            predecessors.push_back(t);
            if (t % 2) predecessors.push_back(t - 1);
        }
        else if (shape == 6)
        {
            predecessors.push_back((t + 1) % len);
            if (t % 2) predecessors.push_back(t - 1);
        }
        else if (!(shape == 2 && t % 7 == 0) && shape != 4)
        {
            predecessors.push_back(t - 1);
            if (shape == 1 && t >= 3 && t % 5 == 0)
            {
                predecessors.push_back(t - 3);
                predecessors.push_back(-1);
            }
            if (shape == 3 && t >= 2)
            {
                predecessors.push_back(t - 2);
                if (t % 3 == 0) predecessors.push_back(t - 2); // preserve duplicates/order
            }
        }
        index.push_back(predecessors.size());
    }
    return da::EdgeMap(nodeIds, index, predecessors);
}

template <bool local, bool penalize>
void inspectBest(const Matrix<penalize>& matrix, Digest& out, bool full)
{
    auto start = matrix.alignBegin();
    int seen = 0;
    while (start != matrix.alignEnd())
    {
        Score score = 123;
        const auto best = matrix.template nextBestAlign<local>(start, score);
        out.add(score);
        out.add(std::distance(matrix.alignBegin(), best));
        if (best == matrix.alignEnd()) break;
        out.add(matrix.queryOffset(best));
        out.add(matrix.targetOffset(best));
        if (!full || ++seen == 12) break;
        start = best + 1;
    }
}

template <bool local, bool penalize>
void inspectTrace(const Aligner<penalize>& aligner, const da::EdgeMap& edges, Digest& out)
{
    std::vector<da::Cigar> cigars;
    Score second = 314;
    try
    {
        const Score score = aligner.template backtrackAllPaths<local>(edges, cigars, second);
        out.add(1);
        out.add(score);
        out.add(second);
        out.add(cigars.size());
        for (const auto& cigar : cigars)
        {
            out.add(cigar.length());
            for (const auto& operation : cigar)
            {
                out.add(operation.code_);
                out.add(operation.value_);
            }
        }
    }
    catch (const std::exception& error)
    {
        // Extreme legacy 16-bit scores can leave no valid traceback. Compare
        // that observable failure too, rather than silently excluding the case.
        out.add(0);
        out.text(error.what());
        out.add(second);
        out.add(cigars.size());
    }
}

template <bool penalize>
void inspectCase(Matrix<penalize>& matrix, Aligner<penalize>& aligner,
    const std::string& query, const std::string& target, const da::EdgeMap& edges,
    Digest& out, bool traceback)
{
    out.add(penalize);
    out.text(query);
    out.text(target);
    matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges);
    // Includes interior V scores, inter-row left-boundary cells, and padding.
    for (auto it = matrix.alignBegin(); it != matrix.alignEnd(); ++it) out.add(*it);
    inspectBest<false>(matrix, out, true);
    inspectBest<true>(matrix, out, true);
    for (int t = -1; t < int(target.size()); ++t)
    {
        for (int q = -1; q < int(query.size()); ++q)
        {
            out.add(matrix.isInsertion(q, t));
            if (t < 0) continue;
            for (auto pred = edges.prevNodesBegin(t); pred != edges.prevNodesEnd(t); ++pred)
            {
                out.add(matrix.isDeletion(q, t, *pred));
                out.add(matrix.isMatch(q, t, *pred));
                out.add(matrix.isMismatch(q, t, *pred));
            }
        }
    }
    if (traceback)
    {
        aligner.align(query.begin(), query.end(), target.begin(), target.end(), edges);
        inspectTrace<false>(aligner, edges, out);
        inspectTrace<true>(aligner, edges, out);
    }
}

template <bool penalize> uint64_t verifyCases(Digest& out, bool exhaustive)
{
    const std::string queryAlphabet = "ACGTNacgt!?";
    const std::string targetAlphabet = "ACGTXBDHKMNRSVWYacgt!?";
    const std::array<std::array<int, 4>, 13> scoring = {{
        {{2, -3, -4, -1}}, {{5, -4, -8, -2}}, {{0, 0, 0, 0}},
        {{32767, -32768, -32768, -1}}, {{-32768, 32767, 32767, 1}},
        {{1, -1, -32768, -32768}}, {{32767, 32767, 32767, 32767}},
        {{-32768, -32768, -32768, -32768}}, {{-1, 1, 1, -1}},
        {{2, -2, 32767, -32768}}, {{2, -2, -32768, 32767}},
        {{2, -2, 32766, 2}}, {{2, -2, -32767, -2}}
    }};
    const int sizes[] = {1, 15, 16, 17, 2, 63, 64, 65, 31, 32, 33, 150, 3, 250, 7};
    uint64_t cases = 0;
    for (const auto& scores : scoring)
    {
        for (int s : scores) out.add(s);
        Matrix<penalize> matrix(Penalty(scores[0], scores[1]), scores[2], scores[3]);
        Aligner<penalize> aligner(Penalty(scores[0], scores[1]), scores[2], scores[3]);
        for (int shape = 0; shape < 5; ++shape)
        {
            for (int pass = 0; pass < 2; ++pass)
            {
                for (std::size_t i = 0; i < sizeof(sizes) / sizeof(*sizes); ++i)
                {
                    const int qLen = sizes[i];
                    const int tLen = 1 + (int(i) * 11 + 7 * pass) % 37;
                    auto query = sequence(qLen, 73 + i + pass * 97, queryAlphabet);
                    auto target = sequence(tLen, 47 + i + pass * 59, targetAlphabet);
                    // Force a repeated one-code target as well as every ambiguity
                    // code, after changing shapes and sizes on the same objects.
                    if (pass == 1 && i % 3 == 0) target.assign(tLen, targetAlphabet[i]);
                    const auto edges = makeEdges(tLen, shape);
                    inspectCase(matrix, aligner, query, target, edges, out, qLen <= 33);
                    ++cases;
                }
            }
        }
        // Exact ties and the public maxRepeats/order behavior with many joins.
        const auto ties = makeEdges(13, 3);
        inspectCase(matrix, aligner, "AAAAAAAAAAAAA", "AAAAAAAAAAAAA", ties, out, true);
        ++cases;
        // All byte encodings, including NUL and bytes >=128, through the actual
        // unsigned-char translation routines rather than a simplified alphabet.
        std::string allBytes;
        for (int i = 0; i < 256; ++i) allBytes.push_back(char(i));
        inspectCase(matrix, aligner, allBytes, allBytes, makeEdges(256, 1), out, false);
        ++cases;
        // Invalid topological order is not a DAG, but exercise candidate3's
        // full-initialization fallback for all addressable self/future edges.
        // Do not traceback cycles: that legacy API need not terminate on them.
        for (int shape = 5; shape <= 6; ++shape)
        {
            for (int qLen : {1, 17, 7, 65, 15})
            {
                inspectCase(matrix, aligner, sequence(qLen, 79, queryAlphabet),
                    "ARNTGCTA", makeEdges(8, shape), out, false);
                ++cases;
            }
        }
    }
    if (exhaustive)
    {
        // Sweep every representable Score independently through match,
        // mismatch, gap-open, and gap-extension. Cross-combinations at the
        // overflow boundaries are covered by the scoring cases above.
        for (int field = 0; field < 4; ++field)
        {
            for (int value = -32768; value <= 32767; ++value)
            {
                int scores[] = {2, -3, -4, -1};
                scores[field] = value;
                Matrix<penalize> matrix(Penalty(scores[0], scores[1]), scores[2], scores[3]);
                Aligner<penalize> aligner(Penalty(scores[0], scores[1]), scores[2], scores[3]);
                const auto edges = makeEdges(3, value % 2 == 0 ? 0 : 3);
                inspectCase(matrix, aligner, "ACN", "ARN", edges, out, false);
                inspectCase(matrix, aligner, "T", "T", makeEdges(1, 0), out, false);
                cases += 2;
            }
        }
    }
    return cases;
}

void benchmark(const std::string& shape, int qLen, int tLen, int repetitions, bool trace)
{
    if (qLen < 1 || tLen < 1 || repetitions < 1) throw std::invalid_argument("positive sizes required");
    const int edgeShape = shape == "join" ? 1 : shape == "empty" ? 2 : shape == "dense" ? 3 : 0;
    const std::string alphabet = shape == "degenerate" ? "ACGTXBDHKMNRSVWY" : "ACGT";
    const auto edges = makeEdges(tLen, edgeShape);
    const auto target = sequence(tLen, 79, alphabet);
    std::vector<std::string> queries;
    for (int i = 0; i < 64; ++i) queries.push_back(sequence(qLen, 97 + i * 17, "ACGTNacgt"));
    Matrix<true> matrix(Penalty(5, -4), -8, -2);
    Aligner<true> aligner(Penalty(5, -4), -8, -2);
    Digest out;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repetitions; ++i)
    {
        const auto& query = queries[i % queries.size()];
        if (trace)
        {
            aligner.align(query.begin(), query.end(), target.begin(), target.end(), edges);
            inspectTrace<false>(aligner, edges, out);
        }
        else
        {
            matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges);
            inspectBest<false>(matrix, out, false);
            out.add(matrix.isInsertion(qLen - 1, tLen - 1));
        }
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "{\"mode\":\"bench\",\"shape\":\"" << shape << "\",\"query\":" << qLen
              << ",\"target\":" << tLen << ",\"repetitions\":" << repetitions
              << ",\"trace\":" << (trace ? "true" : "false")
              << ",\"checksum\":" << out.hash << ",\"words\":" << out.words << "}\n";
    std::cerr << "{\"kernel_seconds\":" << seconds << "}\n";
}

void memorySequence(int scale, int repetitions)
{
    if (scale < 1 || scale > 100000 || repetitions < 1)
        throw std::invalid_argument("positive bounded scale/repetitions required");
    // Adversarial profile-capacity history: skipping an unused profile's resize
    // would let it grow 16->32->64 while baseline grows16->48 and retains48.
    const std::string targets[] = {"ACGTXBDHKMNRSVWY", "A", "CGTXBDHKMNRSVWY", "CGTXBDHKMNRSVWY"};
    const int queryLengths[] = {16, 48, 32, 48};
    Matrix<true> matrix(Penalty(5, -4), -8, -2);
    Digest out;
    for (int repeat = 0; repeat < repetitions; ++repeat)
    {
        for (int stage = 0; stage < 4; ++stage)
        {
            const auto query = sequence(queryLengths[stage] * scale, 79 + stage, "ACGT");
            const auto& target = targets[stage];
            const auto edges = makeEdges(target.size(), 0);
            matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges);
            inspectBest<false>(matrix, out, false);
            out.add(matrix.isInsertion(query.size() - 1, target.size() - 1));
        }
    }
    std::cout << "{\"mode\":\"memory\",\"scale\":" << scale << ",\"repetitions\":" << repetitions
              << ",\"checksum\":" << out.hash << ",\"words\":" << out.words << "}\n";
}

int main(int argc, char** argv)
{
    try
    {
        if (argc >= 2 && std::string(argv[1]) == "verify")
        {
            Digest out;
            bool exhaustive = false;
            for (int i = 2; i < argc; ++i)
            {
                if (std::string(argv[i]) == "--dump") out.dump = true;
                else if (std::string(argv[i]) == "--exhaustive") exhaustive = true;
                else throw std::invalid_argument("unknown verify flag");
            }
            const auto a = verifyCases<true>(out, exhaustive);
            const auto b = verifyCases<false>(out, exhaustive);
            std::cout << "{\"mode\":\"verify\",\"cases\":" << a + b << ",\"checksum\":" << out.hash
                      << ",\"words\":" << out.words << ",\"exhaustive\":" << (exhaustive ? "true" : "false") << "}\n";
            return 0;
        }
        if (argc >= 6 && std::string(argv[1]) == "bench")
        {
            benchmark(argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]),
                argc >= 7 && std::string(argv[6]) == "trace");
            return 0;
        }
        if (argc == 4 && std::string(argv[1]) == "memory")
        {
            memorySequence(std::stoi(argv[2]), std::stoi(argv[3]));
            return 0;
        }
        std::cerr << "usage: dp-probe verify [--dump] [--exhaustive]\n"
                     "       dp-probe bench linear|join|dense|empty|degenerate Q T REP [trace]\n"
                     "       dp-probe memory SCALE REP\n";
        return 2;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
