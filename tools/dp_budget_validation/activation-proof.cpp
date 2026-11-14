// A proof witness must use the sparse algorithm, not only match via fallback.
// Each case checks ordered full production output and its own activation delta.
#define main inheritedPinnedProbeMain
#include "pinned-probe.cpp"
#undef main

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

void runProof(const ProofCase& c, Digest& out)
{
    PinnedDagAligner aligner(c.match, c.mismatch, c.gapOpen, c.gapExt);
    const auto seed = makeSeed(c.graph, c.suffix, 0);
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
    const auto before = graphalign::dagAligner::PinnedBudgetAlignMatrix::stats();
#endif
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
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
    const auto after = graphalign::dagAligner::PinnedBudgetAlignMatrix::stats();
    const uint64_t calls = after.calls - before.calls;
    const uint64_t attempts = after.attempts - before.attempts;
    const uint64_t successes = after.succeeded - before.succeeded;
    const uint64_t fallbacks = after.fallbacks - before.fallbacks;
    if (calls != 1 || successes != (c.requireSparse ? 1u : 0u) || fallbacks != (c.requireSparse ? 0u : 1u)
        || bool(attempts) != c.requireAttempt)
        throw std::runtime_error(c.name + ": expected activation/fallback gate was not exercised");
    std::cerr << "{\"case\":\"" << c.name << "\",\"calls\":" << calls << ",\"attempts\":" << attempts
              << ",\"successes\":" << successes << ",\"fallbacks\":" << fallbacks
              << ",\"dense_cells\":" << after.denseCells - before.denseCells
              << ",\"band_cells\":" << after.bandCells - before.bandCells
              << ",\"inserted\":" << inserted << ",\"deleted\":" << deleted << "}\n";
#endif
}

int main()
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
