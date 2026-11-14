// Differential production-wrapper witnesses for budget-9 SIMD and gap fallback.
// Compile against the original and candidate headers; compare complete stdout.
#define main inheritedPinnedProbeMain
#include "pinned-probe.cpp"
#undef main

void inspectUngappedCase(PinnedDagAligner& aligner, const Graph& graph, const std::string& query,
                        bool suffix, const std::string& name, Digest& output)
{
#if defined(GRAPHTOOLS_UNGAPPED_DP_VALIDATE_SIMD)
    const auto before = graphalign::dagAligner::PinnedBudgetAlignMatrix::stats();
#endif
    inspectExplicit(aligner, graph, query, int(query.size()) + 8, suffix, 0, name, output);
#if defined(GRAPHTOOLS_UNGAPPED_DP_VALIDATE_SIMD)
    const auto after = graphalign::dagAligner::PinnedBudgetAlignMatrix::stats();
    if (after.linearUngappedAttempts == before.linearUngappedAttempts)
        throw std::runtime_error("Linear budget-9 path was not exercised: " + name);
#if defined(GRAPHTOOLS_DP_TEST_FORCE_SCALAR) || defined(GRAPHTOOLS_BUDGET_DP_TEST_FORCE_SCALAR)
    if (after.linearSimdBlocks != before.linearSimdBlocks)
        throw std::runtime_error("Forced scalar unexpectedly executed SIMD");
#elif defined(__aarch64__) && defined(__ARM_NEON)
    if (after.linearSimdBlocks == before.linearSimdBlocks)
        throw std::runtime_error("Native NEON path was not exercised: " + name);
#elif defined(GRAPHTOOLS_DP_TEST_FORCE_AVX2) || defined(GRAPHTOOLS_UNGAPPED_DP_VALIDATE_AVX2)
    if (query.size() >= 32 && after.linearSimdBlocks == before.linearSimdBlocks)
        throw std::runtime_error("AVX2 path was not exercised: " + name);
#endif
#endif
}

int main()
{
    try
    {
        Digest output;
        output.dump = true;
        unsigned cases = 0;
        PinnedDagAligner aligner(5, -4, -8, -2);
        for (int length : { 16, 17, 31, 32, 33, 63, 64, 65, 150, 250 })
        {
            Graph graph(1);
            const std::string target = dna(length + 16, 832773 + length);
            graph.setNodeSeq(0, target);
            for (bool suffix : { false, true })
            {
                const std::string exact = suffix ? target.substr(target.size() - length) : target.substr(0, length);
                inspectUngappedCase(aligner, graph, exact, suffix, "exact", output);
                ++cases;
                // Every mismatch position, including all SIMD lanes and tails.
                for (int position = 0; position < length; ++position)
                {
                    auto query = exact;
                    query[position] = query[position] == 'A' ? 'C' : 'A';
                    inspectUngappedCase(aligner, graph, query, suffix, "one-mismatch", output);
                    ++cases;
                }
#ifndef GRAPHTOOLS_UNGAPPED_DP_SANITIZER_ONLY
                // The sanitizer-only mode isolates the newly added <=1-error
                // paths. The unchanged dense fallback reads logical vector
                // padding and fails libc++ container-overflow instrumentation.
                // Full differential runs retain every fallback witness below.
                for (int position : { 0, length / 2, length - 2 })
                {
                    auto query = exact;
                    query[position] = query[position + 1] = 'N';
                    inspectUngappedCase(aligner, graph, query, suffix, "two-mismatches-fallback", output);
                    ++cases;
                }
                for (int position : { 1, length / 2, length - 2 })
                {
                    auto insertion = exact;
                    insertion.insert(position, "G");
                    inspectUngappedCase(aligner, graph, insertion, suffix, "insertion-fallback", output);
                    ++cases;
                    auto deletion = exact;
                    deletion.erase(position, 1);
                    // Length 15 deliberately takes the unchanged dense gate.
                    if (deletion.size() >= 16)
                        inspectUngappedCase(aligner, graph, deletion, suffix, "deletion-fallback", output);
                    else
                        inspectExplicit(aligner, graph, deletion, int(deletion.size()) + 8, suffix, 0,
                            "short-deletion-fallback", output);
                    ++cases;
                }
#endif
            }
        }
        // Every reference IUPAC mask, both base cases, plus invalid query codes.
        for (char reference : std::string("ACGTBDHKMNRSVWY"))
        {
            char matching = 'A';
            for (char query : std::string("ACGT"))
                if (graphtools::checkIfReferenceBaseMatchesQueryBase(reference, query))
                { matching = query; break; }
            for (int position : { 0, 15, 16, 31, 32, 63, 64 })
            {
                Graph graph(1);
                std::string target(80, matching);
                // Only one degenerate reference base: Graph also builds
                // sequence expansions, so repeated IUPAC nodes explode.
                target[position] = reference;
                graph.setNodeSeq(0, target);
                for (char query : std::string("ACGTNacgtn"))
                {
                    std::string read(65, matching);
                    read[position] = query;
                    inspectUngappedCase(aligner, graph, read, false, "iupac", output);
                    ++cases;
                }
            }
        }
#ifdef GRAPHTOOLS_UNGAPPED_DP_SANITIZER_ONLY
        // Exercise the new nonlinear budget-9 propagation as well, including
        // more than ten tied paths. Query length 128 avoids the unrelated
        // legacy dense-padding issue if a baseline chooses that fallback.
        for (bool suffix : { false, true })
            for (int alternatives : { 2, 20 })
                for (int errors : { 0, 1 })
                {
                    Graph graph(alternatives + 2);
                    graph.setNodeSeq(0, dna(suffix ? 512 : 48, 41377));
                    graph.setNodeSeq(alternatives + 1, dna(suffix ? 48 : 512, 737141));
                    for (int branch = 1; branch <= alternatives; ++branch)
                    {
                        graph.setNodeSeq(branch, dna(48, 11331));
                        graph.addEdge(0, branch);
                        graph.addEdge(branch, alternatives + 1);
                    }
                    const std::string walk = std::string(graph.nodeSeq(0)) + std::string(graph.nodeSeq(1))
                        + std::string(graph.nodeSeq(alternatives + 1));
                    std::string query = suffix ? walk.substr(walk.size() - 128) : walk.substr(0, 128);
                    if (errors)
                        query[20] = query[20] == 'A' ? 'C' : 'A';
                    inspectExplicit(aligner, graph, query, 152, suffix, 0, "sanitized-nonlinear-budget9", output);
                    ++cases;
                }
#endif
        std::cout << "{\"cases\":" << cases << ",\"words\":" << output.words
                  << ",\"checksum\":" << output.hash << "}\n";
#if defined(GRAPHTOOLS_UNGAPPED_DP_VALIDATE_SIMD)
        const auto& stats = graphalign::dagAligner::PinnedBudgetAlignMatrix::stats();
        std::cerr << "{\"ungapped_attempts\":" << stats.ungappedAttempts
                  << ",\"linear_ungapped_attempts\":" << stats.linearUngappedAttempts
                  << ",\"linear_simd_blocks\":" << stats.linearSimdBlocks << "}\n";
#endif
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
