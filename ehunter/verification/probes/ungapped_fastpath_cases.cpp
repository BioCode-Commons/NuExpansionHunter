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

// Verification build: -DEH_VERIFY="FAST_PATHS=OFF;SIMD=SCALAR".
// Differential production-wrapper witnesses for budget-9 SIMD and gap fallback.
// Compare complete stdout between default and verification builds.
#define main pinnedCaseEntry
#include "pinned_align_cases.cpp"
#undef main

void inspectUngappedCase(CheckedPinnedAligner& aligner, const Graph& graph, const std::string& query,
                        bool suffix, const std::string& name, Digest& output)
{
    inspectExplicit(aligner, graph, query, int(query.size()) + 8, suffix, 0, name, output);
}

int main()
{
    try
    {
        Digest output;
        output.dump = true;
        unsigned cases = 0;
        CheckedPinnedAligner aligner(5, -4, -8, -2);
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
        // Exercise the new nonlinear budget-9 propagation as well, including
        // more than ten tied paths. Query length 128 avoids the unrelated
        // reference dense-padding issue if a baseline chooses that fallback.
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
        std::cout << "{\"cases\":" << cases << ",\"words\":" << output.words
                  << ",\"checksum\":" << output.hash << "}\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
