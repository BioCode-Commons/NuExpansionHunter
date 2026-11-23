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

// Verification build: -DEH_VERIFY="BUDGET_SEARCH=ALWAYS;FUSED_ROWS=OFF".
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "graphalign/dagAligner/PinnedBudgetAlignMatrix.hh"

namespace lifecycle
{
// Direct matrix lifecycle and malformed-edge-map differential checks.

using namespace graphalign::dagAligner;

// Report the matrix extent, best score and every tied optimal endpoint in traversal order.
// The sparse matrix may omit endpoints whose costs exceed its successful budget.
void inspect(const char* name, PinnedBudgetAlignMatrix& matrix, const std::string& query,
             const std::string& target, const EdgeMap& edges)
{
    std::cout << name << ' ' << query.size() << ' ' << target.size() << '\n';
    try
    {
        matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges);
        Score score = -1973;
        const auto best = matrix.nextBestAlign<false>(matrix.alignBegin(), score);
        std::cout << "best " << best << ' ' << score << ' ' << matrix.alignEnd() << '\n';
        auto endpoint = best;
        while (endpoint != matrix.alignEnd())
        {
            Score nextScore = -1973;
            const auto next = matrix.nextBestAlign<false>(endpoint + 1, nextScore);
            if (next == matrix.alignEnd() || nextScore != score)
                break;
            std::cout << "optimal-tie " << next << '\n';
            endpoint = next;
        }
    }
    catch (const std::exception& error)
    {
        std::cout << "error " << error.what() << '\n';
    }
}

void runCases()
{
    PinnedBudgetAlignMatrix matrix(BaseMatchingPenaltyMatrix(5, -4), -8, -2);
    const auto chain = [](int length) { return EdgeMap({ { length, length } }, { 0 }); };
    const std::string exact = "ACGTGCTAGCTAGGCTAAGCTGATCGATCGGTACGATGCTAGCTAGCATCGATCGATGATCCGA";
    inspect("compact-exact", matrix, exact.substr(0, 32), exact, chain(exact.size()));
    // Empty init throws before replacing the previous successful matrix state.
    inspect("empty-query", matrix, "", exact, chain(exact.size()));
    Score score = -1973;
    auto best = matrix.nextBestAlign<false>(matrix.alignBegin(), score);
    std::cout << "after-empty-query " << best << ' ' << score << '\n';
    inspect("empty-target", matrix, exact, "", chain(exact.size()));
    score = -1973;
    best = matrix.nextBestAlign<false>(matrix.alignBegin(), score);
    std::cout << "after-empty-target " << best << ' ' << score << '\n';
    auto twoErrors = exact.substr(0, 32);
    twoErrors[0] = twoErrors[31] = 'N';
    inspect("two-error-fallback", matrix, twoErrors, exact, chain(exact.size()));
    inspect("compact-after-fallback", matrix, exact.substr(0, 33), exact, chain(exact.size()));
    inspect("short-dense-gate", matrix, exact.substr(0, 15), exact, chain(exact.size()));
    inspect("edge-map-too-short", matrix, exact.substr(0, 32), exact, chain(16));
    inspect("edge-map-longer-than-target", matrix, exact.substr(0, 16), exact.substr(0, 32), chain(exact.size()));
    // Same valid chain through the raw constructor must retain the generic gate.
    std::vector<int> ids(64, 0), predecessors;
    std::vector<size_t> index;
    for (int i = 0; i < 64; ++i) { predecessors.push_back(i - 1); index.push_back(i); }
    index.push_back(64);
    inspect("raw-chain", matrix, exact.substr(0, 32), exact, EdgeMap(ids, index, predecessors));
    inspect("short-target", matrix, exact.substr(0, 32), exact.substr(0, 16), chain(16));
    const std::string longSequence(1500, 'A');
    inspect("numeric-fallback", matrix, longSequence, longSequence, chain(1500));
    PinnedBudgetAlignMatrix changedScores(BaseMatchingPenaltyMatrix(4, -4), -8, -2);
    inspect("nondefault-score-fallback", changedScores, exact.substr(0, 32), exact, chain(exact.size()));
    inspect("compact-after-everything", matrix, exact.substr(0, 32), exact, chain(exact.size()));
}

}
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "graphalign/dagAligner/PinnedBudgetAlignMatrix.hh"
#include "graphalign/DagAlignerAffine.hh"

namespace endpoints
{
// Endpoint and traceback cases use 150-base queries.

using namespace graphalign::dagAligner;

static EdgeMap branches(int count, int length)
{
    std::vector<int> ids, predecessors;
    std::vector<std::size_t> index;
    for (int t = 0; t < count * length; ++t)
    {
        ids.push_back(t / length);
        index.push_back(predecessors.size());
        predecessors.push_back(t % length == 0 ? -1 : t - 1);
    }
    index.push_back(predecessors.size());
    return EdgeMap(ids, index, predecessors);
}

static void dump(const char* name, const PinnedBudgetAlignMatrix& matrix, const EdgeMap& edges)
{
    std::cout << name << ' ' << matrix.queryLen() << ' ' << matrix.alignEnd() << '\n';
    Score score = -1973;
    const auto endpoint = matrix.nextBestAlign<false>(matrix.alignBegin(), score);
    std::cout << "best " << endpoint << ' ' << score << '\n';
    (void)edges;
}

static void inspect(const char* name, PinnedBudgetAlignMatrix& matrix,
                    const std::string& query, const std::string& target, const EdgeMap& edges)
{
    matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges);
    dump(name, matrix, edges);
    Aligner<PinnedBudgetAlignMatrix, false> aligner(BaseMatchingPenaltyMatrix(5, -4), -8, -2);
    aligner.align(query.begin(), query.end(), target.begin(), target.end(), edges);
    std::vector<Cigar> cigars;
    Score second = -1973;
    const Score best = aligner.backtrackAllPaths<false>(edges, cigars, second);
    using Dense = AffineAlignMatrixVectorized<BaseMatchingPenaltyMatrix, true>;
    Aligner<Dense, false> reference(BaseMatchingPenaltyMatrix(5, -4), -8, -2);
    reference.align(query.begin(), query.end(), target.begin(), target.end(), edges);
    std::vector<Cigar> expected;
    Score ignored = -1973;
    const Score expectedScore = reference.backtrackAllPaths<false>(edges, expected, ignored);
    std::ostringstream actualText, expectedText;
    for (const auto& cigar : cigars) actualText << cigar << '\n';
    for (const auto& cigar : expected) expectedText << cigar << '\n';
    if (best != expectedScore || actualText.str() != expectedText.str())
        throw std::runtime_error("Budget matrix differs from dense traceback");
    std::cout << "ordered-paths " << best << ' ' << cigars.size() << '\n';
    for (const auto& cigar : cigars)
        std::cout << cigar << '\n';
}

void runCases()
{
    const std::string query(150, 'A');
    const auto chain = [](int length) { return EdgeMap({ { length, length } }, { 0 }); };
    PinnedBudgetAlignMatrix matrix(BaseMatchingPenaltyMatrix(5, -4), -8, -2);
    std::string target(600, 'A');
    const auto parallel = branches(4, 150);
    inspect("four-tied-sparse-endpoints", matrix, query, target, parallel);
    // Costs 9, 0, 9, 18: the last row is not finite at the successful budget 9.
    target[74] = 'C'; target[300 + 74] = 'C';
    target[450 + 74] = target[450 + 100] = 'C';
    inspect("ordered-ties-and-suboptimal-endpoints", matrix, query, target, parallel);
    try
    {
        matrix.init(query.begin(), query.begin(), target.begin(), target.end(), parallel);
    }
    catch (const std::exception& e) { std::cout << "empty-query-error " << e.what() << '\n'; }
    dump("sparse-state-after-empty-query", matrix, parallel);
    try
    {
        matrix.init(query.begin(), query.end(), target.begin(), target.begin(), parallel);
    }
    catch (const std::exception& e) { std::cout << "empty-target-error " << e.what() << '\n'; }
    dump("sparse-state-after-empty-target", matrix, parallel);

    target.assign(600, 'A');
    for (int branch = 0; branch < 4; ++branch)
        for (int position : { 20, 60, 100 })
            target[150 * branch + position] = 'C';
    inspect("budget-retry-three-errors", matrix, query, target, parallel);
    inspect("exhausted-to-dense", matrix, query, std::string(600, 'C'), parallel);
    inspect("sparse-after-dense", matrix, query, std::string(600, 'A'), parallel);
    inspect("compact-after-sparse", matrix, query, std::string(1200, 'A'), chain(1200));
    inspect("dense-short-target-after-compact", matrix, query, "AAAA", chain(4));
    // The raw chain uses generic prepareRows/fillUngappedLinear, not the
    // constructor-certified compact single-node route.
    inspect("generic-linear-after-dense", matrix, query, std::string(1200, 'A'), branches(1, 1200));
    inspect("sparse-after-linear", matrix, query, std::string(600, 'A'), parallel);
    PinnedBudgetAlignMatrix changed(BaseMatchingPenaltyMatrix(4, -4), -8, -2);
    changed.init(query.begin(), query.end(), target.begin(), target.end(), parallel);
    dump("non-default-dense", changed, parallel);
}

}

int main()
{
    try { lifecycle::runCases(); endpoints::runCases(); }
    catch (const std::exception& error) { std::cerr << error.what() << "\n"; return 1; }
}
