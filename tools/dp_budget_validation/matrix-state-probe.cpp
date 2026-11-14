// Direct matrix lifecycle and malformed-edge-map differential checks.
#include <iostream>
#include <string>
#include <vector>
#include "graphalign/dagAligner/PinnedBudgetAlignMatrix.hh"

using namespace graphalign::dagAligner;

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
        if (best != matrix.alignEnd())
        {
            Score nextScore = -1973;
            const auto next = matrix.nextBestAlign<false>(best + 1, nextScore);
            std::cout << "next " << next << ' ' << nextScore << '\n';
        }
    }
    catch (const std::exception& error)
    {
        std::cout << "error " << error.what() << '\n';
    }
}

int main()
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
