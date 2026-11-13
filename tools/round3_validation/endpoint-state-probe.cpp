// Compile unchanged against the pre-experiment and prototype headers, then
// compare complete stdout. All query lengths remain 150 bp.
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include "graphalign/dagAligner/PinnedBudgetAlignMatrix.hh"
#include "graphalign/DagAlignerAffine.hh"

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
    // Include the unchanged sentinel score when no endpoint exists, every
    // suffix start and starts beyond the matrix.
    for (int start = 0; start <= matrix.alignEnd() + 1; ++start)
    {
        Score score = -1973;
        const auto endpoint = matrix.nextBestAlign<false>(start, score);
        std::cout << "next " << start << ' ' << endpoint << ' ' << score << '\n';
    }
    uint64_t digest = UINT64_C(1469598103934665603);
    const auto add = [&digest](bool value) {
        digest ^= static_cast<uint64_t>(value);
        digest *= UINT64_C(1099511628211);
    };
    for (int t = 0; t < matrix.alignEnd(); ++t)
        for (int q = 0; q < matrix.queryLen(); ++q)
        {
            add(matrix.isInsertion(q, t));
            for (auto p = edges.prevNodesBegin(t); p != edges.prevNodesEnd(t); ++p)
            {
                add(matrix.isDeletion(q, t, *p));
                add(matrix.isMatch(q, t, *p));
                add(matrix.isMismatch(q, t, *p));
            }
        }
    std::cout << "traceback-flags " << digest << '\n';
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
    std::cout << "ordered-paths " << best << ' ' << second << ' ' << cigars.size() << '\n';
    for (const auto& cigar : cigars)
        std::cout << cigar << '\n';
}

int main()
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
