// Exact cost-bounded affine DP for the default, pinned production alignment.
// The generic full-matrix API remains in AffineAlignMatrixVectorized.hh.
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <ostream>
#include <vector>

#include "AffineAlignMatrixVectorized.hh"
#include "BaseMatchingPenaltyMatrix.hh"

namespace graphalign
{
namespace dagAligner
{
class PinnedBudgetAlignMatrix
{
public:
    typedef BaseMatchingPenaltyMatrix PenaltyMatrix;
    // Only full-query endpoints are exposed to the original ordered traceback.
    typedef int const_iterator;

#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
    struct Stats
    {
        uint64_t calls = 0, attempts = 0, succeeded = 0, fallbacks = 0;
        uint64_t denseCells = 0, bandCells = 0;
        uint64_t successDenseCells = 0, fallbackDenseCells = 0;
        uint64_t edgeRelaxations = 0, initializedCells = 0, insertionCells = 0;
        uint64_t metadataVisits = 0, exhaustedAttempts = 0;
    };
    static Stats& stats()
    {
        static thread_local Stats value;
        return value;
    }
#endif

    PinnedBudgetAlignMatrix(const PenaltyMatrix& penalty, Score gapOpen, Score gapExt)
        : penalty_(penalty)
        , dense_(penalty, gapOpen, gapExt)
        , eligibleScores_(gapOpen == -8 && gapExt == -2 && hasDefaultScores(penalty))
    {
    }

    template <typename QueryIt, typename TargetIt>
    void init(QueryIt queryBegin, QueryIt queryEnd, TargetIt targetBegin, TargetIt targetEnd, const EdgeMap& edges)
    {
        // Preserve the original exceptions and state on empty-input failure.
        if (queryBegin == queryEnd || targetBegin == targetEnd)
        {
            dense_.init(queryBegin, queryEnd, targetBegin, targetEnd, edges);
            return;
        }
        queryLength_ = std::distance(queryBegin, queryEnd);
        targetLength_ = std::distance(targetBegin, targetEnd);
        sparse_ = false;
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
        ++stats().calls;
        stats().denseCells += uint64_t(queryLength_) * targetLength_;
#endif
        // The bound protects score differences used by legacy traceback as well
        // as scores. Non-default scores and unusual graph inputs retain their
        // exact original int16 narrowing behavior through the dense fallback.
        bool eligible = eligibleScores_ && queryLength_ >= 16 && targetLength_ >= 8
            && int64_t(15) * queryLength_ + int64_t(10) * targetLength_ + 10 < 32768;
#ifdef GRAPHTOOLS_BUDGET_DP_DISABLE
        eligible = false;
#endif
        if (eligible && prepareRows(edges))
        {
            uint64_t spent = uint64_t(targetLength_) + edgeCount_;
            const uint64_t limit = uint64_t(queryLength_) * (uint64_t(targetLength_) + edgeCount_) / 2;
            int firstBudget = 9;
            while (firstBudget < minimumEndpointCost_)
                firstBudget *= 2;
            if (spent < limit)
            {
                query_.clear();
                target_.clear();
                penalty_.translateQuery(queryBegin, queryEnd, std::back_inserter(query_));
                penalty_.translateTarget(targetBegin, targetEnd, std::back_inserter(target_));
                for (int budget = firstBudget; spent < limit; budget *= 2)
                {
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
                    ++stats().attempts;
#endif
                    const FillResult result = fill(budget, edges, spent, limit);
                    if (result == FillResult::Found)
                    {
                        sparse_ = true;
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
                        ++stats().succeeded;
                        stats().successDenseCells += uint64_t(queryLength_) * targetLength_;
#endif
                        return;
                    }
                    if (result == FillResult::Exhausted)
                    {
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
                        ++stats().exhaustedAttempts;
#endif
                        break;
                    }
                }
            }
        }
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
        ++stats().fallbacks;
        stats().fallbackDenseCells += uint64_t(queryLength_) * targetLength_;
#endif
        dense_.init(queryBegin, queryEnd, targetBegin, targetEnd, edges);
    }

    const_iterator alignBegin() const { return 0; }
    const_iterator alignEnd() const { return targetLength_; }
    int queryLen() const { return queryLength_; }
    int targetOffset(const_iterator cell) const { return cell; }
    int queryOffset(const_iterator) const { return queryLength_ - 1; }

    template <bool localAlign> const_iterator nextBestAlign(const_iterator start, Score& bestScore) const
    {
        static_assert(!localAlign, "Budget matrix is restricted to pinned full-query endpoints");
        if (start >= targetLength_)
        {
            return alignEnd();
        }
        if (!sparse_)
        {
            const int stride = 1 + (queryLength_ + 15) / 16 * 16;
            const auto cell = dense_.template nextBestAlign<false>(dense_.alignBegin() + start * stride, bestScore);
            return cell == dense_.alignEnd() ? alignEnd() : dense_.targetOffset(cell);
        }
        Cost best = infinity();
        int result = alignEnd();
        for (int t = start; t < targetLength_; ++t)
        {
            const Cost cost = h(queryLength_, t);
            if (cost < best)
            {
                best = cost;
                result = t;
            }
        }
        if (result != alignEnd())
        {
            bestScore = Score(5 * queryLength_ - best);
        }
        return result;
    }

    bool isInsertion(int q, int t) const
    {
        if (!sparse_)
            return dense_.isInsertion(q, t);
        if (q < 0)
            return false;
        const Cell* cell = lookup(q + 1, t);
        return cell && cell->insertion;
    }
    bool isDeletion(int q, int t, int p) const
    {
        if (!sparse_)
            return dense_.isDeletion(q, t, p);
        const int i = q + 1;
        const Cost current = h(i, t);
        return current != infinity()
            && (equalWithIncrement(current, e(i, p), 2) || equalWithIncrement(current, h(i, p), 10));
    }
    bool isMatch(int q, int t, int p) const
    {
        if (!sparse_)
            return dense_.isMatch(q, t, p);
        return q >= 0 && penalty_.isMatch(query_[q], target_[t])
            && equalWithIncrement(h(q + 1, t), h(q, p), 0);
    }
    bool isMismatch(int q, int t, int p) const
    {
        if (!sparse_)
            return dense_.isMismatch(q, t, p);
        return q >= 0 && !penalty_.isMatch(query_[q], target_[t])
            && equalWithIncrement(h(q + 1, t), h(q, p), 9);
    }

private:
    typedef uint16_t Cost;
    static Cost infinity() { return std::numeric_limits<Cost>::max(); }
    struct Row
    {
        int shortest = 0, longest = 0, begin = 0, end = 0;
        int supportBegin = 0, supportEnd = 0;
        std::size_t offset = 0;
    };
    struct Cell
    {
        Cost h, e;
        bool insertion;
    };
    const PenaltyMatrix penalty_;
    AffineAlignMatrixVectorized<PenaltyMatrix, true> dense_;
    const bool eligibleScores_;
    bool sparse_ = false;
    int queryLength_ = 0, targetLength_ = 0, budget_ = 0;
    int minimumEndpointCost_ = 0;
    uint64_t edgeCount_ = 0;
    std::vector<Row> rows_;
    std::vector<Cell> cells_;
    std::vector<PenaltyMatrix::QueryChar> query_;
    std::vector<PenaltyMatrix::TargetChar> target_;

    static bool hasDefaultScores(const PenaltyMatrix& penalty)
    {
        for (int t = 0; t < PenaltyMatrix::ROWS_; ++t)
            for (int q = 0; q < PenaltyMatrix::COLUMNS_; ++q)
                if (penalty(q, t) != (penalty.isMatch(q, t) ? 5 : -4))
                    return false;
        return true;
    }
    bool prepareRows(const EdgeMap& edges)
    {
        rows_.resize(targetLength_);
        edgeCount_ = 0;
        minimumEndpointCost_ = std::numeric_limits<int>::max();
        for (int t = 0; t < targetLength_; ++t)
        {
            const auto first = edges.prevNodesBegin(t), last = edges.prevNodesEnd(t);
            if (first == last)
                return false;
            int shortest = targetLength_ + 1, longest = 0;
            for (auto it = first; it != last; ++it)
            {
                const int p = *it;
                if (p < -1 || p >= t)
                    return false;
                shortest = std::min(shortest, p == -1 ? 1 : rows_[p].shortest + 1);
                longest = std::max(longest, p == -1 ? 1 : rows_[p].longest + 1);
                ++edgeCount_;
            }
            rows_[t].shortest = shortest;
            rows_[t].longest = longest;
            // Every complete-query alignment ending here must pay for the
            // graph/query length difference. Opening penalties are omitted,
            // making this a conservative bound even for legacy boundaries.
            const int endpointBound = 2 * std::max(shortest - queryLength_, 0)
                + 7 * std::max(queryLength_ - longest, 0);
            minimumEndpointCost_ = std::min(minimumEndpointCost_, endpointBound);
        }
        return true;
    }
    static int overlapLength(int firstBegin, int firstEnd, int secondBegin, int secondEnd)
    {
        return std::max(0, std::min(firstEnd, secondEnd) - std::max(firstBegin, secondBegin));
    }
    enum class FillResult { Found, Missing, Exhausted };
    static bool charge(uint64_t amount, uint64_t& spent, uint64_t limit)
    {
        if (amount > limit - spent)
            return false;
        spent += amount;
        return true;
    }
    const Cell* lookup(int i, int t) const
    {
        if (t < 0)
            return nullptr;
        const Row& row = rows_[t];
        return i < row.begin || i >= row.end ? nullptr : &cells_[row.offset + i - row.begin];
    }
    Cost h(int i, int t) const
    {
        if (t == -1)
            return 15 * i <= budget_ ? Cost(15 * i) : infinity();
        const Cell* cell = lookup(i, t);
        return cell ? cell->h : infinity();
    }
    Cost e(int i, int t) const
    {
        if (t == -1)
            return 15 * i <= budget_ ? Cost(15 * i) : infinity();
        const Cell* cell = lookup(i, t);
        return cell ? cell->e : infinity();
    }
    static bool equalWithIncrement(Cost current, Cost previous, int increment)
    {
        return current != infinity() && previous != infinity() && int(current) == int(previous) + increment;
    }
    Cost add(Cost value, int increment) const
    {
        return value != infinity() && int(value) + increment <= budget_ ? Cost(value + increment) : infinity();
    }
    FillResult fill(int budget, const EdgeMap& edges, uint64_t& spent, uint64_t limit)
    {
        budget_ = budget;
        std::size_t usedCells = 0;
        const Cell unreachable = { infinity(), infinity(), false };
        bool endpoint = false;
        const int virtualEnd = std::min(queryLength_, budget / 15) + 1;
        for (int t = 0; t < targetLength_; ++t)
        {
            // Construct a row from the finite predecessor H supports. E is a
            // subset of H after each completed row, so no separate E support
            // is needed. Charge each metadata visit and every cell operation
            // before performing it, including partial attempts that fall back.
            if (!charge(1, spent, limit))
                return FillResult::Exhausted;
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
            ++stats().metadataVisits;
#endif
            Row& row = rows_[t];
            row.offset = usedCells;
            row.supportBegin = row.supportEnd = 0;
            int seedBegin = queryLength_ + 1, seedEnd = 0;
            if (10 * row.shortest <= budget)
            {
                seedBegin = 0;
                seedEnd = 1;
            }
            for (auto p = edges.prevNodesBegin(t); p != edges.prevNodesEnd(t); ++p)
            {
                if (!charge(1, spent, limit))
                    return FillResult::Exhausted;
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
                ++stats().metadataVisits;
#endif
                const int predecessorBegin = *p == -1 ? 0 : rows_[*p].supportBegin;
                const int predecessorEnd = *p == -1 ? virtualEnd : rows_[*p].supportEnd;
                if (predecessorBegin < predecessorEnd)
                {
                    // Deletion consumes the same i; match consumes i+1.
                    seedBegin = std::min(seedBegin, std::max(1, predecessorBegin));
                    seedEnd = std::max(seedEnd, predecessorEnd + 1);
                }
            }
            const int structuralBegin = std::max(0, row.shortest - budget / 2);
            const int structuralEnd = std::min(queryLength_, row.longest + budget / 7) + 1;
            row.begin = std::max(structuralBegin, seedBegin);
            // Every insertion adds at least 7, including the legacy i=0 F
            // boundary. This right expansion therefore contains every finite
            // insertion descendant of a deletion/match/boundary seed.
            row.end = seedEnd == 0 ? row.begin : std::min(structuralEnd, seedEnd + budget / 7);
            if (row.end <= row.begin)
                continue;
            const int width = row.end - row.begin;
            if (!charge(width, spent, limit))
                return FillResult::Exhausted;
            usedCells += width;
            if (cells_.size() < usedCells)
                cells_.resize(usedCells);
            // resize may move storage: obtain all cell pointers only after it.
            Cell* const current = cells_.data() + row.offset;
            std::fill(current, current + width, unreachable);
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
            stats().initializedCells += width;
            stats().bandCells += width;
#endif
            if (row.begin == 0)
                current[0].h = Cost(10 * row.shortest);
            const int begin = std::max(1, row.begin);
            for (auto p = edges.prevNodesBegin(t); p != edges.prevNodesEnd(t); ++p)
            {
                const int predecessorBegin = *p == -1 ? 0 : rows_[*p].supportBegin;
                const int predecessorEnd = *p == -1 ? virtualEnd : rows_[*p].supportEnd;
                const int deletionBegin = std::max(begin, predecessorBegin);
                const int deletionEnd = std::min(row.end, predecessorEnd);
                const int matchBegin = std::max(begin, predecessorBegin + 1);
                const int matchEnd = std::min(row.end, predecessorEnd + 1);
                const int relaxations = std::max(0, deletionEnd - deletionBegin)
                    + std::max(0, matchEnd - matchBegin);
                if (!charge(1 + uint64_t(relaxations), spent, limit))
                    return FillResult::Exhausted;
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
                ++stats().metadataVisits;
                stats().edgeRelaxations += relaxations;
#endif
                if (*p == -1)
                {
                    // The original virtual row initializes H and E alike.
                    // Its deletion extension therefore dominates opening.
                    for (int i = deletionBegin; i < deletionEnd; ++i)
                    {
                        Cell& cell = current[i - row.begin];
                        cell.e = std::min(cell.e, add(Cost(15 * i), 2));
                    }
                    for (int i = matchBegin; i < matchEnd; ++i)
                    {
                        Cell& cell = current[i - row.begin];
                        const int mismatch = penalty_.isMatch(query_[i - 1], target_[t]) ? 0 : 9;
                        cell.h = std::min(cell.h, add(Cost(15 * (i - 1)), mismatch));
                    }
                }
                else if (predecessorBegin < predecessorEnd)
                {
                    const Row& predecessor = rows_[*p];
                    const Cell* const previous = cells_.data() + predecessor.offset;
                    for (int i = deletionBegin; i < deletionEnd; ++i)
                    {
                        Cell& cell = current[i - row.begin];
                        const Cell& source = previous[i - predecessor.begin];
                        cell.e = std::min(cell.e, std::min(add(source.e, 2), add(source.h, 10)));
                    }
                    for (int i = matchBegin; i < matchEnd; ++i)
                    {
                        Cell& cell = current[i - row.begin];
                        const int mismatch = penalty_.isMatch(query_[i - 1], target_[t]) ? 0 : 9;
                        cell.h = std::min(cell.h, add(previous[i - 1 - predecessor.begin].h, mismatch));
                    }
                }
            }
            const int insertionCells = row.end - begin;
            if (!charge(insertionCells, spent, limit))
                return FillResult::Exhausted;
            Cost previousInsertion = row.begin == 0 ? current[0].h : infinity();
            Cost previousBest = previousInsertion;
            if (previousBest != infinity())
            {
                row.supportBegin = 0;
                row.supportEnd = 1;
            }
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
            stats().insertionCells += insertionCells;
#endif
            for (int i = begin; i < row.end; ++i)
            {
                Cell& cell = current[i - row.begin];
                const Cost insertion = std::min(add(previousInsertion, 7), add(previousBest, 15));
                previousInsertion = insertion;
                cell.h = std::min(cell.h, std::min(cell.e, insertion));
                previousBest = cell.h;
                cell.insertion = cell.h != infinity() && cell.h == insertion;
                if (cell.h != infinity())
                {
                    if (row.supportBegin == row.supportEnd)
                        row.supportBegin = i;
                    row.supportEnd = i + 1;
                    endpoint = endpoint || i == queryLength_;
                }
            }
        }
        return endpoint ? FillResult::Found : FillResult::Missing;
    }

    friend std::ostream& operator<<(std::ostream& out, const PinnedBudgetAlignMatrix& matrix)
    {
        return out << "PinnedBudgetAlignMatrix(query=" << matrix.queryLength_ << ", target="
                   << matrix.targetLength_ << ", sparse=" << matrix.sparse_ << ')';
    }
};
} // namespace dagAligner
} // namespace graphalign
