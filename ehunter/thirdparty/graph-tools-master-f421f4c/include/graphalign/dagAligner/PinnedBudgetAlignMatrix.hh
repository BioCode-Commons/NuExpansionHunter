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
// Exact cost-bounded affine DP for the default, pinned production alignment.
// The generic full-matrix API remains in AffineAlignMatrixVectorized.hh.
#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <ostream>
#include <vector>

#include "AffineAlignMatrixVectorized.hh"
#include "BaseMatchingPenaltyMatrix.hh"
#include "graphutils/VerifyPolicy.hh"

namespace graphalign
{
namespace dagAligner
{
// Lets a caller skip the budgeted sparse search of one call and learn how the
// search went. Skipping changes only which exact matrix computes the result:
// for every eligible call the sparse search either returns the dense matrix's
// result or falls back to the dense matrix.
struct BudgetSearch
{
    enum Outcome
    {
        kNotAttempted, // ineligible, empty or compact-linear calls
        kFound,
        kFailed, // the search ran and fell back to the dense matrix
        kSkipped
    };
    bool skip = false;
    Outcome outcome = kNotAttempted;
    // PinnedDagAligner consumes only the best score and optimal tracebacks.
    bool bestEndpointsOnly = false;
};

class PinnedBudgetAlignMatrix
{
public:
    typedef BaseMatchingPenaltyMatrix PenaltyMatrix;
    // Only full-query endpoints are exposed to v5.0.0's ordered traceback.
    typedef int const_iterator;

    PinnedBudgetAlignMatrix(const PenaltyMatrix& penalty, Score gapOpen, Score gapExt)
        : penalty_(penalty)
        , dense_(penalty, gapOpen, gapExt)
        , eligibleScores_(gapOpen == -8 && gapExt == -2 && hasDefaultScores(penalty))
    {
    }

    template <typename QueryIt, typename TargetIt>
    void init(QueryIt queryBegin, QueryIt queryEnd, TargetIt targetBegin, TargetIt targetEnd,
              const EdgeMap& edges, int minimumCost = 0, int preferredInitialBudget = 0,
              BudgetSearch* search = nullptr)
    {
        if (search)
            search->outcome = BudgetSearch::kNotAttempted;
        // Preserve v5.0.0's exceptions and state on empty-input failure.
        if (queryBegin == queryEnd || targetBegin == targetEnd)
        {
            dense_.init(queryBegin, queryEnd, targetBegin, targetEnd, edges);
            return;
        }
        // Empty-input failure above must preserve the previous matrix state.
        finiteEndpoints_.clear();
        queryLength_ = std::distance(queryBegin, queryEnd);
        targetLength_ = std::distance(targetBegin, targetEnd);
        sparse_ = false;
        compactLinear_ = false;
        linearMismatch_ = -2;
        // The bound protects traceback score differences as well as scores.
        // The dense fallback preserves v5.0.0's exact int16 narrowing for
        // non-default scores and unusual graph inputs.
        const bool eligible = eligibleScores_ && queryLength_ >= 16 && targetLength_ >= 8
            && int64_t(15) * queryLength_ + int64_t(10) * targetLength_ + 10 < 32768;
        const int fullTargetLength = targetLength_;
        const bool linear = eligible && edges.hasLinearChain(fullTargetLength);
        const bool compact = graphtools::verify::kFastPaths && linear && minimumCost <= 9
            && queryLength_ <= fullTargetLength;
        const bool boundPrefix = graphtools::verify::kRowStop && !graphtools::verify::kShadow && linear
            && search && search->bestEndpointsOnly && fullTargetLength - queryLength_ > 16;
        // A proven linear chain has only one ungapped full-query endpoint.
        // With at most one mismatch its cost is <=9, below every gap opening.
        // Keep just the mismatch position instead of allocating per-target rows
        // and cells, or translating an unused extension past the query length.
        int translated = 0;
        if (compact || boundPrefix)
        {
            query_.resize(queryLength_);
            target_.resize(queryLength_);
            penalty_.translateQuery(queryBegin, queryEnd, query_.begin());
            penalty_.translateTarget(targetBegin, std::next(targetBegin, queryLength_), target_.begin());
            translated = queryLength_;
            if (compact)
            {
                linearMismatch_ = linearFirstMismatch();
                if (linearMismatch_ >= 0)
                {
                    compactLinear_ = true;
                    return;
                }
            }
            if (boundPrefix)
                targetLength_ = linearTargetPrefixLength(fullTargetLength);
        }
        const bool clipped = targetLength_ != fullTargetLength;
        if (clipped)
            targetEnd = std::next(targetBegin, targetLength_);
        const bool skipped = eligible && search && search->skip && !graphtools::verify::kBudgetSearchAlways;
        if (skipped)
            search->outcome = BudgetSearch::kSkipped;
        if (eligible && !skipped && prepareRows(edges, linear))
        {
            if (search)
                search->outcome = BudgetSearch::kFailed; // until a budget succeeds
            uint64_t spent = uint64_t(targetLength_) + edgeCount_;
            int firstBudget = 9;
            while (firstBudget < std::max(minimumEndpointCost_, minimumCost))
                firstBudget *= 2;
            if (minimumCost > 9 && preferredInitialBudget >= 36 && firstBudget < 36
                && !graphtools::verify::kBudgetSearchAlways)
            {
                // Search a wider complete cost range; this is not a proof
                // that smaller costs are absent. All optimal ties remain in
                // the range, and the work cap and dense fallback apply.
                firstBudget = 36;
            }
            if (attemptLimit(firstBudget, 0, spent) > spent)
            {
                if (!translated)
                {
                    query_.resize(queryLength_);
                    penalty_.translateQuery(queryBegin, queryEnd, query_.begin());
                }
                target_.resize(targetLength_);
                penalty_.translateTarget(std::next(targetBegin, translated), targetEnd, target_.begin() + translated);
                const auto searchStart = std::chrono::steady_clock::now();
                double elapsed = 0; // nanoseconds spent in this search's attempts
                for (int budget = firstBudget;; budget *= 2)
                {
                    const uint64_t limit = attemptLimit(budget, elapsed, spent);
                    if (limit <= spent)
                        break;
                    const uint64_t before = spent;
                    const FillResult result = fill(budget, edges, spent, limit);
                    const double now = std::chrono::duration<double, std::nano>(
                        std::chrono::steady_clock::now() - searchStart).count();
                    recordAttempt(budget, now - elapsed, spent - before);
                    elapsed = now;
                    if (result == FillResult::Found)
                    {
                        sparse_ = true;
                        if (search)
                            search->outcome = BudgetSearch::kFound;
                        return;
                    }
                    if (result == FillResult::Exhausted)
                        break;
                }
            }
        }
        // The row stop needs the default scores (endpoint cost = 5 * query - score)
        // and true scores in the rows it drops: the 16-bit guard for the whole target.
        const bool rowStop = eligibleScores_ && int64_t(15) * queryLength_ + int64_t(10) * fullTargetLength + 10 < 32768;
        dense_.setRowStop(
            rowStop && graphtools::verify::kRowStop ? DenseMatrix::RowStop::kStop : DenseMatrix::RowStop::kOff);
        if (eligible)
        {
            // Only fills that a search could have replaced measure the dense
            // alternative; shorter queries and other scores never compete.
            const auto denseStart = std::chrono::steady_clock::now();
            dense_.init(queryBegin, queryEnd, targetBegin, targetEnd, edges, clipped ? fullTargetLength : 0);
            denseNanos_
                += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - denseStart).count();
            denseCells_ += double(queryLength_) * dense_.keptRows();
        }
        else
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
        if (compactLinear_)
        {
            const int endpoint = queryLength_ - 1;
            if (start > endpoint)
                return alignEnd();
            bestScore = Score(5 * queryLength_ - (linearMismatch_ < queryLength_ ? 9 : 0));
            return endpoint;
        }
        if (!sparse_)
        {
            if (start >= dense_.scanRows())
                return alignEnd();
            const int stride = 1 + (queryLength_ + 15) / 16 * 16;
            const auto cell = dense_.template nextBestAlign<false>(dense_.alignBegin() + start * stride, bestScore);
            return cell == dense_.alignEnd() ? alignEnd() : dense_.targetOffset(cell);
        }
        Cost best = infinity();
        int result = alignEnd();
        // Keep the virtual-row behavior for callers supplying a position
        // before alignBegin(), while normal traceback starts at zero.
        for (int t = start; t < 0; ++t)
        {
            const Cost cost = h(queryLength_, t);
            if (cost < best)
            {
                best = cost;
                result = t;
            }
        }
        // The successful sparse fill records every finite full-query endpoint
        // in increasing target order. Retain all of them, including suboptimal
        // costs, so nextBestAlign(start) and second-best scores stay unchanged.
        for (auto endpoint = std::lower_bound(finiteEndpoints_.begin(), finiteEndpoints_.end(), start);
             endpoint != finiteEndpoints_.end(); ++endpoint)
        {
            const int t = *endpoint;
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
        if (compactLinear_)
            return false;
        if (!sparse_)
            return dense_.isInsertion(q, t);
        if (q < 0)
            return false;
        const Cell* cell = lookup(q + 1, t);
        return cell && cell->insertion;
    }
    bool isDeletion(int q, int t, int p) const
    {
        if (compactLinear_)
            return false;
        if (!sparse_)
            return dense_.isDeletion(q, t, p);
        const int i = q + 1;
        const Cost current = h(i, t);
        return current != infinity()
            && (equalWithIncrement(current, e(i, p), 2) || equalWithIncrement(current, h(i, p), 10));
    }
    bool isMatch(int q, int t, int p) const
    {
        if (compactLinear_)
            return q >= 0 && q < queryLength_ && q == t && p == t - 1 && q != linearMismatch_;
        if (!sparse_)
            return dense_.isMatch(q, t, p);
        return q >= 0 && penalty_.isMatch(query_[q], target_[t])
            && equalWithIncrement(h(q + 1, t), h(q, p), 0);
    }
    bool isMismatch(int q, int t, int p) const
    {
        if (compactLinear_)
            return q >= 0 && q < queryLength_ && q == t && p == t - 1 && q == linearMismatch_;
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
    typedef AffineAlignMatrixVectorized<PenaltyMatrix, true> DenseMatrix;
    DenseMatrix dense_;
    const bool eligibleScores_;
    bool sparse_ = false;
    // Measured on this matrix: wall time and work units of budgeted search
    // attempts, by budget class (larger budgets touch more memory per unit),
    // and wall time and cells (query x rows filled) of the dense fills that
    // eligible calls ended in.
    static constexpr int kBudgetClasses = 6;
    double attemptNanos_[kBudgetClasses] = {}, attemptUnits_[kBudgetClasses] = {};
    double denseNanos_ = 0, denseCells_ = 0;

    // Budget class of a search attempt: 9, 18, 36, 72, 144, and 288 or more.
    static size_t budgetClass(int budget)
    {
        size_t index = 0;
        for (int value = 9; value < budget && index < kBudgetClasses - 1; value *= 2)
            ++index;
        return index;
    }

    void recordAttempt(int budget, double nanos, uint64_t units)
    {
        const size_t index = budgetClass(budget);
        attemptNanos_[index] += nanos;
        attemptUnits_[index] += double(units);
    }

    // The total work (in units, including `spent`) that attempt `budget` may
    // reach: the attempt stops paying once the search's time would exceed the
    // dense alternative's. That fill needs about query + budget / 2 rows (the
    // row stop ends there once an endpoint of cost <= budget is found); it is
    // costed at this matrix's measured time per dense cell, and the remaining
    // time converted to units at the measured time per unit of the attempt's
    // budget class (of all classes until that one is measured). Until both
    // methods have been measured, a cap of half the full dense
    // matrix applies.
    // Contract [K-10]
    // The cap reads the clock, so which exact matrix computes a result differs between runs; it has no functional
    // contract. K-07 and K-09 hold for every skip flag, first budget, limit and proven minimum cost, so every policy
    // outcome gives the reference result. Left to show here: termination and absence of undefined behaviour.
    uint64_t attemptLimit(int budget, double elapsedNanos, uint64_t spent) const
    {
        if (graphtools::verify::kBudgetSearchNever)
            return 0; // every eligible call goes to the dense matrix
        if (graphtools::verify::kBudgetSearchAlways)
            return std::numeric_limits<uint64_t>::max() / 2; // searches are never cut short
        const uint64_t fullCap = uint64_t(double(queryLength_) * (double(targetLength_) + edgeCount_) / 2);
        const size_t index = budgetClass(budget);
        double nanos = attemptNanos_[index], units = attemptUnits_[index];
        if (units <= 0 || nanos <= 0)
        {
            nanos = units = 0;
            for (int i = 0; i < kBudgetClasses; ++i)
                nanos += attemptNanos_[i], units += attemptUnits_[i];
        }
        if (units <= 0 || nanos <= 0 || denseCells_ <= 0)
            return fullCap;
        const double rows = std::min(double(targetLength_), double(queryLength_) + budget / 2.0);
        const double remaining = (denseNanos_ / denseCells_) * queryLength_ * rows - elapsedNanos;
        if (remaining <= 0)
            return spent;
        const double allowed = remaining / (nanos / units);
        const double cap = double(std::numeric_limits<uint64_t>::max() / 4);
        return spent + uint64_t(std::min(allowed, cap));
    }
    bool linearGraph_ = false;
    bool compactLinear_ = false;
    int linearMismatch_ = -2;
    int queryLength_ = 0, targetLength_ = 0, budget_ = 0, filledRows_ = 0;
    int minimumEndpointCost_ = 0;
    uint64_t edgeCount_ = 0;
    std::vector<Row> rows_;
    std::vector<Cell> cells_;
    std::vector<int> finiteEndpoints_;
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
    bool prepareRows(const EdgeMap& edges, bool linear)
    {
        rows_.resize(targetLength_);
        linearGraph_ = true;
        if (linear)
        {
            edgeCount_ = targetLength_;
            minimumEndpointCost_ = 7 * std::max(queryLength_ - targetLength_, 0);
            return true;
        }
        edgeCount_ = 0;
        minimumEndpointCost_ = std::numeric_limits<int>::max();
        for (int t = 0; t < targetLength_; ++t)
        {
            const auto first = edges.prevNodesBegin(t), last = edges.prevNodesEnd(t);
            if (first == last)
                return false;
            linearGraph_ = linearGraph_ && first + 1 == last && *first == t - 1;
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
            // making this a conservative bound under both boundary rules.
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
        if (t < 0 || t >= filledRows_)
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

    int linearTargetPrefixLength(int fullTargetLength) const
    {
        // The ungapped diagonal is feasible at cost 9 per mismatch. On this
        // chain, consuming L target bases and Q query bases needs at least
        // L-Q deletions, each costing at least 2 even at the virtual row.
        // Keep equality to preserve every optimal endpoint and its traceback.
        const int extra = fullTargetLength - queryLength_;
        int cost = 0;
        for (int i = 0; i < queryLength_; ++i)
        {
            if (!penalty_.isMatch(query_[i], target_[i]))
            {
                cost += 9;
                if (cost / 2 >= extra)
                    return fullTargetLength;
            }
        }
        return queryLength_ + cost / 2;
    }

    // At budget 9 every finite path is ungapped: opening a deletion costs 10,
    // opening an insertion costs 15, and the virtual row costs 15 per query
    // base. Its cheaper E extension cannot help: E(0,-1) is not propagated to
    // row i=0, while E(i,-1) for i>0 already exceeds this budget.
    // On a linear graph this leaves exactly the diagonal. Compare its bases in
    // native vectors and retain the first mismatch position for exact traceback.
    int linearFirstMismatch() const
    {
        int firstMismatch = queryLength_;
        int i = 0;
#if defined(__aarch64__) && defined(__ARM_NEON)
        if (!graphtools::verify::kScalarOnly)
            return dp_linear_first_mismatch_neon(query_.data(), target_.data(), queryLength_); // kernel K-09
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
        if (linearAvx2Supported())
            return dp_linear_first_mismatch_avx2(query_.data(), target_.data(), queryLength_); // kernel K-09
#endif
        for (; i < queryLength_; ++i)
        {
            if (!penalty_.isMatch(query_[i], target_[i]))
            {
                if (firstMismatch != queryLength_)
                    return -1;
                firstMismatch = i;
            }
        }
        return firstMismatch;
    }

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    static bool linearAvx2Supported()
    {
        if (graphtools::verify::kNoAvx2)
            return false;
        static const bool supported = __builtin_cpu_supports("avx2");
        return supported;
    }
#endif

    FillResult fillUngappedLinear(uint64_t& spent, uint64_t limit)
    {
        if (!charge(uint64_t(targetLength_) + queryLength_, spent, limit))
            return FillResult::Exhausted;
        if (queryLength_ > targetLength_)
            return FillResult::Missing;
        const int firstMismatch = linearMismatch_ == -2 ? linearFirstMismatch() : linearMismatch_;
        if (firstMismatch < 0)
            return FillResult::Missing;
        if (cells_.size() < static_cast<size_t>(queryLength_))
            cells_.resize(queryLength_);
        for (int t = 0; t < targetLength_; ++t)
        {
            Row& row = rows_[t];
            row.offset = std::min(t, queryLength_);
            row.begin = row.supportBegin = t + 1;
            row.end = row.supportEnd = t < queryLength_ ? t + 2 : t + 1;
            if (t < queryLength_)
                cells_[t] = { Cost(t < firstMismatch ? 0 : 9), infinity(), false };
        }
        finiteEndpoints_.push_back(queryLength_ - 1);
        return FillResult::Found;
    }

    FillResult fillUngapped(const EdgeMap& edges, uint64_t& spent, uint64_t limit)
    {
        if (linearGraph_)
            return fillUngappedLinear(spent, limit);
        std::size_t usedCells = 0;
        bool endpoint = false;
        for (int t = 0; t < targetLength_; ++t)
        {
            if (!charge(1, spent, limit))
                return FillResult::Exhausted;
            Row& row = rows_[t];
            row.offset = usedCells;
            row.supportBegin = row.supportEnd = 0;
            int begin = queryLength_ + 1, end = 0;
            for (auto p = edges.prevNodesBegin(t); p != edges.prevNodesEnd(t); ++p)
            {
                if (!charge(1, spent, limit))
                    return FillResult::Exhausted;
                const int predecessorBegin = *p == -1 ? 0 : rows_[*p].supportBegin;
                const int predecessorEnd = *p == -1 ? 1 : rows_[*p].supportEnd;
                if (predecessorBegin < predecessorEnd)
                {
                    begin = std::min(begin, predecessorBegin + 1);
                    end = std::max(end, predecessorEnd + 1);
                }
            }
            row.begin = begin;
            row.end = std::min(queryLength_ + 1, end);
            if (row.end <= row.begin)
                continue;
            const int width = row.end - row.begin;
            if (!charge(width, spent, limit))
                return FillResult::Exhausted;
            usedCells += width;
            if (cells_.size() < usedCells)
                cells_.resize(usedCells);
            Cell* const current = cells_.data() + row.offset;
            const auto single = edges.prevNodesBegin(t);
            if (graphtools::verify::kFusedRows && single + 1 == edges.prevNodesEnd(t))
            {
                // One predecessor supplies every cell of this ungapped row.
                // Fuse initialization, relaxation and support discovery while
                // charging 1 + 2 * width work units.
                if (!charge(1 + uint64_t(width) * 2, spent, limit))
                    return FillResult::Exhausted;
                const int p = *single;
                const Cell* const previous = p < 0 ? nullptr : cells_.data() + rows_[p].offset;
                for (int i = begin; i < row.end; ++i)
                {
                    const Cost previousBest = p < 0 ? 0 : previous[i - 1 - rows_[p].begin].h;
                    const int mismatch = penalty_.isMatch(query_[i - 1], target_[t]) ? 0 : 9;
                    const Cost best = add(previousBest, mismatch);
                    current[i - begin] = { best, infinity(), false };
                    if (best != infinity())
                    {
                        if (row.supportBegin == row.supportEnd)
                            row.supportBegin = i;
                        row.supportEnd = i + 1;
                        if (i == queryLength_)
                        {
                            endpoint = true;
                            finiteEndpoints_.push_back(t);
                        }
                    }
                }
                continue;
            }
            std::fill(current, current + width, Cell { infinity(), infinity(), false });
            for (auto p = edges.prevNodesBegin(t); p != edges.prevNodesEnd(t); ++p)
            {
                const int predecessorBegin = *p == -1 ? 0 : rows_[*p].supportBegin;
                const int predecessorEnd = *p == -1 ? 1 : rows_[*p].supportEnd;
                const int matchBegin = std::max(begin, predecessorBegin + 1);
                const int matchEnd = std::min(row.end, predecessorEnd + 1);
                if (!charge(1 + uint64_t(std::max(0, matchEnd - matchBegin)), spent, limit))
                    return FillResult::Exhausted;
                for (int i = matchBegin; i < matchEnd; ++i)
                {
                    const Cost previous = *p == -1 ? 0
                        : cells_[rows_[*p].offset + i - 1 - rows_[*p].begin].h;
                    const int mismatch = penalty_.isMatch(query_[i - 1], target_[t]) ? 0 : 9;
                    current[i - begin].h = std::min(current[i - begin].h, add(previous, mismatch));
                }
            }
            if (!charge(width, spent, limit))
                return FillResult::Exhausted;
            for (int i = begin; i < row.end; ++i)
                if (current[i - begin].h != infinity())
                {
                    if (row.supportBegin == row.supportEnd)
                        row.supportBegin = i;
                    row.supportEnd = i + 1;
                    if (i == queryLength_)
                    {
                        endpoint = true;
                        finiteEndpoints_.push_back(t);
                    }
                }
        }
        return endpoint ? FillResult::Found : FillResult::Missing;
    }
    // Contract [K-07]
    // Requires: default scores; 15q + 10t + 10 < 32768; prepareRows returned true.
    // Invariant: after each row, every cell whose true cost is <= budget lies in the stored band with its exact H and E
    //           costs, and every cell outside has cost > budget.
    // Ensures:  Found => best score, optimal endpoints in order and the four predicates on every cell an optimal path
    //           visits equal the dense matrix's. Missing and Exhausted have no observable effect (the caller falls back).
    // Reference: the dense matrix AffineAlignMatrixVectorized<..., true> (EH_VERIFY_BUDGET_SEARCH=NEVER).
    FillResult fill(int budget, const EdgeMap& edges, uint64_t& spent, uint64_t limit)
    {
        budget_ = budget;
        filledRows_ = targetLength_;
        // A failed or exhausted attempt must not leak endpoints into the next
        // budget. Dense and compact modes do not consult this cache.
        finiteEndpoints_.clear();
        if (budget == 9)
            return fillUngapped(edges, spent, limit);
        return linearGraph_ ? fillGapped<true>(budget, edges, spent, limit)
                            : fillGapped<false>(budget, edges, spent, limit);
    }

    template<bool linear>
    // Contract [K-08]
    // Ensures:  for every predecessor row and budget, the fused single-predecessor row writes the same cells, support
    //           range and endpoint flag, and charges the same work, as the general row code below it.
    // Reference: the general row code (EH_VERIFY_FUSED_ROWS=OFF).
    FillResult fillGapped(int budget, const EdgeMap& edges, uint64_t& spent, uint64_t limit)
    {
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
            Row& row = rows_[t];
            row.offset = usedCells;
            row.supportBegin = row.supportEnd = 0;
            int seedBegin = queryLength_ + 1, seedEnd = 0;
            const int shortest = linear ? t + 1 : row.shortest;
            const int longest = linear ? t + 1 : row.longest;
            if (linear && t && 10 * shortest > budget
                && rows_[t - 1].supportEnd == rows_[t - 1].supportBegin)
            {
                // All later rows have one unreachable predecessor and their
                // pinned boundary is over budget. Charge the same two metadata
                // visits per skipped row, including this row's prior charge.
                filledRows_ = t;
                if (!charge(2 * uint64_t(targetLength_ - t) - 1, spent, limit))
                    return FillResult::Exhausted;
                break;
            }
            if (10 * shortest <= budget)
            {
                seedBegin = 0;
                seedEnd = 1;
            }
            // A validated chain has one implicit predecessor. The generic
            // specialization still visits every incoming edge of a join.
            const int preceding = t - 1;
            const int* first = linear ? &preceding : &*edges.prevNodesBegin(t);
            const auto predecessors = linear ? 1 : edges.prevNodesEnd(t) - edges.prevNodesBegin(t);
            const int* last = first + predecessors;
            if (!charge(predecessors, spent, limit))
                return FillResult::Exhausted;
            for (auto p = first; p != last; ++p)
            {
                const int predecessorBegin = *p < 0 ? 0 : rows_[*p].supportBegin;
                const int predecessorEnd = *p < 0 ? virtualEnd : rows_[*p].supportEnd;
                if (predecessorBegin < predecessorEnd)
                {
                    seedBegin = std::min(seedBegin, std::max(1, predecessorBegin));
                    seedEnd = std::max(seedEnd, predecessorEnd + 1);
                }
            }
            const int structuralBegin = std::max(0, shortest - budget / 2);
            const int structuralEnd = std::min(queryLength_, longest + budget / 7) + 1;
            row.begin = std::max(structuralBegin, seedBegin);
            // Every insertion adds at least 7, including the i=0 F
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
            if (graphtools::verify::kFusedRows && predecessors == 1 && width <= 64)
            {
                const int p = *first;
                const int begin = std::max(1, row.begin);
                const int predecessorBegin = p == -1 ? 0 : rows_[p].supportBegin;
                const int predecessorEnd = p == -1 ? virtualEnd : rows_[p].supportEnd;
                const int deletionBegin = std::max(begin, predecessorBegin);
                const int deletionEnd = std::min(row.end, predecessorEnd);
                const int matchBegin = std::max(begin, predecessorBegin + 1);
                const int matchEnd = std::min(row.end, predecessorEnd + 1);
                const int relaxations = std::max(0, deletionEnd - deletionBegin)
                    + std::max(0, matchEnd - matchBegin);
                const int insertionCells = row.end - begin;
                // These are the same predecessor and insertion charges as the
                // general path. Failure may occur earlier within this row, but
                // both paths then discard it and initialize dense DP.
                if (!charge(1 + uint64_t(relaxations) + insertionCells, spent, limit))
                    return FillResult::Exhausted;
                const Cell* const previous = p < 0 ? nullptr : cells_.data() + rows_[p].offset;
                int previousInsertion = infinity();
                if (row.begin == 0)
                {
                    current[0] = { Cost(10 * shortest), infinity(), false };
                    previousInsertion = current[0].h;
                    row.supportBegin = 0;
                    row.supportEnd = 1;
                }
                int previousBest = previousInsertion;
                const int cutoff = std::min(budget, int(infinity()) - 1);
                // Past the last predecessor match/deletion, only insertions
                // remain. Compute that suffix directly and omit unreachable
                // cells from lookups; work-budget accounting stays unchanged.
                const int lastSeed = std::max(begin, std::max(deletionEnd, matchEnd));
                for (int i = begin; i < lastSeed; ++i)
                {
                    int deletion = infinity(), match = infinity();
                    if (i >= deletionBegin && i < deletionEnd)
                    {
                        if (p < 0)
                            deletion = 15 * i + 2;
                        else
                        {
                            const Cell& source = previous[i - rows_[p].begin];
                            deletion = std::min(int(source.e) + 2, int(source.h) + 10);
                        }
                    }
                    if (i >= matchBegin && i < matchEnd)
                    {
                        const Cost source = p < 0 ? Cost(15 * (i - 1)) : previous[i - 1 - rows_[p].begin].h;
                        match = int(source) + (penalty_.isMatch(query_[i - 1], target_[t]) ? 0 : 9);
                    }
                    // Every extension has nonnegative cost. An over-budget
                    // predecessor can never produce a finite cell, so clip at
                    // the two stored values instead of at every relaxation.
                    const int insertion = std::min(previousInsertion + 7, previousBest + 15);
                    const int best = std::min(match, std::min(deletion, insertion));
                    const bool reachable = best <= cutoff;
                    current[i - row.begin] = { reachable ? Cost(best) : infinity(),
                        deletion <= cutoff ? Cost(deletion) : infinity(), reachable && best == insertion };
                    previousInsertion = insertion;
                    previousBest = best;
                    if (reachable)
                    {
                        if (row.supportBegin == row.supportEnd)
                            row.supportBegin = i;
                        row.supportEnd = i + 1;
                    }
                }
                int i = lastSeed;
                int tailCost = std::min(previousInsertion + 7, previousBest + 15);
                if (i < row.end && tailCost <= cutoff)
                {
                    if (row.supportBegin == row.supportEnd)
                        row.supportBegin = i;
                    do
                    {
                        current[i - row.begin] = { Cost(tailCost), infinity(), true };
                        ++i;
                        // Opening from the preceding insertion would cost 15.
                        // Extending the same insertion costs 7.
                        tailCost += 7;
                    } while (i < row.end && tailCost <= cutoff);
                    row.supportEnd = i;
                }
                if (row.supportEnd == queryLength_ + 1)
                {
                    endpoint = true;
                    finiteEndpoints_.push_back(t);
                }
                usedCells -= row.end - i;
                row.end = i;
                continue;
            }
            std::fill(current, current + width, unreachable);
            if (row.begin == 0)
                current[0].h = Cost(10 * shortest);
            const int begin = std::max(1, row.begin);
            for (auto p = first; p != last; ++p)
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
                if (*p == -1)
                {
                    // v5.0.0's virtual row initializes H and E alike.
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
                    if (i == queryLength_)
                    {
                        endpoint = true;
                        finiteEndpoints_.push_back(t);
                    }
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
