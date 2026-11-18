//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from the GraphTools copy shipped with ExpansionHunter v5.0.0: affine matrices use SIMD row kernels, compact insertion traceback, and bounded row stopping.
// Modifications are licensed under the Apache License, Version 2.0.
//

#pragma once

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>

#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#endif

#include "DpKernels.hh"
#include "Details.hh"
#include "graphutils/VerifyPolicy.hh"

namespace graphalign
{

namespace dagAligner
{

    // the 2-d table of scores filled during the alignment
    template <typename PenaltyMatrixT, bool penalizeMove, int step = 16> class AffineAlignMatrixVectorized
    {
    public:
        typedef PenaltyMatrixT PenaltyMatrix;

    private:
        const PenaltyMatrix penaltyMatrix_;
        const Score gapOpen_;
        const Score gapExt_;
        // v5.0.0's traceback compares a promoted gap sum with a narrowed score.
        // Record whether that sum can equal any Score, once per scoring object.
        const bool gapOpeningRepresentable_;
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
        const bool useAvx2Insertion_;
        const bool useAvx512Insertion_;
#endif

        PaddedAlignMatrix<step> v_;
        // Match scores are consumed in the current row; traceback uses v_, e_, and insertionTrace_.
        std::vector<Score> g_;
        // The insertion recurrence is horizontal. Historical f scores are needed
        // only through isInsertion(), so retain that exact predicate as one bit.
        std::vector<Score> insertionLeft_;
        std::vector<uint64_t> insertionTrace_;
        std::size_t insertionWordsPerRow_ = 0;
        PaddedAlignMatrix<step> e_;

        std::vector<typename PenaltyMatrix::QueryChar> query_;
        std::size_t actualQueryLength_ = 0;

    public:
        // Optional early stop of the fill (for callers using the default
        // scores, match 5, on inputs within the 16-bit guard): once no later
        // row can hold a full-query endpoint as good as the best one so far,
        // the remaining rows are not filled and endpoint scans stop there.
        enum class RowStop
        {
            kOff,
            kStop
        };
        void setRowStop(RowStop mode) { rowStop_ = mode; }
        // Rows filled; the target length when no stop applied.
        int keptRows() const { return keptRows_; }
        // Rows that endpoint scans may visit.
        int scanRows() const { return rowStop_ == RowStop::kStop ? keptRows_ : int(target_.size()); }

    private:
        RowStop rowStop_ = RowStop::kOff;
        int keptRows_ = 0;
        // rowStopBound_[t]: lowest endpoint cost any row >= t can have.
        std::vector<int> rowStopBound_, rowShortest_, rowLongest_;
        std::vector<typename PenaltyMatrix::TargetChar> target_;
        std::vector<Score> alignmentPenalties_[PenaltyMatrix::TARGET_CHAR_MAX_ + 1];

    public:
        AffineAlignMatrixVectorized(const PenaltyMatrix& penaltyMatrix, Score gapOpen, Score gapExt)
            : penaltyMatrix_(penaltyMatrix)
            , gapOpen_(gapOpen)
            , gapExt_(gapExt)
            , gapOpeningRepresentable_(Score(int(gapOpen) + int(gapExt)) == int(gapOpen) + int(gapExt))
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
            , useAvx2Insertion_(gapOpen == -8 && gapExt == -2 && avx2InsertionSupported())
            , useAvx512Insertion_(gapOpen == -8 && gapExt == -2 && avx512InsertionSupported())
#endif
        {
        }

        template <typename QueryIt, typename TargetIt>
        void
        init(QueryIt queryBegin, QueryIt queryEnd, TargetIt targetBegin, TargetIt targetEnd, const EdgeMap& edgeMap)
        {
            if (queryEnd == queryBegin)
            {
                throw std::logic_error("Empty query is not allowed.");
            }

            if (targetEnd == targetBegin)
            {
                throw std::logic_error("Empty target is not allowed.");
            }

            const std::size_t queryLength = std::distance(queryBegin, queryEnd);
            // Keep the complete SIMD block inside the vector's size. Resizing
            // from the logical length preserves v5.0.0's padding values
            // and capacity growth across calls.
            query_.resize(actualQueryLength_);
            query_.resize((queryLength + step - 1) / step * step);
            penaltyMatrix_.translateQuery(queryBegin, queryEnd, query_.begin());
            actualQueryLength_ = queryLength;
            target_.clear();
            penaltyMatrix_.translateTarget(targetBegin, targetEnd, std::back_inserter(target_));

            reset(edgeMap);

            fill(edgeMap);
        }

        typedef AlignMatrix::const_iterator const_iterator;
        template <bool localAlign> const_iterator nextBestAlign(const_iterator start, Score& bestScore) const
        {
            return !localAlign ? v_.nextBestAlign(start, queryLen() - 1, bestScore)
                               : v_.nextBestAlign(start, bestScore);
        }
        const_iterator alignBegin() const { return v_.cellOneOne(); }
        const_iterator alignEnd() const { return v_.end(); }
        int targetOffset(const_iterator cell) const { return std::distance(v_.cellOneOne(), cell) / v_.paddedRowLen(); }
        int queryOffset(const_iterator cell) const { return std::distance(v_.cellOneOne(), cell) % v_.paddedRowLen(); }
        int queryLen() const { return actualQueryLength_; }

        bool isInsertion(int q, int t) const
        {
            if (-1 == q)
            {
                return false;
            }

            const std::size_t word = (t + 1) * insertionWordsPerRow_ + q / 64;
            return (insertionTrace_[word] >> (q % 64)) & uint64_t(1);
        }

        bool isDeletion(int q, int t, int p) const
        {
            // q == -1 is ok here, just check the score match as usual
            const Score delExtScore = v_.at(q, t) - e_.at(q, p);
            const Score delOpenScore = v_.at(q, t) - v_.at(q, p);
            return gapExt_ == delExtScore || gapOpen_ + gapExt_ == delOpenScore;
        }

        bool isMatch(int q, int t, int p) const
        {
            if (-1 == q)
            {
                return false;
            }

            typename PenaltyMatrix::QueryChar queryChar = query_[q];
            typename PenaltyMatrix::TargetChar targetChar = target_[t];
            const Score alnScore = v_.at(q, t) - v_.at(q - 1, p);
            return penaltyMatrix_.isMatch(queryChar, targetChar) && penaltyMatrix_(queryChar, targetChar) == alnScore;
        }

        bool isMismatch(int q, int t, int p) const
        {
            if (-1 == q)
            {
                return false;
            }

            typename PenaltyMatrix::QueryChar queryChar = query_[q];
            typename PenaltyMatrix::TargetChar targetChar = target_[t];
            const Score alnScore = v_.at(q, t) - v_.at(q - 1, p);
            return !penaltyMatrix_.isMatch(queryChar, targetChar) && penaltyMatrix_(queryChar, targetChar) == alnScore;
        }

    private:
        // __attribute((noinline))
        // Contract [K-04]
        // Requires: the edge map was built by initializeNodeSpans or the single-node branch (every edge goes from a lower
        //           to a higher row) whenever hasOrderedNodeSpans holds.
        // Ensures:  every cell that fill, the endpoint scan or the traceback reads was written in the same call, so the
        //           observable results do not depend on the matrices' previous contents.
        // Reference: the same reset with every interior cell cleared (fillInterior = true).
        void reset(const EdgeMap& edgeMap)
        {
            const int qLen = actualQueryLength_;
            const int tLen = target_.size();

            // Only profiles referenced by this target can be read in fill().
            // The presence table is alphabet-sized, independent of DP dimensions.
            bool usedTargetChars[PenaltyMatrix::TARGET_CHAR_MAX_ + 1] = {};
            for (const auto tc : target_)
            {
                usedTargetChars[tc] = true;
            }
            for (typename PenaltyMatrix::TargetChar tc = 0; tc <= PenaltyMatrix::TARGET_CHAR_MAX_; ++tc)
            {
                // Resize every profile on each call to preserve v5.0.0's size
                // and capacity growth, including for codes absent from the target.
                alignmentPenalties_[tc].resize((qLen + step - 1) / step * step, 0);
                if (!usedTargetChars[tc])
                {
                    continue;
                }
                for (int q = 0; q < qLen; q += step)
                {
                    computeAlignPenalties(q, tc, &alignmentPenalties_[tc].front() + q);
                }
            }

            const bool ordered = edgeMap.hasOrderedNodeSpans(tLen);
            v_.reset(qLen, tLen, !ordered);
            g_.resize((qLen + step - 1) / step * step);
            insertionLeft_.resize(std::size_t(tLen) + 1);
            insertionLeft_.front() = 0;
            std::fill(insertionLeft_.begin() + 1, insertionLeft_.end(), SCORE_MIN);
            insertionWordsPerRow_ = (std::size_t(qLen) + 63) / 64;
            insertionTrace_.resize(insertionWordsPerRow_ * (std::size_t(tLen) + 1));
            e_.reset(qLen, tLen, !ordered);

            // top left must be 0 and never change
            if (v_.at(-1, -1))
            {
                throw std::logic_error("Incorrectly initialized v_");
            }
            if (insertionLeft_.front())
            {
                throw std::logic_error("Incorrectly initialized f_");
            }
            if (e_.at(-1, -1))
            {
                throw std::logic_error("Incorrectly initialized e_");
            }

            // A single-predecessor row overwrites all interior and padding
            // cells before any successor can read them on a validated DAG.
            // Junctions still accumulate into SCORE_MIN; every row needs its
            // boundary before the first-column recurrence.
            for (int t = 0; t < tLen; ++t)
            {
                if (ordered)
                {
                    v_.at(-1, t) = e_.at(-1, t) = SCORE_MIN;
                    // Rows that take the general (accumulating) row code start from SCORE_MIN.
                    if (edgeMap.prevNodesEnd(t) - edgeMap.prevNodesBegin(t) != 1 || !graphtools::verify::kFusedRows)
                    {
                        std::fill_n(v_.row(0, t), v_.paddedRowLen() - 1, SCORE_MIN);
                        std::fill_n(e_.row(0, t), e_.paddedRowLen() - 1, SCORE_MIN);
                    }
                }
                if (penalizeMove)
                {
                    for (EdgeMap::OffsetEdges::const_iterator prevNodeIndexIt = edgeMap.prevNodesBegin(t);
                         prevNodeIndexIt != edgeMap.prevNodesEnd(t); ++prevNodeIndexIt)
                    {
                        const int p = *prevNodeIndexIt;
                        v_.at(-1, t) = std::max(v_.at(-1, t), Score(v_.at(-1, p) + gapOpen_ + gapExt_));
                        insertionLeft_[t + 1] = std::max(v_.at(-1, t), Score(insertionLeft_[p + 1] + gapOpen_ + gapExt_));
                    }
                }
                else
                {
                    v_.at(-1, t) = 0;
                    insertionLeft_[t + 1] = 0;
                }
            }

            // first row penalizes for insertion
            for (int q = 0; q < queryLen(); ++q)
            {
                v_.at(q, -1) = v_.at(q - 1, -1) + gapOpen_ + gapExt_;
                e_.at(q, -1) = e_.at(q - 1, -1) + gapOpen_ + gapExt_;
            }

            // Preserve the query boundary's observable insertion predicates too.
            Score previousF = insertionLeft_.front();
            for (int start = 0; start < qLen;)
            {
                uint64_t trace = 0;
                const int end = start + std::min(64, qLen - start);
                for (int q = start; q < end; ++q)
                {
                    const Score insExtScore = v_.at(q, -1) - previousF;
                    const Score insOpenScore = v_.at(q, -1) - v_.at(q - 1, -1);
                    const bool insertion = gapExt_ == insExtScore || gapOpen_ + gapExt_ == insOpenScore;
                    trace |= uint64_t(insertion) << (q - start);
                    previousF = SCORE_MIN;
                }
                insertionTrace_[start / 64] = trace;
                start = end;
            }
        }

        // Every full-query alignment ending in a row pays 2 per base by which
        // the shortest path to the row exceeds the query, and 7 per base by
        // which the query exceeds the longest path (opening penalties
        // omitted, so this holds under both boundary rules). Returns
        // false for an edge map that is not a DAG in row order.
        // Contract [K-05]
        // Requires: default scores (5, -4, -8, -2); 15q + 10t + 10 < 32768; edges form a DAG in row order.
        // Ensures:  best score, the optimal endpoints in order and the four traceback predicates on every cell an optimal
        //           path visits equal those of the full fill. The second-best score is not preserved (callers discard it).
        // Lemma:    a full-query alignment ending in row t costs at least
        //           2 * max(shortest(t) - q, 0) + 7 * max(q - longest(t), 0), where cost = 5q - score.
        // Reference: the full fill (RowStop::kOff; EH_VERIFY_ROW_STOP=OFF).
        bool prepareRowStop(const EdgeMap& edgeMap, int qLen, int tLen)
        {
            rowShortest_.resize(tLen);
            rowLongest_.resize(tLen);
            rowStopBound_.assign(std::size_t(tLen) + 1, std::numeric_limits<int>::max());
            for (int t = 0; t < tLen; ++t)
            {
                int shortest = tLen + 1, longest = 0;
                for (auto it = edgeMap.prevNodesBegin(t); it != edgeMap.prevNodesEnd(t); ++it)
                {
                    const int p = *it;
                    if (p < -1 || p >= t)
                        return false;
                    shortest = std::min(shortest, p == -1 ? 1 : rowShortest_[p] + 1);
                    longest = std::max(longest, p == -1 ? 1 : rowLongest_[p] + 1);
                }
                rowShortest_[t] = shortest;
                rowLongest_[t] = longest;
            }
            for (int t = tLen - 1; t >= 0; --t)
            {
                const int bound = 2 * std::max(rowShortest_[t] - qLen, 0) + 7 * std::max(qLen - rowLongest_[t], 0);
                rowStopBound_[t] = std::min(bound, rowStopBound_[t + 1]);
            }
            return true;
        }

        // __attribute((noinline))
        void fill(const EdgeMap& edgeMap)
        {
            const int qLen = actualQueryLength_;
            const int tLen = target_.size();
            keptRows_ = tLen;
            const bool stopping = rowStop_ != RowStop::kOff && prepareRowStop(edgeMap, qLen, tLen);
            // Every row pass covers whole blocks of step cells, as the padded rows allow.
            const int lanes = (qLen + step - 1) / step * step;
            int best = std::numeric_limits<int>::min();

            for (int t = 0; t < tLen; ++t)
            {
                const typename PenaltyMatrix::TargetChar tc = target_[t];
                const Score* const penalties = &alignmentPenalties_[tc].front();
                const auto prevBegin = edgeMap.prevNodesBegin(t);
                const auto prevEnd = edgeMap.prevNodesEnd(t);
                if (graphtools::verify::kFusedRows && prevEnd - prevBegin == 1)
                {
                    // Most target positions have one predecessor. Since the
                    // current E/G/V cells start at SCORE_MIN, their first max
                    // is just the candidate. Fuse the three passes and avoid
                    // materializing G for these rows.
                    fusedRow(t, *prevBegin, penalties, lanes);
                }
                else
                {
                    std::fill(g_.begin(), g_.end(), SCORE_MIN);
                    for (auto prev = prevBegin; prev != prevEnd; ++prev)
                    {
                        deletionRow(t, *prev, lanes);
                        alignRow(*prev, penalties, lanes);
                    }
                    consolidateRow(t, lanes);
                }
                recomputeForInsertion(qLen, t);
                if (stopping && keptRows_ == tLen)
                {
                    best = std::max(best, int(v_.at(qLen - 1, t)));
                    // Strict: a later endpoint with the same cost would be co-optimal.
                    if (rowStopBound_[t + 1] > 5 * qLen - best)
                    {
                        keptRows_ = t + 1;
                        // A shadow build fills every row and checks the dropped ones below.
                        if (rowStop_ == RowStop::kStop && !graphtools::verify::kShadow)
                            break;
                    }
                }
            }
            if (rowStop_ != RowStop::kStop || keptRows_ == tLen)
                return;
            if (graphtools::verify::kShadow)
            {
                // Check the dropped rows against the stop rule, then overwrite them
                // so that any later read of a dropped row changes the output.
                for (int t = keptRows_; t < tLen; ++t)
                    if (int(v_.at(qLen - 1, t)) >= best)
                        throw std::logic_error("Dense row stop dropped a row holding an optimal endpoint");
                for (int t = keptRows_; t < tLen; ++t)
                {
                    std::fill_n(v_.row(-1, t), v_.paddedRowLen(), std::numeric_limits<Score>::max());
                    std::fill_n(e_.row(-1, t), e_.paddedRowLen(), std::numeric_limits<Score>::max());
                }
            }
            v_.limitRows(keptRows_);
            e_.limitRows(keptRows_);
        }

        // __attribute((noinline))
        void computeAlignPenalties(int q, typename PenaltyMatrix::TargetChar tc, Score penalties[step])
        {
            const typename PenaltyMatrix::QueryChar* query = &query_[q];
            for (int i = 0; i < step; ++i)
            {
                penalties[i] = penaltyMatrix_(query[i], tc);
            }
        }

        // Row passes (K-03 rows). The ordered rows of a DAG make the predecessor row p differ from t, so the
        // arrays a pass writes never overlap the ones it reads. A build whose baseline includes AVX2 calls the
        // AVX2 kernels directly, so they can be inlined; a portable x86-64 build selects AVX2 or SSE2 at run
        // time.
#if defined(__aarch64__) && defined(__ARM_NEON)
#define EH_DP_ROW_DISPATCH(NAME, ...)                                                                                  \
    if (graphtools::verify::kScalarOnly)                                                                               \
        dp_##NAME##_row_reference(__VA_ARGS__);                                                                        \
    else                                                                                                               \
        dp_##NAME##_row_neon(__VA_ARGS__);
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && defined(__AVX2__)
#define EH_DP_ROW_DISPATCH(NAME, ...)                                                                                  \
    if (graphtools::verify::kScalarOnly)                                                                               \
        dp_##NAME##_row_reference(__VA_ARGS__);                                                                        \
    else if (graphtools::verify::kNoAvx2)                                                                              \
        dp_##NAME##_row_sse2(__VA_ARGS__);                                                                             \
    else                                                                                                               \
        dp_##NAME##_row_avx2(__VA_ARGS__);
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
        static int rowTier()
        {
            if (graphtools::verify::kScalarOnly)
                return 0;
            static const int tier = !graphtools::verify::kNoAvx2 && __builtin_cpu_supports("avx2") ? 2 : 1;
            return tier;
        }
#define EH_DP_ROW_DISPATCH(NAME, ...)                                                                                  \
    switch (rowTier())                                                                                                 \
    {                                                                                                                  \
    case 2:                                                                                                            \
        dp_##NAME##_row_avx2(__VA_ARGS__);                                                                             \
        break;                                                                                                         \
    case 1:                                                                                                            \
        dp_##NAME##_row_sse2(__VA_ARGS__);                                                                             \
        break;                                                                                                         \
    default:                                                                                                           \
        dp_##NAME##_row_reference(__VA_ARGS__);                                                                        \
    }
#else
#define EH_DP_ROW_DISPATCH(NAME, ...) dp_##NAME##_row_reference(__VA_ARGS__);
#endif

        void fusedRow(int t, int p, const Score* penalties, int lanes)
        {
            EH_DP_ROW_DISPATCH(fused, e_.row(0, p), v_.row(0, p), penalties, gapOpen_, gapExt_, e_.row(0, t),
                               v_.row(0, t), lanes)
        }

        void deletionRow(int t, int p, int lanes)
        {
            EH_DP_ROW_DISPATCH(deletion, e_.row(0, p), v_.row(0, p), gapOpen_, gapExt_, e_.row(0, t), lanes)
        }

        void alignRow(int p, const Score* penalties, int lanes)
        {
            EH_DP_ROW_DISPATCH(align, v_.row(0, p), penalties, g_.data(), lanes)
        }

        void consolidateRow(int t, int lanes)
        {
            EH_DP_ROW_DISPATCH(consolidate, g_.data(), e_.row(0, t), v_.row(0, t), lanes)
        }
#undef EH_DP_ROW_DISPATCH

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
        static bool avx512InsertionSupported()
        {
            if (graphtools::verify::kNoAvx512)
                return false;
            static const bool supported = __builtin_cpu_supports("avx512bw");
            return supported;
        }

        static bool avx2InsertionSupported()
        {
            if (graphtools::verify::kNoAvx2)
                return false;
            static const bool supported = __builtin_cpu_supports("avx2");
            return supported;
        }
#endif

        // __attribute((noinline))
        // Contract [K-02]
        // Requires: nothing (all scores, including sums that wrap in 16 bits; both values of gapOpeningRepresentable_).
        // Ensures:  v equals v5.0.0's row and trace bit q equals v5's isInsertion predicate
        //           gapExt == Score(v[q] - f[q-1]) || gapOpen + gapExt == Score(v[q] - v[q-1]), from the carried F alone.
        // Reference: v5.0.0 recomputeForInsertion and isInsertion with the full f matrix (AffineAlignMatrix.hh).
        void recomputeForInsertion(int qLen, int t)
        {
            Score* v = v_.row(0, t);
            Score previousF = insertionLeft_[t + 1];
            Score previousV = v_.at(-1, t);
            for (int start = 0; start < qLen;)
            {
                uint64_t trace = 0;
                const int end = start + std::min(64, qLen - start);
                int i = start;
                if (gapOpeningRepresentable_)
                {
#if defined(__aarch64__) && defined(__ARM_NEON)
                    // Kernel K-01 (DpKernels.hh), NEON block of 8 cells.
                    if (!graphtools::verify::kScalarOnly && gapOpen_ == -8 && gapExt_ == -2)
                        i = dp_insertion_neon(v, i, end, &previousF, &previousV, &trace, start);
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                    // Kernels K-01 (DpKernels.hh). Only scalar arguments cross the baseline-ISA
                    // boundary; no global -mavx2 or -mavx512bw is required.
                    if (useAvx512Insertion_ && end - i >= 32)
                        i = dp_insertion_avx512(v, i, end, &previousF, &previousV, &trace, start);
                    if (useAvx2Insertion_ && end - i >= 16)
                    {
                        i = dp_insertion_avx2(v, i, end, &previousF, &previousV, &trace, start);
                    }
#endif
                    // Four serial cells share one loop branch and one variable
                    // trace shift. Keep every candidate narrowed to Score before
                    // max; the next cell still consumes this cell's exact F/V.
                    for (; end - i >= 4; i += 4)
                    {
                        const Score extension0 = previousF + gapExt_;
                        const Score opening0 = previousV + gapOpen_ + gapExt_;
                        const Score currentF0 = std::max(extension0, opening0);
                        const Score currentV0 = std::max(v[i], currentF0);
                        v[i] = currentV0;

                        const Score extension1 = currentF0 + gapExt_;
                        const Score opening1 = currentV0 + gapOpen_ + gapExt_;
                        const Score currentF1 = std::max(extension1, opening1);
                        const Score currentV1 = std::max(v[i + 1], currentF1);
                        v[i + 1] = currentV1;

                        const Score extension2 = currentF1 + gapExt_;
                        const Score opening2 = currentV1 + gapOpen_ + gapExt_;
                        const Score currentF2 = std::max(extension2, opening2);
                        const Score currentV2 = std::max(v[i + 2], currentF2);
                        v[i + 2] = currentV2;

                        const Score extension3 = currentF2 + gapExt_;
                        const Score opening3 = currentV2 + gapOpen_ + gapExt_;
                        const Score currentF3 = std::max(extension3, opening3);
                        const Score currentV3 = std::max(v[i + 3], currentF3);
                        v[i + 3] = currentV3;

                        const uint64_t fourBits = uint64_t(currentV0 == currentF0)
                            | (uint64_t(currentV1 == currentF1) << 1)
                            | (uint64_t(currentV2 == currentF2) << 2)
                            | (uint64_t(currentV3 == currentF3) << 3);
                        trace |= fourBits << (i - start);
                        previousF = currentF3;
                        previousV = currentV3;
                    }
                }
                // Use v5.0.0's recurrence for exceptional scoring
                // parameters and for the final zero to three cells of a word.
                for (; i < end; ++i)
                {
                    // f at an unprocessed interior cell is SCORE_MIN. Preserve
                    // each Score conversion from v5.0.0 before comparing the terms.
                    const Score extension = previousF + gapExt_;
                    const Score opening = previousV + gapOpen_ + gapExt_;
                    const Score currentF = std::max(extension, opening);
                    const Score currentV = std::max(v[i], currentF);
                    // Equality of modular score differences is equality with the
                    // narrowed terms. V >= both terms; if the promoted
                    // opening sum is representable, either equality means V == F.
                    const Score traceScore = gapOpeningRepresentable_ ? currentF : extension;
                    const bool insertion = currentV == traceScore;
                    trace |= uint64_t(insertion) << (i - start);
                    v[i] = currentV;
                    previousF = currentF;
                    previousV = currentV;
                }
                insertionTrace_[(t + 1) * insertionWordsPerRow_ + start / 64] = trace;
                start = end;
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const AffineAlignMatrixVectorized& matrix)
        {
            return os << "AffineAlignMatrix(" << matrix.v_ << ")";
        }
    };

} // namespace dagAligner

} // namespace graphalign
