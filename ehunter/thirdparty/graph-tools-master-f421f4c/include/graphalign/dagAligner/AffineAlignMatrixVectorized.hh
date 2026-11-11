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

#pragma once

#include <algorithm>
#include <iostream>

#include "Details.hh"

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
        // The legacy traceback compares a promoted gap sum with a narrowed score.
        // Record whether that sum can equal any Score, once per scoring object.
        const bool gapOpeningRepresentable_;

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
        std::vector<typename PenaltyMatrix::TargetChar> target_;
        std::vector<Score> alignmentPenalties_[PenaltyMatrix::TARGET_CHAR_MAX_ + 1];

    public:
        AffineAlignMatrixVectorized(const PenaltyMatrix& penaltyMatrix, Score gapOpen, Score gapExt)
            : penaltyMatrix_(penaltyMatrix)
            , gapOpen_(gapOpen)
            , gapExt_(gapExt)
            , gapOpeningRepresentable_(Score(int(gapOpen) + int(gapExt)) == int(gapOpen) + int(gapExt))
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

            // avoid "uninitialized read" complaints from valgrind
            query_.resize((std::distance(queryBegin, queryEnd) + step - 1) / step * step);
            query_.clear();
            penaltyMatrix_.translateQuery(queryBegin, queryEnd, std::back_inserter(query_));
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
        int queryLen() const { return query_.size(); }

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
        void reset(const EdgeMap& edgeMap)
        {
            const int qLen = query_.size();
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
                // Preserve every profile's original size/capacity history. Lazy
                // resizing can cause extra geometric growth when an absent code
                // returns after the query length has changed.
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

            v_.reset(qLen, tLen);
            g_.resize((qLen + step - 1) / step * step);
            insertionLeft_.resize(std::size_t(tLen) + 1);
            insertionLeft_.front() = 0;
            std::fill(insertionLeft_.begin() + 1, insertionLeft_.end(), SCORE_MIN);
            insertionWordsPerRow_ = (std::size_t(qLen) + 63) / 64;
            insertionTrace_.resize(insertionWordsPerRow_ * (std::size_t(tLen) + 1));
            e_.reset(qLen, tLen);

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

            // first column penalises for deletion
            for (int t = 0; t < tLen; ++t)
            {
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

        // __attribute((noinline))
        void fill(const EdgeMap& edgeMap)
        {
            const int qLen = query_.size();
            const int tLen = target_.size();

            for (int t = 0; t < tLen; ++t)
            {
                const typename PenaltyMatrix::TargetChar tc = target_[t];
                const Score* const penalties = &alignmentPenalties_[tc].front();
                const auto prevBegin = edgeMap.prevNodesBegin(t);
                const auto prevEnd = edgeMap.prevNodesEnd(t);
                if (prevEnd - prevBegin == 1)
                {
                    // Most target positions have one predecessor. Since the
                    // current E/G/V cells start at SCORE_MIN, their first max
                    // is just the candidate. Fuse the three passes and avoid
                    // materializing G for these rows.
                    for (int q = 0; q < qLen; q += step)
                    {
                        initializeFromPredecessor(q, t, *prevBegin, penalties + q);
                    }
                }
                else
                {
                    std::fill(g_.begin(), g_.end(), SCORE_MIN);
                    for (auto prev = prevBegin; prev != prevEnd; ++prev)
                    {
                        for (int q = 0; q < qLen; q += step)
                        {
                            recomputeForDeletion(q, t, *prev);
                            recomputeForAlign(q, *prev, penalties + q);
                        }
                    }
                    for (int q = 0; q < qLen; q += step)
                    {
                        consolidate(q, t);
                    }
                }
                recomputeForInsertion(qLen, t);
            }
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

        void initializeFromPredecessor(int q, int t, int p, const Score penalties[step])
        {
            const Score* ep = e_.row(q, p);
            const Score* vp = v_.row(q, p);
            Score deletion[step];
            Score match[step];
            for (int i = 0; i < step; ++i)
            {
                // Narrow each candidate before taking max, just as in the
                // general recurrence (including scores near the 16-bit limits).
                const Score extension = ep[i] + gapExt_;
                const Score opening = vp[i] + gapOpen_ + gapExt_;
                deletion[i] = std::max(extension, opening);
                match[i] = vp[i - 1] + penalties[i];
            }
            // Read predecessor blocks before writing the current block, keeping
            // these loops vectorizable without assuming pointers do not alias.
            Score* e = e_.row(q, t);
            for (int i = 0; i < step; ++i)
            {
                e[i] = deletion[i];
            }
            Score* v = v_.row(q, t);
            for (int i = 0; i < step; ++i)
            {
                v[i] = std::max(deletion[i], match[i]);
            }
        }

        // __attribute((noinline))
        void recomputeForAlign(int q, int p, const Score penalties[step])
        {
            Score tmp[step];
            Score* v = v_.row(q - 1, p);
            for (int i = 0; i < step; ++i)
            {
                tmp[i] = v[i] + penalties[i];
            }

            Score* g = g_.data() + q;
            for (int i = 0; i < step; ++i)
            {
                g[i] = g[i] > tmp[i] ? g[i] : tmp[i];
            }
        }

        // __attribute((noinline))
        void recomputeForDeletion(int q, int t, int p)
        {
            Score* ep = e_.row(q, p);
            Score tmpEp[step];
            for (int i = 0; i < step; ++i)
            {
                tmpEp[i] = ep[i] + gapExt_;
            }

            Score* vp = v_.row(q, p);
            Score tmpVp[step];
            for (int i = 0; i < step; ++i)
            {
                tmpVp[i] = vp[i] + gapOpen_ + gapExt_;
            }

            Score tmp[step];
            for (int i = 0; i < step; ++i)
            {
                tmp[i] = tmpEp[i] > tmpVp[i] ? tmpEp[i] : tmpVp[i];
            }

            Score* e = e_.row(q, t);
            for (int i = 0; i < step; ++i)
            {
                e[i] = e[i] > tmp[i] ? e[i] : tmp[i];
            }
        }

        // __attribute((noinline))
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
                // Preserve the original recurrence for exceptional scoring
                // parameters and for the final zero to three cells of a word.
                for (; i < end; ++i)
                {
                    // f at an unprocessed interior cell is SCORE_MIN. Preserve
                    // each original Score conversion before comparing the terms.
                    const Score extension = previousF + gapExt_;
                    const Score opening = previousV + gapOpen_ + gapExt_;
                    const Score currentF = std::max(extension, opening);
                    const Score currentV = std::max(v[i], currentF);
                    // Equality of modular score differences is equality with the
                    // original narrowed terms. V >= both terms; if the promoted
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

        // __attribute((noinline))
        void consolidate(int q, int t)
        {
            Score tmp[step];
            Score* e = e_.row(q, t);
            Score* g = g_.data() + q;
            for (int i = 0; i < step; ++i)
            {
                tmp[i] = g[i] > e[i] ? g[i] : e[i];
            }

            Score* v = v_.row(q, t);
            for (int i = 0; i < step; ++i)
            {
                v[i] = v[i] > tmp[i] ? v[i] : tmp[i];
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const AffineAlignMatrixVectorized& matrix)
        {
            return os << "AffineAlignMatrix(" << matrix.v_ << ")";
        }
    };

} // namespace dagAligner

} // namespace graphalign
