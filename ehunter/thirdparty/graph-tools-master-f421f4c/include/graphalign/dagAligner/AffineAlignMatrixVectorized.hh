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

#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#endif

#include "Details.hh"

#if defined(GRAPHTOOLS_DP_TEST_FORCE_AVX2) && defined(GRAPHTOOLS_DP_TEST_FORCE_SCALAR)
#error "Conflicting DP test dispatch overrides"
#endif

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
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
        const bool useAvx2Insertion_;
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
            // Keep the complete SIMD block inside the vector's size. First
            // restore the previous logical size so resize preserves exactly
            // the legacy padding values and capacity-growth history.
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
            const int qLen = actualQueryLength_;
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

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
        static bool avx2InsertionSupported()
        {
#if defined(GRAPHTOOLS_DP_TEST_FORCE_SCALAR)
            return false;
#elif defined(GRAPHTOOLS_DP_TEST_FORCE_AVX2)
            // Test-only override for translators whose CPUID hides AVX2.
            // Production builds always use the runtime capability check.
            return true;
#else
            static const bool supported = __builtin_cpu_supports("avx2");
            return supported;
#endif
        }

#ifdef GRAPHTOOLS_DP_TEST_INSTRUMENT_AVX2
    public:
        static uint64_t& avx2InsertionTestBlockCounter()
        {
            static thread_local uint64_t counter = 0;
            return counter;
        }

    private:
#endif

        // Only scalar arguments cross the baseline-ISA/AVX2 boundary. GCC and
        // Clang emit this helper separately, without requiring global -mavx2.
        static int __attribute__((target("avx2"), noinline)) recomputeForInsertionAvx2(
            Score* v, int i, int end, Score& previousF, Score& previousV, uint64_t& trace, int start)
        {
            // Across sixteen cells the earliest opening can lose 10 + 15*2
            // points. This guard also protects every dominated scalar opening,
            // so legacy int16 narrowing cannot change which candidate wins.
            const Score safeMinimum = SCORE_MIN + 40;
            const __m256i safeMinimumVector = _mm256_set1_epi16(safeMinimum);
            const __m256i bottom = _mm256_set1_epi16(SCORE_MIN);
            for (; end - i >= 16; i += 16)
            {
                const __m256i base = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(v + i));
                if (previousF < safeMinimum || previousV < safeMinimum
                    || _mm256_movemask_epi8(_mm256_cmpgt_epi16(safeMinimumVector, base)))
                {
                    break;
                }
#ifdef GRAPHTOOLS_DP_TEST_INSTRUMENT_AVX2
                ++avx2InsertionTestBlockCounter();
#endif
                // permute + alignr shifts across the 128-bit lane boundary.
                // The vacated positions contain negative infinity, not zero.
                __m256i insertion = _mm256_add_epi16(
                    _mm256_insert_epi16(
                        _mm256_alignr_epi8(base, _mm256_permute2x128_si256(bottom, base, 0x20), 14),
                        previousV, 0),
                    _mm256_set1_epi16(-10));
                insertion = _mm256_insert_epi16(
                    insertion, std::max(Score(previousF - 2), Score(previousV - 10)), 0);

                __m256i shifted = _mm256_adds_epi16(insertion, _mm256_set1_epi16(-2));
                insertion = _mm256_max_epi16(insertion,
                    _mm256_alignr_epi8(shifted, _mm256_permute2x128_si256(bottom, shifted, 0x20), 14));
                shifted = _mm256_adds_epi16(insertion, _mm256_set1_epi16(-4));
                insertion = _mm256_max_epi16(insertion,
                    _mm256_alignr_epi8(shifted, _mm256_permute2x128_si256(bottom, shifted, 0x20), 12));
                shifted = _mm256_adds_epi16(insertion, _mm256_set1_epi16(-8));
                insertion = _mm256_max_epi16(insertion,
                    _mm256_alignr_epi8(shifted, _mm256_permute2x128_si256(bottom, shifted, 0x20), 8));
                shifted = _mm256_adds_epi16(insertion, _mm256_set1_epi16(-16));
                insertion = _mm256_max_epi16(insertion, _mm256_permute2x128_si256(bottom, shifted, 0x20));

                const __m256i current = _mm256_max_epi16(base, insertion);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(v + i), current);
                const __m256i equal = _mm256_cmpeq_epi16(current, insertion);
                const __m128i traceBytes
                    = _mm_packs_epi16(_mm256_castsi256_si128(equal), _mm256_extracti128_si256(equal, 1));
                const uint64_t sixteenBits = _mm_movemask_epi8(traceBytes);
                trace |= sixteenBits << (i - start);
                previousF = Score(_mm256_extract_epi16(insertion, 15));
                previousV = Score(_mm256_extract_epi16(current, 15));
            }
            return i;
        }
#endif

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
#if defined(__aarch64__) && defined(__ARM_NEON) && !defined(GRAPHTOOLS_DP_TEST_FORCE_SCALAR)
                    if (gapOpen_ == -8 && gapExt_ == -2)
                    {
                        // With nonpositive opening scores and no narrowing
                        // overflow, opening from the preceding F cannot beat
                        // extending it. The eight-cell recurrence is therefore
                        // a max-plus prefix scan over the unmodified V cells.
                        // A 24-point margin covers every scalar candidate in
                        // this block, including the dominated gap openings.
                        const Score safeMinimum = SCORE_MIN + 24;
                        const int16x8_t bottom = vdupq_n_s16(SCORE_MIN);
                        const uint8x8_t traceWeights = { 1, 2, 4, 8, 16, 32, 64, 128 };
                        for (; end - i >= 8; i += 8)
                        {
                            const int16x8_t base = vld1q_s16(v + i);
                            if (previousF < safeMinimum || previousV < safeMinimum
                                || vminvq_s16(base) < safeMinimum)
                            {
                                // Keep the original modular int16 behavior
                                // when the prefix-scan algebra is not valid.
                                break;
                            }

                            int16x8_t insertion = vaddq_s16(
                                vextq_s16(vdupq_n_s16(previousV), base, 7), vdupq_n_s16(-10));
                            insertion = vsetq_lane_s16(
                                std::max(Score(previousF - 2), Score(previousV - 10)), insertion, 0);
                            // Saturation applies only to negative-infinity
                            // padding lanes; the range guard protects scores.
                            insertion = vmaxq_s16(insertion,
                                vextq_s16(bottom, vqaddq_s16(insertion, vdupq_n_s16(-2)), 7));
                            insertion = vmaxq_s16(insertion,
                                vextq_s16(bottom, vqaddq_s16(insertion, vdupq_n_s16(-4)), 6));
                            insertion = vmaxq_s16(insertion,
                                vextq_s16(bottom, vqaddq_s16(insertion, vdupq_n_s16(-8)), 4));

                            const int16x8_t current = vmaxq_s16(base, insertion);
                            vst1q_s16(v + i, current);
                            const uint8x8_t traceLanes = vmovn_u16(vceqq_s16(current, insertion));
                            const uint64_t eightBits = vaddv_u8(vand_u8(traceLanes, traceWeights));
                            trace |= eightBits << (i - start);
                            previousF = vgetq_lane_s16(insertion, 7);
                            previousV = vgetq_lane_s16(current, 7);
                        }
                    }
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                    if (useAvx2Insertion_ && end - i >= 16)
                    {
                        i = recomputeForInsertionAvx2(v, i, end, previousF, previousV, trace, start);
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
