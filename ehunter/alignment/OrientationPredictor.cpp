//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: orientation prediction uses encoded k-mer scans with threshold bounds and SIMD packing.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "alignment/OrientationPredictor.hh"

#include <algorithm>
#include <cstring>
#include <limits>
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <tmmintrin.h>
#define EH_PACKED_STRANDS 1
#define EH_PACKED_STRANDS_TARGET __attribute__((target("ssse3")))
#elif defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#define EH_PACKED_STRANDS 1
#define EH_PACKED_STRANDS_TARGET
#endif
#include <iostream>
#include <list>
#include <vector>

#include "graphcore/Path.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "core/VerifyPolicy.hh"

using graphtools::KmerIndex;
using graphtools::Path;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

// Fast codes for ASCII A/C/G/T in either case. Other bytes retain the original
// toupper semantics on forward scans and are invalid on reverse scans.
template <bool reverseComplement>
static uint32_t encodeOrientationBase(int raw, bool& valid)
{
    static const uint8_t codes[32] = {
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };
    const auto byte = static_cast<uint8_t>(raw);
    const uint32_t code = codes[byte & 31];
    if ((byte & 0xc0) == 0x40 && code < 4)
    {
        valid = true;
        return reverseComplement ? 3 - code : code;
    }
    if (reverseComplement)
    {
        valid = false;
        return 0;
    }
    const int base = ::toupper(raw);
    valid = base == 'A' || base == 'C' || base == 'G' || base == 'T';
    return (base == 'C') + 2 * (base == 'G') + 3 * (base == 'T');
}

// A failing scan may return a partial count below minimumMatches. A passing
// forward scan remains exact because its count determines which strand wins.
template <bool reverseComplement = false, bool stopAtThreshold = false>
static int scanNonoverlappingKmerMatches(const string& query, const KmerIndex& kmerIndex, int minimumMatches)
{
    const std::size_t kmerLength = kmerIndex.kmerLength();
    // At position p, at most floor((length - p) / k) further nonoverlapping
    // matches can fit. Update the equivalent position bound only on a hit,
    // avoiding a division or extra threshold check at every rejected k-mer.
    size_t lastPossibleStart = query.length();
    if (kmerLength > 0)
    {
        const size_t requiredMatches = static_cast<size_t>(std::max(minimumMatches, 1));
        if (requiredMatches > query.length() / kmerLength)
        {
            return 0;
        }
        lastPossibleStart -= requiredMatches * kmerLength;
    }
    int matchCount = 0;
    size_t position = 0;
    size_t scanned = 0;
    size_t validSuffixLength = 0;
    uint32_t key = 0;
    // KmerIndex uses at most 32 bits; widening the mask shift also covers k=16.
    const uint32_t mask = static_cast<uint32_t>((uint64_t(1) << (2 * kmerLength)) - 1);
    while (position <= lastPossibleStart)
    {
        // Keep the scan inside the current window, including when the window is empty.
        scanned = std::max(scanned, position);
        // Carry the exact encoding of the last k bases; each input base is read once.
        while (scanned < position + kmerLength)
        {
            bool valid;
            const int raw = reverseComplement ? query[query.size() - 1 - scanned++] : query[scanned++];
            const uint32_t code = encodeOrientationBase<reverseComplement>(raw, valid);
            key = ((key << 2) | code) & mask;
            if (valid)
                ++validSuffixLength;
            else
                validSuffixLength = 0;
        }

        if (validSuffixLength >= kmerLength && kmerIndex.containsEncodedKmer(key))
        {
            ++matchCount;
            position += kmerLength;
            if (matchCount < minimumMatches)
            {
                lastPossibleStart += kmerLength;
            }
            else if (stopAtThreshold)
            {
                return matchCount;
            }
        }
        else
        {
            ++position;
        }
    }
    return matchCount;
}

#ifdef EH_PACKED_STRANDS
// Store 16 canonical bases per word in sequence order, with a zero guard
// word for windows crossing the final boundary. Both strands share one pass.
struct PackedStrands
{
    static constexpr size_t maxBases = 256;
    uint32_t forward[maxBases / 16 + 1], reverse[maxBases / 16 + 1];
    size_t reverseOffset;

    EH_PACKED_STRANDS_TARGET bool encode(const string& query)
    {
        const size_t blocks = (query.size() + 15) / 16;
        reverseOffset = 16 * blocks - query.size();
        forward[blocks] = reverse[blocks] = 0;
        for (size_t i = 0; i < blocks; ++i)
        {
            const char* chars = query.data() + 16 * i;
            char tail[16];
            if (i + 1 == blocks && reverseOffset)
            {
                std::memset(tail, 'A', sizeof(tail));
                std::memcpy(tail, chars, 16 - reverseOffset);
                chars = tail;
            }
#if defined(__x86_64__)
            const auto bases = _mm_or_si128(_mm_loadu_si128(reinterpret_cast<const __m128i*>(chars)), _mm_set1_epi8(32));
            const auto a = _mm_cmpeq_epi8(bases, _mm_set1_epi8('a'));
            const auto c = _mm_cmpeq_epi8(bases, _mm_set1_epi8('c'));
            const auto g = _mm_cmpeq_epi8(bases, _mm_set1_epi8('g'));
            const auto t = _mm_cmpeq_epi8(bases, _mm_set1_epi8('t'));
            if (_mm_movemask_epi8(_mm_or_si128(_mm_or_si128(a, c), _mm_or_si128(g, t))) != 0xffff)
                return false;
            const auto codes = _mm_or_si128(_mm_and_si128(_mm_or_si128(c, t), _mm_set1_epi8(1)),
                                            _mm_and_si128(_mm_or_si128(g, t), _mm_set1_epi8(2)));
            const auto pairs = _mm_maddubs_epi16(codes, _mm_set1_epi16(0x0104));
            const auto quads = _mm_madd_epi16(pairs, _mm_set1_epi32(0x00010010));
            const uint32_t word = _mm_cvtsi128_si32(_mm_shuffle_epi8(quads, _mm_setr_epi8(
                12, 8, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)));
#else
            const uint8x16_t bases = vorrq_u8(vld1q_u8(reinterpret_cast<const uint8_t*>(chars)), vdupq_n_u8(32));
            const uint8x16_t a = vceqq_u8(bases, vdupq_n_u8('a'));
            const uint8x16_t c = vceqq_u8(bases, vdupq_n_u8('c'));
            const uint8x16_t g = vceqq_u8(bases, vdupq_n_u8('g'));
            const uint8x16_t t = vceqq_u8(bases, vdupq_n_u8('t'));
            if (vminvq_u8(vorrq_u8(vorrq_u8(a, c), vorrq_u8(g, t))) != 255)
                return false;
            const uint8x16_t codes = vorrq_u8(vandq_u8(vorrq_u8(c, t), vdupq_n_u8(1)),
                                              vandq_u8(vorrq_u8(g, t), vdupq_n_u8(2)));
            // Base j of each group of four gets weight 4^(3 - j); the four groups then become the bytes of
            // the word, the first group in the most significant byte.
            static const uint8_t weightValues[16] = { 64, 16, 4, 1, 64, 16, 4, 1, 64, 16, 4, 1, 64, 16, 4, 1 };
            const uint32x4_t quads = vpaddlq_u16(vpaddlq_u8(vmulq_u8(codes, vld1q_u8(weightValues))));
            const uint8x8_t quadBytes = vmovn_u16(vcombine_u16(vmovn_u32(quads), vdup_n_u16(0)));
            const uint32_t word = __builtin_bswap32(vget_lane_u32(vreinterpret_u32_u8(quadBytes), 0));
#endif
            forward[i] = word;
            uint32_t reversed = __builtin_bswap32(~word);
            reversed = ((reversed >> 4) & 0x0f0f0f0fu) | ((reversed & 0x0f0f0f0fu) << 4);
            reverse[blocks - i - 1] = ((reversed >> 2) & 0x33333333u) | ((reversed & 0x33333333u) << 2);
        }
        return true;
    }
};

static int countPackedMatches(const uint32_t* words, size_t offset, size_t length, const KmerIndex& index,
                              int minimum, bool stopAtThreshold)
{
    const size_t k = index.kmerLength();
    const size_t required = static_cast<size_t>(std::max(1, minimum));
    if (required > length / k)
        return 0;
    size_t last = offset + length - required * k;
    const auto mask = std::numeric_limits<uint32_t>::max() >> (32 - 2 * k);
    int count = 0;
    for (size_t pos = offset; pos <= last;)
    {
        const size_t block = pos / 16;
        const uint64_t bases = (uint64_t(words[block]) << 32) | words[block + 1];
        const uint32_t key = static_cast<uint32_t>(bases >> (64 - 2 * ((pos % 16) + k))) & mask;
        if (!index.containsEncodedKmer(key))
        {
            ++pos;
            continue;
        }
        ++count;
        pos += k;
        if (count < minimum)
            last += k;
        else if (stopAtThreshold)
            break;
    }
    return count;
}
#endif

// Contract [K-23]
// Requires: 1 <= k <= 16; the index answers membership as v5.0.0's does (K-15).
// Ensures:  the three-way decision equals v5.0.0's countNonoverlappingKmerMatches and predict (the counts may stop
//           early); the SSSE3 or NEON packed scan equals the scalar scan (EH_VERIFY_SIMD=SCALAR).
OrientationPrediction OrientationPredictor::predict(const std::string& query) const
{
#ifdef EH_PACKED_STRANDS
#if defined(__x86_64__)
    static const bool hasSsse3 = !verify::kScalarOnly && __builtin_cpu_supports("ssse3");
#else
    const bool hasSsse3 = !verify::kScalarOnly; // NEON is part of AArch64
#endif
    if (hasSsse3 && query.size() >= 16 && query.size() <= PackedStrands::maxBases && kmerLength_ > 0)
    {
        PackedStrands packed;
        if (packed.encode(query))
        {
            const int forward = countPackedMatches(packed.forward, 0, query.size(), kmerIndex_, minKmerMatchesToPass_, false);
            const bool passes = forward >= minKmerMatchesToPass_;
            if (passes && static_cast<size_t>(forward) == query.size() / kmerLength_)
                return OrientationPrediction::kAlignsInOriginalOrientation;
            const int threshold = passes ? forward + 1 : minKmerMatchesToPass_;
            const int reverse = countPackedMatches(packed.reverse, packed.reverseOffset, query.size(), kmerIndex_, threshold, true);
            if (reverse >= threshold)
                return OrientationPrediction::kAlignsInReverseComplementOrientation;
            return passes ? OrientationPrediction::kAlignsInOriginalOrientation : OrientationPrediction::kDoesNotAlign;
        }
    }
#endif
    const int numForwardMatches = scanNonoverlappingKmerMatches(query, kmerIndex_, minKmerMatchesToPass_);
    // Neither orientation can contain more than floor(query length / k)
    // nonoverlapping matches. Forward wins ties, so a passing forward score
    // at that upper bound cannot be improved by scanning the reverse strand.
    if (kmerIndex_.kmerLength() > 0 && numForwardMatches >= minKmerMatchesToPass_
        && static_cast<size_t>(numForwardMatches) == query.size() / kmerIndex_.kmerLength())
    {
        return OrientationPrediction::kAlignsInOriginalOrientation;
    }
    const bool forwardPasses = numForwardMatches >= minKmerMatchesToPass_;
    // Passing forward counts are exact. A failed forward scan need not be:
    // reverse only has to pass the threshold in that case. Forward wins ties.
    const int reverseMatchesToWin = forwardPasses ? numForwardMatches + 1 : minKmerMatchesToPass_;
    const int numReverseComplementMatches
        = scanNonoverlappingKmerMatches<true, true>(query, kmerIndex_, reverseMatchesToWin);
    if (numReverseComplementMatches >= reverseMatchesToWin)
    {
        return OrientationPrediction::kAlignsInReverseComplementOrientation;
    }
    return forwardPasses ? OrientationPrediction::kAlignsInOriginalOrientation : OrientationPrediction::kDoesNotAlign;
}

std::ostream& operator<<(std::ostream& out, OrientationPrediction orientationPrediction)
{
    switch (orientationPrediction)
    {
    case OrientationPrediction::kAlignsInOriginalOrientation:
        out << "kAlignsInOriginalOrientation";
        break;
    case OrientationPrediction::kAlignsInReverseComplementOrientation:
        out << "kAlignsInReverseComplementOrientation";
        break;
    case OrientationPrediction::kDoesNotAlign:
        out << "kDoesNotAlign";
    }

    return out;
}

}
