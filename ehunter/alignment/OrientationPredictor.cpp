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

#include "alignment/OrientationPredictor.hh"

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>

#include "graphcore/Path.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/SequenceOperations.hh"

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

OrientationPrediction OrientationPredictor::predict(const std::string& query) const
{
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
