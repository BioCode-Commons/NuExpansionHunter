//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Konrad Scheffler <kscheffler@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include "genotyping/TwoAlleleStrGenotyper.hh"

#include <stack>

#include "core/LogSum.hh"

using std::unordered_set;

namespace ehunter
{
namespace strgt
{

struct Ci
{
    Ci(int begin, int end)
        : begin(begin)
        , end(end)
    {
    }
    int begin;
    int end;
};

struct CiAndLoglik
{
    CiAndLoglik(int startSize, int endSize, double loglik)
        : startSize(startSize)
        , endSize(endSize)
        , loglik(loglik)
    {
    }
    int startSize;
    int endSize;
    double loglik;
};

using LikelihoodFn = double (TwoAlleleGenotyper::*)(int, int);
static Ci getCiAlongX(int bestX, int bestY, TwoAlleleGenotyper* genotyper, LikelihoodFn likelihood)
{
    int xFrom = bestX;
    int xTo = bestX;
    int yFrom = bestY;
    int yTo = bestY;

    std::stack<CiAndLoglik> ciCandidates;
    double topGtLoglik = (genotyper->*likelihood)(bestX, bestY);
    double totalLoglik = topGtLoglik;
    ciCandidates.emplace(xFrom, xTo, totalLoglik);

    const int kMaxIntervalWidth = 750;
    double likelihoodRatio = 1;
    double llShiftLeft0 = 0;
    double llShiftLeft1 = 0;
    double llShiftLeft2 = 0;
    double llShiftRight0 = 0;
    double llShiftRight1 = 0;
    double llShiftRight2 = 0;
    bool refreshLeft = true;
    bool refreshRight = true;
    while (likelihoodRatio >= 0.01 && xTo - xFrom <= kMaxIntervalWidth)
    {
        // The unmoved boundary still has exactly the same three candidates.
        if (refreshLeft)
        {
            llShiftLeft0 = (genotyper->*likelihood)(xFrom - 1, yFrom - 1);
            llShiftLeft1 = (genotyper->*likelihood)(xFrom - 1, yFrom);
            llShiftLeft2 = (genotyper->*likelihood)(xFrom - 1, yFrom + 1);
            refreshLeft = false;
        }
        double llShiftLeft = std::max(std::max(llShiftLeft0, llShiftLeft1), llShiftLeft2);

        if (refreshRight)
        {
            llShiftRight0 = (genotyper->*likelihood)(xTo + 1, yTo + 1);
            llShiftRight1 = (genotyper->*likelihood)(xTo + 1, yTo);
            llShiftRight2 = (genotyper->*likelihood)(xTo + 1, yTo - 1);
            refreshRight = false;
        }
        double llShiftRight = std::max(std::max(llShiftRight0, llShiftRight1), llShiftRight2);

        double gtLoglik;
        if (llShiftLeft >= llShiftRight)
        {
            xFrom -= 1;
            // Preserve the original strict comparisons: a tie between 0 and 2
            // keeps y unchanged even when candidate 1 is not the maximum.
            gtLoglik = llShiftLeft1;
            if (llShiftLeft0 > llShiftLeft1 && llShiftLeft0 > llShiftLeft2)
            {
                yFrom -= 1;
                gtLoglik = llShiftLeft0;
            }
            else if (llShiftLeft2 > llShiftLeft0 && llShiftLeft2 > llShiftLeft1)
            {
                yFrom += 1;
                gtLoglik = llShiftLeft2;
            }

            refreshLeft = true;
        }
        else
        {
            xTo += 1;
            gtLoglik = llShiftRight1;
            if (llShiftRight0 > llShiftRight1 && llShiftRight0 > llShiftRight2)
            {
                yTo += 1;
                gtLoglik = llShiftRight0;
            }
            else if (llShiftRight2 > llShiftRight0 && llShiftRight2 > llShiftRight1)
            {
                yTo -= 1;
                gtLoglik = llShiftRight2;
            }

            refreshRight = true;
        }

        totalLoglik = getLogSum(totalLoglik, gtLoglik);
        ciCandidates.emplace(xFrom, xTo, totalLoglik);
        likelihoodRatio = std::exp(gtLoglik - topGtLoglik);
    }

    auto ciCandidate = ciCandidates.top();
    while (ciCandidates.size() > 1)
    {
        ciCandidates.pop();
        auto nextCiCandidate = ciCandidates.top();
        const double ciProbability = std::exp(nextCiCandidate.loglik - totalLoglik);
        if (ciProbability >= 0.95)
        {
            ciCandidate = nextCiCandidate;
        }
        else
        {
            break;
        }
    }

    return { ciCandidate.startSize, ciCandidate.endSize };
}

RepeatGenotype TwoAlleleGenotyper::genotype(const unordered_set<int>& alleleSizeCandidates)
{
    RepeatGenotype gt = getMostLikelyGenotype(alleleSizeCandidates);
    const int bestShortSize = gt.shortAlleleSizeInUnits();
    const int bestLongSize = gt.longAlleleSizeInUnits();
    Ci shortStrCi = getCiAlongX(bestShortSize, bestLongSize, this, &TwoAlleleGenotyper::getShortAndLongAlleleLoglik);
    Ci longStrCi = getCiAlongX(bestLongSize, bestShortSize, this, &TwoAlleleGenotyper::getLongAndShortAlleleLoglik);

    gt.setShortAlleleSizeInUnitsCi(shortStrCi.begin, shortStrCi.end);
    gt.setLongAlleleSizeInUnitsCi(longStrCi.begin, longStrCi.end);

    return gt;
}

RepeatGenotype TwoAlleleGenotyper::getMostLikelyGenotype(const unordered_set<int>& alleleSizeCandidates)
{
    double maxGtLoglik = std::numeric_limits<double>::lowest();
    int bestShortAlleleSize = 0;
    int bestLongAlleleSize = 0;

    for (int shortAlleleSize : alleleSizeCandidates)
    {
        for (int longAlleleSize : alleleSizeCandidates)
        {
            if (shortAlleleSize > longAlleleSize)
            {
                continue;
            }

            const double gtLoglik = getShortAndLongAlleleLoglik(shortAlleleSize, longAlleleSize);

            if (maxGtLoglik < gtLoglik)
            {
                maxGtLoglik = gtLoglik;
                bestShortAlleleSize = shortAlleleSize;
                bestLongAlleleSize = longAlleleSize;
            }
        }
    }

    return RepeatGenotype(motifLen_, { bestShortAlleleSize, bestLongAlleleSize });
}

double TwoAlleleGenotyper::getShortAndLongAlleleLoglik(int shortAlleleSize, int longAlleleSize)
{
    if (shortAlleleSize < 0 || longAlleleSize < 0 || shortAlleleSize > longAlleleSize)
    {
        return std::numeric_limits<double>::lowest();
    }

    const int shortAlleleLen = shortAlleleSize * motifLen_ + fragLen_ + 1;
    const int longAlleleLen = longAlleleSize * motifLen_ + fragLen_ + 1;
    const double shortAlleleFrac = static_cast<double>(shortAlleleLen) / (shortAlleleLen + longAlleleLen);
    double genotypeLoglik = 0;
    const int numFrags = fragLogliks_.numFrags();
    if (numFrags == 0)
    {
        return genotypeLoglik;
    }

    const double shortAllelePrior = std::log(shortAlleleFrac);
    const double longAllelePrior = std::log(1.0 - shortAlleleFrac);
    const double mismapPrior = std::log(0.001);
    const double rightmapPrior = std::log(1.0 - 0.001);
    const bool isHomozygous = shortAlleleSize == longAlleleSize;

    for (int fragIndex = 0; fragIndex != numFrags; ++fragIndex)
    {
        double fragLoglikForShortAllele = fragLogliks_.getLoglik(fragIndex, shortAlleleSize);
        double fragLoglikForLongAllele
            = isHomozygous ? fragLoglikForShortAllele : fragLogliks_.getLoglik(fragIndex, longAlleleSize);

        const double shortAlleleTerm = shortAllelePrior + fragLoglikForShortAllele;
        const double longAlleleTerm = longAllelePrior + fragLoglikForLongAllele;

        const double loglikGivenRightmap = getLogSum(shortAlleleTerm, longAlleleTerm);
        const double loglikGivenMismap = topFragLogliks_[fragIndex];
        genotypeLoglik += getLogSum(mismapPrior + loglikGivenMismap, rightmapPrior + loglikGivenRightmap);
    }

    return genotypeLoglik;
}

double TwoAlleleGenotyper::getLongAndShortAlleleLoglik(int longAlleleSize, int shortAlleleSize)
{
    return getShortAndLongAlleleLoglik(shortAlleleSize, longAlleleSize);
}

}
}
