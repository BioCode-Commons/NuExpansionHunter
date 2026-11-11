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

#include "genotyping/OneAlleleStrGenotyper.hh"

#include <limits>
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

using LikelihoodFn = double (OneAlleleGenotyper::*)(int);
static Ci getCiAlongX(int& bestX, OneAlleleGenotyper* genotyper, LikelihoodFn likelihood)
{
    int xFrom = bestX;
    int xTo = bestX;

    std::stack<CiAndLoglik> ciCandidates;
    double topGtLoglik = (genotyper->*likelihood)(bestX);
    double totalLoglik = topGtLoglik;
    ciCandidates.emplace(xFrom, xTo, totalLoglik);

    const int kMaxIntervalWidth = 750;
    double likelihoodRatio = 1;
    double llShiftLeft = 0;
    double llShiftRight = 0;
    bool refreshLeft = true;
    bool refreshRight = true;
    while (likelihoodRatio >= 0.01 && xTo - xFrom <= kMaxIntervalWidth)
    {
        // Only the moved boundary has a new neighboring genotype.
        if (refreshLeft)
        {
            llShiftLeft = (genotyper->*likelihood)(xFrom - 1);
            refreshLeft = false;
        }
        if (refreshRight)
        {
            llShiftRight = (genotyper->*likelihood)(xTo + 1);
            refreshRight = false;
        }

        double gtLoglik;
        if (llShiftLeft >= llShiftRight)
        {
            xFrom -= 1;
            gtLoglik = llShiftLeft;
            refreshLeft = true;

            if (gtLoglik > topGtLoglik)
            {
                topGtLoglik = gtLoglik;
                bestX = xFrom;
            }
        }
        else
        {
            xTo += 1;
            gtLoglik = llShiftRight;
            refreshRight = true;

            if (gtLoglik > topGtLoglik)
            {
                topGtLoglik = gtLoglik;
                bestX = xTo;
            }
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

RepeatGenotype OneAlleleGenotyper::genotype(const unordered_set<int>& alleleSizeCandidates)
{
    RepeatGenotype gt = getMostLikelyGenotype(alleleSizeCandidates);
    int bestSize = gt.shortAlleleSizeInUnits();

    Ci strCi = getCiAlongX(bestSize, this, &OneAlleleGenotyper::getAlleleLoglik);
    gt = RepeatGenotype(gt.repeatUnitLen(), { bestSize });
    gt.setShortAlleleSizeInUnitsCi(strCi.begin, strCi.end);

    return gt;
}

RepeatGenotype OneAlleleGenotyper::getMostLikelyGenotype(const unordered_set<int>& alleleSizeCandidates)
{
    double maxGtLoglik = std::numeric_limits<double>::lowest();
    int bestMotifCount = 0;

    for (int motifCount : alleleSizeCandidates)
    {
        const double gtLoglik = getAlleleLoglik(motifCount);

        if (maxGtLoglik < gtLoglik)
        {
            maxGtLoglik = gtLoglik;
            bestMotifCount = motifCount;
        }
    }

    return RepeatGenotype(motifLen_, { bestMotifCount });
}

double OneAlleleGenotyper::getAlleleLoglik(int motifCount)
{
    if (motifCount < 0)
    {
        return std::numeric_limits<double>::lowest();
    }

    double genotypeLoglik = 0;
    const double mismapPrior = std::log(0.001);
    const double rightmapPrior = std::log(1.0 - 0.001);

    for (int fragIndex = 0; fragIndex != fragLogliks_.numFrags(); ++fragIndex)
    {
        const double fragLoglik = fragLogliks_.getLoglik(fragIndex, motifCount);
        const double loglikGivenRightmap = fragLoglik;
        const double loglikGivenMismap = topFragLogliks_[fragIndex];
        genotypeLoglik += getLogSum(mismapPrior + loglikGivenMismap, rightmapPrior + loglikGivenRightmap);
    }

    return genotypeLoglik;
}

}
}
