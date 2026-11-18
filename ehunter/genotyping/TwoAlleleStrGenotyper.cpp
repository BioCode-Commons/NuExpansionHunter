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
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: candidate grids support parallel scoring with shared fragment contributions and ordered likelihood reduction.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "genotyping/TwoAlleleStrGenotyper.hh"

#include <cstring>
#include <stack>
#include <unordered_map>

#include "core/LogSum.hh"
#include "core/ThreadPool.hh"

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

RepeatGenotype TwoAlleleGenotyper::genotype(const unordered_set<int>& alleleSizeCandidates, ctpl::thread_pool* gridWorkers)
{
    const size_t n = alleleSizeCandidates.size();
    const bool largeGrid = gridWorkers && n > 1 && double(n) * (n + 1) * fragLogliks_.numFrags() >= 200000;
    RepeatGenotype gt = largeGrid ? scoreGrid(alleleSizeCandidates, *gridWorkers) : getMostLikelyGenotype(alleleSizeCandidates);
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

// A large candidate grid repeatedly evaluates the same fragment likelihoods.
// Two fragments can share a result only when every input to the calculation
// is bit-identical for every candidate. Sum the shared results in the original
// fragment order so likelihood and tie comparisons retain scalar rounding.
struct TwoAlleleGenotyper::FragmentGroups
{
    std::vector<int> representatives;
    std::vector<size_t> groupByFragment;

    FragmentGroups(const std::vector<double>& top, const std::vector<const std::vector<double>*>& columns)
    {
        std::vector<size_t> hashes(top.size(), 0);
        auto hashColumn = [&](const std::vector<double>& column) {
            for (size_t i = 0; i < top.size(); ++i)
            {
                uint64_t bits;
                std::memcpy(&bits, &column[i], sizeof(bits));
                hashes[i] = (hashes[i] ^ static_cast<size_t>(bits ^ (bits >> 32))) * 16777619u;
            }
        };
        hashColumn(top);
        for (const auto* column : columns)
            if (column) hashColumn(*column);
        auto equal = [&](size_t a, size_t b) {
            if (std::memcmp(&top[a], &top[b], sizeof(double))) return false;
            for (const auto* column : columns)
                if (column && std::memcmp(&(*column)[a], &(*column)[b], sizeof(double))) return false;
            return true;
        };
        auto hash = [&](size_t i) { return hashes[i]; };
        std::unordered_map<size_t, size_t, decltype(hash), decltype(equal)> lookup(top.size(), hash, equal);
        groupByFragment.reserve(top.size());
        for (size_t i = 0; i < top.size(); ++i)
        {
            const auto group = lookup.emplace(i, representatives.size());
            if (group.second) representatives.push_back(static_cast<int>(i));
            groupByFragment.push_back(group.first->second);
        }
    }
};

// Contract [K-26]
// Requires: the same compiler and math library, -ffp-contract=off and no -ffast-math.
// Ensures:  every likelihood is the same double as getMostLikelyGenotype's, and the grid returns the first maximum
//           in (i, j) order.
// Reference: getMostLikelyGenotype (the serial search) and the ungrouped likelihoods.
RepeatGenotype TwoAlleleGenotyper::scoreGrid(const unordered_set<int>& candidates, ctpl::thread_pool& workers)
{
    struct Best
    {
        double likelihood = std::numeric_limits<double>::lowest();
        int shortSize = 0, longSize = 0;
    };
    const std::vector<int> sizes(candidates.begin(), candidates.end());
    std::vector<const std::vector<double>*> columns(sizes.size(), nullptr);
    for (size_t i = 0; i < sizes.size(); ++i)
        if (sizes[i] >= 0)
            columns[i] = &fragLogliks_.forAllele(sizes[i]);
    const FragmentGroups groups(topFragLogliks_, columns);
    std::vector<Best> rows(sizes.size());
    std::vector<std::future<void>> work;
    work.reserve(sizes.size());
    try
    {
        for (size_t i = 0; i < sizes.size(); ++i)
            work.push_back(workers.push([&, i](int) {
                if (!columns[i]) return;
                std::vector<double> values(groups.representatives.size());
                for (size_t j = 0; j < sizes.size(); ++j)
                {
                    if (sizes[i] > sizes[j]) continue;
                    const double value = score(sizes[i], sizes[j], *columns[i], *columns[j], &groups, &values);
                    if (rows[i].likelihood < value)
                    {
                        rows[i].likelihood = value;
                        rows[i].shortSize = sizes[i];
                        rows[i].longSize = sizes[j];
                    }
                }
            }));
        for (auto& job : work) job.get();
    }
    catch (...)
    {
        for (auto& job : work)
            if (job.valid()) job.wait();
        throw;
    }
    // Reduce in the scalar loop's original order, including tied likelihoods.
    Best best;
    for (const auto& row : rows)
        if (best.likelihood < row.likelihood)
            best = row;
    return RepeatGenotype(motifLen_, {best.shortSize, best.longSize});
}

double TwoAlleleGenotyper::getShortAndLongAlleleLoglik(int shortAlleleSize, int longAlleleSize)
{
    if (shortAlleleSize < 0 || longAlleleSize < 0 || shortAlleleSize > longAlleleSize)
    {
        return std::numeric_limits<double>::lowest();
    }

    if (!fragLogliks_.numFrags()) return 0;
    const auto& shortValues = fragLogliks_.forAllele(shortAlleleSize);
    const auto& longValues = shortAlleleSize == longAlleleSize ? shortValues : fragLogliks_.forAllele(longAlleleSize);
    return score(shortAlleleSize, longAlleleSize, shortValues, longValues);
}

double TwoAlleleGenotyper::score(int shortAlleleSize, int longAlleleSize, const std::vector<double>& shortValues,
                                 const std::vector<double>& longValues, const FragmentGroups* groups,
                                 std::vector<double>* values) const
{
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
    auto contribution = [&](int fragIndex) {
        const double fragLoglikForShortAllele = shortValues[fragIndex];
        const double fragLoglikForLongAllele = longValues[fragIndex];

        const double shortAlleleTerm = shortAllelePrior + fragLoglikForShortAllele;
        const double longAlleleTerm = longAllelePrior + fragLoglikForLongAllele;

        const double loglikGivenRightmap = getLogSum(shortAlleleTerm, longAlleleTerm);
        const double loglikGivenMismap = topFragLogliks_[fragIndex];
        return getLogSum(mismapPrior + loglikGivenMismap, rightmapPrior + loglikGivenRightmap);
    };
    if (groups && groups->representatives.size() < groups->groupByFragment.size())
    {
        for (size_t i = 0; i < groups->representatives.size(); ++i)
            (*values)[i] = contribution(groups->representatives[i]);
        for (const auto id : groups->groupByFragment)
            genotypeLoglik += (*values)[id];
    }
    else
    {
        for (int i = 0; i < numFrags; ++i)
            genotypeLoglik += contribution(i);
    }

    return genotypeLoglik;
}

double TwoAlleleGenotyper::getLongAndShortAlleleLoglik(int longAlleleSize, int shortAlleleSize)
{
    return getShortAndLongAlleleLoglik(shortAlleleSize, longAlleleSize);
}

}
}
