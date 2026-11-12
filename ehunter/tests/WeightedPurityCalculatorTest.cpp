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

#include "core/WeightedPurityCalculator.hh"

#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <map>
#include <random>

#include "gmock/gmock.h"

using std::string;
using std::vector;
using testing::DoubleNear;

using namespace ehunter;

TEST(CalculatingWeightedPurityScore, PerfectRepeat_Calculated)
{
    const string sequence = "GGCCCCGGCCCC";
    WeightedPurityCalculator wpCalculator("GGCCGG");
    EXPECT_THAT(wpCalculator.score(sequence), DoubleNear(1.0, 0.005));
}

TEST(CalculatingWeightedPurityScore, ImperfectRepeat_Calculated)
{
    WeightedPurityCalculator wpCalculator("AACCCC");
    EXPECT_THAT(wpCalculator.score("ACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA"), DoubleNear(1.0, 0.005));
    EXPECT_THAT(wpCalculator.score("tCCCCttCCCCttCCCCttCCCCtTCCCCttCCCCT"), DoubleNear(0.75, 0.005));
}

namespace
{
double referencePurity(const vector<string>& permutations, const string& query)
{
    const std::map<char, string> expansions = { { 'A', "A" }, { 'C', "C" }, { 'G', "G" }, { 'T', "T" },
        { 'B', "CGT" }, { 'D', "AGT" }, { 'H', "ACT" }, { 'K', "GT" }, { 'M', "AC" }, { 'N', "ACGT" },
        { 'R', "AG" }, { 'S', "CG" }, { 'V', "ACG" }, { 'W', "AT" }, { 'Y', "CT" } };
    vector<double> scores;
    for (const auto& repeat : permutations)
    {
        double score = 0;
        for (size_t position = 0; position < query.size(); ++position)
        {
            const char base = query[position];
            const auto found = expansions.find(repeat[position % repeat.size()]);
            const bool match = found != expansions.end()
                && found->second.find(static_cast<char>(std::toupper(base))) != string::npos;
            score += match ? 1.0 : (string("acgt").find(base) != string::npos ? 0.5 : -1.0);
        }
        scores.push_back(score);
    }
    return *std::max_element(scores.begin(), scores.end()) / static_cast<double>(query.size());
}

graphtools::TextView permutationAt(const void* context, size_t index, const void*, size_t)
{
    return static_cast<const vector<string>*>(context)->at(index);
}
}

TEST(CalculatingWeightedPurityScore, OwnedAndPreparedMatchModuloReference)
{
    std::mt19937 random(80351);
    const string alphabet = "ACGTacgtNX";
    for (const string repeat : { "A", "CG", "CAG", "AACCCC", "BDHKMNRSVWY", "X" })
    {
        const WeightedPurityCalculator owned(repeat);
        const auto& permutations = owned.permutations();
        const WeightedPurityCalculator prepared(graphtools::CatalogRange<graphtools::TextView>(
            &permutations, permutations.size(), permutationAt));
        EXPECT_TRUE(std::isnan(owned.score("")));
        EXPECT_TRUE(std::isnan(prepared.score("")));
        for (int trial = 0; trial < 200; ++trial)
        {
            string query(1 + random() % 200, 'A');
            for (char& base : query)
                base = alphabet[random() % alphabet.size()];
            const double expected = referencePurity(permutations, query);
            EXPECT_EQ(expected, owned.score(query)) << repeat << " " << query;
            EXPECT_EQ(expected, prepared.score(query)) << repeat << " " << query;
        }
    }
}

TEST(CalculatingWeightedPurityScore, CutoffMatchesFullScoreAtBoundaries)
{
    std::mt19937 random(19351);
    const string alphabet = "ACGTacgtNX!?";
    for (const string repeat : { "A", "CG", "CAG", "AACCCC", "BDHKMNRSVWY", "X", "acg" })
    {
        const WeightedPurityCalculator owned(repeat);
        const auto& permutations = owned.permutations();
        const WeightedPurityCalculator prepared(graphtools::CatalogRange<graphtools::TextView>(
            &permutations, permutations.size(), permutationAt));
        for (int trial = 0; trial < 150; ++trial)
        {
            string query(trial % 101, 'A');
            for (char& base : query) base = alphabet[random() % alphabet.size()];
            const double score = owned.score(query);
            for (double cutoff : { -std::numeric_limits<double>::infinity(), -1.0, 0.0, 0.9, 1.0,
                     std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN(), score,
                     std::nextafter(score, -std::numeric_limits<double>::infinity()),
                     std::nextafter(score, std::numeric_limits<double>::infinity()) })
            {
                EXPECT_EQ(score >= cutoff, owned.scoreAtLeast(query, cutoff));
                EXPECT_EQ(prepared.score(query) >= cutoff, prepared.scoreAtLeast(query, cutoff));
            }
        }
    }
}
