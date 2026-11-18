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
// This file differs from ExpansionHunter v5.0.0: tests compare direct token parsing with the five reference regular expressions.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "io/GraphBlueprint.hh"

#include <functional>
#include <regex>

#include <boost/optional.hpp>

#include "gtest/gtest.h"

using std::string;
using std::vector;

using namespace ehunter;

TEST(SplittingStringsIntoTokens, ValidStrings_Split)
{
    const string regex = "ATGC(CAG)+GTCG(AAA|TTT)(AGTC)?(CAG)*";
    auto tokens = tokenizeRegex(regex);

    vector<string> expectedTokens = { "ATGC", "(CAG)+", "GTCG", "(AAA|TTT)", "(AGTC)?", "(CAG)*" };
    ASSERT_EQ(expectedTokens, tokens);
}

TEST(ParsingTokens, TypicalTokens_Parsed)
{
    TokenParser parser;
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kInsertionOrDeletion, { "AGTC" } };
        EXPECT_EQ(expectedResult, parser.parse("(AGTC)?"));
    }
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kSkippableRepeat, { "CAG" } };
        EXPECT_EQ(expectedResult, parser.parse("(CAG)*"));
    }
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kUnskippableRepeat, { "CAG" } };
        EXPECT_EQ(expectedResult, parser.parse("(CAG)+"));
    }
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kInterruption, { "GTCG" } };
        EXPECT_EQ(expectedResult, parser.parse("GTCG"));
    }
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kSwap, { "AAA", "TTT" } };
        EXPECT_EQ(expectedResult, parser.parse("(AAA|TTT)"));
    }
}

TEST(ParsingGraphBlueprintTokens, RejectsMalformedDelimitersAndAmbiguityCodes)
{
    TokenParser parser;
    for (const std::string& token : std::vector<std::string>{
             "", "()?", "(|A)", "(A|)", "(A|C|T)", "(A|C)*", "(A)*+", "(A)!", "A:", "a",
             std::string("A\0C", 3), "(ACX)*", "(A|X)", "(R|K)\n" })
        EXPECT_THROW(parser.parse(token), std::logic_error) << token;
    FeatureTypeAndSequences expected = { GraphBlueprintFeatureType::kSwap, { "ACGTBDHKMNSRVWY", "N" } };
    EXPECT_EQ(expected, parser.parse("(ACGTBDHKMNSRVWY|N)"));
}

// Kernel K-29: the token parser accepts exactly the language of ExpansionHunter v5.0.0's five regular
// expressions and returns the same feature type and sequences. Exhaustive over every token of up to five
// symbols from an alphabet that covers each character class of the expressions.
TEST(ParsingTokens, ParserMatchesTheFiveRegularExpressionsOnEveryShortToken)
{
    const std::regex insertionOrDeletion("^\\([ACGTBDHKMNSRVWY]+\\)\\?$");
    const std::regex skippableRepeat("^\\([ACGTBDHKMNSRVWY]+\\)\\*$");
    const std::regex unskippableRepeat("^\\([ACGTBDHKMNSRVWY]+\\)\\+$");
    const std::regex swap("^\\([ACGTBDHKMNSRVWY]+\\|[ACGTBDHKMNSRVWY]+\\)$");
    const std::regex interruption("^[ACGTBDHKMNSRVWY]+$");
    const auto reference = [&](const std::string& token) -> boost::optional<FeatureTypeAndSequences> {
        if (std::regex_match(token, insertionOrDeletion))
            return FeatureTypeAndSequences { GraphBlueprintFeatureType::kInsertionOrDeletion, { token.substr(1, token.size() - 3) } };
        if (std::regex_match(token, skippableRepeat))
            return FeatureTypeAndSequences { GraphBlueprintFeatureType::kSkippableRepeat, { token.substr(1, token.size() - 3) } };
        if (std::regex_match(token, unskippableRepeat))
            return FeatureTypeAndSequences { GraphBlueprintFeatureType::kUnskippableRepeat, { token.substr(1, token.size() - 3) } };
        if (std::regex_match(token, swap))
        {
            const auto sequence = token.substr(1, token.size() - 2);
            const auto pipe = sequence.find('|');
            return FeatureTypeAndSequences { GraphBlueprintFeatureType::kSwap, { sequence.substr(0, pipe), sequence.substr(pipe + 1) } };
        }
        if (std::regex_match(token, interruption))
            return FeatureTypeAndSequences { GraphBlueprintFeatureType::kInterruption, { token } };
        return boost::none;
    };
    const std::string alphabet = "AN()|*+?a\n";
    TokenParser parser;
    std::string token;
    size_t checked = 0;
    std::function<void()> visit = [&]() {
        const auto expected = reference(token);
        if (expected)
            EXPECT_EQ(*expected, parser.parse(token)) << token;
        else
            EXPECT_THROW(parser.parse(token), std::logic_error) << token;
        ++checked;
        if (token.size() == 5)
            return;
        for (const char c : alphabet)
        {
            token.push_back(c);
            visit();
            token.pop_back();
        }
    };
    visit();
    EXPECT_EQ(checked, size_t(1 + 10 + 100 + 1000 + 10000 + 100000));
}
