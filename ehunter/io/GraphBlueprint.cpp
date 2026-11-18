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
// This file differs from ExpansionHunter v5.0.0: graph tokens are parsed by direct delimiter and alphabet checks.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "io/GraphBlueprint.hh"

#include <algorithm>
#include <sstream>

using graphtools::NodeId;
using std::string;
using std::vector;

namespace ehunter
{

class TokenizationHelper
{
public:
    TokenizationHelper(const string& regex)
        : regex_(regex)
        , currentSymbolIter_(regex_.begin())
    {
    }

    bool reachedEnd() const { return currentSymbolIter_ == regex_.end(); }
    void advance() { ++currentSymbolIter_; }

    char currentSymbol() const { return *currentSymbolIter_; }
    bool pointingAtBase() const { return kBaseSymbols.find(*currentSymbolIter_) != string::npos; }

    bool pointingAtTokenTerminator() const
    {
        // Last character is always a token terminator
        if (currentSymbolIter_ + 1 == regex_.end())
        {
            return true;
        }

        if (isCountQuantifier(*currentSymbolIter_))
        {
            return true;
        }

        char nextSymbol = *(currentSymbolIter_ + 1);
        if (*currentSymbolIter_ == ')' && !isCountQuantifier(nextSymbol))
        {
            return true;
        }

        if (nextSymbol == '(')
        {
            return true;
        }

        return false;
    }

private:
    bool isCountQuantifier(char symbol) const { return kCountQuantifiers.find(symbol) != string::npos; }

    static const string kBaseSymbols;
    static const string kCountQuantifiers;

    const string& regex_;
    string::const_iterator currentSymbolIter_;
};

const string TokenizationHelper::kBaseSymbols("ACGTBDHKMNSRVWY");
const string TokenizationHelper::kCountQuantifiers("*+?");

vector<string> tokenizeRegex(const string& regex)
{
    vector<string> tokens;
    string token;

    TokenizationHelper tokenizationHelper(regex);
    while (!tokenizationHelper.reachedEnd())
    {
        token += tokenizationHelper.currentSymbol();
        if (tokenizationHelper.pointingAtTokenTerminator())
        {
            tokens.push_back(token);
            token.clear();
        }

        tokenizationHelper.advance();
    }

    return tokens;
}

// Contract [K-29]
// Ensures:  the parser accepts exactly the language of v5.0.0's five regular expressions (interruption, skippable
//           and unskippable repeat, insertion or deletion, swap) and yields the same type and sequences.
// Reference: the five expressions, kept in GraphBlueprintTest as the executable check.
FeatureTypeAndSequences TokenParser::parse(const string& token) const
{
    static const char bases[] = "ACGTBDHKMNSRVWY";
    const size_t n = token.size();
    if (n && token.find_first_not_of(bases) == string::npos)
        return { GraphBlueprintFeatureType::kInterruption, { token } };
    if (n >= 4 && token.front() == '(')
    {
        const auto delimiter = token.find_first_not_of(bases, 1);
        if (delimiter == n - 2 && token[delimiter] == ')')
        {
            const auto sequence = token.substr(1, n - 3);
            switch (token.back())
            {
            case '?': return { GraphBlueprintFeatureType::kInsertionOrDeletion, { sequence } };
            case '*': return { GraphBlueprintFeatureType::kSkippableRepeat, { sequence } };
            case '+': return { GraphBlueprintFeatureType::kUnskippableRepeat, { sequence } };
            }
        }
        if (delimiter > 1 && delimiter < n - 2 && token[delimiter] == '|' && token.back() == ')'
            && token.find_first_not_of(bases, delimiter + 1) == n - 1)
            return { GraphBlueprintFeatureType::kSwap,
                     { token.substr(1, delimiter - 1), token.substr(delimiter + 1, n - delimiter - 2) } };
    }
    throw std::logic_error("Could not parse the token " + token);
}

std::ostream& operator<<(std::ostream& out, GraphBlueprintFeatureType tokenType)
{
    switch (tokenType)
    {
    case GraphBlueprintFeatureType::kLeftFlank:
        out << "GraphBlueprintFeatureType::kLeftFlank";
        break;
    case GraphBlueprintFeatureType::kRightFlank:
        out << "GraphBlueprintFeatureType::kRightFlank";
        break;
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
        out << "GraphBlueprintFeatureType::kInsertionOrDeletion";
        break;
    case GraphBlueprintFeatureType::kInterruption:
        out << "GraphBlueprintFeatureType::kInterruption";
        break;
    case GraphBlueprintFeatureType::kSkippableRepeat:
        out << "GraphBlueprintFeatureType::kSkippableRepeat";
        break;
    case GraphBlueprintFeatureType::kUnskippableRepeat:
        out << "GraphBlueprintFeatureType::kUnskippableRepeat";
        break;
    case GraphBlueprintFeatureType::kSwap:
        out << "GraphBlueprintFeatureType::kSwap";
        break;
    default:
        throw std::logic_error("Encountered unknown token type");
    }

    return out;
}

bool isSkippable(GraphBlueprintFeatureType featureType)
{
    switch (featureType)
    {
    case GraphBlueprintFeatureType::kLeftFlank:
    case GraphBlueprintFeatureType::kRightFlank:
    case GraphBlueprintFeatureType::kInterruption:
    case GraphBlueprintFeatureType::kUnskippableRepeat:
    case GraphBlueprintFeatureType::kSwap:
        return false;
    case GraphBlueprintFeatureType::kSkippableRepeat:
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
        return true;
    default:
        throw std::logic_error("Encountered unrecognized feature blueprint type");
    }
}

GraphBlueprint decodeFeaturesFromRegex(const string& regex)
{
    GraphBlueprint blueprint;

    const auto tokens = tokenizeRegex(regex);
    TokenParser parser;

    NodeId firstUnusedNodeId = 0;
    for (int index = 0; index != static_cast<int>(tokens.size()); ++index)
    {
        const auto& token = tokens[index];
        auto featureTypeAndSequences = parser.parse(token);

        auto featureType = featureTypeAndSequences.first;
        const auto& sequences = featureTypeAndSequences.second;

        if (index == 0)
        {
            if (featureType == GraphBlueprintFeatureType::kInterruption)
            {
                featureType = GraphBlueprintFeatureType::kLeftFlank;
            }

            blueprint.push_back(GraphBlueprintFeature(featureType, sequences, { firstUnusedNodeId }));
            ++firstUnusedNodeId;
        }
        else if (index == static_cast<int>(tokens.size()) - 1)
        {
            if (featureType == GraphBlueprintFeatureType::kInterruption)
            {
                featureType = GraphBlueprintFeatureType::kRightFlank;
            }

            blueprint.push_back(GraphBlueprintFeature(featureType, sequences, { firstUnusedNodeId }));
            ++firstUnusedNodeId;
        }
        else
        {
            vector<NodeId> nodeIds;
            for (int nodeIndex = 0; nodeIndex != static_cast<int>(sequences.size()); ++nodeIndex)
            {
                nodeIds.push_back(firstUnusedNodeId);
                ++firstUnusedNodeId;
            }

            blueprint.push_back(GraphBlueprintFeature(featureType, sequences, nodeIds));
        }
    }

    return blueprint;
}

bool doesFeatureDefineVariant(GraphBlueprintFeatureType featureType)
{
    switch (featureType)
    {
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
    case GraphBlueprintFeatureType::kSkippableRepeat:
    case GraphBlueprintFeatureType::kUnskippableRepeat:
    case GraphBlueprintFeatureType::kSwap:
        return true;

    case GraphBlueprintFeatureType::kLeftFlank:
    case GraphBlueprintFeatureType::kRightFlank:
    case GraphBlueprintFeatureType::kInterruption:
        return false;

    default:
        std::stringstream encoding;
        encoding << featureType;
        throw std::logic_error("Unrecognized feature type: " + encoding.str());
    }
}

}
