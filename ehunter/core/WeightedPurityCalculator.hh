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

#pragma once
#include "graphutils/MappedData.hh"
#include <string>
#include <vector>

namespace ehunter
{

class WeightedPurityCalculator
{
public:
    WeightedPurityCalculator(const std::string& repeatUnit);
    explicit WeightedPurityCalculator(graphtools::CatalogRange<graphtools::TextView> prepared)
        : prepared_(std::move(prepared))
    {
        if (prepared_.empty())
            throw std::invalid_argument("Missing prepared repeat permutations");
    }
    const std::vector<std::string>& permutations() const { return repeatUnits_; }
    double score(const std::string& querySequence) const;
    bool scoreAtLeast(const std::string& querySequence, double cutoff) const;

private:
    double score(graphtools::TextView repeatUnit, const std::string& querySequence) const;
    std::vector<std::string> computeCircularPermutations(std::string sequence) const;
    std::vector<std::string> repeatUnits_;
    graphtools::CatalogRange<graphtools::TextView> prepared_;
};

}
