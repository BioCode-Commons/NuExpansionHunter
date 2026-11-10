//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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

#include <string>

#include "core/WeightedPurityCalculator.hh"

namespace ehunter
{
namespace locus
{

class IrrPairFinder
{
public:
    explicit IrrPairFinder(std::string motif);
    IrrPairFinder(graphtools::TextView motif, graphtools::CatalogRange<graphtools::TextView> permutations)
        : borrowedMotif_(motif)
        , borrowed_(true)
        , purityCalculator_(std::move(permutations))
    {
    }

    graphtools::TextView targetMotif() const { return borrowed_ ? borrowedMotif_ : graphtools::TextView(targetMotif_); }
    bool check(const std::string& read, const std::string& mate) const;

private:
    std::string targetMotif_;
    graphtools::TextView borrowedMotif_;
    bool borrowed_ = false;
    WeightedPurityCalculator purityCalculator_;
    double purityCutoff_ = 0.90;
};

}
}
