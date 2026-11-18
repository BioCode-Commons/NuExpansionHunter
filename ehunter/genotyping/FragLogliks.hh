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
// This file differs from ExpansionHunter v5.0.0: fragment likelihoods are cached as vectors indexed by allele size.
// Modifications are licensed under the Apache License, Version 2.0.
//

#pragma once

#include <unordered_map>
#include <utility>
#include <vector>

#include "genotyping/AlignMatrix.hh"

namespace ehunter
{
namespace strgt
{

class FragLogliks
{
public:
    FragLogliks(int motifLen, int readLen, int fragLen, const AlignMatrix* alignMatrixPtr)
        : motifLen_(motifLen)
        , readLen_(readLen)
        , fragLen_(fragLen)
        , alignMatrix_(*alignMatrixPtr)
    {
        assert(alignMatrix_.numReads() % 2 == 0);
    }
    int numFrags() const { return alignMatrix_.numReads() / 2; }
    const std::vector<double>& forAllele(int alleleMotifCount);
    double getLoglik(int fragIndex, int alleleMotifCount) { return forAllele(alleleMotifCount).at(fragIndex); }

private:
    double computeLoglik(const StrAlign& readAlign, const StrAlign& mateAlign, int alleleMotifCount,
                         double logPossibleOrigins) const;

    int motifLen_;
    int readLen_;
    int fragLen_;
    const AlignMatrix& alignMatrix_;
    std::unordered_map<int, std::vector<double>> fragLogliksBySize_;
};

}
}
