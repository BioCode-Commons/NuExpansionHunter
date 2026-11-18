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
// This file differs from ExpansionHunter v5.0.0: query collections construct indexes from catalogs and reuse prepared masks.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "sample/GenomeQueryCollection.hh"

using ehunter::locus::LocusAnalyzer;
using std::unique_ptr;
using std::vector;

namespace ehunter
{

namespace
{
void addRegions(GenomeMask& mask, const LocusSpecification& locus)
{
    for (const auto& region : locus.targetReadExtractionRegions())
        mask.addRegion(region.contigIndex(), region.start(), region.end());
    for (const auto& region : locus.offtargetReadExtractionRegions())
        mask.addRegion(region.contigIndex(), region.start(), region.end());
}
}

GenomeQueryCollection::GenomeQueryCollection(const RegionCatalog& catalog, unsigned workers)
    : analyzerFinder(catalog, workers, &targetRegionMask)
{
}

GenomeQueryCollection::GenomeQueryCollection(
    vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers, const RegionCatalog* catalog)
    : analyzerFinder(locusAnalyzers, catalog)
{
    if (catalog && catalog->mappedBytes().data)
        targetRegionMask.usePrepared(catalog->mappedBytes(), catalog->preparedIntervals());
    else
        for (const auto& analyzer : locusAnalyzers)
            addRegions(targetRegionMask, analyzer->locusSpec());
}

}
