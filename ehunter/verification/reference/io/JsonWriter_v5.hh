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
// This file differs from ExpansionHunter v5.0.0: re-hosted in namespace `reference` for equivalence checking.
// Modifications are licensed under the Apache License, Version 2.0.
//

#pragma once
#include "io/JsonWriter.hh"
#include <iomanip>
#include <sstream>
namespace ehunter
{
namespace reference
{
using Json = nlohmann::json;
using std::string;
inline std::string writeReference(
    const SampleParameters& sampleParams_, const ReferenceContigInfo& contigInfo_,
    const RegionCatalog& regionCatalog_, const SampleFindings& sampleFindings_)
{
    std::ostringstream out;
    Json sampleParametersRecord;
    sampleParametersRecord["SampleId"] = sampleParams_.id();
    sampleParametersRecord["Sex"] = streamToString(sampleParams_.sex());

    Json resultsRecord;
    const unsigned locusCount(sampleFindings_.size());
    for (unsigned locusIndex(0); locusIndex < locusCount; ++locusIndex)
    {
        const auto locusSpec = regionCatalog_[locusIndex];
        const LocusFindings& locusFindings = sampleFindings_[locusIndex];
        const std::string locusId(locusSpec.locusId());

        Json locusRecord;
        locusRecord["LocusId"] = locusId;
        locusRecord["Coverage"] = locusFindings.stats.depth();
        locusRecord["ReadLength"] = locusFindings.stats.meanReadLength();
        locusRecord["FragmentLength"] = locusFindings.stats.medianFragLength();
        locusRecord["AlleleCount"] = static_cast<int>(locusFindings.stats.alleleCount());

        Json variantRecords;
        for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
        {
            const string& variantId = variantIdAndFindings.first;
            const auto variantSpec = locusSpec.getVariantSpecById(variantId);

            VariantJsonWriter variantWriter(contigInfo_, locusSpec, variantSpec);
            variantIdAndFindings.second->accept(&variantWriter);
            variantRecords[variantId] = variantWriter.record();
        }

        if (!variantRecords.empty())
        {
            locusRecord["Variants"] = variantRecords;
        }
        resultsRecord[locusId] = locusRecord;
    }

    Json sampleRecords;
    if (!resultsRecord.empty())
    {
        sampleRecords["LocusResults"] = resultsRecord;
    }
    sampleRecords["SampleParameters"] = sampleParametersRecord;

    out << std::setw(2) << sampleRecords << std::endl;
    return out.str();
}

}
}
