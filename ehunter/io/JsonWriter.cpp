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
// This file differs from ExpansionHunter v5.0.0: JSON output serializes ordered locus blocks with optional shared workers.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "io/JsonWriter.hh"
#include "core/OrderedBlocks.hh"

#include <algorithm>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/optional.hpp>

#include "core/Common.hh"
#include "core/ReadSupportCalculator.hh"

namespace ehunter
{

using std::map;
using std::string;
using Json = nlohmann::json;
using boost::optional;
using std::to_string;
using std::vector;

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter)
{
    jsonWriter.write(out);
    return out;
}

JsonWriter::JsonWriter(
    const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
    const SampleFindings& sampleFindings, unsigned threads, ctpl::thread_pool* shared)
    : sampleParams_(sampleParams)
    , contigInfo_(contigInfo)
    , regionCatalog_(regionCatalog)
    , sampleFindings_(sampleFindings)
    , threads_(std::max(1u, threads))
    , shared_(shared)
{
}

// Contract [K-33]
// Ensures:  for every catalog and set of findings, the bytes equal v5.0.0's whole-document writer
//           (verification/reference/io/JsonWriter_v5.hh; tests/JsonWriterTest.cpp).
void JsonWriter::write(std::ostream& out)
{
    Json sampleParametersRecord;
    sampleParametersRecord["SampleId"] = sampleParams_.id();
    sampleParametersRecord["Sex"] = streamToString(sampleParams_.sex());

    // Keep only borrowed IDs and indexes for ordering. The catalog owns these
    // strings for the entire write, including when it is memory mapped.
    using LocusOutputEntry = std::pair<graphtools::TextView, size_t>;
    vector<LocusOutputEntry> loci;
    loci.reserve(sampleFindings_.size());
    for (size_t i = 0; i < sampleFindings_.size(); ++i)
        loci.emplace_back(regionCatalog_[i].locusId(), i);
    std::sort(
        loci.begin(), loci.end(), [](const LocusOutputEntry& left, const LocusOutputEntry& right)
        {
            const int order = left.first.compare(right.first);
            return order < 0 || (order == 0 && left.second > right.second);
        });

    // Match v5.0.0's whole-document bytes for strings, floating point values
    // and nested records by using its serializer and adjusting indentation.
    auto writeRecord = [](std::ostream& out, const Json& record, size_t indentation)
    {
        const string encoded = record.dump(2);
        const string padding(indentation, ' ');
        size_t start = 0;
        for (size_t end = encoded.find('\n'); end != string::npos; end = encoded.find('\n', start))
        {
            out.write(encoded.data() + start, end - start + 1);
            out.write(padding.data(), padding.size());
            start = end + 1;
        }
        out.write(encoded.data() + start, encoded.size() - start);
    };

    out.width(0);
    out << "{\n";
    if (!loci.empty())
        out << "  \"LocusResults\": {\n";
    auto format = [&](int, size_t begin, size_t end) {
        std::ostringstream buffer;
        for (size_t position = begin; position < end; ++position)
        {
            // Match v5.0.0's final-entry rule for duplicate catalog IDs.
            // Sorting equal IDs by descending index selects that entry.
            if (position && loci[position].first == loci[position - 1].first)
                continue;
            const size_t locusIndex = loci[position].second;
            const LocusSpecification& locusSpec = regionCatalog_[locusIndex];
            const LocusFindings& locusFindings = sampleFindings_[locusIndex];
            const std::string& locusId(locusSpec.locusId());

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
                const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

                VariantJsonWriter variantWriter(contigInfo_, locusSpec, variantSpec);
                variantIdAndFindings.second->accept(&variantWriter);
                variantRecords[variantId] = variantWriter.takeRecord();
            }

            if (!variantRecords.empty())
            {
                locusRecord["Variants"] = std::move(variantRecords);
            }
            if (position != 0)
                buffer << ",\n";
            buffer << "    " << Json(locusId).dump() << ": ";
            writeRecord(buffer, locusRecord, 4);
        }
        return buffer.str();
    };
    orderedBlocks(loci.size(), 256, threads_, format,
        [&](const string& block) { out.write(block.data(), block.size()); }, shared_);

    if (!loci.empty())
        out << "\n  },\n";
    out << "  \"SampleParameters\": ";
    writeRecord(out, sampleParametersRecord, 2);
    out << "\n}" << std::endl;
}

static string encodeGenotype(const RepeatGenotype& genotype)
{
    string encoding = std::to_string(genotype.shortAlleleSizeInUnits());

    if (genotype.numAlleles() == 2)
    {
        encoding = encoding + "/" + std::to_string(genotype.longAlleleSizeInUnits());
    }

    return encoding;
}

void VariantJsonWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    assert(variantSpec_.classification().type == VariantType::kRepeat);

    const RepeatFindings& repeatFindings = *repeatFindingsPtr;

    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["ReferenceRegion"] = encode(contigInfo_, variantSpec_.referenceLocus());
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);

    const auto repeatNodeId = variantSpec_.nodes().front();
    const auto& repeatUnit = locusSpec_.regionGraph().nodeSeq(repeatNodeId);
    record_["RepeatUnit"] = repeatUnit;

    record_["CountsOfSpanningReads"] = streamToString(repeatFindings.countsOfSpanningReads());
    record_["CountsOfFlankingReads"] = streamToString(repeatFindings.countsOfFlankingReads());
    record_["CountsOfInrepeatReads"] = streamToString(repeatFindings.countsOfInrepeatReads());

    if (repeatFindings.optionalGenotype())
    {
        record_["Genotype"] = encodeGenotype(*repeatFindings.optionalGenotype());
        record_["GenotypeConfidenceInterval"] = streamToString(*repeatFindings.optionalGenotype());
    }
}

void VariantJsonWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    const SmallVariantFindings& findings = *smallVariantFindingsPtr;
    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);
    record_["ReferenceRegion"] = encode(contigInfo_, variantSpec_.referenceLocus());
    record_["CountOfRefReads"] = findings.numRefReads();
    record_["CountOfAltReads"] = findings.numAltReads();
    record_["StatusOfRefAllele"] = streamToString(findings.refAllelePresenceStatus().status);
    record_["LogLikelihoodRefAllelePresent"] = streamToString(findings.refAllelePresenceStatus().logLikelihoodRatio);
    record_["StatusOfAltAllele"] = streamToString(findings.altAllelePresenceStatus().status);
    record_["LogLikelihoodAltAllelePresent"] = streamToString(findings.altAllelePresenceStatus().logLikelihoodRatio);
    if (findings.optionalGenotype())
    {
        record_["Genotype"] = streamToString(*findings.optionalGenotype());
    }
}

}
