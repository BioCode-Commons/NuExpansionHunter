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
// This file differs from ExpansionHunter v5.0.0: VCF records are sorted and serialized in ordered parallel blocks.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "io/VcfWriter.hh"
#include "core/OrderedBlocks.hh"

#include <algorithm>
#include <deque>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/algorithm/string/join.hpp>

#include "core/ReadSupportCalculator.hh"
#include "io/VcfHeader.hh"
#include "io/VcfWriterHelpers.hh"

using boost::optional;
using std::deque;
using std::map;
using std::ostream;
using std::pair;
using std::set;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

static string computeFilterSymbol(GenotypeFilter filter)
{
    vector<string> encoding;
    if (static_cast<bool>(filter & GenotypeFilter::kLowDepth))
    {
        encoding.push_back("LowDepth");
    }

    if (encoding.empty())
    {
        return "PASS";
    }
    else
    {
        return boost::algorithm::join(encoding, ";");
    }
}

void writeBodyHeader(const string& sampleName, ostream& out)
{
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
}

std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter)
{
    outputVcfHeader(vcfWriter.regionCatalog_, vcfWriter.sampleFindings_, out, vcfWriter.threads_, vcfWriter.shared_);
    writeBodyHeader(vcfWriter.sampleId_, out);
    vcfWriter.writeBody(out);
    return out;
}

VcfWriter::VcfWriter(
    std::string sampleId, Reference& reference, const RegionCatalog& regionCatalog,
    const SampleFindings& sampleFindings, unsigned threads, ctpl::thread_pool* shared)
    : sampleId_(std::move(sampleId))
    , reference_(reference)
    , regionCatalog_(regionCatalog)
    , sampleFindings_(sampleFindings)
    , threads_(std::max(1u, threads))
    , shared_(shared && shared->size() ? shared : nullptr)
{
}

// Contract [K-33]
// Ensures:  for every catalog and set of findings, the bytes equal v5.0.0's serial VCF body and VcfHeader::dumpTo.
void VcfWriter::writeBody(ostream& out)
{
    const auto variants = sortedVariants();

    const auto* fasta = dynamic_cast<const FastaReference*>(&reference_);
    const size_t grain = 256;
    const size_t workers = fasta ? std::min<size_t>(threads_, (variants.size() + grain - 1) / grain) : 1;
    std::vector<std::unique_ptr<FastaReference>> handles(shared_ ? static_cast<size_t>(shared_->size()) : workers);
    auto format = [&](int worker, size_t begin, size_t end) {
        if (workers > 1 && !handles[worker])
            handles[worker].reset(new FastaReference(fasta->path(), reference_.contigInfo()));
        Reference& reference = handles[worker] ? *handles[worker] : reference_;
        std::ostringstream buffer;
        for (size_t i = begin; i < end; ++i)
        {
            const auto& variant = variants[i];
            const auto locusSpec = regionCatalog_[variant.locus];
            const auto depth = sampleFindings_[variant.locus].stats.depth();
            VariantVcfWriter variantWriter(reference, locusSpec, depth, variant.specification, buffer);
            variant.findings->accept(&variantWriter);
        }
        return buffer.str();
    };
    orderedBlocks(variants.size(), grain, workers, format,
        [&](const std::string& block) { out.write(block.data(), block.size()); }, shared_);
}

std::vector<VcfWriter::OutputVariant> VcfWriter::sortedVariants() const
{
    std::vector<OutputVariant> variants;
    variants.reserve(sampleFindings_.size());
    auto collect = [&](int, size_t begin, size_t end) {
        std::vector<OutputVariant> block;
        block.reserve(end - begin);
        for (size_t locus = begin; locus < end; ++locus)
        {
            const auto spec = regionCatalog_[locus];
            for (const auto& entry : sampleFindings_[locus].findingsForEachVariant)
                block.push_back({locus, spec.getVariantSpecById(entry.first), entry.second.get()});
        }
        return block;
    };
    orderedBlocks(sampleFindings_.size(), 1024, threads_, collect, [&](std::vector<OutputVariant> block) {
        variants.insert(variants.end(), std::make_move_iterator(block.begin()), std::make_move_iterator(block.end()));
    }, shared_);
    const auto less = [](const OutputVariant& left, const OutputVariant& right) {
        const auto& a = left.specification.referenceLocus();
        const auto& b = right.specification.referenceLocus();
        if (a != b) return a < b;
        if (left.locus != right.locus) return left.locus < right.locus;
        return left.specification.id().compare(right.specification.id()) < 0;
    };
    if (!std::is_sorted(variants.begin(), variants.end(), less))
        std::sort(variants.begin(), variants.end(), less);
    return variants;
}

static string createRepeatAlleleSymbol(int repeatSize) { return "<STR" + std::to_string(repeatSize) + ">"; }

static string computeAltSymbol(const optional<RepeatGenotype>& optionalGenotype, int referenceSizeInUnits)
{
    if (!optionalGenotype)
    {
        return ".";
    }
    const RepeatGenotype& genotype = *optionalGenotype;

    vector<string> alleleEncodings;

    if (genotype.shortAlleleSizeInUnits() != referenceSizeInUnits)
    {
        alleleEncodings.push_back(createRepeatAlleleSymbol(genotype.shortAlleleSizeInUnits()));
    }

    if (genotype.longAlleleSizeInUnits() != referenceSizeInUnits
        && genotype.shortAlleleSizeInUnits() != genotype.longAlleleSizeInUnits())
    {
        alleleEncodings.push_back(createRepeatAlleleSymbol(genotype.longAlleleSizeInUnits()));
    }

    if (alleleEncodings.empty())
    {
        return ".";
    }

    return boost::algorithm::join(alleleEncodings, ",");
}

static string computeInfoFields(const VariantSpecification& variantSpec, const string& repeatUnit)
{
    const auto& referenceLocus = variantSpec.referenceLocus();
    const int referenceSizeInBp = referenceLocus.length();
    const int referenceSizeInUnits = referenceSizeInBp / repeatUnit.length();

    vector<string> fields;
    fields.push_back("END=" + std::to_string(referenceLocus.end()));
    fields.push_back("REF=" + std::to_string(referenceSizeInUnits));
    fields.push_back("RL=" + std::to_string(referenceSizeInBp));
    fields.push_back("RU=" + repeatUnit);
    fields.push_back("VARID=" + variantSpec.id());
    fields.push_back("REPID=" + variantSpec.id());

    return boost::algorithm::join(fields, ";");
}

static ReadType determineSupportType(const CountTable& spanningCounts, const CountTable& flankingCounts, int repeatSize)
{
    if (spanningCounts.countOf(repeatSize) != 0)
    {
        return ReadType::kSpanning;
    }
    else if (flankingCounts.countOf(repeatSize) != 0)
    {
        return ReadType::kFlanking;
    }

    return ReadType::kRepeat;
}

static string computeAlleleFields(
    const VariantSpecification& variantSpec, const string& repeatUnit, const RepeatFindings& repeatFindings)
{
    if (!repeatFindings.optionalGenotype())
    {
        if (repeatFindings.alleleCount() == AlleleCount::kOne)
        {
            // genotype sources alleleSizes confidenceIntervals spanningReadCounts flankingReadCounts repeatReadCounts
            return ".:.:.:.:.:.:.";
        }
        else
        {
            return "./.:./.:./.:./.:./.:./.:.";
        }
    }

    const RepeatGenotype& genotype = *repeatFindings.optionalGenotype();
    const auto referenceLocus = variantSpec.referenceLocus();
    const int referenceSizeInBp = referenceLocus.length();
    const int referenceSizeInUnits = referenceSizeInBp / repeatUnit.length();

    ReadSupportCalculator readSupportCalculator(
        repeatFindings.countsOfSpanningReads(), repeatFindings.countsOfFlankingReads(),
        repeatFindings.countsOfInrepeatReads());

    VcfAlleleFields alleleFields(referenceSizeInUnits);

    const int shortAlleleSize = genotype.shortAlleleSizeInUnits();
    ReadType shortAlleleSupportType = determineSupportType(
        repeatFindings.countsOfSpanningReads(), repeatFindings.countsOfFlankingReads(), shortAlleleSize);

    alleleFields.addAlleleInfo(
        shortAlleleSize, shortAlleleSupportType, genotype.shortAlleleSizeInUnitsCi(),
        readSupportCalculator.getCountOfConsistentSpanningReads(shortAlleleSize),
        readSupportCalculator.getCountOfConsistentFlankingReads(shortAlleleSize),
        readSupportCalculator.getCountOfConsistentRepeatReads(shortAlleleSize));

    if (genotype.numAlleles() == 2)
    {
        const int longAlleleSize = genotype.longAlleleSizeInUnits();
        ReadType longAlleleSupportType = determineSupportType(
            repeatFindings.countsOfSpanningReads(), repeatFindings.countsOfFlankingReads(), longAlleleSize);

        alleleFields.addAlleleInfo(
            genotype.longAlleleSizeInUnits(), longAlleleSupportType, genotype.longAlleleSizeInUnitsCi(),
            readSupportCalculator.getCountOfConsistentSpanningReads(longAlleleSize),
            readSupportCalculator.getCountOfConsistentFlankingReads(longAlleleSize),
            readSupportCalculator.getCountOfConsistentRepeatReads(longAlleleSize));
    }

    return alleleFields.encode();
}

void VariantVcfWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    const auto& referenceLocus = variantSpec_.referenceLocus();
    const auto repeatNodeId = variantSpec_.nodes().front();
    const string& repeatUnit = locusSpec_.regionGraph().nodeSeq(repeatNodeId);
    const int referenceSizeInUnits = referenceLocus.length() / repeatUnit.length();
    const string infoFields = computeInfoFields(variantSpec_, repeatUnit);

    const int posPreceedingRepeat1based = referenceLocus.start();
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    const string leftFlankingBase
        = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

    const string altSymbol = computeAltSymbol(repeatFindingsPtr->optionalGenotype(), referenceSizeInUnits);
    const string alleleFields = computeAlleleFields(variantSpec_, repeatUnit, *repeatFindingsPtr);
    const string sampleFields = alleleFields + ":" + std::to_string(locusDepth_);

    string genotypeFilter = computeFilterSymbol(repeatFindingsPtr->genotypeFilter());

    vector<string> vcfRecordElements = { contigName,
                                         to_string(posPreceedingRepeat1based),
                                         ".",
                                         leftFlankingBase,
                                         altSymbol,
                                         ".",
                                         genotypeFilter,
                                         infoFields,
                                         "GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC",
                                         sampleFields };

    out_ << boost::algorithm::join(vcfRecordElements, "\t") << std::endl;
}

void VariantVcfWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    const auto& referenceLocus = variantSpec_.referenceLocus();
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    string refSequence;
    string altSequence;
    int64_t startPosition = -1;

    if ((variantSpec_.classification().subtype == VariantSubtype::kSwap)
        || (variantSpec_.classification().subtype == VariantSubtype::kSMN))
    {
        assert(variantSpec_.optionalRefNode());
        const auto refNode = *variantSpec_.optionalRefNode();
        const int refNodeIndex = refNode == variantSpec_.nodes().front() ? 0 : 1;
        const int altNodeIndex = refNode == variantSpec_.nodes().front() ? 1 : 0;

        const auto refNodeId = variantSpec_.nodes()[refNodeIndex];
        const auto altNodeId = variantSpec_.nodes()[altNodeIndex];

        refSequence = locusSpec_.regionGraph().nodeSeq(refNodeId);
        altSequence = locusSpec_.regionGraph().nodeSeq(altNodeId);
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start() + 1;
    }
    else if (variantSpec_.classification().subtype == VariantSubtype::kDeletion)
    {
        const string refFlankingBase
            = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

        const int refNodeId = variantSpec_.nodes().front();
        refSequence = refFlankingBase + locusSpec_.regionGraph().nodeSeq(refNodeId);
        altSequence = refFlankingBase;
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start();
    }
    else if (variantSpec_.classification().subtype == VariantSubtype::kInsertion)
    {
        const string refFlankingBase
            = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

        const int altNodeId = variantSpec_.nodes().front();
        refSequence = refFlankingBase;
        altSequence = refFlankingBase + locusSpec_.regionGraph().nodeSeq(altNodeId);
        // Conversion from 0-based to 1-based coordinates
        startPosition = referenceLocus.start();
    }
    else
    {
        std::ostringstream encoding;
        encoding << variantSpec_.classification().type << "/" << variantSpec_.classification().subtype;
        throw std::logic_error("Unable to generate VCF record for " + encoding.str());
    }

    const string infoFields = "VARID=" + variantSpec_.id();

    vector<string> sampleFields;
    vector<string> sampleValues;

    sampleFields.emplace_back("GT");
    const auto& optionalGenotype = smallVariantFindingsPtr->optionalGenotype();
    if (optionalGenotype)
    {
        sampleValues.push_back(streamToString(*optionalGenotype));
    }
    else
    {
        sampleValues.emplace_back(smallVariantFindingsPtr->alleleCount() == AlleleCount::kOne ? "." : "./.");
    }

    sampleFields.emplace_back("AD");
    std::ostringstream adEncoding;
    adEncoding << smallVariantFindingsPtr->numRefReads() << "," << smallVariantFindingsPtr->numAltReads();
    sampleValues.push_back(adEncoding.str());

    if (variantSpec_.classification().subtype == VariantSubtype::kSMN)
    {
        string dst;
        switch (smallVariantFindingsPtr->refAllelePresenceStatus().status)
        {
        case AlleleStatus::kAbsent:
            dst = "+";
            break;
        case AlleleStatus::kPresent:
            dst = "-";
            break;
        case AlleleStatus::kUncertain:
            dst = "?";
            break;
        }
        sampleFields.emplace_back("DST");
        sampleValues.push_back(dst);
        sampleFields.emplace_back("RPL");
        sampleValues.push_back(streamToString(smallVariantFindingsPtr->refAllelePresenceStatus().logLikelihoodRatio));
    }

    sampleFields.emplace_back("LC");
    sampleValues.push_back(std::to_string(locusDepth_));

    const string sampleField = boost::algorithm::join(sampleFields, ":");
    const string sampleValue = boost::algorithm::join(sampleValues, ":");

    string genotypeFilter = computeFilterSymbol(smallVariantFindingsPtr->genotypeFilter());

    vector<string> line { contigName,
                          streamToString(startPosition),
                          ".",
                          refSequence,
                          altSequence,
                          ".",
                          genotypeFilter,
                          infoFields,
                          sampleField,
                          sampleValue };
    out_ << boost::algorithm::join(line, "\t") << std::endl;
}

}
