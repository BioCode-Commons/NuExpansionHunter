//
// NuExpansionHunter
// Copyright 2026 NuExpansionHunter contributors
//
// This file is part of NuExpansionHunter, a modified version of
// Illumina ExpansionHunter v5.0.0 (Copyright 2016-2019 Illumina, Inc.).
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
#include "verification/reference/io/JsonWriter_v5.hh"
#include "io/JsonWriter.hh"
#include "io/PreparedCatalog.hh"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unistd.h>

#include "graphcore/GraphBuilders.hh"
#include "gtest/gtest.h"

namespace
{
using namespace ehunter;
using Json = nlohmann::json;

using ehunter::reference::writeReference;

struct JsonFixture
{
    SampleParameters parameters { "sample\"\\\t\n\r UTF-8: \xc3\xa9 \xe4\xb8\xad", Sex::kFemale };
    ReferenceContigInfo contigs { { { "chr1", 10000 } } };
    RegionCatalog catalog;
    SampleFindings findings;

    void addLocus(const std::string& id, double coverage)
    {
        catalog.emplace_back(LocusSpecification(
            id, ChromType::kAutosome, { GenomicRegion(0, 100, 200) },
            graphtools::makeDeletionGraph("ACGT", "CAG", "TTTT"), {}, GenotyperParameters()));
        findings.emplace_back(LocusStats(AlleleCount::kTwo, 150, 350, coverage));
    }

    std::string write()
    {
        JsonWriter writer(parameters, contigs, catalog, findings);
        std::ostringstream out;
        // Existing callers may set width before using the stream insertion API.
        out << std::setw(19) << writer;
        EXPECT_EQ(0, out.width());
        return out.str();
    }

    void expectReferenceBytes() { EXPECT_EQ(writeReference(parameters, contigs, catalog, findings), write()); }
};

struct TemporaryFile
{
    TemporaryFile()
    {
        char name[] = "/tmp/eh-json-XXXXXX";
        const int fd = mkstemp(name);
        if (fd < 0)
            throw std::runtime_error("Cannot reserve JSON test file");
        close(fd);
        path = name;
    }
    ~TemporaryFile() { std::remove(path.c_str()); }
    std::string path;
};
}

TEST(JsonWriterStreaming, EmptyCatalogRetainsOnlySampleParameters)
{
    JsonFixture data;
    data.expectReferenceBytes();
    EXPECT_EQ(0u, Json::parse(data.write()).count("LocusResults"));
}

TEST(JsonWriterStreaming, SortsUnescapedIdsAndKeepsLastDuplicate)
{
    JsonFixture data;
    data.addLocus("z", 1.0);
    data.addLocus("duplicate", 2.0);
    data.addLocus("\xe4\xb8\xad", 3.0);
    data.addLocus("a", 4.0);
    data.addLocus("\"quoted\\name\t\n\r", 5.0);
    data.addLocus(std::string("embedded\0null", 13), 6.0);
    data.addLocus("", 7.0);
    data.addLocus("duplicate", 8.0);
    data.addLocus("\xc3\xa9", 9.0);
    data.addLocus("duplicate", 10.0);
    data.expectReferenceBytes();
    const auto document = Json::parse(data.write());
    EXPECT_EQ(8u, document["LocusResults"].size());
    EXPECT_EQ(10.0, document["LocusResults"]["duplicate"]["Coverage"].get<double>());
}

TEST(JsonWriterStreaming, PreservesFloatingPointFormattingAndNonfiniteValues)
{
    JsonFixture data;
    const std::vector<double> values { 0.0,
                                       -0.0,
                                       1.0 / 3.0,
                                       1.0e-20,
                                       1.0e20,
                                       std::numeric_limits<double>::denorm_min(),
                                       std::numeric_limits<double>::max(),
                                       std::numeric_limits<double>::infinity(),
                                       -std::numeric_limits<double>::infinity(),
                                       std::numeric_limits<double>::quiet_NaN() };
    for (size_t i = 0; i < values.size(); ++i)
        data.addLocus(std::to_string(values.size() - i), values[i]);
    data.expectReferenceBytes();
}

TEST(JsonWriterStreaming, PreservesNestedVariantRecordsAndOptionalGenotypes)
{
    JsonFixture data;
    data.addLocus("variants", 37.1234567890123);
    auto locus = data.catalog[0];
    locus.addVariantSpecification(
        "z_repeat", { VariantType::kRepeat, VariantSubtype::kCommonRepeat }, GenomicRegion(0, 110, 119), { 1 },
        boost::none);
    locus.addVariantSpecification(
        "a\t\"repeat", { VariantType::kRepeat, VariantSubtype::kRareRepeat }, GenomicRegion(0, 110, 119), { 1 },
        boost::none);
    locus.addVariantSpecification(
        "small", { VariantType::kSmallVariant, VariantSubtype::kDeletion }, GenomicRegion(0, 110, 113), { 1 }, 1);
    data.catalog = RegionCatalog { locus };

    auto& variants = data.findings[0].findingsForEachVariant;
    RepeatGenotype genotype(3, { 5, 12 });
    genotype.setShortAlleleSizeInUnitsCi(4, 6);
    genotype.setLongAlleleSizeInUnitsCi(10, 14);
    variants.emplace(
        "z_repeat", std::unique_ptr<VariantFindings>(new RepeatFindings(
                        CountTable({ { 5, 7 }, { 12, 8 } }), CountTable(std::map<int32_t, int32_t> { { 8, 3 } }),
                        CountTable(std::map<int32_t, int32_t> { { 9, 2 } }),
                        AlleleCount::kTwo, genotype, GenotypeFilter())));
    variants.emplace(
        "a\t\"repeat", std::unique_ptr<VariantFindings>(new RepeatFindings(
                            {}, {}, {}, AlleleCount::kTwo, boost::none, GenotypeFilter::kLowDepth)));
    variants.emplace(
        "small", std::unique_ptr<VariantFindings>(new SmallVariantFindings(
                     17, 13, { AlleleStatus::kPresent, 12.3456 }, { AlleleStatus::kUncertain, -0.123456 },
                     AlleleCount::kTwo, SmallVariantGenotype(AlleleType::kRef, AlleleType::kAlt), GenotypeFilter())));
    data.expectReferenceBytes();
}

TEST(JsonWriterStreaming, MatchesReferenceForLargeUnsortedCatalog)
{
    JsonFixture data;
    for (size_t i = 0; i < 2048; ++i)
        data.addLocus("locus_" + std::to_string((i * 97) % 1021), double(i) / 7.0);
    data.expectReferenceBytes();
}

TEST(JsonWriterStreaming, PreservesBytesWithMappedCatalogIds)
{
    JsonFixture data;
    data.addLocus("z", 12.5);
    data.addLocus("duplicate", 2.0);
    data.addLocus("a\"\n\xc3\xa9", 0.0);
    data.addLocus("duplicate", 7.25);
    const auto expected = writeReference(data.parameters, data.contigs, data.catalog, data.findings);

    TemporaryFile archive, reference, catalog;
    std::remove(archive.path.c_str());
    std::ofstream(reference.path) << ">chr1\nACGT\n";
    std::ofstream(catalog.path) << "JSON test catalog identity";
    InputPaths inputs("unused.bam", reference.path, catalog.path);
    HeuristicParameters params(1000, 20, true, graphtools::AlignerType::DAG_ALIGNER, 2, 1, 1, 2, 1);
    prepareCatalogCache(archive.path, data.catalog, inputs, params, data.contigs);
    data.catalog = loadPreparedCatalog(archive.path, inputs, params, data.contigs);
    EXPECT_EQ(expected, data.write());
}
