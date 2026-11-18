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
#include <fstream>
#include "io/BamletWriter.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphBuilders.hh"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <future>
#include <random>
#include <unistd.h>

namespace
{
struct TemporaryBam
{
    TemporaryBam()
    {
        char name[] = "/tmp/eh-writer-XXXXXX";
        const int fd = mkstemp(name);
        if (fd < 0)
            throw std::runtime_error("Cannot create test BAM");
        close(fd);
        path = name;
    }
    ~TemporaryBam() { std::remove(path.c_str()); }
    std::string path;
};

struct WriterFixture
{
    WriterFixture()
        : graph(graphtools::makeDeletionGraph("ACGTAC", "C", "TGTA"))
        , contigs({ { "chr1", 1000 } })
    {
        graph.setGraphId("test-locus");
        catalog.push_back(ehunter::LocusSpecification(
            "test-locus", ehunter::ChromType::kAutosome, { ehunter::GenomicRegion(0, 0, 100) }, graph,
            { { 0, ehunter::GenomicRegion(0, 10, 16) } },
            ehunter::GenotyperParameters()));
    }

    void writeRecords(const std::string& path, std::size_t budget, int compressionThreads)
    {
        ehunter::BamletWriter writer(path, contigs, catalog, compressionThreads, budget);
        const auto alignment = graphtools::decodeGraphAlignment(0, "0[5M]", &graph);
        std::vector<std::future<void>> producers;
        for (int producer = 0; producer < 4; ++producer)
        {
            producers.emplace_back(std::async(std::launch::async, [&, producer] {
                for (int record = 0; record < 500; ++record)
                {
                    // Repeated names deliberately exercise duplicate multiplicity.
                    const auto name = "read-" + std::to_string(producer) + "-" + std::to_string(record % 10);
                    writer.write("test-locus", name, "aCGTa", record % 2 == 0, record % 3 == 0,
                        record % 5 == 0, alignment);
                }
            }));
        }
        for (auto& producer : producers)
            producer.get();
        writer.finish();
        writer.finish(); // Finishing twice must not enqueue into a stopped consumer.
    }

    graphtools::Graph graph;
    ehunter::ReferenceContigInfo contigs;
    ehunter::RegionCatalog catalog;
};

std::pair<std::string, std::vector<std::string>> readBam(const std::string& path)
{
    std::unique_ptr<htsFile, decltype(&hts_close)> file(hts_open(path.c_str(), "rb"), hts_close);
    if (!file)
        throw std::runtime_error("Cannot open test BAM");
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(file.get()), bam_hdr_destroy);
    if (!header)
        throw std::runtime_error("Cannot read test BAM header");
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    std::vector<std::string> records;
    int status;
    while ((status = sam_read1(file.get(), header.get(), record.get())) >= 0)
    {
        kstring_t text = { 0, 0, nullptr };
        const int formatted = sam_format1(header.get(), record.get(), &text);
        if (formatted < 0)
        {
            std::free(text.s);
            throw std::runtime_error("Cannot format test BAM record");
        }
        records.emplace_back(text.s, text.l);
        std::free(text.s);
    }
    EXPECT_EQ(status, -1);
    std::sort(records.begin(), records.end());
    return { std::string(header->text, header->l_text), std::move(records) };
}
}

TEST(BamletWriterQueue, PreservesAllFieldsAndDuplicatesWithTinyBudget)
{
    WriterFixture fixture;
    TemporaryBam unbounded, bounded;
    fixture.writeRecords(unbounded.path, 0, 2);
    fixture.writeRecords(bounded.path, 1, 2); // Every record is larger than the budget.
    const auto expected = readBam(unbounded.path);
    ASSERT_EQ(expected.second.size(), 2000u);
    EXPECT_EQ(readBam(bounded.path), expected);
}

TEST(BamletWriterQueue, SynchronousCompressionPreservesOutput)
{
    WriterFixture fixture;
    TemporaryBam threaded, synchronous;
    fixture.writeRecords(threaded.path, 1024, 2);
    fixture.writeRecords(synchronous.path, 1024, 0);
    EXPECT_EQ(readBam(synchronous.path), readBam(threaded.path));
}

TEST(BamletWriterQueue, DrainsWorkerBuffersAndInterleavedSingleRecordsInOrder)
{
    WriterFixture fixture;
    TemporaryBam output;
    const auto graphAlignment = graphtools::decodeGraphAlignment(0, "0[5M]", &fixture.graph);
    std::vector<std::string> names;
    {
        ehunter::BamletWriter writer(output.path, fixture.contigs, fixture.catalog, 3, 1024, true);
        ehunter::BamletWriter::BatchBuffer buffer(&writer);
        auto write = [&](const std::string& name) {
            names.push_back(name);
            writer.write("test-locus", name, std::string(151, 'a'), true, false, true, graphAlignment);
        };
        write("before");
        for (int group = 0; group < 12; ++group)
        {
            ehunter::BamletWriter::Batch scope(&writer, &buffer);
            for (int i = 0; i < 170; ++i)
                write("duplicate-" + std::to_string(i % 31));
            scope.finish();
        }
        buffer.flush();
        write("after");
        writer.finish();
    }
    std::unique_ptr<htsFile, decltype(&hts_close)> file(hts_open(output.path.c_str(), "rb"), hts_close);
    ASSERT_TRUE(file);
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(file.get()), bam_hdr_destroy);
    ASSERT_TRUE(header);
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    for (const auto& name : names)
    {
        ASSERT_GE(sam_read1(file.get(), header.get(), record.get()), 0);
        EXPECT_STREQ(bam_get_qname(record.get()), name.c_str());
        EXPECT_EQ(record->core.pos, 10);
        EXPECT_EQ(record->core.l_qseq, 151);
    }
    EXPECT_EQ(sam_read1(file.get(), header.get(), record.get()), -1);
}

TEST(BamletWriterQueue, PositionalWritesStayOrderedAcrossOrdinaryRecordsAndFinalEof)
{
    WriterFixture fixture;
    TemporaryBam output;
    auto alignment = graphtools::decodeGraphAlignment(0, "0[5M]", &fixture.graph);
    std::string query(151, 'A');
    uint32_t random = 17;
    {
        ehunter::BamletWriter writer(output.path, fixture.contigs, fixture.catalog, 4, 16 * 1024 * 1024, true);
        ehunter::BamletWriter::BatchBuffer buffer(&writer);
        for (int group = 0; group < 3; ++group)
        {
            ehunter::BamletWriter::Batch batch(&writer, &buffer);
            for (int i = 0; i < 24000; ++i)
            {
                for (char& base : query)
                {
                    random = random * 1664525u + 1013904223u;
                    base = "ACGTacgt"[random >> 29];
                }
                writer.write("test-locus", std::to_string(group * 24000 + i), query, true, false, false, alignment);
            }
            batch.finish();
            buffer.flush();
            writer.write("test-locus","boundary",query,true,false,false,alignment);
        }
        writer.finish();
    }
    std::ifstream bytes(output.path, std::ios::binary | std::ios::ate);
    ASSERT_GT(bytes.tellg(), std::streampos(3 * 1024 * 1024));
    std::unique_ptr<htsFile, decltype(&hts_close)> file(hts_open(output.path.c_str(), "rb"), hts_close);
    ASSERT_TRUE(file);
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(file.get()),bam_hdr_destroy);
    ASSERT_TRUE(header);
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(),bam_destroy1);
    for (int group=0; group<3; ++group)
    {
        for (int i=0; i<24000; ++i)
        {
            ASSERT_GE(sam_read1(file.get(),header.get(),record.get()),0);
            ASSERT_EQ(std::string(bam_get_qname(record.get())),std::to_string(group*24000+i));
        }
        ASSERT_GE(sam_read1(file.get(),header.get(),record.get()),0);
        ASSERT_STREQ(bam_get_qname(record.get()),"boundary");
    }
    ASSERT_EQ(sam_read1(file.get(),header.get(),record.get()),-1);
}

TEST(BamletWriterEncoding, PreservesAllAsciiBaseMappingsQualitiesAndOddPadding)
{
    WriterFixture fixture;
    TemporaryBam output;
    {
        ehunter::BamletWriter writer(output.path, fixture.contigs, fixture.catalog, 0);
        const auto alignment = graphtools::decodeGraphAlignment(0, "0[5M]", &fixture.graph);
        for (unsigned length = 0; length <= 130; ++length)
        {
            std::string query;
            for (unsigned index = 0; index < length; ++index)
                query += static_cast<char>(index % 128);
            // A long qname also exercises allocation across rounding boundaries.
            const std::string name = std::string(length, 'q') + std::to_string(length);
            writer.write("test-locus", name, query, true, false, true, alignment);
        }
        writer.finish();
    }
    std::unique_ptr<htsFile, decltype(&hts_close)> file(hts_open(output.path.c_str(), "rb"), hts_close);
    ASSERT_TRUE(file);
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(file.get()), bam_hdr_destroy);
    ASSERT_TRUE(header);
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    for (unsigned length = 0; length <= 130; ++length)
    {
        ASSERT_GE(sam_read1(file.get(), header.get(), record.get()), 0);
        ASSERT_EQ(record->core.l_qseq, length);
        EXPECT_STREQ(bam_get_qname(record.get()), (std::string(length, 'q') + std::to_string(length)).c_str());
        EXPECT_EQ(record->core.tid, 0);
        EXPECT_EQ(record->core.pos, 10);
        EXPECT_EQ(record->core.flag, BAM_FUNMAP | BAM_FPAIRED | BAM_FMUNMAP | BAM_FMREVERSE | BAM_FREAD1);
        for (unsigned index = 0; index < length; ++index)
        {
            const unsigned character = index % 128;
            EXPECT_EQ(bam_seqi(bam_get_seq(record.get()), index), seq_nt16_table[character]);
            EXPECT_EQ(bam_get_qual(record.get())[index], (character >= 'A' && character <= 'Z') ? 40 : 0);
        }
        if (length % 2)
        {
            EXPECT_EQ(bam_get_seq(record.get())[length / 2] & 15, 0);
        }
        const auto* tag = bam_aux_get(record.get(), "XG");
        ASSERT_NE(tag, nullptr);
        EXPECT_STREQ(bam_aux2Z(tag), "test-locus,0,0[5M]");
    }
    EXPECT_EQ(sam_read1(file.get(), header.get(), record.get()), -1);
}

TEST(BamletWriterEncoding, PreservesVectorTailsMixedCaseLettersAndScalarFallback)
{
    WriterFixture fixture;
    TemporaryBam output;
    std::vector<std::string> queries;
    const std::string alphabet = "ACGTNacgtn";
    std::mt19937 random(196527);
    for (int read = 0; read != 200; ++read)
    {
        std::string query(150, 'A');
        for (char& base : query)
            base = alphabet[random() % alphabet.size()];
        queries.push_back(std::move(query));
    }
    for (int character = 0; character != 128; ++character)
    {
        for (int offset : { 0, 15, 16, 31, 32, 127, 128, 149 })
        {
            auto query = queries[character];
            query[offset] = static_cast<char>(character);
            queries.push_back(std::move(query));
        }
    }
    for (unsigned length = 0; length <= 160; ++length)
    {
        std::string query(length, 'A');
        for (unsigned i = 0; i < length; ++i)
            query[i] = (i % 3 == 0 ? 'a' : 'A') + i % 26;
        queries.push_back(std::move(query));
    }
    {
        ehunter::BamletWriter writer(output.path, fixture.contigs, fixture.catalog, 0, 4096);
        const auto alignment = graphtools::decodeGraphAlignment(0, "0[5M]", &fixture.graph);
        for (std::size_t index = 0; index != queries.size(); ++index)
            writer.write("test-locus", "read-" + std::to_string(index), queries[index], true, false, true, alignment);
        writer.finish();
    }
    std::unique_ptr<htsFile, decltype(&hts_close)> file(hts_open(output.path.c_str(), "rb"), hts_close);
    ASSERT_TRUE(file);
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(file.get()), bam_hdr_destroy);
    ASSERT_TRUE(header);
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    for (std::size_t index = 0; index != queries.size(); ++index)
    {
        ASSERT_GE(sam_read1(file.get(), header.get(), record.get()), 0);
        EXPECT_STREQ(bam_get_qname(record.get()), ("read-" + std::to_string(index)).c_str());
        ASSERT_EQ(record->core.l_qseq, queries[index].size());
        EXPECT_EQ(record->core.tid, 0);
        EXPECT_EQ(record->core.pos, 10);
        if (queries[index].size() % 2)
        {
            EXPECT_EQ(bam_get_seq(record.get())[queries[index].size() / 2] & 15, 0);
        }
        for (size_t base = 0; base != queries[index].size(); ++base)
        {
            const unsigned char character = queries[index][base];
            EXPECT_EQ(bam_seqi(bam_get_seq(record.get()), base), seq_nt16_table[character]);
            EXPECT_EQ(bam_get_qual(record.get())[base], isupper(queries[index][base]) ? 40 : 0);
        }
    }
    EXPECT_EQ(sam_read1(file.get(), header.get(), record.get()), -1);
}

TEST(BamletWriterEncoding, PreservesGraphTagForAllOperationsAndRepeatedNodes)
{
    auto graph = graphtools::makeStrGraph(std::string(200, 'A'), "CAG", std::string(200, 'T'));
    graph.setGraphId("test,long:graph-id");
    ehunter::ReferenceContigInfo contigs({ { "chr1", 1000 } });
    ehunter::RegionCatalog catalog;
    catalog.push_back(ehunter::LocusSpecification(
        "test-locus", ehunter::ChromType::kAutosome, { ehunter::GenomicRegion(0, 0, 500) }, graph,
        { { 0, ehunter::GenomicRegion(0, 0, 200) } }, ehunter::GenotyperParameters()));
    std::vector<graphtools::GraphAlignment> alignments;
    alignments.push_back(graphtools::decodeGraphAlignment(10, "0[150M]", &graph));
    alignments.push_back(graphtools::decodeGraphAlignment(10, "0[0S150M0I0M0D0N0X]", &graph));
    alignments.push_back(graphtools::decodeGraphAlignment(10, "0[2S50M3N1X1I2D93M]", &graph));
    std::string repeatCigar = "0[30M]";
    for (int index = 0; index != 40; ++index)
        repeatCigar += "1[3M]";
    alignments.push_back(graphtools::decodeGraphAlignment(170, repeatCigar, &graph));
    alignments.push_back(graphtools::decodeGraphAlignment(0, "2[150M]", &graph));
    TemporaryBam output;
    {
        ehunter::BamletWriter writer(output.path, contigs, catalog, 0);
        for (const auto& alignment : alignments)
        {
            ASSERT_EQ(alignment.queryLength(), 150);
            writer.write("test-locus", "read", std::string(150, 'A'), true, false, true, alignment);
        }
        writer.finish();
    }
    std::unique_ptr<htsFile, decltype(&hts_close)> file(hts_open(output.path.c_str(), "rb"), hts_close);
    ASSERT_TRUE(file);
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(file.get()), bam_hdr_destroy);
    ASSERT_TRUE(header);
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    for (const auto& alignment : alignments)
    {
        ASSERT_GE(sam_read1(file.get(), header.get(), record.get()), 0);
        const auto expected = std::string(graph.graphId()) + "," + std::to_string(alignment.path().startPosition())
            + "," + alignment.generateCigar();
        const auto* tag = bam_aux_get(record.get(), "XG");
        ASSERT_NE(tag, nullptr);
        EXPECT_STREQ(bam_aux2Z(tag), expected.c_str());
    }
    EXPECT_EQ(sam_read1(file.get(), header.get(), record.get()), -1);
}

TEST(BamletWriterEncoding, ProjectsTheFirstMappedNodeAfterUnmappedNodes)
{
    WriterFixture fixture;
    // Only the final node is mapped in this catalog.
    ehunter::RegionCatalog catalog;
    catalog.push_back(ehunter::LocusSpecification(
        "test-locus", ehunter::ChromType::kAutosome, { ehunter::GenomicRegion(0, 0, 100) }, fixture.graph,
        { { 2, ehunter::GenomicRegion(0, 30, 34) } }, ehunter::GenotyperParameters()));
    TemporaryBam output;
    {
        ehunter::BamletWriter writer(output.path, fixture.contigs, catalog, 0);
        const auto alignment = graphtools::decodeGraphAlignment(0, "0[6M]2[4M]", &fixture.graph);
        writer.write("test-locus", "pair", "ACGTACTGTA", true, false, false, alignment);
        writer.finish();
    }
    const auto records = readBam(output.path).second;
    ASSERT_EQ(records.size(), 1u);
    EXPECT_EQ(records.front(), "pair\t77\tchr1\t31\t0\t*\t*\t0\t0\tACGTACTGTA\tIIIIIIIIII\tXG:Z:test-locus,0,0[6M]2[4M]");
}

TEST(BamletWriterEncoding, RetainsValidationOfNodesAfterTheFirstMappedNode)
{
    WriterFixture fixture;
    TemporaryBam output;
    ehunter::BamletWriter writer(output.path, fixture.contigs, fixture.catalog, 0);
    const auto alignment = graphtools::decodeGraphAlignment(0, "0[6M]2[4M]", &fixture.graph);
    // A caller can mutate the backing graph after the Path was constructed.
    // Subpath bounds must be checked even when node 0 is mapped.
    fixture.graph.setNodeSeq(2, "T");
    EXPECT_THROW(writer.write("test-locus", "pair", "ACGTACTGTA", true, false, false, alignment), std::logic_error);
    writer.finish();
}

#ifdef __linux__
TEST(BamletWriterQueue, WriteFailureDrainsBoundedQueueAndPropagates)
{
    WriterFixture fixture;
    EXPECT_THROW(fixture.writeRecords("/dev/full", 1, 0), std::exception);
}
#endif
