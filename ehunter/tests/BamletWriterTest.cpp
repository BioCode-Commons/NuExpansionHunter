#include "io/BamletWriter.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphBuilders.hh"
#include "gtest/gtest.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <future>
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
            EXPECT_EQ(bam_get_seq(record.get())[length / 2] & 15, 0);
        const auto* tag = bam_aux_get(record.get(), "XG");
        ASSERT_NE(tag, nullptr);
        EXPECT_STREQ(bam_aux2Z(tag), "test-locus,0,0[5M]");
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
    // The old subpath construction rejected this even though node 0 is mapped.
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
