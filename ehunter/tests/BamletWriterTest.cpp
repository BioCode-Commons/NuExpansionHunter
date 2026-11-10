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

#ifdef __linux__
TEST(BamletWriterQueue, WriteFailureDrainsBoundedQueueAndPropagates)
{
    WriterFixture fixture;
    EXPECT_THROW(fixture.writeRecords("/dev/full", 1, 0), std::exception);
}
#endif
