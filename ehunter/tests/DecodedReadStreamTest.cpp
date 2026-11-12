#include "gtest/gtest.h"

#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <boost/filesystem.hpp>

#include "sample/DecodedReadStream.hh"

namespace
{
struct RecordSnapshot
{
    int32_t contig, mateContig;
    int64_t position, matePosition;
    std::string name, sequence;
    bool firstMate, reversed;

    bool operator==(const RecordSnapshot& other) const
    {
        return contig == other.contig && mateContig == other.mateContig
            && position == other.position && matePosition == other.matePosition
            && name == other.name && sequence == other.sequence
            && firstMate == other.firstMate && reversed == other.reversed;
    }
};

RecordSnapshot snapshot(const ehunter::DecodedStreamRecord& record)
{
    return { record.contig, record.mateContig, record.position, record.matePosition,
        record.read.fragmentId(), record.read.sequence(), record.read.isFirstMate(), record.read.isReversed() };
}

struct StreamRun
{
    std::vector<RecordSnapshot> records;
    std::size_t peakCache = 0, largeRecords = 0;
    std::string error;
};

class DecodedReadStreamTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        directory_ = boost::filesystem::temp_directory_path()
            / boost::filesystem::unique_path("nuexpansionhunter-decode-%%%%-%%%%-%%%%");
        boost::filesystem::create_directory(directory_);
        input_ = directory_ / "reads.sam";
        reference_ = directory_ / "reference.fa";
        std::ofstream reference(reference_.string());
        reference << ">chr1\nACGT\n>chr2\nTGCA\n";
        reference.close();
        output_.open(input_.string());
        ASSERT_TRUE(output_.good());
        output_ << "@HD\tVN:1.6\tSO:unknown\n"
                   "@SQ\tSN:chr1\tLN:1000000000\n"
                   "@SQ\tSN:chr2\tLN:1000000000\n";
    }

    void TearDown() override
    {
        output_.close();
        boost::system::error_code error;
        boost::filesystem::remove_all(directory_, error);
    }

    void add(const std::string& name, uint32_t flags, int32_t contig, int64_t position,
        int32_t mateContig, int64_t matePosition, const std::string& sequence,
        const std::string& quality = "*")
    {
        const auto chromosome = [](int32_t value)
        { return value < 0 ? std::string("*") : "chr" + std::to_string(value + 1); };
        const auto cigar = (contig < 0 || sequence == "*")
            ? std::string("*") : std::to_string(sequence.size()) + "M";
        output_ << name << '\t' << flags << '\t' << chromosome(contig) << '\t'
                << (contig < 0 ? 0 : position + 1) << "\t60\t" << cigar << '\t'
                << chromosome(mateContig) << '\t' << (mateContig < 0 ? 0 : matePosition + 1)
                << "\t0\t" << sequence << '\t' << quality << '\n';
    }

    void finish() { output_.close(); }

    static ehunter::GenomeMask broadMask()
    {
        ehunter::GenomeMask mask;
        mask.addRegion(0, 0, 10000000);
        mask.addRegion(1, 0, 10000000);
        return mask;
    }

    StreamRun run(unsigned workers, const ehunter::GenomeMask& mask, bool catchError = false)
    {
        ehunter::htshelpers::HtsFileStreamer source(input_.string(), reference_.string(), 0);
        ehunter::DecodedReadStream stream(source, mask, workers);
        StreamRun result;
        try
        {
            while (auto record = stream.next())
                result.records.push_back(snapshot(*record));
            // Repeated EOF must not restart a batch or duplicate the final item.
            EXPECT_FALSE(stream.next());
            result.peakCache = stream.peakWorkerCacheBytes();
            result.largeRecords = stream.synchronousLargeRecords();
        }
        catch (const std::exception& error)
        {
            if (!catchError)
                throw;
            result.error = error.what();
        }
        return result;
    }

    boost::filesystem::path directory_, input_, reference_;
    std::ofstream output_;
};

TEST_F(DecodedReadStreamTest, PreservesQualityThresholdAmbiguousBasesFlagsAndCoordinates)
{
    add("same-name", BAM_FPAIRED | BAM_FREAD1 | BAM_FREVERSE, 0, 17, 1, 31,
        "ACGTNRY", "!56I!5I");
    add("same-name", BAM_FPAIRED | BAM_FREAD2, 1, 31, 0, 17,
        "=ACMGRSVTWYHKDBN"); // Absent qualities decode as 255, preserving uppercase.
    finish();
    const std::vector<RecordSnapshot> expected = {
        { 0, 1, 17, 31, "same-name", "acGTnrY", true, true },
        { 1, 0, 31, 17, "same-name", "=ACMGRSVTWYHKDBN", false, false }
    };
    for (const unsigned workers : { 0u, 1u, 2u, 4u })
    {
        SCOPED_TRACE(workers);
        const auto actual = run(workers, broadMask());
        EXPECT_EQ(actual.records, expected);
        EXPECT_LE(actual.peakCache, ehunter::DecodedReadStream::workerCacheLimit);
    }
}

TEST_F(DecodedReadStreamTest, PreservesInputOrderAcrossCountBoundariesAndRepeatedNames)
{
    const std::string bases(151, 'A');
    std::vector<RecordSnapshot> expected;
    for (int i = 0; i != 5123; ++i)
    {
        const std::string name = "repeated-" + std::to_string(i % 7);
        const bool first = i % 2 == 0;
        const bool reversed = i % 3 == 0;
        add(name, BAM_FPAIRED | (first ? BAM_FREAD1 : BAM_FREAD2)
                | (reversed ? BAM_FREVERSE : 0),
            i % 2, i * 7, (i + 1) % 2, i * 11, bases);
        expected.push_back({ i % 2, (i + 1) % 2, i * 7, i * 11,
            name, bases, first, reversed });
    }
    finish();
    for (const unsigned workers : { 0u, 1u, 2u, 4u })
    {
        SCOPED_TRACE(workers);
        const auto actual = run(workers, broadMask());
        EXPECT_EQ(actual.records, expected);
        EXPECT_EQ(actual.largeRecords, 0u);
        EXPECT_LE(actual.peakCache, ehunter::DecodedReadStream::workerCacheLimit);
    }
}

TEST_F(DecodedReadStreamTest, PreservesInputOrderAcrossRawByteBoundaries)
{
    const std::string bases(1023, 'C');
    const std::string qualities(1023, '5'); // All exactly at the low-quality cutoff.
    std::vector<RecordSnapshot> expected;
    for (int i = 0; i != 1301; ++i)
    {
        const std::string name = "raw-boundary-" + std::to_string(i);
        add(name, BAM_FPAIRED | BAM_FREAD1, 0, i, 1, i + 9, bases, qualities);
        expected.push_back({ 0, 1, i, i + 9, name, std::string(1023, 'c'), true, false });
    }
    finish();
    for (const unsigned workers : { 1u, 2u, 4u })
    {
        SCOPED_TRACE(workers);
        const auto actual = run(workers, broadMask());
        EXPECT_EQ(actual.records, expected);
        EXPECT_EQ(actual.largeRecords, 0u);
        EXPECT_LE(actual.peakCache, ehunter::DecodedReadStream::workerCacheLimit);
    }
}

TEST_F(DecodedReadStreamTest, AppliesOriginalFiltersAndStopsAtFirstUnmappedPrimary)
{
    ehunter::GenomeMask mask;
    mask.addRegion(0, 0, 100);
    add("selected", BAM_FPAIRED | BAM_FREAD1, 0, 9, 0, 9000, "ACGT");
    add("outside", BAM_FPAIRED | BAM_FREAD1, 0, 5000, 0, 9000, "ACGT");
    add("unpaired", BAM_FREAD1, 0, 10, 0, 10, "ACGT");
    add("secondary", BAM_FPAIRED | BAM_FSECONDARY, 0, 11, 0, 10, "ACGT");
    add("supplementary", BAM_FPAIRED | BAM_FSUPPLEMENTARY, 0, 12, 0, 10, "ACGT");
    add("mate-selected", BAM_FPAIRED | BAM_FREAD2, 0, 5000, 0, 13, "TGCA");
    add("unmapped", BAM_FPAIRED | BAM_FUNMAP, -1, -1, 0, 14, "ACGT");
    add("after-unmapped", BAM_FPAIRED | BAM_FREAD1, 0, 15, 0, 15, "ACGT");
    finish();
    const std::vector<RecordSnapshot> expected = {
        { 0, 0, 9, 9000, "selected", "ACGT", true, false },
        { 0, 0, 5000, 13, "mate-selected", "TGCA", false, false }
    };
    for (const unsigned workers : { 0u, 1u, 2u, 4u })
    {
        SCOPED_TRACE(workers);
        EXPECT_EQ(run(workers, mask).records, expected);
    }
}

TEST_F(DecodedReadStreamTest, EmptyAndFullyFilteredStreamsFinish)
{
    add("outside", BAM_FPAIRED | BAM_FREAD1, 0, 5000, 0, 9000, "ACGT");
    finish();
    const ehunter::GenomeMask emptyMask;
    for (const unsigned workers : { 0u, 1u, 4u })
        EXPECT_TRUE(run(workers, emptyMask).records.empty());

    std::ofstream input(input_.string(), std::ios::trunc);
    input << "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000000\n";
    input.close();
    for (const unsigned workers : { 0u, 1u, 4u })
        EXPECT_TRUE(run(workers, broadMask()).records.empty());
}

TEST_F(DecodedReadStreamTest, LargeRecordFallbackDrainsEarlierBatchesAndResumes)
{
    // Both a first-record and a mid-stream fallback exercise empty/nonempty rings.
    const std::string large(400001, 'G');
    add("large-first", BAM_FPAIRED | BAM_FREAD1, 0, 0, 1, 1, large);
    for (int i = 0; i != 2051; ++i)
        add("before-" + std::to_string(i), BAM_FPAIRED | BAM_FREAD1,
            0, i + 1, 1, i + 2, "ACGT");
    add("large-middle", BAM_FPAIRED | BAM_FREAD2, 0, 3000, 1, 4000, large);
    for (int i = 0; i != 1025; ++i)
        add("after-" + std::to_string(i), BAM_FPAIRED | BAM_FREAD2,
            0, i + 4000, 1, i + 5000, "TGCA");
    add("large-last", BAM_FPAIRED | BAM_FREAD1, 0, 7000, 1, 8000, large);
    finish();
    const auto reference = run(0, broadMask());
    ASSERT_EQ(reference.records.size(), 3079u);
    for (const unsigned workers : { 1u, 2u, 4u })
    {
        SCOPED_TRACE(workers);
        const auto actual = run(workers, broadMask());
        EXPECT_EQ(actual.records, reference.records);
        EXPECT_EQ(actual.largeRecords, 3u);
        EXPECT_LE(actual.peakCache, ehunter::DecodedReadStream::workerCacheLimit);
    }
}

TEST_F(DecodedReadStreamTest, InvalidSequencePreservesEarlierDecodedPrefixAndError)
{
    for (int i = 0; i != 2057; ++i)
        add("valid-" + std::to_string(i), BAM_FPAIRED | BAM_FREAD1, 0, i, 1, i, "ACGT");
    add("empty-sequence", BAM_FPAIRED | BAM_FREAD1, 0, 3000, 1, 3000, "*");
    add("after-error", BAM_FPAIRED | BAM_FREAD1, 0, 3001, 1, 3001, "ACGT");
    finish();
    const auto reference = run(0, broadMask(), true);
    ASSERT_FALSE(reference.error.empty());
    ASSERT_EQ(reference.records.size(), 2057u);
    for (const unsigned workers : { 1u, 2u, 4u })
    {
        SCOPED_TRACE(workers);
        const auto actual = run(workers, broadMask(), true);
        EXPECT_EQ(actual.records, reference.records);
        EXPECT_EQ(actual.error, reference.error);
    }
}

TEST_F(DecodedReadStreamTest, InputParseFailurePreservesEarlierDecodedPrefixAndError)
{
    for (int i = 0; i != 2057; ++i)
        add("valid-" + std::to_string(i), BAM_FPAIRED | BAM_FREAD1, 0, i, 1, i, "ACGT");
    output_ << "not-a-valid-SAM-record\n";
    finish();
    const auto reference = run(0, broadMask(), true);
    ASSERT_FALSE(reference.error.empty());
    ASSERT_EQ(reference.records.size(), 2057u);
    for (const unsigned workers : { 1u, 2u, 4u })
    {
        SCOPED_TRACE(workers);
        const auto actual = run(workers, broadMask(), true);
        EXPECT_EQ(actual.records, reference.records);
        EXPECT_EQ(actual.error, reference.error);
    }
}

TEST_F(DecodedReadStreamTest, EarlyDestructionJoinsOutstandingWorkers)
{
    for (int i = 0; i != 5000; ++i)
        add("valid-" + std::to_string(i), BAM_FPAIRED | BAM_FREAD1,
            0, i, 1, i, std::string(151, 'A'));
    finish();
    const auto mask = broadMask();
    for (const unsigned workers : { 1u, 2u, 4u })
    {
        ehunter::htshelpers::HtsFileStreamer source(input_.string(), reference_.string(), 0);
        ehunter::DecodedReadStream stream(source, mask, workers);
        const auto first = stream.next();
        ASSERT_TRUE(first);
        EXPECT_EQ(first->read.fragmentId(), "valid-0");
        // The remaining submitted batches still own storage when shutdown joins.
    }
}
}
