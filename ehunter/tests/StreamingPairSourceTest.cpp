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
#include "sample/StreamingPairSource.hh"
#include "sample/ParallelPairSource.hh"

#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>
#include "gtest/gtest.h"

#include "core/TaskExecutor.hh"
#include "graphcore/GraphBuilders.hh"
#include "sample/GenomeMask.hh"
#include "sample/HtsFileStreamer.hh"

namespace
{

using namespace ehunter;

// Commits every pair of a batch, in input order, on the publishing thread.
PrepareStreamingBatch publishEach(std::function<void(StreamingInputPair)> publish)
{
    return [publish](const std::vector<StreamingInputPair*>& pairs) -> std::function<void()> {
        return [publish, &pairs] {
            for (auto* pair : pairs)
                if (pair)
                    publish(std::move(*pair));
        };
    };
}

class StreamingPairSourceTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        directory_ = boost::filesystem::temp_directory_path()
            / boost::filesystem::unique_path("eh-pair-source-%%%%-%%%%-%%%%");
        boost::filesystem::create_directory(directory_);
        inputPath_ = directory_ / "reads.sam";
        referencePath_ = directory_ / "reference.fa";
        std::ofstream reference(referencePath_.string());
        reference << ">chr1\nACGT\n>chr2\nTGCA\n";
        output_.open(inputPath_.string());
        ASSERT_TRUE(output_.good());
        output_ << "@HD\tVN:1.6\tSO:unknown\n"
                   "@SQ\tSN:chr1\tLN:100000\n"
                   "@SQ\tSN:chr2\tLN:100000\n";

        LocusSpecification spec(
            "target", ChromType::kAutosome, { GenomicRegion(0, 100, 200) },
            graphtools::makeDeletionGraph("ACGTACGT", "C", "TGCATGCA"), {}, GenotyperParameters());
        HeuristicParameters params(1000, 20, false, graphtools::AlignerType::DAG_ALIGNER, 2, 1, 1, 2, 1);
        analyzers_.emplace_back(
            new locus::LocusAnalyzer(spec, params, std::make_shared<graphtools::BlankAlignmentWriter>()));
        finder_.reset(new AnalyzerFinder(analyzers_));
        mask_.addRegion(0, 100, 200);
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
        const auto chromosome = [](int32_t value) {
            return value < 0 ? std::string("*") : "chr" + std::to_string(value + 1);
        };
        const auto cigar = (contig < 0 || sequence == "*")
            ? std::string("*") : std::to_string(sequence.size()) + "M";
        output_ << name << '\t' << flags << '\t' << chromosome(contig) << '\t'
                << (contig < 0 ? 0 : position + 1) << "\t60\t" << cigar << '\t'
                << chromosome(mateContig) << '\t' << (mateContig < 0 ? 0 : matePosition + 1)
                << "\t0\t" << sequence << '\t' << quality << '\n';
    }

    void addLocal(const std::string& name, int64_t position, const std::string& sequence,
        uint32_t flags = BAM_FPAIRED | BAM_FREAD1)
    {
        add(name, flags, 0, position, 0, 110, sequence);
    }

    std::unique_ptr<htshelpers::HtsFileStreamer> openInput()
    {
        output_.close();
        return std::unique_ptr<htshelpers::HtsFileStreamer>(
            new htshelpers::HtsFileStreamer(inputPath_.string(), referencePath_.string(), 0));
    }

    boost::filesystem::path directory_, inputPath_, referencePath_;
    std::ofstream output_;
    GenomeMask mask_;
    std::vector<std::unique_ptr<locus::LocusAnalyzer>> analyzers_;
    std::unique_ptr<AnalyzerFinder> finder_;
};

TEST_F(StreamingPairSourceTest, ParallelLanesPreservePairOrderAcrossReusedInputBatches)
{
    const unsigned names = 17001;
    for (unsigned round = 0; round < 3; ++round)
    {
        for (unsigned i = 0; i < names; ++i)
            addLocal(std::to_string(i), 110, "AACG");
        for (unsigned i = names; i > 0; --i)
            addLocal(std::to_string(i - 1), 110, "GATT", BAM_FPAIRED | BAM_FREAD2);
    }
    auto input = openInput();
    unsigned completed = 0;
    streamParallelPairs(*input, mask_, *finder_, 4, [] {}, 0, publishEach([&](StreamingInputPair pair) {
        EXPECT_EQ(pair.read.fragmentId(), std::to_string(names - 1 - completed % names));
        EXPECT_EQ(pair.read.fragmentId(), pair.mate.fragmentId());
        EXPECT_EQ(pair.mate.decode().sequence(), "AACG");
        EXPECT_EQ(pair.read.decode().sequence(), "GATT");
        ASSERT_EQ(pair.bundles.size(), 1u);
        EXPECT_EQ(pair.bundles.front().inputType, AnalyzerInputType::kBothReads);
        ++completed;
    }));
    EXPECT_EQ(completed, names * 3);
}

TEST_F(StreamingPairSourceTest, ParallelReaderJoinsOnScanPublicationAndStagingFailures)
{
    for (unsigned i = 0; i < 1024; ++i)
    {
        const auto name = std::to_string(i);
        addLocal(name, 110, "AACG");
        addLocal(name, 110, "GATT", BAM_FPAIRED | BAM_FREAD2);
    }
    const auto discard = [](StreamingInputPair) {};
    const auto scanError = [] { throw std::runtime_error("injected scan failure"); };
    auto input = openInput();
    EXPECT_THROW(streamParallelPairs(*input, mask_, *finder_, 2, scanError, 0, publishEach(discard)), std::runtime_error);
    input = openInput();
    const auto publishError = [](StreamingInputPair) { throw std::runtime_error("injected publication failure"); };
    EXPECT_THROW(streamParallelPairs(*input, mask_, *finder_, 2, [] {}, 0, publishEach(publishError)), std::runtime_error);
    input = openInput();
    const auto stageError = [](const std::vector<StreamingInputPair*>&) -> std::function<void()> {
        throw std::runtime_error("injected staging failure");
    };
    EXPECT_THROW(streamParallelPairs(*input, mask_, *finder_, 2, [] {}, 2, stageError), std::runtime_error);
}

TEST_F(StreamingPairSourceTest, IndexedPhysicalRangesPreserveDistantMatesAndSkipSupplementaryReads)
{
    constexpr unsigned count = 30000;
    for (unsigned i = 0; i < count; ++i)
        add(std::to_string(i), BAM_FPAIRED | BAM_FREAD1, 0, 110, 0, 65536, "AACG");
    // Store uncompressed BGZF to force the interval boundary past the scanner's
    // minimum distance, even on this small synthetic input.
    for (unsigned i = 0; i < count; ++i)
        add(std::to_string(i), BAM_FPAIRED | BAM_FSUPPLEMENTARY, 0, 140, 0, 110, "CCCC");
    output_ << "cross-block\t2049\tchr1\t141\t60\t4M99996S\t=\t111\t0\t"
            << std::string(100000, 'C') << '\t' << std::string(100000, 'I') << '\n';
    output_ << "complex\t2049\tchr1\t141\t60\t2M1I2M\t=\t111\t0\tAACGN\tII!II\n";
    for (unsigned i = count; i > 0; --i)
        add(std::to_string(i - 1), BAM_FPAIRED | BAM_FREAD2, 0, 65536, 0, 110, "GATT");
    output_.close();
    const auto file = (directory_ / "reads.bam").string();
    auto* source = sam_open(inputPath_.string().c_str(), "r");
    ASSERT_TRUE(source);
    auto* target = sam_open(file.c_str(), "wb0");
    ASSERT_TRUE(target);
    auto* header = sam_hdr_read(source);
    ASSERT_TRUE(header);
    ASSERT_EQ(sam_hdr_write(target, header), 0);
    auto* record = bam_init1();
    int status;
    while ((status = sam_read1(source, header, record)) >= 0)
        ASSERT_GE(sam_write1(target, header, record), 0);
    EXPECT_EQ(status, -1);
    bam_destroy1(record);
    bam_hdr_destroy(header);
    EXPECT_EQ(sam_close(source), 0);
    ASSERT_EQ(sam_close(target), 0);
    ASSERT_EQ(sam_index_build(file.c_str(), 0), 0);
    {
        htshelpers::HtsFileStreamer ordinary(file, referencePath_.string());
        htshelpers::HtsFileStreamer streaming(file, referencePath_.string());
        htshelpers::HtsFileStreamer::StreamingRecord view;
        while (ordinary.nextPhysicalRecord())
        {
            ASSERT_TRUE(streaming.nextStreamingRecord(view));
            const auto* expected = ordinary.currentAlignment();
            ASSERT_EQ(ordinary.virtualOffset(), streaming.virtualOffset());
            ASSERT_EQ(expected->core.tid, view.tid);
            ASSERT_EQ(expected->core.mtid, view.mtid);
            ASSERT_EQ(expected->core.pos, view.pos);
            ASSERT_EQ(expected->core.mpos, view.mpos);
            ASSERT_EQ(expected->core.flag, view.flag);
            ASSERT_EQ(expected->core.l_qseq, view.length);
            ASSERT_STREQ(bam_get_qname(expected), view.name);
            ASSERT_EQ(0, std::memcmp(bam_get_seq(expected), view.sequence, (view.length + 1) / 2));
            ASSERT_EQ(0, std::memcmp(bam_get_qual(expected), view.quality, view.length));
        }
        ASSERT_FALSE(streaming.nextStreamingRecord(view));
    }
    htshelpers::HtsFileStreamer input(file, referencePath_.string(), 4);
    const auto indexed = input.bamRangeStarts();
    ASSERT_GT(indexed.size(), 1u);
    EXPECT_GT(input.bamRangeStarts(2).size(), indexed.size());
    unsigned completed = 0;
    streamParallelPairs(input, mask_, *finder_, 2, [] {}, 2, publishEach([&](StreamingInputPair pair) {
        EXPECT_EQ(pair.read.fragmentId(), std::to_string(count - 1 - completed));
        EXPECT_EQ(pair.mate.decode().sequence(), "AACG");
        EXPECT_EQ(pair.read.decode().sequence(), "GATT");
        ++completed;
    }));
    EXPECT_EQ(completed, count);
}

TEST_F(StreamingPairSourceTest, EmitsCompletionOrderAndOwnsReturnedReads)
{
    addLocal("a", 110, "AAAA");
    addLocal("b", 111, "CCCC");
    addLocal("b", 112, "GGGG", BAM_FPAIRED | BAM_FREAD2);
    addLocal("a", 113, "TTTT", BAM_FPAIRED | BAM_FREAD2);
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    auto first = source.next();
    ASSERT_TRUE(first);
    EXPECT_EQ(first->read.fragmentId(), "b");
    EXPECT_EQ(first->read.decode().sequence(), "GGGG");
    EXPECT_EQ(first->mate.decode().sequence(), "CCCC");
    EXPECT_EQ(source.pendingReadCount(), 1u);
    auto second = source.next();
    ASSERT_TRUE(second);
    EXPECT_EQ(second->read.fragmentId(), "a");
    EXPECT_EQ(second->read.decode().sequence(), "TTTT");
    EXPECT_EQ(second->mate.decode().sequence(), "AAAA");
    EXPECT_EQ(source.pendingReadCount(), 0u);
    EXPECT_FALSE(source.next());
    EXPECT_FALSE(source.next());
    EXPECT_EQ(first->read.decode().sequence(), "GGGG");
    ASSERT_EQ(first->bundles.size(), 1u);
    EXPECT_EQ(first->bundles.front().locusIndex, 0u);
    EXPECT_EQ(first->bundles.front().inputType, AnalyzerInputType::kBothReads);
}

TEST_F(StreamingPairSourceTest, RepeatedNamesUseConsecutiveAdmittedEncountersWithoutMateFlagMatching)
{
    // All occurrences deliberately have the same mate number.
    addLocal("repeat", 110, "AAAA");
    addLocal("repeat", 111, "CCCC");
    addLocal("repeat", 112, "GGGG");
    addLocal("repeat", 113, "TTTT");
    addLocal("repeat", 114, "NNNN");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    auto first = source.next();
    auto second = source.next();
    ASSERT_TRUE(first);
    ASSERT_TRUE(second);
    EXPECT_EQ(first->read.decode().sequence(), "CCCC");
    EXPECT_EQ(first->mate.decode().sequence(), "AAAA");
    EXPECT_EQ(second->read.decode().sequence(), "TTTT");
    EXPECT_EQ(second->mate.decode().sequence(), "GGGG");
    EXPECT_TRUE(first->read.decode().isFirstMate());
    EXPECT_TRUE(first->mate.decode().isFirstMate());
    EXPECT_FALSE(source.next());
    EXPECT_EQ(source.pendingReadCount(), 0u);
}

TEST_F(StreamingPairSourceTest, FilteredOccurrencesDoNotChangePairing)
{
    addLocal("same", 110, "AAAA");
    addLocal("same", 111, "CCCC", BAM_FREAD1);
    addLocal("same", 112, "GGGG", BAM_FPAIRED | BAM_FSECONDARY);
    addLocal("same", 113, "TTTT", BAM_FPAIRED | BAM_FSUPPLEMENTARY);
    add("same", BAM_FPAIRED, 0, 5000, 0, 6000, "NNNN");
    addLocal("same", 114, "ACGT", BAM_FPAIRED | BAM_FREAD2);
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_EQ(pair->read.decode().sequence(), "ACGT");
    EXPECT_EQ(pair->mate.decode().sequence(), "AAAA");
    EXPECT_FALSE(source.next());
}

TEST_F(StreamingPairSourceTest, RouteMissConsumesBothOccurrencesBeforeNameReuse)
{
    add("reuse", BAM_FPAIRED, 0, 500, 0, 510, "AAAA");
    add("reuse", BAM_FPAIRED, 0, 510, 0, 500, "CCCC");
    addLocal("reuse", 110, "GGGG");
    addLocal("reuse", 120, "TTTT");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_EQ(pair->read.decode().sequence(), "TTTT");
    EXPECT_EQ(pair->mate.decode().sequence(), "GGGG");
    EXPECT_FALSE(source.next());
}

TEST_F(StreamingPairSourceTest, UsesArrivingMateCoordinatesAndStoredMateLength)
{
    // The stored record is far outside the locus. Its own coordinates must not
    // replace the arriving record's mate coordinates during routing.
    add("remote", BAM_FPAIRED | BAM_FREAD1, 0, 5000, 0, 110, "ACGTAC");
    add("remote", BAM_FPAIRED | BAM_FREAD2, 1, 6000, 0, 190, "TG");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_EQ(pair->read.decode().sequence(), "TG");
    EXPECT_EQ(pair->mate.decode().sequence(), "ACGTAC");
    ASSERT_EQ(pair->bundles.size(), 1u);
    EXPECT_EQ(pair->bundles.front().locusIndex, 0u);
    EXPECT_EQ(pair->bundles.front().inputType, AnalyzerInputType::kBothReads);
    EXPECT_FALSE(source.next());
}

TEST_F(StreamingPairSourceTest, StoredLengthControlsMateContainmentAtIntervalBoundary)
{
    // A six-base mate at 195 extends beyond [100,200]; a four-base mate fits.
    // The arriving record is two bases in both pairs.
    add("length", BAM_FPAIRED, 1, 5000, 0, 110, "AAAAAA");
    add("length", BAM_FPAIRED, 1, 5001, 0, 195, "CC");
    add("length", BAM_FPAIRED, 1, 5002, 0, 110, "TTTT");
    add("length", BAM_FPAIRED, 1, 5003, 0, 195, "GG");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_EQ(pair->read.decode().sequence(), "GG");
    EXPECT_EQ(pair->mate.decode().sequence(), "TTTT");
    EXPECT_FALSE(source.next());
}

TEST_F(StreamingPairSourceTest, PreservesQualityCasingMissingQualityAndOrientation)
{
    add("quality", BAM_FPAIRED | BAM_FREAD1 | BAM_FREVERSE, 0, 110, 0, 120, "ACGTN", "456I!");
    add("quality", BAM_FPAIRED | BAM_FREAD2, 0, 120, 0, 110, "TGCA");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_EQ(pair->mate.decode().sequence(), "acGTn");
    EXPECT_EQ(pair->read.decode().sequence(), "TGCA");
    EXPECT_TRUE(pair->mate.decode().isReversed());
    EXPECT_TRUE(pair->mate.decode().isFirstMate());
    EXPECT_FALSE(pair->read.decode().isReversed());
    EXPECT_TRUE(pair->read.decode().isSecondMate());
}

TEST_F(StreamingPairSourceTest, StopsAtFirstPrimaryWithoutContigAndReleasesMissingMates)
{
    addLocal("missing", 110, "AAAA");
    add("unmapped", BAM_FPAIRED | BAM_FUNMAP, -1, -1, 0, 110, "CCCC");
    addLocal("missing", 120, "GGGG");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    EXPECT_FALSE(source.next());
    EXPECT_EQ(source.pendingReadCount(), 0u);
    EXPECT_FALSE(source.next());
}

TEST_F(StreamingPairSourceTest, UnmappedFlagWithCoordinatesDoesNotChangeAdmission)
{
    addLocal("placed", 110, "AAAA", BAM_FPAIRED | BAM_FUNMAP | BAM_FREAD1);
    addLocal("placed", 120, "CCCC", BAM_FPAIRED | BAM_FUNMAP | BAM_FREAD2);
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_EQ(pair->read.decode().sequence(), "CCCC");
    EXPECT_EQ(pair->mate.decode().sequence(), "AAAA");
}

TEST_F(StreamingPairSourceTest, EmptyQueryStillFailsOnRouteMiss)
{
    add("empty", BAM_FPAIRED | BAM_FREAD1, 0, 500, 0, 510, "AAAA");
    add("empty", BAM_FPAIRED | BAM_FREAD2, 0, 510, 0, 500, "*");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    try
    {
        (void)source.next();
        FAIL() << "Empty admitted query must fail even if the completed pair misses all loci";
    }
    catch (const std::logic_error& error)
    {
        EXPECT_EQ(std::string(error.what()), "Encountered empty query for empty/2");
    }
}

TEST_F(StreamingPairSourceTest, EmptyInputIsStableAtEof)
{
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    EXPECT_FALSE(source.next());
    EXPECT_FALSE(source.next());
    EXPECT_EQ(source.pendingReadCount(), 0u);
}

TEST_F(StreamingPairSourceTest, OwnMaskHitShortCircuitsInvalidMateMask)
{
    graphtools::BlobBuilder builder;
    std::vector<PreparedContigIntervals> roots(2);
    roots[0].mask = builder.append(uint64_t(1));
    roots[1].mask = { 4096, 1 };
    const auto rootSlice = builder.append(roots);
    GenomeMask mask;
    mask.usePrepared({ builder.bytes.data(), builder.bytes.size() }, rootSlice);
    add("selected", BAM_FPAIRED, 0, 110, 1, 0, "AAAA");
    add("selected", BAM_FPAIRED, 0, 120, 1, 0, "CCCC");
    auto input = openInput();
    StreamingPairSource source(*input, mask, *finder_);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_EQ(pair->read.decode().sequence(), "CCCC");
}

TEST_F(StreamingPairSourceTest, MaskErrorsPrecedeUnpairedFiltering)
{
    graphtools::BlobBuilder builder;
    std::vector<PreparedContigIntervals> roots(2);
    roots[0].mask = builder.append(uint64_t(1));
    roots[1].mask = { 4096, 1 };
    const auto rootSlice = builder.append(roots);
    GenomeMask mask;
    mask.usePrepared({ builder.bytes.data(), builder.bytes.size() }, rootSlice);
    add("unpaired", 0, 0, 2048, 1, 0, "AAAA");
    auto input = openInput();
    StreamingPairSource source(*input, mask, *finder_);
    EXPECT_THROW(source.next(), std::runtime_error);
}

TEST_F(StreamingPairSourceTest, PackedReadOwnsOddLengthIupacBasesAndRawQualities)
{
    // Include every BAM nibble, quality cutoff boundaries and two vector decode
    // blocks followed by an odd tail. The next record replaces the input buffer.
    const std::string alphabet = "=ACMGRSVTWYHKDBN";
    const std::string qualityAlphabet = "!456I";
    std::string sequence, quality;
    for (size_t index = 0; index != 65; ++index)
    {
        sequence += alphabet[index % alphabet.size()];
        quality += qualityAlphabet[index % qualityAlphabet.size()];
    }
    add("packed-owned", BAM_FPAIRED | BAM_FREAD1 | BAM_FREVERSE, 0, 110, 0, 120, sequence, quality);
    addLocal("overwrite", 130, std::string(257, 'T'));
    auto input = openInput();
    ASSERT_TRUE(input->trySeekingToNextPrimaryAlignment());
    const Read expected = input->decodeRead();
    PackedRead packed(*input->currentAlignment());
    PackedRead copied(packed);
    PackedRead moved(std::move(copied));
    ASSERT_TRUE(input->trySeekingToNextPrimaryAlignment());
    input.reset();

    for (const PackedRead* owned : { &packed, &moved })
    {
        const Read decoded = owned->decode();
        EXPECT_EQ(65u, owned->length());
        EXPECT_EQ("packed-owned", owned->fragmentId());
        EXPECT_EQ(expected.sequence(), decoded.sequence());
        EXPECT_EQ(expected.fragmentId(), decoded.fragmentId());
        EXPECT_EQ(expected.isFirstMate(), decoded.isFirstMate());
        EXPECT_EQ(expected.isReversed(), decoded.isReversed());
        EXPECT_GE(owned->retainedBytes(), sizeof(PackedRead) + 33u + 65u);
    }
}

TEST_F(StreamingPairSourceTest, PackedReadPreservesMissingQualitiesAndSecondMateFlag)
{
    const std::string sequence = std::string(64, 'N') + "R";
    add("missing-qualities", BAM_FPAIRED | BAM_FREAD2, 0, 110, 0, 120, sequence);
    auto input = openInput();
    ASSERT_TRUE(input->trySeekingToNextPrimaryAlignment());
    PackedRead packed(*input->currentAlignment());
    input.reset();
    const Read decoded = packed.decode();
    EXPECT_EQ(sequence, decoded.sequence());
    EXPECT_EQ(sequence.size(), packed.length());
    EXPECT_TRUE(decoded.isSecondMate());
    EXPECT_FALSE(decoded.isReversed());
}

TEST_F(StreamingPairSourceTest, PackedReadPreservesEmptyNameAndQueryValidationOrder)
{
    addLocal("invalid", 110, "ACG");
    auto input = openInput();
    ASSERT_TRUE(input->trySeekingToNextPrimaryAlignment());
    bam1_t invalid = *input->currentAlignment();
    std::vector<uint8_t> bytes(invalid.data, invalid.data + invalid.l_data);
    invalid.data = bytes.data();
    invalid.core.l_qseq = 0;
    try
    {
        PackedRead packed(invalid);
        FAIL() << "Empty query must fail during admission";
    }
    catch (const std::logic_error& error)
    {
        EXPECT_EQ("Encountered empty query for invalid/1", std::string(error.what()));
    }
    bytes.front() = '\0';
    try
    {
        PackedRead packed(invalid);
        FAIL() << "Empty name must take precedence over empty query";
    }
    catch (const std::logic_error& error)
    {
        EXPECT_EQ("Encountered an empty fragment id", std::string(error.what()));
    }
}

TEST_F(StreamingPairSourceTest, EmptyFirstEncounterFailsBeforeEnteringPendingTable)
{
    addLocal("empty-first", 110, "*");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_);
    EXPECT_THROW(source.next(), std::logic_error);
    EXPECT_EQ(0u, source.pendingReadCount());
}

TEST_F(StreamingPairSourceTest, ColdAndRawStorageReleaseExactlyTheSameRead)
{
    const std::string alphabet = "=ACMGRSVTWYHKDBN";
    const std::string qualities = "!456I";
    std::string sequence, quality;
    for (size_t i = 0; i < 257; ++i)
    {
        sequence += alphabet[i % alphabet.size()];
        quality += qualities[i % qualities.size()];
    }
    add("cold-raw-long-name-outside-inline-storage", BAM_FPAIRED | BAM_FREAD2 | BAM_FREVERSE,
        0, 110, 0, 120, sequence, quality);
    auto input = openInput();
    ASSERT_TRUE(input->trySeekingToNextPrimaryAlignment());
    const auto expected = input->decodeRead();
    PackedRead raw(*input->currentAlignment());
    PackedRead cold(*input->currentAlignment(), true);
    input.reset();
    EXPECT_FALSE(raw.isDecoded());
    EXPECT_TRUE(cold.isDecoded());
    EXPECT_GT(raw.rawPayloadBytes(), 0u);
    EXPECT_EQ(cold.rawPayloadBytes(), 0u);
    EXPECT_LE(sizeof(PackedRead), sizeof(Read));
    for (auto* packed : { &raw, &cold })
    {
        EXPECT_EQ(packed->decode().sequence(), expected.sequence());
        const auto released = packed->releaseRead();
        EXPECT_EQ(released.fragmentId(), expected.fragmentId());
        EXPECT_EQ(released.sequence(), expected.sequence());
        EXPECT_EQ(released.isFirstMate(), expected.isFirstMate());
        EXPECT_EQ(released.isReversed(), expected.isReversed());
        EXPECT_EQ(packed->rawPayloadBytes(), 0u);
    }
}

TEST_F(StreamingPairSourceTest, RawCreditReturnsWhenAnEmittedMateLeavesPendingStorage)
{
    addLocal("a", 110, "AAAA");
    add("b", BAM_FPAIRED | BAM_FREAD2 | BAM_FREVERSE, 0, 111, 0, 110, "CCCC", "!56I");
    addLocal("a", 112, "GGGG");
    addLocal("c", 113, "TTTT");
    addLocal("b", 114, "ACGT");
    addLocal("c", 115, "TGCA");
    auto probe = openInput();
    ASSERT_TRUE(probe->trySeekingToNextPrimaryAlignment());
    const size_t oneRaw = PackedRead(*probe->currentAlignment()).rawPayloadBytes();
    probe.reset();
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_, oneRaw);
    auto first = source.next();
    ASSERT_TRUE(first);
    EXPECT_EQ(first->read.fragmentId(), "a");
    EXPECT_FALSE(first->mate.isDecoded());
    EXPECT_EQ(source.pendingRawBytes(), 0u);
    EXPECT_EQ(source.pendingColdReadCount(), 1u);
    auto second = source.next();
    ASSERT_TRUE(second);
    EXPECT_EQ(second->read.fragmentId(), "b");
    EXPECT_TRUE(second->mate.isDecoded());
    const auto coldMate = second->mate.releaseRead();
    EXPECT_EQ(coldMate.sequence(), "ccCC");
    EXPECT_TRUE(coldMate.isSecondMate());
    EXPECT_TRUE(coldMate.isReversed());
    EXPECT_EQ(source.pendingRawBytes(), oneRaw);
    EXPECT_EQ(source.pendingColdReadCount(), 0u);
    auto third = source.next();
    ASSERT_TRUE(third);
    EXPECT_FALSE(third->mate.isDecoded());
    EXPECT_EQ(third->mate.decode().sequence(), "TTTT");
    EXPECT_FALSE(source.next());
    EXPECT_EQ(source.pendingRawBytes(), 0u);
}

TEST_F(StreamingPairSourceTest, RouteMissReturnsRawCreditBeforeAnotherFirstEncounter)
{
    add("miss", BAM_FPAIRED, 0, 500, 0, 510, "AAAA");
    addLocal("cold", 110, "CCCC");
    add("miss", BAM_FPAIRED, 0, 510, 0, 500, "GGGG");
    addLocal("raw-again", 111, "TTTT");
    addLocal("cold", 112, "ACGT");
    addLocal("raw-again", 113, "TGCA");
    auto probe = openInput();
    ASSERT_TRUE(probe->trySeekingToNextPrimaryAlignment());
    const size_t oneRaw = PackedRead(*probe->currentAlignment()).rawPayloadBytes();
    probe.reset();
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_, oneRaw);
    auto first = source.next();
    ASSERT_TRUE(first);
    EXPECT_EQ(first->read.fragmentId(), "cold");
    EXPECT_TRUE(first->mate.isDecoded());
    EXPECT_EQ(source.pendingRawBytes(), oneRaw);
    auto second = source.next();
    ASSERT_TRUE(second);
    EXPECT_EQ(second->read.fragmentId(), "raw-again");
    EXPECT_FALSE(second->mate.isDecoded());
    EXPECT_FALSE(source.next());
    EXPECT_EQ(source.pendingRawBytes(), 0u);
}

TEST_F(StreamingPairSourceTest, ZeroRawBudgetPreservesQualityNameReuseAndMissingMates)
{
    add("repeat", BAM_FPAIRED, 0, 500, 0, 510, "AAAA");
    add("repeat", BAM_FPAIRED, 0, 510, 0, 500, "CCCC");
    add("repeat", BAM_FPAIRED | BAM_FREAD1, 0, 110, 0, 120, "ACGT", "!56I");
    addLocal("repeat", 120, "TGCA");
    addLocal("missing", 130, "NNNN");
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_, 0);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_TRUE(pair->mate.isDecoded());
    EXPECT_EQ(pair->mate.releaseRead().sequence(), "acGT");
    EXPECT_EQ(pair->read.releaseRead().sequence(), "TGCA");
    EXPECT_FALSE(source.next());
    EXPECT_EQ(source.pendingRawBytes(), 0u);
    EXPECT_EQ(source.pendingColdReadCount(), 0u);
}

TEST_F(StreamingPairSourceTest, ManyMissingMatesRespectActualRawCapacityBudget)
{
    for (int index = 0; index < 512; ++index)
        addLocal("orphan-" + std::to_string(index), 110 + index % 40, std::string(16, 'A'));
    addLocal("orphan-0", 170, "GGGG");
    auto probe = openInput();
    ASSERT_TRUE(probe->trySeekingToNextPrimaryAlignment());
    const size_t oneRaw = PackedRead(*probe->currentAlignment()).rawPayloadBytes();
    probe.reset();
    const size_t budget = 2 * oneRaw + oneRaw / 2;
    auto input = openInput();
    StreamingPairSource source(*input, mask_, *finder_, budget);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_EQ(pair->mate.releaseRead().sequence(), std::string(16, 'A'));
    EXPECT_EQ(source.pendingReadCount(), 511u);
    EXPECT_EQ(source.pendingRawBytes(), oneRaw);
    EXPECT_FALSE(source.next());
    EXPECT_EQ(source.pendingReadCount(), 0u);
    EXPECT_EQ(source.pendingRawBytes(), 0u);
    EXPECT_EQ(source.pendingColdReadCount(), 0u);
}

TEST_F(StreamingPairSourceTest, RawCapacityBeyondLogicalPayloadUsesColdStorage)
{
    add("capacity", BAM_FPAIRED | BAM_FREAD1, 0, 110, 0, 130, std::string(16, 'A'), std::string(16, '5'));
    addLocal("capacity", 130, "GGGG");
    auto input = openInput();
    // The payload has 24 logical bytes, but its terminator alone requires at
    // least one more retained byte, independently of allocator growth policy.
    StreamingPairSource source(*input, mask_, *finder_, 24);
    auto pair = source.next();
    ASSERT_TRUE(pair);
    EXPECT_TRUE(pair->mate.isDecoded());
    EXPECT_EQ(pair->mate.releaseRead().sequence(), std::string(16, 'a'));
    EXPECT_FALSE(source.next());
}

TEST_F(StreamingPairSourceTest, CancellationPollingPersistsAcrossPairsAndStopsFilteredScan)
{
    for (int pair = 0; pair < 3; ++pair)
    {
        const auto name = "ready-" + std::to_string(pair);
        addLocal(name, 110, "AAAA");
        addLocal(name, 120, "CCCC", BAM_FPAIRED | BAM_FREAD2);
    }
    for (int record = 0; record < 2050; ++record)
        add("filtered-" + std::to_string(record), BAM_FPAIRED, 0, 5000 + record, 0, 8000, "GGGG");
    auto input = openInput();
    TaskExecutor executor(1);
    unsigned polls = 0;
    StreamingPairSource source(*input, mask_, *finder_, StreamingPairSource::defaultPendingRawBytes,
        [&] {
            ++polls;
            if (executor.isCancelled())
                executor.rethrowIfFailed();
        });
    for (int pair = 0; pair < 3; ++pair)
    {
        const auto result = source.next();
        ASSERT_TRUE(result);
        EXPECT_EQ(result->read.fragmentId(), "ready-" + std::to_string(pair));
        EXPECT_EQ(polls, 0u); // No per-pair callback.
    }
    executor.cancel(std::make_exception_ptr(std::runtime_error("recorded worker failure")));
    try
    {
        (void)source.next();
        FAIL() << "Filtered input must observe the recorded worker failure before EOF";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_EQ(std::string(error.what()), "recorded worker failure");
    }
    EXPECT_EQ(polls, 1u);
    // Six records were consumed by earlier next() calls; the 1024th primary
    // record is therefore filtered record 1017, before its mask is examined.
    EXPECT_EQ(input->currentReadPosition(), 6017);
    EXPECT_EQ(source.pendingReadCount(), 0u);
}

TEST_F(StreamingPairSourceTest, CancellationPollingStopsUnmatchedGrowthBeforeNextAdmission)
{
    for (int record = 0; record < 2050; ++record)
        addLocal("orphan-" + std::to_string(record), 110, "AAAA");
    auto input = openInput();
    unsigned polls = 0;
    StreamingPairSource source(*input, mask_, *finder_, StreamingPairSource::defaultPendingRawBytes,
        [&] {
            ++polls;
            throw std::runtime_error("cancelled pending scan");
        });
    EXPECT_THROW(source.next(), std::runtime_error);
    EXPECT_EQ(polls, 1u);
    EXPECT_EQ(source.pendingReadCount(), 1023u);
}

}
