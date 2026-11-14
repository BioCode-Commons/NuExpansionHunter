#include "locus/LocusAnalyzer.hh"

#include <cstdint>
#include <exception>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

namespace
{
using namespace ehunter;
using namespace ehunter::locus;

class RecordingWriter : public graphtools::AlignmentWriter
{
public:
    std::vector<std::string> records;
    int throwOnCall = -1;

    void write(
        graphtools::TextView locusId, const std::string& fragment, const std::string& query,
        bool first, bool reversed, bool mateReversed, const graphtools::GraphAlignment& alignment) override
    {
        if (static_cast<int>(records.size()) == throwOnCall)
            throw std::runtime_error("writer failure");
        std::ostringstream record;
        record << std::string(locusId) << ' ' << fragment << ' ' << query << ' '
               << first << ' ' << reversed << ' ' << mateReversed << ' ' << alignment;
        records.push_back(record.str());
    }
};

HeuristicParameters shortReadParameters()
{
    return { 1000, 20, true, graphtools::AlignerType::DAG_ALIGNER, 4, 0, 0, 4, 1 };
}

LocusSpecification makeSpec(const std::string& structure, bool rare)
{
    auto graph = makeRegionGraph(decodeFeaturesFromRegex(structure));
    GenotyperParameters parameters;
    parameters.minLocusCoverage = 0;
    parameters.minBreakpointSpanningReads = 0;
    LocusSpecification spec(
        "kernel", ChromType::kAutosome, { GenomicRegion(0, 0, 1000) }, graph, {}, parameters);
    spec.addVariantSpecification(
        "repeat", VariantClassification(VariantType::kRepeat,
            rare ? VariantSubtype::kRareRepeat : VariantSubtype::kCommonRepeat),
        GenomicRegion(0, 100, 103), { 1 }, 1);
    return spec;
}

void expectAlignedPairsEqual(const LocusAligner::AlignedPair& expected, const LocusAligner::AlignedPair& observed)
{
    ASSERT_EQ(bool(expected.first), bool(observed.first));
    ASSERT_EQ(bool(expected.second), bool(observed.second));
    if (expected.first)
        EXPECT_EQ(*expected.first, *observed.first);
    if (expected.second)
        EXPECT_EQ(*expected.second, *observed.second);
}

void expectReadEqual(const Read& expected, const Read& observed)
{
    EXPECT_EQ(expected.readId(), observed.readId());
    EXPECT_EQ(expected.sequence(), observed.sequence());
    EXPECT_EQ(expected.isReversed(), observed.isReversed());
}

std::string dna(size_t length, uint32_t state)
{
    std::string sequence;
    for (size_t i = 0; i < length; ++i)
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        sequence += "ACGT"[state & 3];
    }
    return sequence;
}

std::string repeat(size_t copies)
{
    std::string sequence;
    while (copies--)
        sequence += "CAG";
    return sequence;
}

struct ReadTask
{
    Read read;
    boost::optional<Read> mate;
    RegionType region;
};

std::vector<ReadTask> makeTasks(const std::string& left, const std::string& right)
{
    std::vector<ReadTask> tasks;
    for (size_t index = 0; index < 48; ++index)
    {
        const auto haplotype = left + repeat(4 + index % 3) + right;
        std::string first = haplotype.substr(left.size() - 60, 80);
        std::string second = haplotype.substr(left.size() - 20, 80);
        bool hasMate = true;
        RegionType region = RegionType::kTarget;
        switch (index % 12)
        {
        case 0: // Both start on one flank, providing a fragment-length observation.
            first = left.substr(20, 80);
            second = left.substr(120, 80);
            break;
        case 3: // Qualifying costly DP calls mixed with easy reads train workspace hints.
            for (size_t offset = 32; offset < 78; offset += 5)
            {
                first[offset] = first[offset] == 'A' ? 'C' : 'A';
                second[offset] = second[offset] == 'G' ? 'T' : 'G';
            }
            break;
        case 4:
            first.insert(37, "GTA");
            first.erase(58, 3);
            break;
        case 5:
            first = left.substr(80, 100);
            hasMate = false;
            break;
        case 6:
            first = second = std::string(80, 'N');
            break;
        case 7:
            first = left.substr(30, 80);
            second = std::string(80, 'N');
            break;
        case 8: // A long aligned IRR supplies the template for off-target evidence.
            first = repeat(27);
            second = left.substr(40, 81);
            break;
        case 9:
            first = second = repeat(27);
            region = RegionType::kOfftarget;
            break;
        case 10: // No flank support: target rejection must take the IRR fallback.
            first = second = repeat(27);
            break;
        case 11:
            first = second = std::string(80, 'N');
            region = RegionType::kOfftarget;
            break;
        default:
            break;
        }
        const std::string name = "fragment-" + std::to_string(index);
        ReadTask task { Read(ReadId(name, MateNumber::kFirstMate), first, false), {}, region };
        if (hasMate)
            task.mate.emplace(ReadId(name, MateNumber::kSecondMate), second, true);
        if (index % 12 == 2)
        {
            task.read.reverseComplement();
            task.mate->reverseComplement();
        }
        tasks.push_back(std::move(task));
    }
    return tasks;
}

TEST(LocusParallelKernel, AlignmentComputationDefersOutputAndPreservesOrientationArguments)
{
    auto graph = makeRegionGraph(decodeFeaturesFromRegex("ATATTA(C)*GGCGGC"));
    auto writer = std::make_shared<RecordingWriter>();
    LocusAligner aligner("kernel", &graph, shortReadParameters(), writer);
    Read read(ReadId("fragment", MateNumber::kFirstMate), "GGTAAT", true);
    Read mate(ReadId("fragment", MateNumber::kSecondMate), "GCCGCC", false);
    Read composedRead = read;
    Read composedMate = mate;
    graphtools::AlignerSelector selector(graphtools::AlignerType::DAG_ALIGNER);
    const auto result = aligner.alignMates(read, &mate, selector);

    ASSERT_TRUE(result.first && result.second);
    EXPECT_EQ(graphtools::decodeGraphAlignment(2, "0[4M]1[1M]1[1M]", &graph), *result.first);
    EXPECT_EQ(graphtools::decodeGraphAlignment(0, "2[6M]", &graph), *result.second);
    EXPECT_EQ("ATTACC", read.sequence());
    EXPECT_EQ("GGCGGC", mate.sequence());
    EXPECT_FALSE(read.isReversed());
    EXPECT_TRUE(mate.isReversed());
    EXPECT_TRUE(writer->records.empty());

    aligner.writeAlignedPair(read, &mate, result);
    ASSERT_EQ(2u, writer->records.size());
    EXPECT_NE(std::string::npos, writer->records[0].find("fragment ATTACC 1 0 0 "));
    EXPECT_NE(std::string::npos, writer->records[1].find("fragment GGCGGC 0 1 1 "));
    const auto expectedRecords = writer->records;
    writer->records.clear();
    graphtools::AlignerSelector fresh(graphtools::AlignerType::DAG_ALIGNER);
    expectAlignedPairsEqual(result, aligner.align(composedRead, &composedMate, fresh));
    EXPECT_EQ(expectedRecords, writer->records);
}

TEST(LocusParallelKernel, ConcurrentSameLocusComputationAndOrderedReplayMatchDirectProcessing)
{
    // Deliberately prevent a perfect CAG-only read from borrowing flank bases.
    const auto left = dna(317, 17) + "TTT";
    const auto right = std::string("TTT") + dna(317, 923);
    auto spec = makeSpec(left + "(CAG)*" + right, true);
    HeuristicParameters params(1000, 20, true, graphtools::AlignerType::DAG_ALIGNER);
    auto directWriter = std::make_shared<RecordingWriter>();
    auto replayWriter = std::make_shared<RecordingWriter>();
    LocusAnalyzer direct(spec, params, directWriter);
    LocusAnalyzer replay(spec, params, replayWriter);
    auto directTasks = makeTasks(left, right);
    auto replayTasks = directTasks;
    auto freshTasks = directTasks;
    graphtools::AlignerSelector directSelector(params.alignerType());
    for (auto& task : directTasks)
        direct.processMates(task.read, task.mate ? &*task.mate : nullptr, task.region, directSelector);

    std::vector<LocusAnalyzer::AlignedPair> results(replayTasks.size());
    std::vector<std::exception_ptr> failures(4);
    std::vector<std::thread> workers;
    for (size_t worker = 0; worker < failures.size(); ++worker)
    {
        workers.emplace_back([&, worker]() {
            try
            {
                // Construct independently: copying AlignerSelector would share scratch.
                graphtools::AlignerSelector selector(params.alignerType());
                for (size_t offset = worker; offset < replayTasks.size(); offset += failures.size())
                {
                    const size_t index = replayTasks.size() - 1 - offset;
                    auto& task = replayTasks[index];
                    results[index] = replay.computeMates(
                        task.read, task.mate ? &*task.mate : nullptr, task.region, selector);
                }
            }
            catch (...)
            {
                failures[worker] = std::current_exception();
            }
        });
    }
    for (auto& worker : workers)
        worker.join();
    for (const auto& failure : failures)
        if (failure)
            std::rethrow_exception(failure);
    EXPECT_TRUE(replayWriter->records.empty());

    for (size_t index = 0; index < replayTasks.size(); ++index)
    {
        SCOPED_TRACE(index);
        auto& task = replayTasks[index];
        auto& freshTask = freshTasks[index];
        graphtools::AlignerSelector fresh(params.alignerType());
        const auto freshResult = replay.computeMates(
            freshTask.read, freshTask.mate ? &*freshTask.mate : nullptr, freshTask.region, fresh);
        expectAlignedPairsEqual(freshResult, results[index]);
        expectReadEqual(directTasks[index].read, task.read);
        expectReadEqual(freshTask.read, task.read);
        if (task.mate)
        {
            expectReadEqual(*directTasks[index].mate, *task.mate);
            expectReadEqual(*freshTask.mate, *task.mate);
        }
        replay.consumeMates(task.read, task.mate ? &*task.mate : nullptr, task.region, results[index]);
    }
    EXPECT_EQ(directWriter->records, replayWriter->records);
    const auto expected = direct.analyze(Sex::kFemale, boost::none);
    const auto observed = replay.analyze(Sex::kFemale, boost::none);
    EXPECT_EQ(expected.stats, observed.stats);
    ASSERT_GT(expected.stats.medianFragLength(), 0);
    const auto* before = dynamic_cast<const RepeatFindings*>(expected.findingsForEachVariant.at("repeat").get());
    const auto* after = dynamic_cast<const RepeatFindings*>(observed.findingsForEachVariant.at("repeat").get());
    ASSERT_NE(nullptr, before);
    ASSERT_NE(nullptr, after);
    ASSERT_TRUE(before->optionalGenotype());
    EXPECT_GT(before->countsOfInrepeatReads().countOf(27), 4);
    EXPECT_EQ(*before, *after);
    EXPECT_EQ(before->alleleCount(), after->alleleCount());
    EXPECT_EQ(before->genotypeFilter(), after->genotypeFilter());
}

TEST(LocusParallelKernel, OfftargetComputationIsNoopAndErrorsRemainInConsumption)
{
    auto writer = std::make_shared<RecordingWriter>();
    LocusAnalyzer analyzer(makeSpec("ATATTA(C)*GGCGGC", false), shortReadParameters(), writer);
    Read read(ReadId("fragment", MateNumber::kFirstMate), "GGTAAT", true);
    Read mate(ReadId("fragment", MateNumber::kSecondMate), "GCCGCC", false);
    const Read originalRead = read;
    const Read originalMate = mate;
    graphtools::AlignerSelector selector(graphtools::AlignerType::DAG_ALIGNER);
    const auto result = analyzer.computeMates(read, &mate, RegionType::kOfftarget, selector);
    EXPECT_FALSE(result.first);
    EXPECT_FALSE(result.second);
    expectReadEqual(originalRead, read);
    expectReadEqual(originalMate, mate);
    EXPECT_TRUE(writer->records.empty());
    EXPECT_THROW(analyzer.consumeMates(read, &mate, RegionType::kOfftarget, result), std::logic_error);
    EXPECT_NO_THROW(analyzer.consumeMates(read, nullptr, RegionType::kOfftarget, result));
    EXPECT_EQ(0, analyzer.analyze(Sex::kFemale, boost::none).stats.meanReadLength());
}

TEST(LocusParallelKernel, WriterFailurePrecedesEvidenceMutation)
{
    auto writer = std::make_shared<RecordingWriter>();
    writer->throwOnCall = 1;
    LocusAnalyzer analyzer(makeSpec("ATATTA(C)*GGCGGC", false), shortReadParameters(), writer);
    Read read(ReadId("fragment", MateNumber::kFirstMate), "GGTAAT", true);
    Read mate(ReadId("fragment", MateNumber::kSecondMate), "GCCGCC", false);
    graphtools::AlignerSelector selector(graphtools::AlignerType::DAG_ALIGNER);
    const auto result = analyzer.computeMates(read, &mate, RegionType::kTarget, selector);
    ASSERT_TRUE(result.first && result.second);
    EXPECT_TRUE(writer->records.empty());
    EXPECT_THROW(analyzer.consumeMates(read, &mate, RegionType::kTarget, result), std::runtime_error);
    EXPECT_EQ(1u, writer->records.size());
    EXPECT_EQ(0, analyzer.analyze(Sex::kFemale, boost::none).stats.meanReadLength());
}
}
