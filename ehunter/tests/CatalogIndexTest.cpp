#include "core/WeightedPurityCalculator.hh"
#include "graphalign/KmerIndex.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/PathOperations.hh"
#include "graphio/GraphJson.hh"
#include "graphutils/SequenceOperations.hh"
#include "io/PreparedCatalog.hh"
#include "sample/AnalyzerFinder.hh"
#include "sample/GenomeQueryCollection.hh"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <map>
#include <unistd.h>

namespace
{
struct ArchiveFile
{
    ArchiveFile()
    {
        char name[] = "/tmp/eh-kmer-XXXXXX";
        int fd = mkstemp(name);
        if (fd < 0)
            throw std::runtime_error("Cannot reserve test archive");
        close(fd);
        path = name;
        std::remove(path.c_str());
    }
    ~ArchiveFile() { std::remove(path.c_str()); }
    std::string path;
};

std::map<std::string, std::vector<graphtools::Path>> enumerate(const graphtools::Graph& graph, int k)
{
    std::map<std::string, std::vector<graphtools::Path>> expected;
    for (graphtools::NodeId node = 0; node < graph.numNodes(); ++node)
        for (size_t position = 0; position < graph.nodeSeq(node).size(); ++position)
            for (const auto& path :
                 graphtools::extendPathEnd(graphtools::Path(&graph, position, { node }, position), k))
            {
                std::vector<std::string> sequences;
                graphtools::expandReferenceSequence(path.seq(), sequences);
                for (const auto& sequence : sequences)
                    expected[sequence].push_back(path);
            }
    return expected;
}
}

TEST(CompactCatalogIndex, MatchesEnumeratedPathsAndMembership)
{
    const auto graph = graphtools::makeDeletionGraph("ACNACG", "GG", "CATACGT");
    for (int k = 1; k <= 6; ++k)
    {
        const auto expected = enumerate(graph, k);
        const graphtools::KmerIndex index(graph, k), membership(graph, k, true);
        const char* bases = "ACGT";
        for (uint32_t encoded = 0; encoded < (uint32_t(1) << (2 * k)); ++encoded)
        {
            std::string key(k, 'A');
            uint32_t word = encoded;
            for (int j = k - 1; j >= 0; --j)
            {
                key[j] = bases[word & 3];
                word >>= 2;
            }
            const auto it = expected.find(key);
            const bool present = it != expected.end();
            ASSERT_EQ(present, index.contains(key));
            ASSERT_EQ(present, membership.contains(key));
            if (present)
            {
                ASSERT_EQ(it->second, index.getPaths(key));
            }
        }
        EXPECT_FALSE(index.contains("N"));
        EXPECT_EQ(index, graphtools::KmerIndex(index));
    }
}

TEST(CompactCatalogIndex, MappedArchiveRoundTripAndValidation)
{
    auto graph = graphtools::makeDeletionGraph("ACGTACGTACGTACGT", "GG", "CATACGTACGTACGT");
    ArchiveFile file;
    graphtools::KmerIndex::prepareArchive(
        file.path, 2, [&](std::ostream& out, uint64_t slot)
        { graphtools::KmerIndex(graph, slot ? 10 : 14, slot != 0).writeRecord(out); });
    auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
    graphtools::KmerIndex loaded(graph, 14, false, archive, 0);
    graphtools::KmerIndex membership(graph, 10, true, archive, 1);
    EXPECT_EQ(loaded, graphtools::KmerIndex(graph, 14));
    EXPECT_EQ(membership, graphtools::KmerIndex(graph, 10, true));
    EXPECT_THROW(graphtools::KmerIndex(graph, 13, false, archive, 0), std::runtime_error);
    EXPECT_THROW(graphtools::KmerIndex(graph, 14, true, archive, 0), std::runtime_error);
    graph.setNodeSeq(0, "TCGTACGTACGTACGT");
    EXPECT_THROW(graphtools::KmerIndex(graph, 14, false, archive, 0), std::runtime_error);
    EXPECT_THROW(
        graphtools::KmerIndex::prepareArchive(file.path, 0, [](std::ostream&, uint64_t) { }), std::runtime_error);
    archive.reset();
    std::fstream damaged(file.path, std::ios::in | std::ios::out | std::ios::binary);
    damaged.seekg(-1, std::ios::end);
    char byte;
    damaged.get(byte);
    damaged.seekp(-1, std::ios::end);
    damaged.put(byte ^ 1);
    damaged.close();
    auto corrupt = std::make_shared<graphtools::KmerIndexArchive>(file.path);
    EXPECT_THROW(graphtools::KmerIndex(graph, 10, true, corrupt, 1), std::runtime_error);
    ArchiveFile incomplete;
    EXPECT_THROW(
        graphtools::KmerIndex::prepareArchive(
            incomplete.path, 2,
            [&](std::ostream& out, uint64_t slot)
            {
                if (slot == 1)
                    throw std::runtime_error("Simulated preparation failure");
                graphtools::KmerIndex(graph, 10, true).writeRecord(out);
            }),
        std::runtime_error);
    EXPECT_FALSE(std::ifstream(incomplete.path).good());
}

TEST(CatalogCoordinateSweep, PreservesRemoteMatesOfftargetsDuplicatesAndOrder)
{
    using namespace ehunter;
    RegionCatalog catalog;
    std::vector<std::unique_ptr<locus::LocusAnalyzer>> analyzers;
    catalog.reserve(180);
    HeuristicParameters params(1000, 20, true, graphtools::AlignerType::DAG_ALIGNER, 2, 1, 1, 2, 1);
    auto writer = std::make_shared<graphtools::BlankAlignmentWriter>();
    for (int i = 0; i < 180; ++i)
    {
        const int start = (i * 37) % 2000;
        LocusSpecification spec(
            std::to_string(i), ChromType::kAutosome,
            { GenomicRegion(i % 3, start, start + 50 + i % 90), GenomicRegion(i % 3, start, start + 50 + i % 90) },
            graphtools::makeDeletionGraph("ACGT", "C", "ACGT"), {}, GenotyperParameters());
        spec.setOfftargetReadExtractionRegions({ GenomicRegion((i + 1) % 3, 3000 + start, 3150 + start) });
        catalog.push_back(std::move(spec));
    }
    for (const auto& spec : catalog)
        analyzers.emplace_back(new locus::LocusAnalyzer(spec, params, writer));
    AnalyzerFinder finder(analyzers);
    ArchiveFile archive, reference, catalogFile;
    std::ofstream(reference.path) << ">chr1\nACGT\n";
    std::ofstream(catalogFile.path) << "coordinate-oracle-catalog";
    ReferenceContigInfo contigs({ { "chr1", 10000 }, { "chr2", 10000 }, { "chr3", 10000 } });
    InputPaths inputs("unused.bam", reference.path, catalogFile.path);
    prepareCatalogCache(archive.path, catalog, inputs, params, contigs);
    auto mapped = loadPreparedCatalog(archive.path, inputs, params, contigs);
    auto mappedParams = params;
    mappedParams.catalogIndexes = mapped.archive();
    std::vector<std::unique_ptr<locus::LocusAnalyzer>> mappedAnalyzers;
    for (size_t i = 0; i < mapped.size(); ++i)
        mappedAnalyzers.emplace_back(new locus::LocusAnalyzer(mapped[i], mappedParams, writer, i * 2 + 1));
    AnalyzerFinder mappedFinder(mappedAnalyzers, &mapped);
    for (int contig = 0; contig < 3; ++contig)
        for (int start = 0; start < 5500; start += 11)
            for (int mateContig = -1; mateContig < 3; ++mateContig)
            {
                const int mateStart = (start * 113) % 5500;
                const auto expected
                    = finder.query(contig, start, start + start % 120, mateContig, mateStart, mateStart + 30);
                const auto observed
                    = finder.queryStreaming(contig, start, start + start % 120, mateContig, mateStart, mateStart + 30);
                const auto mappedRandom
                    = mappedFinder.query(contig, start, start + start % 120, mateContig, mateStart, mateStart + 30);
                const auto mappedSweep = mappedFinder.queryStreaming(
                    contig, start, start + start % 120, mateContig, mateStart, mateStart + 30);
                ASSERT_EQ(expected.size(), observed.size());
                ASSERT_EQ(expected.size(), mappedRandom.size());
                ASSERT_EQ(expected.size(), mappedSweep.size());
                for (size_t i = 0; i < expected.size(); ++i)
                {
                    EXPECT_EQ(expected[i].locusIndex, observed[i].locusIndex);
                    EXPECT_EQ(expected[i].regionType, observed[i].regionType);
                    EXPECT_EQ(expected[i].inputType, observed[i].inputType);
                    EXPECT_EQ(expected[i].locusIndex, mappedRandom[i].locusIndex);
                    EXPECT_EQ(expected[i].regionType, mappedRandom[i].regionType);
                    EXPECT_EQ(expected[i].inputType, mappedRandom[i].inputType);
                    EXPECT_EQ(expected[i].locusIndex, mappedSweep[i].locusIndex);
                    EXPECT_EQ(expected[i].regionType, mappedSweep[i].regionType);
                    EXPECT_EQ(expected[i].inputType, mappedSweep[i].inputType);
                }
            }
    // A backward current-read query must still use the exact random-access path.
    const auto expected = finder.query(0, 37, 45, 2, 4000, 4030);
    const auto observed = finder.queryStreaming(0, 37, 45, 2, 4000, 4030);
    ASSERT_EQ(expected.size(), observed.size());
    for (size_t i = 0; i < expected.size(); ++i)
        EXPECT_EQ(expected[i].locusIndex, observed[i].locusIndex);
}

TEST(PreparedCatalog, WholeCatalogRoundTripAndIndependentHandleLifetime)
{
    using namespace ehunter;
    ArchiveFile archive, reference, catalogFile;
    std::ofstream(reference.path) << ">chr1\nACGTACGGCATAC\n";
    std::ofstream(catalogFile.path) << "unit-test-catalog-v1";
    ReferenceContigInfo contigs({ { "chr1", 10000 }, { "chr2", 10000 } });
    InputPaths inputs("unused.bam", reference.path, catalogFile.path);
    HeuristicParameters params(1000, 20, true, graphtools::AlignerType::DAG_ALIGNER, 2, 1, 1, 2, 1);
    auto graph = graphtools::makeDeletionGraph("ACGTAC", "GG", "CATAC");
    graph.setGraphId("mapped-graph");
    graph.setNodeName(0, "left");
    graph.addLabelToEdge(0, 1, "a");
    graph.addLabelToEdge(0, 1, "b");
    LocusSpecification spec(
        "locus", ChromType::kAutosome, { GenomicRegion(0, 90, 140) }, graph,
        { { 0, GenomicRegion(0, 100, 106) }, { 2, GenomicRegion(0, 108, 113) } }, GenotyperParameters());
    spec.setOfftargetReadExtractionRegions({ GenomicRegion(1, 200, 250) });
    spec.addVariantSpecification(
        "repeat", { VariantType::kRepeat, VariantSubtype::kRareRepeat }, GenomicRegion(0, 106, 108), { 1 }, 1);
    RegionCatalog original { spec };
    prepareCatalogCache(archive.path, original, inputs, params, contigs);
    auto mapped = loadPreparedCatalog(archive.path, inputs, params, contigs);
    ASSERT_EQ(1u, mapped.size());
    EXPECT_EQ(graphtools::graphToJson(graph), graphtools::graphToJson(mapped[0].regionGraph()));
    const auto loadedLocus = mapped[0];
    EXPECT_TRUE(loadedLocus.isPrepared());
    EXPECT_EQ(
        std::vector<GenomicRegion>(spec.targetReadExtractionRegions()), loadedLocus.targetReadExtractionRegions());
    EXPECT_EQ(
        std::vector<GenomicRegion>(spec.offtargetReadExtractionRegions()),
        loadedLocus.offtargetReadExtractionRegions());
    auto bytes = mapped.mappedBytes();
    const auto sequence = loadedLocus.regionGraph().nodeSeq(0);
    EXPECT_GE(sequence.data(), bytes.data);
    EXPECT_LT(sequence.data() + sequence.size(), bytes.data + bytes.size);
    const auto heldVariant = [&]()
    {
        auto temporary = loadPreparedCatalog(archive.path, inputs, params, contigs);
        return temporary[0].getVariantSpecById("repeat");
    }();
    auto heldGraph = [&]()
    {
        auto temporary = loadPreparedCatalog(archive.path, inputs, params, contigs);
        return temporary[0].regionGraph();
    }();
    EXPECT_EQ("ACGTAC", heldGraph.nodeSeq(0));
    EXPECT_THROW(heldGraph.setNodeSeq(0, "AAAAAA"), std::logic_error);
    EXPECT_EQ("repeat", heldVariant.id());
    WeightedPurityCalculator prepared(heldVariant.purityPermutations()), plain(std::string("GG"));
    for (const auto& read : std::vector<std::string> { "GGGGGG", "ggGgGG", "ACGTAC", "CCGGCC" })
        EXPECT_DOUBLE_EQ(plain.score(read), prepared.score(read));
    std::vector<std::unique_ptr<locus::LocusAnalyzer>> ownedAnalyzers, mappedAnalyzers;
    auto writer = std::make_shared<graphtools::BlankAlignmentWriter>();
    ownedAnalyzers.emplace_back(new locus::LocusAnalyzer(spec, params, writer));
    auto cachedParams = params;
    cachedParams.catalogIndexes = mapped.archive();
    mappedAnalyzers.emplace_back(new locus::LocusAnalyzer(mapped[0], cachedParams, writer, 1));
    GenomeQueryCollection ownedQueries(ownedAnalyzers), mappedQueries(mappedAnalyzers, &mapped);
    for (int contig = -1; contig < 3; ++contig)
        for (int start = 0; start < 500; ++start)
        {
            EXPECT_EQ(
                ownedQueries.targetRegionMask.query(contig, start),
                mappedQueries.targetRegionMask.query(contig, start));
            const auto expected = ownedQueries.analyzerFinder.query(contig, start, start + 20, 1, 215, 235);
            const auto observed = mappedQueries.analyzerFinder.queryStreaming(contig, start, start + 20, 1, 215, 235);
            ASSERT_EQ(expected.size(), observed.size());
            for (size_t i = 0; i < expected.size(); ++i)
            {
                EXPECT_EQ(expected[i].locusIndex, observed[i].locusIndex);
                EXPECT_EQ(expected[i].regionType, observed[i].regionType);
                EXPECT_EQ(expected[i].inputType, observed[i].inputType);
            }
        }
    auto mismatched = params;
    mismatched = HeuristicParameters(1200, 20, true, graphtools::AlignerType::DAG_ALIGNER, 2, 1, 1, 2, 1);
    EXPECT_THROW(loadPreparedCatalog(archive.path, inputs, mismatched, contigs), std::runtime_error);
    std::ofstream(catalogFile.path) << "changed catalog";
    EXPECT_THROW(loadPreparedCatalog(archive.path, inputs, params, contigs), std::runtime_error);
}
