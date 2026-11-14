#include "core/WeightedPurityCalculator.hh"
#include "graphalign/KmerIndex.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/PathOperations.hh"
#include "graphio/GraphJson.hh"
#include "graphutils/KmerEncoding.hh"
#include "graphutils/SequenceOperations.hh"
#include "io/PreparedCatalog.hh"
#include "sample/AnalyzerFinder.hh"
#include "sample/GenomeQueryCollection.hh"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <random>
#include <sstream>
#include <utility>
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

void expectSameBundles(
    const std::vector<ehunter::AnalyzerBundle>& expected, const std::vector<ehunter::AnalyzerBundle>& observed)
{
    ASSERT_EQ(expected.size(), observed.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(expected[i].locusIndex, observed[i].locusIndex);
        EXPECT_EQ(expected[i].regionType, observed[i].regionType);
        EXPECT_EQ(expected[i].inputType, observed[i].inputType);
    }
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

TEST(CompactCatalogIndex, ViewsPreservePathOrderMultiplicityAndRepeatedNodes)
{
    const auto graph = graphtools::makeStrGraph("ACN", "A", "GTA");
    for (int k = 1; k <= 6; ++k)
    {
        SCOPED_TRACE(k);
        const auto expected = enumerate(graph, k);
        const graphtools::KmerIndex index(graph, k), membership(graph, k, true);
        const graphtools::TwoBitKmerEncoder encoder(k);
        std::vector<size_t> uniqueNodeCounts(graph.numNodes(), 0);
        std::map<std::pair<graphtools::NodeId, graphtools::NodeId>, size_t> uniqueEdgeCounts;
        for (const auto& item : expected)
        {
            const auto encoded = encoder.encode(item.first);
            const auto view = index.findEncodedPaths(encoded);
            ASSERT_EQ(item.second.size(), view.size());
            EXPECT_EQ(item.second.size(), index.numPaths(item.first));
            EXPECT_EQ(item.second, index.getPaths(item.first));
            for (size_t i = 0; i < item.second.size(); ++i)
                EXPECT_EQ(item.second[i], view[i]);
            EXPECT_THROW(view[view.size()], std::out_of_range);
            EXPECT_TRUE(membership.containsEncodedKmer(encoded));
            EXPECT_THROW(membership.findEncodedPaths(encoded).size(), std::logic_error);
            if (item.second.size() == 1)
            {
                const auto& nodes = item.second.front().nodeIds();
                for (auto node : nodes)
                    ++uniqueNodeCounts[node];
                for (size_t i = 1; i < nodes.size(); ++i)
                    ++uniqueEdgeCounts[{ nodes[i - 1], nodes[i] }];
            }
        }
        for (graphtools::NodeId node = 0; node < graph.numNodes(); ++node)
        {
            EXPECT_EQ(uniqueNodeCounts[node], index.numUniqueKmersOverlappingNode(node));
            for (auto next : graph.successors(node))
                EXPECT_EQ(uniqueEdgeCounts[std::make_pair(node, next)], index.numUniqueKmersOverlappingEdge(node, next));
        }
    }
}

TEST(CompactCatalogIndex, EmptyIndexesAndMembershipSurviveMapping)
{
    graphtools::Graph graph(1);
    graph.setNodeSeq(0, "ACG");
    ArchiveFile file;
    graphtools::KmerIndex::prepareArchive(file.path, 2, [&](std::ostream& out, uint64_t slot) {
        graphtools::KmerIndex(graph, 16, slot != 0).writeRecord(out);
    });
    auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
    for (bool membership : { false, true })
    {
        const graphtools::KmerIndex index(graph, 16, membership), mapped(graph, 16, membership, archive, membership);
        for (const auto* candidate : { &index, &mapped })
        {
            EXPECT_TRUE(candidate->kmers().empty());
            EXPECT_TRUE(candidate->encode().empty());
            EXPECT_FALSE(candidate->containsEncodedKmer(0));
            EXPECT_FALSE(candidate->containsEncodedKmer(0xffffffffU));
            EXPECT_EQ(0u, candidate->findEncodedPaths(0).size());
            EXPECT_THROW(candidate->findEncodedPaths(0)[0], std::out_of_range);
            EXPECT_THROW(candidate->getPaths(std::string(16, 'A')), std::out_of_range);
        }
    }
}

TEST(CompactCatalogIndex, SixteenBaseKeysKeepAllBitsForMiniAndMultiplePaths)
{
    graphtools::Graph graph(4);
    graph.setNodeSeq(0, std::string(16, 'A'));
    graph.setNodeSeq(1, std::string(16, 'T'));
    graph.setNodeSeq(2, std::string(16, 'T'));
    graph.setNodeSeq(3, std::string(8, 'C') + std::string(8, 'G'));
    const graphtools::KmerIndex index(graph, 16), membership(graph, 16, true);
    ArchiveFile file;
    graphtools::KmerIndex::prepareArchive(file.path, 2, [&](std::ostream& out, uint64_t slot) {
        (slot ? membership : index).writeRecord(out);
    });
    auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
    const graphtools::KmerIndex mapped(graph, 16, false, archive, 0), mappedMembership(graph, 16, true, archive, 1);
    const std::vector<graphtools::Path> aPaths { graphtools::Path(&graph, 0, { 0 }, 16) };
    const std::vector<graphtools::Path> tPaths {
        graphtools::Path(&graph, 0, { 1 }, 16), graphtools::Path(&graph, 0, { 2 }, 16)
    };
    for (const auto* candidate : { &index, &mapped })
    {
        EXPECT_EQ(aPaths, candidate->getPaths(std::string(16, 'A')));
        EXPECT_EQ(tPaths, candidate->getPaths(std::string(16, 'T')));
        EXPECT_EQ(1u, candidate->findEncodedPaths(0).size());
        EXPECT_EQ(2u, candidate->findEncodedPaths(0xffffffffU).size());
    }
    for (const auto* candidate : { &index, &membership, &mapped, &mappedMembership })
    {
        EXPECT_TRUE(candidate->containsEncodedKmer(0));
        EXPECT_TRUE(candidate->containsEncodedKmer(0xffffffffU));
        EXPECT_FALSE(candidate->containsEncodedKmer(0x55555555U)); // C repeated 16 times
        EXPECT_FALSE(candidate->contains(std::string(16, 'N')));
    }
}

TEST(CompactCatalogIndex, MiniStartBoundaryAndFullPathFallbackAreLossless)
{
    const std::string tail = "ATGCGTACGTTAGCAATCGT";
    graphtools::Graph graph(1);
    graph.setNodeSeq(0, std::string(65535, 'C') + tail);
    const graphtools::KmerIndex index(graph, 16);
    ArchiveFile file;
    graphtools::KmerIndex::prepareArchive(file.path, 1, [&](std::ostream& out, uint64_t) { index.writeRecord(out); });
    auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
    const graphtools::KmerIndex mapped(graph, 16, false, archive, 0);
    ASSERT_LT(index.kmers().size(), 65535u); // Narrow buckets, with both mini and full entries.
    for (int position : { 65535, 65536 })
    {
        SCOPED_TRACE(position);
        const std::string key = tail.substr(position - 65535, 16);
        const std::vector<graphtools::Path> expected { graphtools::Path(&graph, position, { 0 }, position + 16) };
        EXPECT_EQ(expected, index.getPaths(key));
        EXPECT_EQ(expected, mapped.getPaths(key));
    }
    EXPECT_EQ(index, mapped);
}

TEST(CompactCatalogIndex, MiniNodeIdBoundaryFallsBackWithoutTruncation)
{
    graphtools::Graph graph(65537);
    for (graphtools::NodeId node = 0; node < 65535; ++node)
        graph.setNodeSeq(node, "C"); // These isolated nodes cannot provide a 16-mer.
    graph.setNodeSeq(65535, std::string(16, 'A'));
    graph.setNodeSeq(65536, std::string(16, 'T'));
    const graphtools::KmerIndex index(graph, 16);
    ArchiveFile file;
    graphtools::KmerIndex::prepareArchive(file.path, 1, [&](std::ostream& out, uint64_t) { index.writeRecord(out); });
    auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
    const graphtools::KmerIndex mapped(graph, 16, false, archive, 0);
    const std::vector<graphtools::Path> mini { graphtools::Path(&graph, 0, { 65535 }, 16) };
    const std::vector<graphtools::Path> full { graphtools::Path(&graph, 0, { 65536 }, 16) };
    for (const auto* candidate : { &index, &mapped })
    {
        EXPECT_EQ(mini, candidate->getPaths(std::string(16, 'A')));
        EXPECT_EQ(full, candidate->getPaths(std::string(16, 'T')));
        EXPECT_EQ(1u, candidate->numUniqueKmersOverlappingNode(65535));
        EXPECT_EQ(1u, candidate->numUniqueKmersOverlappingNode(65536));
    }
}

TEST(CompactCatalogIndex, LargeRandomGraphExercisesWideBuckets)
{
    std::string sequence(70000, 'A');
    uint32_t state = 173;
    for (char& base : sequence)
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        base = "ACGT"[state & 3];
    }
    graphtools::Graph graph(1);
    graph.setNodeSeq(0, sequence);
    const graphtools::KmerIndex index(graph, 16);
    const auto expected = enumerate(graph, 16);
    ASSERT_GT(expected.size(), 65535u);
    EXPECT_EQ(expected.size(), index.kmers().size());
    ArchiveFile file;
    graphtools::KmerIndex::prepareArchive(file.path, 1, [&](std::ostream& out, uint64_t) { index.writeRecord(out); });
    auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
    const graphtools::KmerIndex mapped(graph, 16, false, archive, 0);
    const graphtools::TwoBitKmerEncoder encoder(16);
    for (const auto& item : expected)
    {
        const auto key = encoder.encode(item.first);
        ASSERT_TRUE(index.containsEncodedKmer(key));
        ASSERT_TRUE(mapped.containsEncodedKmer(key));
        ASSERT_EQ(item.second, index.getPaths(item.first));
        ASSERT_EQ(item.second, mapped.getPaths(item.first));
    }
    EXPECT_FALSE(index.containsEncodedKmer(0));
    EXPECT_FALSE(mapped.containsEncodedKmer(0));
}

TEST(CompactCatalogIndex, BucketIndexWidthBoundaryPreservesAllKeys)
{
    // B(4,8) has all 65,536 eight-base words exactly once. Trimming one
    // character gives precisely 65,535 entries, the last uint16_t index.
    std::string sequence;
    std::vector<int> digits(33, 0);
    std::function<void(int, int)> appendCycle = [&](int depth, int period) {
        if (depth > 8)
        {
            if (8 % period == 0)
                for (int i = 1; i <= period; ++i)
                    sequence += "ACGT"[digits[i]];
            return;
        }
        digits[depth] = digits[depth - period];
        appendCycle(depth + 1, period);
        for (int value = digits[depth - period] + 1; value < 4; ++value)
        {
            digits[depth] = value;
            appendCycle(depth + 1, depth);
        }
    };
    appendCycle(1, 1);
    ASSERT_EQ(65536u, sequence.size());
    sequence += sequence.substr(0, 7);
    const graphtools::TwoBitKmerEncoder encoder(8);
    for (size_t entryCount : { 65535u, 65536u })
    {
        SCOPED_TRACE(entryCount);
        graphtools::Graph graph(1);
        graph.setNodeSeq(0, sequence.substr(0, entryCount + 7));
        const graphtools::KmerIndex index(graph, 8);
        ASSERT_EQ(entryCount, index.kmers().size());
        ArchiveFile file;
        graphtools::KmerIndex::prepareArchive(file.path, 1, [&](std::ostream& out, uint64_t) { index.writeRecord(out); });
        auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
        const graphtools::KmerIndex mapped(graph, 8, false, archive, 0);
        for (size_t position = 0; position < entryCount; ++position)
        {
            const auto key = encoder.encode(sequence.substr(position, 8));
            ASSERT_TRUE(index.containsEncodedKmer(key));
            ASSERT_TRUE(mapped.containsEncodedKmer(key));
        }
        const auto lastKey = sequence.substr(entryCount - 1, 8);
        const auto position = static_cast<int32_t>(entryCount - 1);
        const std::vector<graphtools::Path> expected { graphtools::Path(&graph, position, { 0 }, position + 8) };
        EXPECT_EQ(expected, index.getPaths(lastKey));
        EXPECT_EQ(expected, mapped.getPaths(lastKey));
        if (entryCount == 65535)
        {
            const auto missingKey = encoder.encode(sequence.substr(65535, 8));
            EXPECT_FALSE(index.containsEncodedKmer(missingKey));
            EXPECT_FALSE(mapped.containsEncodedKmer(missingKey));
        }
    }
}

TEST(CompactCatalogIndex, MappedCopyAndMoveKeepArchiveAlive)
{
    const auto graph = graphtools::makeStrGraph("ACGTACGTACGTACGT", "A", "TTTTACGTACGTACGT");
    const graphtools::KmerIndex expected(graph, 14);
    auto retained = [&]() {
        ArchiveFile file;
        graphtools::KmerIndex::prepareArchive(file.path, 1, [&](std::ostream& out, uint64_t) { expected.writeRecord(out); });
        auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
        const graphtools::KmerIndex mapped(graph, 14, false, archive, 0);
        graphtools::KmerIndex copy(mapped);
        graphtools::KmerIndex assigned(graph, 14);
        assigned = copy;
        archive.reset();
        return assigned;
    }(); // The source indexes, external archive handle, and file name are gone.
    graphtools::KmerIndex moved(std::move(retained));
    graphtools::KmerIndex assigned(graph, 14);
    assigned = std::move(moved);
    EXPECT_EQ(expected, assigned);
    EXPECT_EQ(expected.encode(), assigned.encode());
}

namespace
{
std::string randomReference(size_t length, uint32_t seed)
{
    std::string sequence(length, 'A');
    for (auto& base : sequence)
    {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        base = "ACGT"[seed & 3];
    }
    return sequence;
}

bool usesPrefixIndex(const graphtools::KmerIndex& index)
{
    std::ostringstream stream;
    index.writeRecord(stream);
    const auto record = stream.str();
    uint64_t identitySize;
    std::memcpy(&identitySize, record.data() + 24, sizeof(identitySize));
    uint32_t buckets;
    std::memcpy(&buckets, record.data() + 24 + 16 + identitySize + 4, sizeof(buckets));
    return buckets == 0;
}
}

TEST(CompactCatalogIndex, PrefixBucketsPreserveAllKeysOrderedPathsAndArchiveLifetime)
{
    auto left = randomReference(1000, 7813), right = randomReference(1000, 5627);
    // Include expanded bases, a repeated k-mer, and paths crossing a loop.
    left[150] = 'N';
    right.replace(170, 16, left.substr(500, 16));
    const auto graph = graphtools::makeStrGraph(left, "CAG", right);
    for (int k : { 10, 14 })
    {
        const bool memberOnly = k == 10;
        const auto expected = enumerate(graph, k);
        const graphtools::KmerIndex original(graph, k, memberOnly);
        ASSERT_TRUE(usesPrefixIndex(original));
        auto retained = [&]() {
            ArchiveFile file;
            graphtools::KmerIndex::prepareArchive(file.path, 1, [&](std::ostream& out, uint64_t) {
                original.writeRecord(out);
            });
            auto archive = std::make_shared<graphtools::KmerIndexArchive>(file.path);
            graphtools::KmerIndex mapped(graph, k, memberOnly, archive, 0);
            graphtools::KmerIndex copy(mapped);
            archive.reset();
            return copy;
        }();
        graphtools::KmerIndex moved(std::move(retained));
        EXPECT_EQ(original.encode(), moved.encode());
        EXPECT_EQ(original.kmers(), moved.kmers());
        const graphtools::TwoBitKmerEncoder encoder(k);
        for (const auto& item : expected)
        {
            const auto key = encoder.encode(item.first);
            ASSERT_TRUE(original.containsEncodedKmer(key));
            ASSERT_TRUE(moved.containsEncodedKmer(key));
            if (!memberOnly)
            {
                ASSERT_EQ(item.second, original.getPaths(item.first));
                ASSERT_EQ(item.second, moved.getPaths(item.first));
                const auto view = moved.findEncodedPaths(key);
                ASSERT_EQ(item.second.size(), view.size());
                for (size_t i = 0; i < view.size(); ++i)
                    ASSERT_EQ(item.second[i], view[i]);
            }
        }
        uint32_t random = 817239;
        for (int trial = 0; trial < 100000; ++trial)
        {
            random ^= random << 13;
            random ^= random >> 17;
            random ^= random << 5;
            const auto key = random & ((uint32_t(1) << (2 * k)) - 1);
            const bool present = expected.count(encoder.decode(key));
            ASSERT_EQ(present, original.containsEncodedKmer(key));
            ASSERT_EQ(present, moved.containsEncodedKmer(key));
        }
        EXPECT_FALSE(original.containsEncodedKmer(0xffffffffU));
    }
}

TEST(CompactCatalogIndex, PrefixPackingFallbackPreservesLongAndSkewedGraphs)
{
    for (bool longGraph : { false, true })
    {
        auto sequence = randomReference(longGraph ? 20000 : 1500, 71673);
        if (!longGraph)
            for (size_t i = 0; i + 16 < sequence.size(); i += 16)
                sequence.replace(i, 6, "AAAAAA");
        graphtools::Graph graph(1);
        graph.setNodeSeq(0, sequence);
        const graphtools::KmerIndex index(graph, 14);
        ASSERT_FALSE(usesPrefixIndex(index));
        for (const auto& item : enumerate(graph, 14))
            ASSERT_EQ(item.second, index.getPaths(item.first));
    }
}

TEST(CompactCatalogIndex, PreviousArchiveVersionIsRejected)
{
    const auto graph = graphtools::makeDeletionGraph("ACGT", "A", "TGCA");
    ArchiveFile file;
    graphtools::KmerIndex::prepareArchive(file.path, 1, [&](std::ostream& out, uint64_t) {
        graphtools::KmerIndex(graph, 4).writeRecord(out);
    });
    std::fstream oldFormat(file.path, std::ios::in | std::ios::out | std::ios::binary);
    const uint64_t previousMagic = 0x0000000345484b49ULL;
    oldFormat.write(reinterpret_cast<const char*>(&previousMagic), sizeof(previousMagic));
    oldFormat.close();
    EXPECT_THROW(graphtools::KmerIndexArchive(file.path), std::runtime_error);
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
    graph.setNodeSeq(0, "ACGTACGTACGTACGT");
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
    try
    {
        graphtools::KmerIndex unexpected(graph, 10, true, corrupt, 1);
        FAIL() << "Expected a checksum failure";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_NE(std::string::npos, std::string(error.what()).find("checksum"));
    }
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

    // Random read lengths and remote mate coordinates exercise both advancing
    // sweeps and their backward-query fallback, using the unchanged query as oracle.
    std::mt19937 random(41723);
    for (int round = 0; round < 4; ++round)
    {
        AnalyzerFinder streaming(analyzers), mappedStreaming(mappedAnalyzers, &mapped);
        for (int position = 0; position < 5500; position += 3)
        {
            SCOPED_TRACE(position);
            const int contig = int(random() % 3);
            const int start = position % 2 ? int(random() % (position + 1)) : position;
            const int end = start + random() % 180;
            const int mateContig = int(random() % 5) - 1;
            const int mateStart = random() % 5500;
            const int mateEnd = mateStart + random() % 180;
            const auto legacy = finder.query(contig, start, end, mateContig, mateStart, mateEnd);
            expectSameBundles(legacy, streaming.queryStreaming(contig, start, end, mateContig, mateStart, mateEnd));
            expectSameBundles(
                legacy, mappedStreaming.queryStreaming(contig, start, end, mateContig, mateStart, mateEnd));
        }
    }
}

TEST(CatalogCoordinateSweep, EmptyAndSingletonCoalescingMatchesLegacyQueries)
{
    using namespace ehunter;
    HeuristicParameters params(1000, 20, true, graphtools::AlignerType::DAG_ALIGNER, 2, 1, 1, 2, 1);
    auto writer = std::make_shared<graphtools::BlankAlignmentWriter>();
    std::vector<std::unique_ptr<locus::LocusAnalyzer>> analyzers;
    for (int contig = 0; contig < 2; ++contig)
    {
        LocusSpecification spec(
            std::to_string(contig), ChromType::kAutosome, { GenomicRegion(contig, 100, 200) },
            graphtools::makeDeletionGraph("ACGT", "C", "ACGT"), {}, GenotyperParameters());
        spec.setOfftargetReadExtractionRegions(
            { GenomicRegion(contig, 300, 400), GenomicRegion(contig, 2000, 2100) });
        analyzers.emplace_back(new locus::LocusAnalyzer(spec, params, writer));
    }
    AnalyzerFinder legacy(analyzers), streaming(analyzers);
    const std::vector<int> positions { 0, 99, 100, 150, 200, 250, 300, 350, 400, 1000, 2000, 2050, 2200 };
    for (int contig = -1; contig < 3; ++contig)
        for (int start : positions)
            for (int mateContig = -1; mateContig < 3; ++mateContig)
                for (int mateStart : positions)
                    for (int length : { 0, 20, 120 })
                    {
                        SCOPED_TRACE(start);
                        SCOPED_TRACE(mateStart);
                        expectSameBundles(
                            legacy.query(contig, start, start + length, mateContig, mateStart, mateStart + length),
                            streaming.queryStreaming(
                                contig, start, start + length, mateContig, mateStart, mateStart + length));
                    }
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
