#include "io/PreparedCatalog.hh"
#include "sample/PreparedIntervals.hh"
#include "thirdparty/intervaltree/IntervalTree.h"
#include <array>
#include <boost/uuid/detail/sha1.hpp>
#include <fstream>
#include <map>
#include <sstream>

namespace ehunter
{
#ifndef EH_PREPARED_SCHEMA_ID
#error "Prepared catalog builds require the CMake source fingerprint"
#endif
namespace
{
using namespace graphtools;
using Digest = std::array<unsigned int, 5>;
struct Signature
{
    Digest catalog, reference, contigs;
    std::array<char, 64> builder;
    int32_t extension, orientationK, alignmentK;
};
struct Root
{
    uint64_t magic, bytes;
    Slice loci, ids, intervals;
    Signature signature;
    uint64_t locusWidth, variantWidth, regionWidth, graphWidth;
};
constexpr uint64_t catalogMagic = 0x000000014548434dULL;

Digest hashFile(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("Cannot fingerprint cache input: " + path);
    boost::uuids::detail::sha1 hash;
    std::array<char, 1024 * 1024> buffer;
    while (in)
    {
        in.read(buffer.data(), buffer.size());
        hash.process_bytes(buffer.data(), in.gcount());
    }
    if (!in.eof())
        throw std::runtime_error("Cannot read cache input: " + path);
    Digest result;
    hash.get_digest(*reinterpret_cast<unsigned int(*)[5]>(result.data()));
    return result;
}
Signature signature(const InputPaths& inputs, const HeuristicParameters& params, const ReferenceContigInfo& contigs)
{
    std::ostringstream text;
    text << contigs;
    boost::uuids::detail::sha1 hash;
    const auto bytes = text.str();
    hash.process_bytes(bytes.data(), bytes.size());
    Digest contigHash;
    hash.get_digest(*reinterpret_cast<unsigned int(*)[5]>(contigHash.data()));
    std::array<char, 64> builder;
    const std::string schema = EH_PREPARED_SCHEMA_ID;
    if (schema.size() != builder.size())
        throw std::logic_error("Invalid prepared schema fingerprint");
    std::copy(schema.begin(), schema.end(), builder.begin());
    return { hashFile(inputs.catalog()),
             hashFile(inputs.reference()),
             contigHash,
             builder,
             params.regionExtensionLength(),
             params.orientationPredictorKmerLen(),
             params.kmerLenForAlignment() };
}
bool sameSignature(const Signature& a, const Signature& b)
{
    return a.catalog == b.catalog && a.reference == b.reference && a.contigs == b.contigs && a.builder == b.builder
        && a.extension == b.extension && a.orientationK == b.orientationK && a.alignmentK == b.alignmentK;
}

using Tree = IntervalTree<size_t, PreparedInterval>;
Slice freezeTree(const Tree* tree, BlobBuilder& out, uint64_t& rank, std::vector<PreparedInterval>& sorted)
{
    if (!tree)
        return { 0, 0 };
    const auto record = out.append(PreparedIntervalNode {});
    std::vector<PreparedInterval> intervals;
    for (const auto& item : tree->localIntervals())
    {
        auto interval = item.value;
        interval.rank = rank++;
        intervals.push_back(interval);
        sorted.push_back(interval);
    }
    const auto values = out.append(intervals);
    const auto left = freezeTree(tree->leftChild(), out, rank, sorted);
    const auto right = freezeTree(tree->rightChild(), out, rank, sorted);
    out.replace(record.offset, PreparedIntervalNode { tree->centerCoordinate(), values, left, right });
    return record;
}

Slice freezeIntervals(const RegionCatalog& catalog, BlobBuilder& out)
{
    std::map<int32_t, std::vector<Tree::interval>> byContig;
    for (size_t index = 0; index < catalog.size(); ++index)
    {
        const auto locus = catalog[index];
        auto add = [&](graphtools::CatalogRange<GenomicRegion> regions, uint32_t regionType)
        {
            for (const auto& region : regions)
            {
                PreparedInterval interval { uint64_t(region.start()), uint64_t(region.end()), index, 0, regionType, 0 };
                byContig[region.contigIndex()].emplace_back(region.start(), region.end(), interval);
            }
        };
        add(locus.targetReadExtractionRegions(), 0);
        add(locus.offtargetReadExtractionRegions(), 1);
    }
    std::vector<PreparedContigIntervals> roots;
    if (!byContig.empty())
        roots.resize(size_t(byContig.rbegin()->first) + 1);
    for (auto& item : byContig)
    {
        uint64_t maxBin = 0;
        for (const auto& interval : item.second)
            maxBin = std::max(maxBin, uint64_t(interval.stop) >> 10);
        std::vector<uint64_t> mask(maxBin / 64 + 1, 0);
        for (const auto& interval : item.second)
            for (uint64_t bin = uint64_t(interval.start) >> 10; bin <= (uint64_t(interval.stop) >> 10); ++bin)
                mask[bin / 64] |= uint64_t(1) << (bin % 64);
        Tree tree(std::move(item.second));
        std::vector<PreparedInterval> sorted;
        uint64_t rank = 0;
        const auto treeSlice = freezeTree(&tree, out, rank, sorted);
        std::sort(
            sorted.begin(), sorted.end(),
            [](const PreparedInterval& a, const PreparedInterval& b) { return a.start < b.start; });
        roots[item.first] = { item.first, 0, treeSlice, out.append(sorted), out.append(mask) };
    }
    return out.append(roots);
}

std::string freezeCatalog(const RegionCatalog& catalog, const Signature& identity)
{
    BlobBuilder out;
    Root root {};
    out.append(root);
    std::vector<MappedLocus> loci;
    loci.reserve(catalog.size());
    for (const auto& locus : catalog)
        loci.push_back(locus.freeze(out));
    root.loci = out.append(loci);
    size_t buckets = 1;
    while (buckets < catalog.size() * 2)
        buckets *= 2;
    std::vector<uint64_t> ids(buckets, 0);
    for (size_t i = 0; i < catalog.size(); ++i)
    {
        const auto id = catalog[i].locusId();
        size_t bucket = hashText(id) & (buckets - 1);
        while (ids[bucket] && catalog[ids[bucket] - 1].locusId() != id)
            bucket = (bucket + 1) & (buckets - 1);
        if (!ids[bucket])
            ids[bucket] = i + 1;
    }
    root.ids = out.append(ids);
    root.intervals = freezeIntervals(catalog, out);
    root.magic = catalogMagic;
    root.bytes = out.bytes.size();
    root.signature = identity;
    root.locusWidth = sizeof(MappedLocus);
    root.variantWidth = sizeof(MappedVariant);
    root.regionWidth = sizeof(GenomicRegion);
    root.graphWidth = sizeof(MappedGraph);
    out.replace(0, root);
    return std::move(out.bytes);
}
}

void prepareCatalogCache(
    const std::string& path, const RegionCatalog& catalog, const InputPaths& inputs, const HeuristicParameters& params,
    const ReferenceContigInfo& contigs)
{
    const auto identity = signature(inputs, params, contigs);
    graphtools::KmerIndex::prepareArchive(
        path, uint64_t(catalog.size()) * 2 + 1,
        [&](std::ostream& out, uint64_t slot)
        {
            if (slot == 0)
            {
                graphtools::KmerIndex::writeBlobRecord(out, freezeCatalog(catalog, identity));
                return;
            }
            const uint64_t indexSlot = slot - 1;
            const bool membershipOnly = (indexSlot % 2) != 0;
            const int k = membershipOnly ? params.orientationPredictorKmerLen() : params.kmerLenForAlignment();
            const auto locus = catalog[indexSlot / 2];
            graphtools::KmerIndex index(locus.regionGraph(), k, membershipOnly);
            index.writeRecord(out);
        });
}

RegionCatalog loadPreparedCatalog(
    const std::string& path, const InputPaths& inputs, const HeuristicParameters& params,
    const ReferenceContigInfo& contigs)
{
    auto archive = std::make_shared<graphtools::KmerIndexArchive>(path);
    size_t length;
    const char* data = archive->record(0, length);
    graphtools::MappedBytes bytes(data, length, archive);
    const auto& root = *bytes.array<Root>({ 0, 1 });
    if (root.magic != catalogMagic || root.bytes != length || root.locusWidth != sizeof(MappedLocus)
        || root.variantWidth != sizeof(MappedVariant) || root.regionWidth != sizeof(GenomicRegion)
        || root.graphWidth != sizeof(graphtools::MappedGraph) || archive->slotCount() != root.loci.count * 2 + 1)
        throw std::runtime_error("Unsupported prepared catalog format");
    if (!sameSignature(root.signature, signature(inputs, params, contigs)))
        throw std::runtime_error("Prepared catalog input/reference/contig/parameter signature mismatch");
    bytes.array<uint64_t>(root.ids);
    bytes.array<PreparedContigIntervals>(root.intervals);
    return RegionCatalog(bytes, root.loci, std::move(archive), root.ids, root.intervals);
}
}
