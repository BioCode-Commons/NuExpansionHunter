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
#include "io/PreparedCatalog.hh"
#include "core/OrderedBlocks.hh"
#include "sample/PreparedIntervals.hh"
#include "thirdparty/intervaltree/IntervalTree.h"
#include <array>
#include <cstring>
#include <memory>
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

using ContigIntervals = std::map<int32_t, std::vector<Tree::interval>>;

// Adds the target and off-target regions of one locus, in catalog order.
void addIntervals(ContigIntervals& byContig, size_t index, const GenomicRegion* targets, uint64_t targetCount,
                  const GenomicRegion* offtargets, uint64_t offtargetCount)
{
    auto add = [&](const GenomicRegion* regions, uint64_t count, uint32_t regionType)
    {
        for (uint64_t i = 0; i < count; ++i)
        {
            const auto& region = regions[i];
            PreparedInterval interval { uint64_t(region.start()), uint64_t(region.end()), index, 0, regionType, 0 };
            byContig[region.contigIndex()].emplace_back(region.start(), region.end(), interval);
        }
    };
    add(targets, targetCount, 0);
    add(offtargets, offtargetCount, 1);
}

uint64_t roundUp16(uint64_t value) { return (value + 15) / 16 * 16; }

// Contigs are independent: their trees are built and frozen on workers, each
// at a final offset measured in a first pass. Parts start at multiples of 16
// bytes, so padding, and therefore every part's size, is the same in both
// passes. Appends to out, which may itself start at a base offset.
Slice freezeIntervals(ContigIntervals& byContig, BlobBuilder& out, unsigned workers)
{
    std::vector<PreparedContigIntervals> roots;
    if (!byContig.empty())
        roots.resize(size_t(byContig.rbegin()->first) + 1);
    struct Contig
    {
        int32_t id;
        std::vector<Tree::interval>* intervals;
        std::unique_ptr<Tree> tree;
        std::vector<uint64_t> mask;
    };
    std::vector<Contig> contigs;
    for (auto& item : byContig)
    {
        uint64_t maxBin = 0;
        for (const auto& interval : item.second)
            maxBin = std::max(maxBin, uint64_t(interval.stop) >> 10);
        std::vector<uint64_t> mask(maxBin / 64 + 1, 0);
        for (const auto& interval : item.second)
            for (uint64_t bin = uint64_t(interval.start) >> 10; bin <= (uint64_t(interval.stop) >> 10); ++bin)
                mask[bin / 64] |= uint64_t(1) << (bin % 64);
        contigs.push_back({ item.first, &item.second, nullptr, std::move(mask) });
    }
    const auto freezeContig = [&](size_t index, BlobBuilder& part) {
        auto& contig = contigs[index];
        std::vector<PreparedInterval> sorted;
        uint64_t rank = 0;
        const auto treeSlice = freezeTree(contig.tree.get(), part, rank, sorted);
        std::sort(
            sorted.begin(), sorted.end(),
            [](const PreparedInterval& a, const PreparedInterval& b) { return a.start < b.start; });
        const auto sortedSlice = part.append(sorted);
        return PreparedContigIntervals { contig.id, 0, treeSlice, sortedSlice, part.append(contig.mask) };
    };
    std::vector<uint64_t> bases(contigs.size() + 1, 0);
    orderedBlocks(contigs.size(), 1, workers, [&](int, size_t index, size_t) {
        contigs[index].tree.reset(new Tree(std::move(*contigs[index].intervals)));
        BlobBuilder scratch;
        freezeContig(index, scratch);
        bases[index + 1] = roundUp16(scratch.bytes.size());
        return true;
    }, [](bool) {});
    bases[0] = roundUp16(out.base + out.bytes.size());
    for (size_t index = 0; index < contigs.size(); ++index)
        bases[index + 1] += bases[index];
    out.bytes.resize(bases[contigs.size()] - out.base, '\0');
    orderedBlocks(contigs.size(), 1, workers, [&](int, size_t index, size_t) {
        BlobBuilder part(bases[index]);
        roots[contigs[index].id] = freezeContig(index, part);
        if (part.bytes.size() > bases[index + 1] - bases[index])
            throw std::logic_error("Frozen interval block changed size");
        std::memcpy(&out.bytes[bases[index] - out.base], part.bytes.data(), part.bytes.size());
        return true;
    }, [](bool) {});
    return out.append(roots);
}

// Locus ID hash table (open addressing over 2x the loci, rounded to a power of
// two), from IDs already frozen in a blob.
std::vector<uint64_t> locusIdTable(const std::vector<MappedLocus>& loci, const graphtools::MappedBytes& blob)
{
    size_t buckets = 1;
    while (buckets < loci.size() * 2)
        buckets *= 2;
    std::vector<uint64_t> ids(buckets, 0);
    for (size_t i = 0; i < loci.size(); ++i)
    {
        const auto id = blob.text(loci[i].id);
        size_t bucket = hashText(id) & (buckets - 1);
        while (ids[bucket] && blob.text(loci[ids[bucket] - 1].id) != id)
            bucket = (bucket + 1) & (buckets - 1);
        if (!ids[bucket])
            ids[bucket] = i + 1;
    }
    return ids;
}

// Serial layout for prepared archives.
void freezeLayout(const RegionCatalog& catalog, BlobBuilder& out, Root& root)
{
    std::vector<MappedLocus> loci;
    loci.reserve(catalog.size());
    for (const auto& locus : catalog)
        loci.push_back(locus.freeze(out));
    root.loci = out.append(loci);
    const graphtools::MappedBytes blob(out.bytes.data(), out.bytes.size());
    root.ids = out.append(locusIdTable(loci, blob));
    ContigIntervals byContig;
    for (size_t index = 0; index < loci.size(); ++index)
    {
        const graphtools::MappedBytes view(out.bytes.data(), out.bytes.size());
        addIntervals(byContig, index, view.array<GenomicRegion>(loci[index].targets), loci[index].targets.count,
            view.array<GenomicRegion>(loci[index].offtargets), loci[index].offtargets.count);
    }
    root.intervals = freezeIntervals(byContig, out, 1);
}

std::string freezeCatalog(const RegionCatalog& catalog, const Signature& identity)
{
    BlobBuilder out;
    Root root {};
    out.append(root);
    freezeLayout(catalog, out, root);
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

// Contract [K-28]
// Ensures:  every accessor used by alignment, routing, genotyping and the writers returns the same value on the
//           compact catalog as on the owned one; findLocusIndex agrees; a cache accepted by the signature check is
//           byte-equal to the layout frozen from the parsed catalog.
// Reference: the owned catalog (EH_VERIFY_LAYOUT=REFERENCE keeps it in streaming).
RegionCatalog compactCatalog(const RegionCatalog& catalog, unsigned workers)
{
    // Pass 1, on workers: measure each block of frozen loci and collect its
    // intervals in catalog order.
    const size_t grain = 1024, count = catalog.size(), blocks = (count + grain - 1) / grain;
    std::vector<uint64_t> bases(blocks + 1, 0);
    std::vector<ContigIntervals> blockIntervals(blocks);
    orderedBlocks(count, grain, workers, [&](int, size_t begin, size_t end) {
        BlobBuilder scratch;
        auto& intervals = blockIntervals[begin / grain];
        for (size_t i = begin; i < end; ++i)
        {
            const auto locus = catalog[i];
            locus.freeze(scratch);
            const std::vector<GenomicRegion> targets = locus.targetReadExtractionRegions();
            const std::vector<GenomicRegion> offtargets = locus.offtargetReadExtractionRegions();
            addIntervals(intervals, i, targets.data(), targets.size(), offtargets.data(), offtargets.size());
        }
        bases[begin / grain + 1] = roundUp16(scratch.bytes.size());
        return true;
    }, [](bool) {});
    bases[0] = roundUp16(sizeof(Root));
    for (size_t block = 0; block < blocks; ++block)
        bases[block + 1] += bases[block];
    // Tail at known offsets: the locus records, the ID table, then the interval trees.
    size_t buckets = 1;
    while (buckets < count * 2)
        buckets *= 2;
    const Slice lociSlice { bases[blocks], count };
    const Slice idsSlice { roundUp16(lociSlice.offset + count * sizeof(MappedLocus)), buckets };
    BlobBuilder tail(roundUp16(idsSlice.offset + buckets * sizeof(uint64_t)));
    ContigIntervals byContig;
    for (auto& intervals : blockIntervals)
    {
        for (auto& item : intervals)
        {
            auto& target = byContig[item.first];
            target.insert(target.end(), std::make_move_iterator(item.second.begin()),
                std::make_move_iterator(item.second.end()));
        }
        ContigIntervals().swap(intervals);
    }
    const Slice intervalsSlice = freezeIntervals(byContig, tail, workers);
    // One uninitialized allocation; workers write the locus blocks in place.
    // Padding between parts is never read.
    const size_t size = tail.base + tail.bytes.size();
    std::shared_ptr<char> storage(new char[size], std::default_delete<char[]>());
    std::memset(storage.get(), 0, bases[0]);
    std::vector<MappedLocus> loci(count);
    orderedBlocks(count, grain, workers, [&](int, size_t begin, size_t end) {
        const size_t block = begin / grain;
        BlobBuilder part(bases[block]);
        for (size_t i = begin; i < end; ++i)
            loci[i] = catalog[i].freeze(part);
        if (part.bytes.size() > bases[block + 1] - bases[block])
            throw std::logic_error("Frozen catalog block changed size");
        std::memcpy(storage.get() + bases[block], part.bytes.data(), part.bytes.size());
        return true;
    }, [](bool) {});
    std::memcpy(storage.get() + lociSlice.offset, loci.data(), count * sizeof(MappedLocus));
    const graphtools::MappedBytes blob(storage.get(), size);
    const auto ids = locusIdTable(loci, blob);
    std::memcpy(storage.get() + idsSlice.offset, ids.data(), buckets * sizeof(uint64_t));
    std::memcpy(storage.get() + tail.base, tail.bytes.data(), tail.bytes.size());
    return RegionCatalog(std::shared_ptr<const char>(std::move(storage)), size, lociSlice, idsSlice, intervalsSlice);
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
