//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
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

#include "graphalign/KmerIndex.hh"

#include <algorithm>
#include <boost/config.hpp>
#include <boost/uuid/detail/sha1.hpp>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <sparsepp/spp.h>

#include "graphcore/PathOperations.hh"
#include "graphutils/KmerEncoding.hh"
#include "graphutils/PairHashing.hh"
#include "graphutils/SequenceOperations.hh"

using std::list;
using std::string;
using std::unordered_set;
using std::vector;

namespace graphtools
{

namespace
{
template <typename T> T readValue(const char* data)
{
    T value;
    std::memcpy(&value, data, sizeof(value));
    return value;
}
template <typename T> void appendValue(std::string& data, const T& value)
{
    data.append(reinterpret_cast<const char*>(&value), sizeof(value));
}
template <typename T> void writeValue(char*& cursor, const T& value)
{
    std::memcpy(cursor, &value, sizeof(value));
    cursor += sizeof(value);
}
std::string graphIdentity(const Graph& graph)
{
    std::string result;
    appendValue(result, uint64_t(graph.numNodes()));
    for (NodeId node = 0; node < graph.numNodes(); ++node)
    {
        appendValue(result, uint64_t(graph.nodeSeq(node).size()));
        result += graph.nodeSeq(node);
        appendValue(result, uint64_t(graph.successors(node).size()));
        for (auto successor : graph.successors(node))
            appendValue(result, successor);
    }
    return result;
}
std::string checksum(const char* data, size_t size)
{
    boost::uuids::detail::sha1 sha;
    sha.process_bytes(data, size);
    unsigned int digest[5];
    sha.get_digest(digest);
    return std::string(reinterpret_cast<const char*>(digest), sizeof(digest));
}
struct Entry
{
    uint32_t key, payload, count;
};
struct StoredPath
{
    int32_t start, end;
    uint32_t nodes, count;
};
uint32_t narrow(size_t value)
{
    if (value > std::numeric_limits<uint32_t>::max())
        throw std::overflow_error("K-mer index exceeds 32-bit offsets");
    return static_cast<uint32_t>(value);
}
uint32_t hashKey(uint32_t key)
{
    key ^= key >> 16;
    key *= 0x7feb352dU;
    key ^= key >> 15;
    key *= 0x846ca68bU;
    return key ^ (key >> 16);
}
// Version 3 stores narrow hash buckets and omits counts for the minipath prefix.
// Reject earlier archives rather than interpreting their layout as this one.
constexpr uint64_t archiveMagic = 0x0000000345484b49ULL;
constexpr size_t indexHeaderSize = 20;
}

KmerIndexArchive::KmerIndexArchive(const std::string& path)
{
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0)
        throw std::runtime_error("Cannot open k-mer archive: " + path);
    struct stat statbuf;
    if (fstat(fd, &statbuf) != 0 || statbuf.st_size < 16)
    {
        close(fd);
        throw std::runtime_error("Invalid k-mer archive size");
    }
    size_ = static_cast<size_t>(statbuf.st_size);
    void* mapping = mmap(nullptr, size_, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (mapping == MAP_FAILED)
        throw std::runtime_error("Cannot map k-mer archive");
    data_ = static_cast<const char*>(mapping);
    if (readValue<uint64_t>(data_) != archiveMagic || slotCount() > (size_ - 16) / 16)
    {
        munmap(const_cast<char*>(data_), size_);
        data_ = nullptr;
        throw std::runtime_error("Unsupported or truncated k-mer archive");
    }
}
KmerIndexArchive::~KmerIndexArchive()
{
    if (data_)
        munmap(const_cast<char*>(data_), size_);
}
uint64_t KmerIndexArchive::slotCount() const { return readValue<uint64_t>(data_ + 8); }
const char* KmerIndexArchive::record(uint64_t slot, size_t& size) const
{
    if (slot >= slotCount())
        throw std::runtime_error("K-mer archive catalog size mismatch");
    const uint64_t offset = readValue<uint64_t>(data_ + 16 + slot * 16);
    const uint64_t bytes = readValue<uint64_t>(data_ + 24 + slot * 16);
    if (offset < 16 + slotCount() * 16 || offset > size_ || bytes > size_ - offset || bytes < 24 || offset % 8)
        throw std::runtime_error("Invalid k-mer archive record");
    size = static_cast<size_t>(bytes - 24);
    const char* payload = data_ + offset + 24;
    if (checksum(payload, size) != std::string(data_ + offset, 20))
        throw std::runtime_error("K-mer archive checksum mismatch");
    return payload;
}

void KmerIndex::prepareArchive(
    const std::string& path, uint64_t slots, const std::function<void(std::ostream&, uint64_t)>& writeSlot)
{
    // Publish only a complete file. A sibling temporary keeps link() on the
    // same filesystem and atomically rejects an already-published destination.
    std::string pattern = path + ".tmp.XXXXXX";
    std::vector<char> temporary(pattern.begin(), pattern.end());
    temporary.push_back('\0');
    const int fd = mkstemp(temporary.data());
    if (fd < 0)
        throw std::runtime_error("Cannot create archive temporary: " + path);
    close(fd);
    try
    {
        std::ofstream out(temporary.data(), std::ios::binary | std::ios::trunc);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        const uint64_t zero = 0;
        out.write(reinterpret_cast<const char*>(&archiveMagic), 8);
        out.write(reinterpret_cast<const char*>(&slots), 8);
        for (uint64_t slot = 0; slot < slots; ++slot)
        {
            out.write(reinterpret_cast<const char*>(&zero), 8);
            out.write(reinterpret_cast<const char*>(&zero), 8);
        }
        for (uint64_t slot = 0; slot < slots; ++slot)
        {
            const auto position = static_cast<uint64_t>(out.tellp());
            const std::string padding((8 - position % 8) % 8, '\0');
            out.write(padding.data(), padding.size());
            const uint64_t start = static_cast<uint64_t>(out.tellp());
            writeSlot(out, slot);
            const uint64_t end = static_cast<uint64_t>(out.tellp());
            const uint64_t bytes = end - start;
            out.seekp(16 + slot * 16);
            out.write(reinterpret_cast<const char*>(&start), 8);
            out.write(reinterpret_cast<const char*>(&bytes), 8);
            out.seekp(end);
        }
        out.close();
        if (link(temporary.data(), path.c_str()) != 0)
            throw std::runtime_error("Cannot publish archive (destination may already exist): " + path);
        std::remove(temporary.data());
    }
    catch (...)
    {
        std::remove(temporary.data());
        throw;
    }
}

struct MiniPath
{
    bool operator==(const MiniPath& other) const
    {
        return (start_position == other.start_position) and (node_id == other.node_id);
    }

    using pos_t = uint16_t;
    using node_t = uint16_t;
    pos_t start_position;
    node_t node_id;
};

using KmerKey_t = uint32_t;

typedef spp::sparse_hash_map<KmerKey_t, MiniPath> StringToMiniPathsMap;
typedef spp::sparse_hash_map<KmerKey_t, std::vector<Path>> StringToPathsMap;

struct KmerIndex::Impl
{
    explicit Impl(
        const Graph& graph, size_t kmer_len_, bool membershipOnly, std::shared_ptr<const KmerIndexArchive> archive,
        uint64_t slot);
    Impl(const Impl& other)
        : graph_(other.graph_)
        , kmer_len(other.kmer_len)
        , kmer_coder(other.kmer_len)
        , membershipOnly(other.membershipOnly)
        , storage(other.storage)
        , archive(other.archive)
    {
        bind(archive ? other.data : storage.data(), other.dataSize);
    }
    void addKmerPathsStartingAtNode(NodeId node_id);
    void addKmerPaths(const std::list<Path>& kmer_paths);
    void freeze();
    void bind(const char* bytes, size_t size);
    uint32_t key(uint32_t index) const
    {
        return readValue<uint32_t>(data + entriesOffset + size_t(index) * sizeof(uint32_t));
    }
    Entry entry(uint32_t index) const
    {
        if (membershipOnly)
            return { key(index), 0, 0 };
        const auto payload = readValue<uint32_t>(data + entriesOffset + (size_t(entryCount) + index) * sizeof(uint32_t));
        const uint32_t count = index < miniCount
            ? 0
            : readValue<uint32_t>(data + countsOffset + size_t(index - miniCount) * sizeof(uint32_t));
        return { key(index), payload, count };
    }
    template <typename Bucket> BOOST_FORCEINLINE uint32_t findInTable(uint32_t key) const;
    BOOST_FORCEINLINE uint32_t find(uint32_t key) const;
    Path path(uint32_t entryIndex, size_t pathIndex) const;
    Path MiniPathToPath(const MiniPath& miniPath) const;

public:
    const Graph& graph_;

public:
    size_t kmer_len;
    TwoBitKmerEncoder kmer_coder;
    bool membershipOnly;
    std::string storage;
    std::shared_ptr<const KmerIndexArchive> archive;
    const char* data = nullptr;
    size_t dataSize = 0, entriesOffset = 0, countsOffset = 0, pathsOffset = 0, nodesOffset = 0;
    uint32_t entryCount = 0, bucketCount = 0, pathCount = 0, nodeCount = 0, miniCount = 0;
    struct Builder
    {
        StringToMiniPathsMap minis;
        StringToPathsMap paths;
        std::vector<uint32_t> keys;
    };
    std::unique_ptr<Builder> builder;
};

KmerIndex::Impl::Impl(
    const Graph& graph, size_t kmer_len_, bool membership, std::shared_ptr<const KmerIndexArchive> mappedArchive,
    uint64_t slot)
    : graph_(graph)
    , kmer_len(kmer_len_)
    , kmer_coder(kmer_len)
    , membershipOnly(membership)
    , archive(std::move(mappedArchive))
{
    if (archive)
    {
        size_t size;
        const char* bytes = archive->record(slot, size);
        const std::string identity = graphIdentity(graph);
        if (size < 16 || readValue<uint64_t>(bytes) != identity.size() || identity.size() > size - 16
            || std::memcmp(bytes + 16, identity.data(), identity.size()) != 0
            || readValue<uint32_t>(bytes + 8) != kmer_len
            || readValue<uint32_t>(bytes + 12) != uint32_t(membershipOnly))
            throw std::runtime_error("K-mer archive graph, catalog order, or parameters mismatch");
        bind(bytes + 16 + identity.size(), size - 16 - identity.size());
        return;
    }
    builder.reset(new Builder);
    for (NodeId node_id = 0; node_id != graph.numNodes(); ++node_id)
    {
        addKmerPathsStartingAtNode(node_id);
    }
    freeze();
}

void KmerIndex::Impl::addKmerPathsStartingAtNode(NodeId node_id)
{
    const string node_seq = graph_.nodeSeq(node_id);
    vector<NodeId> node_list;
    node_list.push_back(node_id);
    for (size_t pos = 0; pos != node_seq.length(); ++pos)
    {
        Path path(&graph_, static_cast<int32_t>(pos), node_list, static_cast<int32_t>(pos));
        addKmerPaths(extendPathEnd(path, static_cast<int32_t>(kmer_len)));
    }
}

Path KmerIndex::Impl::MiniPathToPath(const MiniPath& miniPath) const
{
    return Path(&graph_, miniPath.start_position, { miniPath.node_id }, miniPath.start_position + kmer_len);
}

void KmerIndex::Impl::addKmerPaths(const list<Path>& kmer_paths)
{
    auto& kmer_to_minipaths_map = builder->minis;
    auto& kmer_to_paths_map = builder->paths;
    auto& membershipKeys = builder->keys;
    for (const Path& kmer_path : kmer_paths)
    {
        vector<string> expanded_sequences;
        expandReferenceSequence(kmer_path.seq(), expanded_sequences);
        for (const auto& expanded_kmer_seq : expanded_sequences)
        {
            const auto expanded_kmer_key = kmer_coder.encode(expanded_kmer_seq);
            if (membershipOnly)
            {
                membershipKeys.push_back(expanded_kmer_key);
                continue;
            }
            auto pathIter(kmer_to_paths_map.find(expanded_kmer_key));
            if (pathIter != kmer_to_paths_map.end())
            {
                pathIter->second.push_back(kmer_path);
                continue;
            }

            auto miniPathIter(kmer_to_minipaths_map.find(expanded_kmer_key));
            if (miniPathIter != kmer_to_minipaths_map.end())
            {
                kmer_to_paths_map[expanded_kmer_key] = { MiniPathToPath(miniPathIter->second), kmer_path };
                kmer_to_minipaths_map.erase(miniPathIter);
            }
            else
            {
                if ((kmer_path.numNodes() == 1)
                    and ((kmer_path.endPosition() - kmer_path.startPosition()) == static_cast<int>(kmer_len))
                    and (kmer_path.startPosition() >= 0)
                    and (kmer_path.startPosition() <= std::numeric_limits<MiniPath::pos_t>::max())
                    and (kmer_path.nodeIds().front() <= std::numeric_limits<MiniPath::node_t>::max()))
                {
                    MiniPath miniPath { static_cast<MiniPath::pos_t>(kmer_path.startPosition()),
                                        static_cast<MiniPath::node_t>(kmer_path.nodeIds().front()) };
                    kmer_to_minipaths_map.emplace(expanded_kmer_key, miniPath);
                    continue;
                }
                else
                {
                    kmer_to_paths_map[expanded_kmer_key] = { kmer_path };
                }
            }
        }
    }
}

void KmerIndex::Impl::freeze()
{
    auto& kmer_to_minipaths_map = builder->minis;
    auto& kmer_to_paths_map = builder->paths;
    auto& membershipKeys = builder->keys;
    std::vector<Entry> entries;
    std::vector<StoredPath> paths;
    std::vector<NodeId> nodes;
    if (membershipOnly)
    {
        std::sort(membershipKeys.begin(), membershipKeys.end());
        membershipKeys.erase(std::unique(membershipKeys.begin(), membershipKeys.end()), membershipKeys.end());
        for (auto key : membershipKeys)
            entries.push_back({ key, 0, 0 });
    }
    else
    {
        for (const auto& item : kmer_to_minipaths_map)
            entries.push_back(
                { item.first, uint32_t(item.second.start_position) | (uint32_t(item.second.node_id) << 16), 0 });
        for (const auto& item : kmer_to_paths_map)
        {
            entries.push_back({ item.first, narrow(paths.size()), narrow(item.second.size()) });
            for (const auto& path : item.second)
            {
                paths.push_back(
                    { path.startPosition(), path.endPosition(), narrow(nodes.size()), narrow(path.nodeIds().size()) });
                nodes.insert(nodes.end(), path.nodeIds().begin(), path.nodeIds().end());
            }
        }
        // Keep the previous minipath-map then full-path-map traversal order,
        // including the observable diagnostic encode() order.
    }
    uint32_t buckets = 1;
    while (uint64_t(buckets) * 3 < uint64_t(entries.size()) * 4)
    {
        if (buckets > std::numeric_limits<uint32_t>::max() / 2)
            throw std::overflow_error("K-mer hash table too large");
        buckets *= 2;
    }
    std::vector<uint32_t> table(buckets, 0);
    for (size_t i = 0; i < entries.size(); ++i)
    {
        uint32_t bucket = hashKey(entries[i].key) & (buckets - 1);
        while (table[bucket])
            bucket = (bucket + 1) & (buckets - 1);
        table[bucket] = narrow(i + 1);
    }
    const size_t bucketWidth = entries.size() <= std::numeric_limits<uint16_t>::max() ? 2 : 4;
    const size_t countBytes = membershipOnly ? 0 : (entries.size() - kmer_to_minipaths_map.size()) * sizeof(uint32_t);
    storage.resize(
        indexHeaderSize + table.size() * bucketWidth + entries.size() * (membershipOnly ? 4 : 8)
        + countBytes + paths.size() * sizeof(StoredPath) + nodes.size() * sizeof(NodeId));
    char* cursor = &storage[0];
    writeValue(cursor, narrow(entries.size()));
    writeValue(cursor, buckets);
    writeValue(cursor, narrow(paths.size()));
    writeValue(cursor, narrow(nodes.size()));
    writeValue(cursor, membershipOnly ? uint32_t(0) : narrow(kmer_to_minipaths_map.size()));
    // Bucket values are one-based entry indexes; zero is the empty sentinel.
    // Preserve the full 32-bit representation when any value would not fit.
    if (entries.size() <= std::numeric_limits<uint16_t>::max())
        for (auto index : table)
            writeValue(cursor, static_cast<uint16_t>(index));
    else
    {
        std::memcpy(cursor, table.data(), table.size() * sizeof(uint32_t));
        cursor += table.size() * sizeof(uint32_t);
    }
    // Keep keys contiguous in both index types: probing never loads path
    // payloads and uses the same fixed stride for membership and alignment.
    for (const auto& item : entries)
        writeValue(cursor, item.key);
    if (!membershipOnly)
    {
        // freeze() already places minipaths first. Their count is implicitly
        // one path (logical count zero), so only full entries need a count.
        for (const auto& item : entries)
            writeValue(cursor, item.payload);
        for (size_t i = kmer_to_minipaths_map.size(); i < entries.size(); ++i)
            writeValue(cursor, entries[i].count);
    }
    if (!paths.empty())
    {
        std::memcpy(cursor, paths.data(), paths.size() * sizeof(StoredPath));
        cursor += paths.size() * sizeof(StoredPath);
    }
    if (!nodes.empty())
        std::memcpy(cursor, nodes.data(), nodes.size() * sizeof(NodeId));
    builder.reset();
    storage.shrink_to_fit();
    bind(storage.data(), storage.size());
}

void KmerIndex::Impl::bind(const char* bytes, size_t size)
{
    if (size < indexHeaderSize)
        throw std::runtime_error("Truncated compact index");
    data = bytes;
    dataSize = size;
    entryCount = readValue<uint32_t>(bytes);
    bucketCount = readValue<uint32_t>(bytes + 4);
    pathCount = readValue<uint32_t>(bytes + 8);
    nodeCount = readValue<uint32_t>(bytes + 12);
    miniCount = readValue<uint32_t>(bytes + 16);
    if (miniCount > entryCount || (membershipOnly && miniCount != 0))
        throw std::runtime_error("Invalid compact index minipath count");
    const size_t bucketWidth = entryCount <= std::numeric_limits<uint16_t>::max() ? 2 : 4;
    entriesOffset = indexHeaderSize + uint64_t(bucketCount) * bucketWidth;
    countsOffset = entriesOffset + uint64_t(entryCount) * (membershipOnly ? 4 : 8);
    pathsOffset = countsOffset + (membershipOnly ? 0 : uint64_t(entryCount - miniCount) * sizeof(uint32_t));
    nodesOffset = pathsOffset + uint64_t(pathCount) * sizeof(StoredPath);
    if (!bucketCount || (bucketCount & (bucketCount - 1)) || entryCount >= bucketCount
        || nodesOffset + uint64_t(nodeCount) * sizeof(NodeId) != size)
        throw std::runtime_error("Invalid compact index layout");
}

BOOST_FORCEINLINE uint32_t KmerIndex::Impl::find(uint32_t key) const
{
    // Dispatch outside the probe loop; misses and collisions only read keys,
    // never the separate full-path counts.
    return entryCount <= std::numeric_limits<uint16_t>::max() ? findInTable<uint16_t>(key) : findInTable<uint32_t>(key);
}

template <typename Bucket> BOOST_FORCEINLINE uint32_t KmerIndex::Impl::findInTable(uint32_t query) const
{
    uint32_t bucket = hashKey(query) & (bucketCount - 1);
    const char* const keys = data + entriesOffset;
    for (uint32_t probe = 0; probe < bucketCount; ++probe)
    {
        const uint32_t index = readValue<Bucket>(data + indexHeaderSize + size_t(bucket) * sizeof(Bucket));
        if (!index)
            return entryCount;
        if (index > entryCount)
            throw std::runtime_error("Invalid compact index bucket");
        if (readValue<uint32_t>(keys + size_t(index - 1) * sizeof(uint32_t)) == query)
            return index - 1;
        bucket = (bucket + 1) & (bucketCount - 1);
    }
    throw std::runtime_error("Compact index has no empty bucket");
}

Path KmerIndex::Impl::path(uint32_t entryIndex, size_t pathIndex) const
{
    if (membershipOnly)
        throw std::logic_error("Membership index does not store paths");
    const auto item = entry(entryIndex);
    if (!item.count)
        return MiniPathToPath({ uint16_t(item.payload), uint16_t(item.payload >> 16) });
    if (pathIndex >= item.count || uint64_t(item.payload) + item.count > pathCount)
        throw std::runtime_error("Invalid compact path range");
    const auto stored = readValue<StoredPath>(data + pathsOffset + (item.payload + pathIndex) * sizeof(StoredPath));
    if (uint64_t(stored.nodes) + stored.count > nodeCount)
        throw std::runtime_error("Invalid compact node range");
    std::vector<NodeId> pathNodes(stored.count);
    std::memcpy(pathNodes.data(), data + nodesOffset + stored.nodes * sizeof(NodeId), stored.count * sizeof(NodeId));
    return Path(&graph_, stored.start, std::move(pathNodes), stored.end);
}

KmerIndex::KmerIndex(
    const Graph& graph, int32_t kmer_len, bool membershipOnly, std::shared_ptr<const KmerIndexArchive> archive,
    uint64_t archiveSlot)
    : pimpl_(new Impl(graph, static_cast<size_t>(kmer_len), membershipOnly, std::move(archive), archiveSlot))
{
}

KmerIndex::KmerIndex(const KmerIndex& other)
    : pimpl_(new Impl(*other.pimpl_))
{
    if (!pimpl_->archive)
        pimpl_->bind(pimpl_->storage.data(), pimpl_->storage.size());
}

KmerIndex::KmerIndex(KmerIndex&& other) noexcept
    : pimpl_(std::move(other.pimpl_))
{
}

KmerIndex& KmerIndex::operator=(const KmerIndex& other)
{
    if (this != &other)
    {
        pimpl_.reset(new Impl(*other.pimpl_));
        if (!pimpl_->archive)
            pimpl_->bind(pimpl_->storage.data(), pimpl_->storage.size());
    }
    return *this;
}

KmerIndex& KmerIndex::operator=(KmerIndex&& other) noexcept
{
    pimpl_ = std::move(other.pimpl_);
    return *this;
}

KmerIndex::~KmerIndex() = default;

bool KmerIndex::operator==(const KmerIndex& other) const
{
    if (pimpl_->kmer_len != other.pimpl_->kmer_len || pimpl_->membershipOnly != other.pimpl_->membershipOnly
        || kmers() != other.kmers())
        return false;
    if (pimpl_->membershipOnly)
        return true;
    for (const auto& kmer : kmers())
        if (getPaths(kmer) != other.getPaths(kmer))
            return false;
    return true;
}

static string encodePaths(const std::vector<Path>& paths)
{
    std::vector<string> path_encodings;
    for (const auto& path : paths)
    {
        path_encodings.push_back(path.encode());
    }
    return boost::algorithm::join(path_encodings, ",");
}

size_t KmerIndex::kmerLength() const { return pimpl_->kmer_len; }

string KmerIndex::encode() const
{
    std::vector<string> kv_encodings;
    for (uint32_t i = 0; i < pimpl_->entryCount; ++i)
    {
        const string kmer = pimpl_->kmer_coder.decode(pimpl_->entry(i).key);
        kv_encodings.emplace_back(
            "{" + kmer + "->" + (pimpl_->membershipOnly ? "present" : encodePaths(getPaths(kmer))) + "}");
    }
    return boost::algorithm::join(kv_encodings, ",");
}

bool KmerIndex::contains(const std::string& kmer) const
{
    if ((kmer.size() != pimpl_->kmer_len) or (kmer.find_last_not_of("ACGT") != std::string::npos))
    {
        return false;
    }

    return containsEncodedKmer(pimpl_->kmer_coder.encode(kmer));
}

bool KmerIndex::containsEncodedKmer(uint32_t key) const { return pimpl_->find(key) != pimpl_->entryCount; }

size_t KmerIndex::numPaths(const std::string& kmer) const
{
    if ((kmer.size() != pimpl_->kmer_len) or (kmer.find_last_not_of("ACGT") != std::string::npos))
    {
        return 0;
    }

    return findEncodedPaths(pimpl_->kmer_coder.encode(kmer)).size();
}

std::vector<Path> KmerIndex::getPaths(const std::string& kmer) const
{
    const auto view = findEncodedPaths(pimpl_->kmer_coder.encode(kmer));
    if (!view.size())
        throw std::out_of_range("K-mer is absent");
    std::vector<Path> paths;
    paths.reserve(view.size());
    for (size_t i = 0; i < view.size(); ++i)
        paths.push_back(view[i]);
    return paths;
}

size_t KmerIndex::PathView::size() const
{
    if (!found_)
        return 0;
    if (index_->pimpl_->membershipOnly)
        throw std::logic_error("Membership index does not count paths");
    return std::max(uint32_t(1), index_->pimpl_->entry(entry_).count);
}

Path KmerIndex::PathView::operator[](size_t path_index) const
{
    if (path_index >= size())
        throw std::out_of_range("K-mer path index");
    return index_->pimpl_->path(entry_, path_index);
}

KmerIndex::PathView KmerIndex::findEncodedPaths(uint32_t kmer_key) const
{
    PathView view;
    view.index_ = this;
    view.entry_ = pimpl_->find(kmer_key);
    view.found_ = view.entry_ != pimpl_->entryCount;
    return view;
}

unordered_set<string> KmerIndex::kmers() const
{
    unordered_set<string> kmers;
    for (uint32_t i = 0; i < pimpl_->entryCount; ++i)
        kmers.insert(pimpl_->kmer_coder.decode(pimpl_->entry(i).key));
    return kmers;
}

size_t KmerIndex::numUniqueKmersOverlappingNode(NodeId node_id) const
{
    size_t count = 0;
    for (uint32_t i = 0; i < pimpl_->entryCount; ++i)
    {
        if (pimpl_->entry(i).count > 1)
            continue;
        const auto path = pimpl_->path(i, 0);
        for (auto node : path.nodeIds())
            if (node == node_id)
                ++count;
    }
    return count;
}

size_t KmerIndex::numUniqueKmersOverlappingEdge(NodeId from, NodeId to) const
{
    size_t count = 0;
    for (uint32_t i = 0; i < pimpl_->entryCount; ++i)
    {
        if (pimpl_->entry(i).count > 1)
            continue;
        const auto path = pimpl_->path(i, 0);
        const auto& nodes = path.nodeIds();
        for (size_t j = 1; j < nodes.size(); ++j)
            if (nodes[j - 1] == from && nodes[j] == to)
                ++count;
    }
    return count;
}

void KmerIndex::writeRecord(std::ostream& out) const
{
    const auto identity = graphIdentity(pimpl_->graph_);
    std::string record;
    appendValue(record, uint64_t(identity.size()));
    appendValue(record, uint32_t(pimpl_->kmer_len));
    appendValue(record, uint32_t(pimpl_->membershipOnly));
    record += identity;
    record.append(pimpl_->data, pimpl_->dataSize);
    writeBlobRecord(out, record);
}

void KmerIndex::writeBlobRecord(std::ostream& out, const std::string& record)
{
    const auto digest = checksum(record.data(), record.size());
    out.write(digest.data(), digest.size());
    const uint32_t padding = 0;
    out.write(reinterpret_cast<const char*>(&padding), sizeof(padding));
    out.write(record.data(), record.size());
}

std::ostream& operator<<(std::ostream& os, const KmerIndex& kmer_index)
{
    os << kmer_index.encode();
    return os;
}
}
