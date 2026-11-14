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

#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <unordered_set>

#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

namespace graphtools
{

// A single read-only mapped archive shared by all locus indexes in a process.
// The preparation process writes slots in catalog order before publishing it.
class KmerIndexArchive
{
public:
    explicit KmerIndexArchive(const std::string& path);
    ~KmerIndexArchive();
    KmerIndexArchive(const KmerIndexArchive&) = delete;
    KmerIndexArchive& operator=(const KmerIndexArchive&) = delete;
    const char* record(uint64_t slot, size_t& size) const;
    uint64_t slotCount() const;

private:
    const char* data_ = nullptr;
    size_t size_ = 0;
};

// Kmer index holds paths that correspond to each kmer that appears in the graph and supports a few standard operations.
class KmerIndex
{
public:
    explicit KmerIndex(
        const Graph& graph, int32_t kmer_len = 12, bool membershipOnly = false,
        std::shared_ptr<const KmerIndexArchive> archive = {}, uint64_t archiveSlot = 0);
    void writeRecord(std::ostream& out) const;
    static void writeBlobRecord(std::ostream& out, const std::string& bytes);
    static void prepareArchive(
        const std::string& path, uint64_t slots, const std::function<void(std::ostream&, uint64_t)>& writeSlot);
    KmerIndex(const KmerIndex& other);
    KmerIndex(KmerIndex&& other) noexcept;
    KmerIndex& operator=(const KmerIndex& other);
    KmerIndex& operator=(KmerIndex&& other) noexcept;
    ~KmerIndex();
    bool operator==(const KmerIndex& other) const;
    std::string encode() const;

    /// \brief Return all paths for kmer
    ///
    /// Calling this with kmer that isn't contained in the index triggers a runtime error
    std::vector<Path> getPaths(const std::string& kmer) const;

    // A non-owning view of paths for one encoded kmer. The index must outlive the
    // view and remain unchanged. Each requested Path is materialized separately,
    // so a multi-path lookup does not copy the whole path vector.
    class PathView
    {
    public:
        size_t size() const;
        Path operator[](size_t path_index) const;

    private:
        friend class KmerIndex;
        const KmerIndex* index_ = nullptr;
        uint32_t entry_ = 0;
        bool found_ = false;
        bool mini_ = false;
    };

    // The caller supplies an A/C/G/T key encoded at this index's kmer length.
    PathView findEncodedPaths(uint32_t kmer_key) const;

    bool contains(const std::string& kmer) const;
    /// Query a valid k-base A/C/G/T key encoded with TwoBitKmerEncoder.
    bool containsEncodedKmer(uint32_t key) const;
    size_t numPaths(const std::string& kmer) const;
    std::unordered_set<std::string> kmers() const;
    size_t kmerLength() const;

    size_t numUniqueKmersOverlappingNode(NodeId node_id) const;
    size_t numUniqueKmersOverlappingEdge(NodeId from, NodeId to) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

std::ostream& operator<<(std::ostream& os, const KmerIndex& kmer_index);
}
