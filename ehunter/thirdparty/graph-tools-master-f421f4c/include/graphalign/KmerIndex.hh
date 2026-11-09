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

#include <memory>
#include <cstdint>
#include <string>
#include <unordered_set>

#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

namespace graphtools
{

// Kmer index holds paths that correspond to each kmer that appears in the graph and supports a few standard operations.
class KmerIndex
{
public:
    explicit KmerIndex(const Graph& graph, int32_t kmer_len = 12);
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
        const std::vector<Path>* paths_ = nullptr;
        uint16_t mini_start_position_ = 0;
        uint16_t mini_node_id_ = 0;
        bool is_minipath_ = false;
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
