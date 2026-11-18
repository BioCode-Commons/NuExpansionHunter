//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>,
//         Chris Saunders <csaunders@illumina.com>
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
//
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: BAM serialization uses producer batches, asynchronous sinks, optional codecs, and explicit completion.
// Modifications are licensed under the Apache License, Version 2.0.
//

#pragma once

#include <array>
#include <exception>
#include <memory>
#include <string>
#include <unordered_map>

#include "boost/noncopyable.hpp"
// cppcheck-suppress missingInclude
#include "htslib/hts.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"

#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"
#include "io/BamRecordSink.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphcore/GraphReferenceMapping.hh"
#include "graphio/AlignmentWriter.hh"
#include "locus/LocusSpecification.hh"

namespace ctpl { class thread_pool; }

namespace ehunter
{

#ifdef __linux__
class PositionalBamOutput;
#endif

/// Thread-safe BAM serialization with optional producer-local batch scopes.
/// A separate owned sink performs output independently of compute workers.
///
class BamletWriter : public graphtools::AlignmentWriter, private boost::noncopyable
{
public:
    BamletWriter(
        const std::string& bamletPath, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
        int compressionThreads, std::size_t maxQueueBytes = BamRecordSink::defaultBudgetBytes,
        bool fastCompression = false, unsigned validationThreads = 1);
    ~BamletWriter() override;

    // Called after all producers finish; propagate asynchronous write/close failure.
    void finish();

    class Batch;
    // Owned by one compute worker. After producers have joined, flush its
    // final partial batch before closing the sink.
    class BatchBuffer
    {
    public:
        explicit BatchBuffer(BamletWriter* owner) : owner_(owner) {}
        BatchBuffer(const BatchBuffer&) = delete;
        BatchBuffer& operator=(const BatchBuffer&) = delete;
        BatchBuffer(BatchBuffer&&) = default;
        void flush();

    private:
        friend class Batch;
        friend class BamletWriter;
        BamletWriter* owner_;
        BamRecordSink::Bytes data;
        size_t count = 0;
        std::string lastId, alignmentEncoding;
        size_t lastIndex = 0, nodeCount = 0;
        bool hasIndex = false;
        uint8_t cachedNodes = 0;
        std::array<boost::optional<GenomicRegion>, 8> projections;
    };

    // Collect writes made by this thread to this writer. Explicit finish()
    // submits the remaining records and can throw. Destruction discards only
    // unsubmitted records and restores the previous scope without throwing.
    // nullptr is a no-op. Scopes must be finished/destroyed in stack order and
    // the writer must outlive every scope. Large scopes submit bounded pieces.
    class Batch
    {
    public:
        explicit Batch(BamletWriter* writer, BatchBuffer* buffer = nullptr);
        ~Batch() noexcept;
        Batch(const Batch&) = delete;
        Batch& operator=(const Batch&) = delete;
        void finish();

    private:
        friend class BamletWriter;
        uint8_t* append(size_t recordBytes);
        size_t locusIndex(graphtools::TextView id);
        BamletWriter* writer_;
        Batch* previous_ = nullptr;
        BatchBuffer owned_;
        BatchBuffer* buffer_;
        boost::optional<GenomicRegion> projection(graphtools::NodeId node) const;
        size_t nodeCount() const { return buffer_->nodeCount; }
        bool finished_ = false;
    };

    /// Thread safe
    void write(
        graphtools::TextView locusId, const std::string& fragmentName, const std::string& query, bool isFirstMate,
        bool isReversed, bool isMateReversed, const graphtools::GraphAlignment& alignment) override;

private:
    void writeHeader();

    /// Thread safe
    void write(
        int32_t contigId, int64_t position, const std::string& fragmentName, const std::string& query,
        bool isFirstMate, bool isReversed, bool isMateReversed, const graphtools::GraphAlignment& alignment);

    void submitRecord(BamRecordSink::Record record);
    static thread_local Batch* activeBatch_;

    int fileDescriptor_ = -1; // Borrowed from filePtr_ until the checked close.
    std::unique_ptr<htsFile, decltype(&hts_close)> filePtr_;
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> bamHeader_;
    ReferenceContigInfo contigInfo_;
    std::vector<int32_t> bamContigIds_;

    const RegionCatalog& catalog_;

    std::unique_ptr<ctpl::thread_pool> encoderPool_;
#ifdef __linux__
    std::unique_ptr<PositionalBamOutput> positional_;
#endif
    std::unique_ptr<BamRecordSink> sink_;
    std::exception_ptr finishError_;
};

}
