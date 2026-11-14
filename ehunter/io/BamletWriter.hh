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

#pragma once

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

namespace ehunter
{

/// Thread-safe BAM serialization with optional producer-local batch scopes.
/// A separate owned sink performs output independently of compute workers.
///
class BamletWriter : public graphtools::AlignmentWriter, private boost::noncopyable
{
public:
    BamletWriter(
        const std::string& bamletPath, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
        int compressionThreads, std::size_t maxQueueBytes = BamRecordSink::defaultBudgetBytes);
    ~BamletWriter() override;

    // Called after all producers finish; propagate asynchronous write/close failure.
    void finish();

    // Collect writes made by this thread to this writer. Explicit finish()
    // submits the remaining records and can throw. Destruction discards only
    // unsubmitted records and restores the previous scope without throwing.
    // nullptr is a no-op. Scopes must be finished/destroyed in stack order and
    // the writer must outlive every scope. Large scopes submit bounded pieces.
    class Batch
    {
    public:
        explicit Batch(BamletWriter* writer);
        ~Batch() noexcept;
        Batch(const Batch&) = delete;
        Batch& operator=(const Batch&) = delete;
        void finish();

    private:
        friend class BamletWriter;
        void append(BamRecordSink::Record record);
        void submit();
        BamletWriter* writer_;
        Batch* previous_ = nullptr;
        BamRecordSink::Records records_;
        std::size_t bytes_ = 0;
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

    std::unique_ptr<htsFile, decltype(&hts_close)> filePtr_;
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> bamHeader_;
    ReferenceContigInfo contigInfo_;
    std::vector<int32_t> bamContigIds_;

    const RegionCatalog& catalog_;

    std::unique_ptr<BamRecordSink> sink_;
    std::exception_ptr finishError_;
    bool sinkStatsLogged_ = false;
};

}
