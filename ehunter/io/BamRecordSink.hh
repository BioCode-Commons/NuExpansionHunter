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
// Owned BAM batches with a dedicated output consumer and byte backpressure.
#pragma once

#include <boost/container/vector.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <future>
#include <memory>
#include <vector>

#include "htslib/sam.h"

namespace ehunter
{

// Contract [K-36]
// Ensures:  batches are written first-in first-out; resident bytes never exceed the budget plus one batch.
class BamRecordSink
{
public:
    struct DestroyRecord
    {
        void operator()(bam1_t* record) const noexcept { bam_destroy1(record); }
    };
    using Record = std::unique_ptr<bam1_t, DestroyRecord>;
    using Records = std::vector<Record>;
    using WriteRecord = std::function<void(const bam1_t*)>;
    using WriteBytes = std::function<void(const uint8_t*, std::size_t)>;
    using Bytes = boost::container::vector<uint8_t, std::allocator<uint8_t>>;
    using SetPreparedCapacity = std::function<void(size_t)>;
    using PrepareAccounted = std::function<std::future<Bytes>(Bytes, SetPreparedCapacity)>;
    enum class WriteErrorPolicy
    {
        kThrowImmediately,
        kDeferUntilFinish
    };
    static constexpr std::size_t defaultBudgetBytes = 8 * 1024 * 1024;

    // Callback runs on one dedicated consumer thread. The callback and any
    // borrowed file/header must remain valid until finish() returns or throws.
    // Zero explicitly selects an unbounded sink. A batch larger than a nonzero
    // budget is allowed only when no other batch is resident.
    explicit BamRecordSink(WriteRecord writeRecord, std::size_t maxBytes = defaultBudgetBytes,
        WriteBytes writeBytes = WriteBytes(), PrepareAccounted prepareAccounted = PrepareAccounted());
    ~BamRecordSink();
    BamRecordSink(const BamRecordSink&) = delete;
    BamRecordSink& operator=(const BamRecordSink&) = delete;

    // Ownership is consumed on success and failure. Producers can block for
    // credits, but the consumer never requires a producer's compute executor.
    // Deferred mode preserves seeking compatibility: after a consumer write
    // failure, release submissions without throwing until explicit finish().
    // Validation/allocation errors and submissions after finish still throw.
    void submit(Records records, WriteErrorPolicy policy = WriteErrorPolicy::kThrowImmediately);
    void submitEncoded(Bytes data, std::size_t recordCount, size_t workspaceBytes = 0);
    // Close admission, drain accepted batches, join, and propagate write errors.
    // Concurrent/repeated finish calls are supported. Producers should normally
    // be joined first; a producer waiting during finish is awakened and rejected.
    void finish();

    static std::size_t recordBytes(const bam1_t& record);
    // Includes queued and consumer-owned batches, excluding producer-local
    // records, HTS compression buffers, and allocator/container overhead.
    std::size_t residentBytes() const;

private:
    void enqueue(Records records, Bytes data, std::size_t count, std::size_t bytes, WriteErrorPolicy policy);
    struct State;
    std::shared_ptr<State> state_;
};

}
