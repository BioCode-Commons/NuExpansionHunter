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
#include "io/BamRecordSink.hh"
#include "io/BamletWriter.hh"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <future>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <unistd.h>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphBuilders.hh"
#include "gtest/gtest.h"

namespace
{

using ehunter::BamRecordSink;
using ehunter::BamletWriter;

BamRecordSink::Records records(std::initializer_list<int> positions, std::size_t payload = 0)
{
    BamRecordSink::Records result;
    result.reserve(positions.size());
    for (const int position : positions)
    {
        BamRecordSink::Record record(bam_init1());
        if (!record)
            throw std::bad_alloc();
        record->core.pos = position;
        if (payload)
        {
            record->data = static_cast<uint8_t*>(std::malloc(payload));
            if (!record->data)
                throw std::bad_alloc();
            record->m_data = payload;
        }
        result.push_back(std::move(record));
    }
    return result;
}

class WriteGate
{
public:
    void enter()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        entered_ = true;
        condition_.notify_all();
        condition_.wait(lock, [&] { return released_; });
    }
    bool waitUntilEntered()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        return condition_.wait_for(lock, std::chrono::seconds(5), [&] { return entered_; });
    }
    void release()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        released_ = true;
        condition_.notify_all();
    }

private:
    std::mutex mutex_;
    std::condition_variable condition_;
    bool entered_ = false, released_ = false;
};

TEST(BamRecordSink, SynchronousPreparationCanReportItsCapacityBeforeReturningTheFuture)
{
    std::vector<uint8_t> output;
    BamRecordSink::PrepareAccounted prepare = [](BamRecordSink::Bytes raw, BamRecordSink::SetPreparedCapacity release) {
        BamRecordSink::Bytes encoded {raw.front()};
        BamRecordSink::Bytes().swap(raw);
        release(encoded.capacity());
        std::promise<BamRecordSink::Bytes> promise;
        promise.set_value(std::move(encoded));
        return promise.get_future();
    };
    BamRecordSink sink([](const bam1_t*) {}, 65, [&](const uint8_t* data, size_t n) {
        output.insert(output.end(), data, data+n);
    }, prepare);
    sink.submitEncoded(BamRecordSink::Bytes(64, 42), 1, 1);
    sink.finish();
    EXPECT_EQ(output, (std::vector<uint8_t> {42}));
    EXPECT_EQ(sink.residentBytes(), 0u);
}

TEST(BamRecordSink, CompletedCompressionReturnsUnusedCapacityWhileFileWriterIsBlocked)
{
    std::promise<void> started, resume;
    auto gate = resume.get_future().share();
    std::vector<uint8_t> output;
    BamRecordSink::PrepareAccounted prepare = [](BamRecordSink::Bytes raw, BamRecordSink::SetPreparedCapacity release) {
        return std::async(std::launch::async, [](BamRecordSink::Bytes raw, BamRecordSink::SetPreparedCapacity release) {
            BamRecordSink::Bytes encoded {raw.front()};
            BamRecordSink::Bytes().swap(raw);
            release(encoded.capacity());
            return encoded;
        }, std::move(raw), std::move(release));
    };
    BamRecordSink sink([](const bam1_t*) {}, 66, [&](const uint8_t* data, size_t) {
        if (output.empty()) { started.set_value(); gate.wait(); }
        output.push_back(*data);
    }, prepare);
    sink.submitEncoded(BamRecordSink::Bytes(64, 1), 1, 1);
    started.get_future().wait();
    auto second = std::async(std::launch::async, [&] { sink.submitEncoded(BamRecordSink::Bytes(64, 2), 1, 1); });
    const auto status = second.wait_for(std::chrono::seconds(5));
    resume.set_value();
    second.get();
    sink.finish();
    EXPECT_EQ(status, std::future_status::ready);
    EXPECT_EQ(output, (std::vector<uint8_t> {1,2}));
}

TEST(BamRecordSink, PreservesBatchAndRecordOrderAndFinishesTwice)
{
    std::vector<int> observed;
    BamRecordSink sink([&](const bam1_t* record) { observed.push_back(record->core.pos); }, 1024);
    sink.submit(records({ 1, 2, 2 }));
    sink.submit(records({ 3, 4 }));
    sink.finish();
    sink.finish();
    EXPECT_EQ(observed, (std::vector<int> { 1, 2, 2, 3, 4 }));
    EXPECT_EQ(sink.residentBytes(), 0u);
    EXPECT_THROW(sink.submit(records({ 5 })), std::logic_error);
}

TEST(BamRecordSink, EmptyBatchesAndEmptyStreamFinish)
{
    std::atomic<int> called(0);
    BamRecordSink sink([&](const bam1_t*) { ++called; });
    sink.submit({});
    sink.finish();
    sink.finish();
    EXPECT_EQ(called.load(), 0);
    EXPECT_EQ(sink.residentBytes(), 0u);
}

TEST(BamRecordSink, KeepsConsumerOwnedPayloadChargedUntilWriteCompletes)
{
    WriteGate gate;
    std::atomic<int> called(0);
    auto first = records({ 1 }, 512);
    const auto bytes = BamRecordSink::recordBytes(*first.front()) + first.capacity() * sizeof(BamRecordSink::Record);
    BamRecordSink sink([&](const bam1_t*) { ++called; gate.enter(); }, bytes);
    sink.submit(std::move(first));
    EXPECT_TRUE(gate.waitUntilEntered());
    EXPECT_EQ(sink.residentBytes(), bytes);
    auto second = std::async(std::launch::async, [&] { sink.submit(records({ 2 }, 512)); });
    EXPECT_EQ(second.wait_for(std::chrono::milliseconds(20)), std::future_status::timeout);
    gate.release();
    second.get();
    sink.finish();
    EXPECT_EQ(called.load(), 2);
}

TEST(BamRecordSink, OversizedBatchRunsExclusivelyWithoutDeadlock)
{
    WriteGate gate;
    std::vector<int> observed;
    BamRecordSink sink([&](const bam1_t* record) {
        gate.enter();
        observed.push_back(record->core.pos);
    }, 1);
    sink.submit(records({ 1, 2 }, 4096));
    EXPECT_TRUE(gate.waitUntilEntered());
    EXPECT_GT(sink.residentBytes(), 1u);
    auto second = std::async(std::launch::async, [&] { sink.submit(records({ 3 }, 16)); });
    EXPECT_EQ(second.wait_for(std::chrono::milliseconds(20)), std::future_status::timeout);
    gate.release();
    second.get();
    sink.finish();
    EXPECT_EQ(observed, (std::vector<int> { 1, 2, 3 }));
}

TEST(BamRecordSink, ExplicitZeroBudgetDoesNotBlockOnResidentPayload)
{
    WriteGate gate;
    BamRecordSink sink([&](const bam1_t*) { gate.enter(); }, 0);
    sink.submit(records({ 1 }, 4096));
    EXPECT_TRUE(gate.waitUntilEntered());
    auto second = std::async(std::launch::async, [&] { sink.submit(records({ 2 }, 4096)); });
    const auto status = second.wait_for(std::chrono::seconds(2));
    gate.release();
    second.get();
    sink.finish();
    EXPECT_EQ(status, std::future_status::ready);
}

TEST(BamRecordSink, ConsumerFailureWakesBlockedProducersAndPersistsThroughFinish)
{
    WriteGate gate;
    std::atomic<int> called(0);
    BamRecordSink sink([&](const bam1_t*) {
        ++called;
        gate.enter();
        throw std::runtime_error("synthetic BAM write failure");
    }, 1);
    sink.submit(records({ 1, 2 }, 512));
    EXPECT_TRUE(gate.waitUntilEntered());
    auto producer = std::async(std::launch::async, [&] {
        try
        {
            sink.submit(records({ 3 }, 512));
            return std::string("unexpected success");
        }
        catch (const std::exception& error)
        {
            return std::string(error.what());
        }
    });
    EXPECT_EQ(producer.wait_for(std::chrono::milliseconds(20)), std::future_status::timeout);
    gate.release();
    EXPECT_EQ(producer.get(), "synthetic BAM write failure");
    EXPECT_THROW(sink.finish(), std::runtime_error);
    EXPECT_THROW(sink.finish(), std::runtime_error);
    EXPECT_THROW(sink.submit(records({ 4 })), std::runtime_error);
    EXPECT_EQ(sink.residentBytes(), 0u);
    EXPECT_EQ(called.load(), 1);
}

TEST(BamRecordSink, DeferredWriteFailureReleasesWaitersButFinishStillThrows)
{
    WriteGate gate;
    std::atomic<int> called(0);
    BamRecordSink sink([&](const bam1_t*) {
        ++called;
        gate.enter();
        throw std::runtime_error("synthetic deferred BAM failure");
    }, 1);
    const auto deferred = BamRecordSink::WriteErrorPolicy::kDeferUntilFinish;
    sink.submit(records({ 1, 2 }, 512), deferred);
    EXPECT_TRUE(gate.waitUntilEntered());
    auto producer = std::async(std::launch::async, [&] {
        sink.submit(records({ 3 }, 512), deferred);
    });
    EXPECT_EQ(producer.wait_for(std::chrono::milliseconds(20)), std::future_status::timeout);
    gate.release();
    EXPECT_NO_THROW(producer.get());
    EXPECT_NO_THROW(sink.submit(records({ 4 }), deferred));
    EXPECT_THROW(sink.submit(records({ 5 })), std::runtime_error);

    auto invalid = records({ 6 });
    invalid.emplace_back(nullptr);
    EXPECT_THROW(sink.submit(std::move(invalid), deferred), std::invalid_argument);
    try
    {
        sink.finish();
        FAIL() << "finish must report the deferred write failure";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_STREQ(error.what(), "synthetic deferred BAM failure");
    }
    EXPECT_THROW(sink.finish(), std::runtime_error);
    EXPECT_THROW(sink.submit(records({ 7 }), deferred), std::runtime_error);
    EXPECT_EQ(sink.residentBytes(), 0u);
    EXPECT_EQ(called.load(), 1);
}

TEST(BamRecordSink, DeferredSubmissionPreservesOutputAndRejectsAfterFinish)
{
    std::vector<int> observed;
    BamRecordSink sink([&](const bam1_t* record) { observed.push_back(record->core.pos); }, 1024);
    const auto deferred = BamRecordSink::WriteErrorPolicy::kDeferUntilFinish;
    sink.submit(records({ 1, 2, 2 }), deferred);
    sink.submit(records({ 3 }));
    sink.finish();
    EXPECT_EQ(observed, (std::vector<int> { 1, 2, 2, 3 }));
    EXPECT_THROW(sink.submit(records({ 4 }), deferred), std::logic_error);
}

TEST(BamRecordSink, FinishWakesCapacityWaitersWhileDrainingAcceptedWork)
{
    WriteGate gate;
    BamRecordSink sink([&](const bam1_t*) { gate.enter(); }, 1);
    sink.submit(records({ 1 }, 512));
    EXPECT_TRUE(gate.waitUntilEntered());
    auto producer = std::async(std::launch::async, [&] {
        try
        {
            sink.submit(records({ 2 }, 512));
            return false;
        }
        catch (const std::logic_error&)
        {
            return true;
        }
    });
    auto finisher = std::async(std::launch::async, [&] { sink.finish(); });
    const auto status = producer.wait_for(std::chrono::seconds(2));
    gate.release();
    EXPECT_TRUE(producer.get());
    finisher.get();
    EXPECT_EQ(status, std::future_status::ready);
    EXPECT_EQ(sink.residentBytes(), 0u);
}

TEST(BamRecordSink, ConcurrentProducersPreserveEveryRecordAndDuplicate)
{
    std::vector<int> observed;
    BamRecordSink sink([&](const bam1_t* record) { observed.push_back(record->core.pos); }, 1024);
    std::vector<std::future<void>> producers;
    for (int producer = 0; producer < 4; ++producer)
        producers.push_back(std::async(std::launch::async, [&, producer] {
            for (int batch = 0; batch < 100; ++batch)
                sink.submit(records({ producer, producer, batch + 10 }));
        }));
    for (auto& producer : producers)
        producer.get();
    auto finishA = std::async(std::launch::async, [&] { sink.finish(); });
    auto finishB = std::async(std::launch::async, [&] { sink.finish(); });
    finishA.get();
    finishB.get();
    ASSERT_EQ(observed.size(), 1200u);
    for (int producer = 0; producer < 4; ++producer)
        EXPECT_EQ(std::count(observed.begin(), observed.end(), producer), 200);
    for (int batch = 0; batch < 100; ++batch)
        EXPECT_EQ(std::count(observed.begin(), observed.end(), batch + 10), 4);
}

TEST(BamRecordSink, ManySingleRecordProducersStayWithinBudgetAndFinish)
{
    const std::size_t payload = 256;
    const std::size_t batchBytes = BamRecordSink::recordBytes(*records({ 0 }, payload).front()) + sizeof(BamRecordSink::Record);
    const std::size_t budget = 16 * batchBytes;
    std::size_t written = 0, peak = 0;
    BamRecordSink* sinkPointer = nullptr;
    BamRecordSink sink([&](const bam1_t*) {
        ++written;
        peak = std::max(peak, sinkPointer->residentBytes());
    }, budget);
    sinkPointer = &sink;
    std::vector<std::future<void>> producers;
    for (int producer = 0; producer < 16; ++producer)
        producers.push_back(std::async(std::launch::async, [&, producer] {
            for (int record = 0; record < 500; ++record)
                sink.submit(records({ producer }, payload), BamRecordSink::WriteErrorPolicy::kDeferUntilFinish);
        }));
    for (auto& producer : producers)
        producer.get();
    sink.finish();
    EXPECT_EQ(written, 16u * 500u);
    EXPECT_LE(peak, budget + batchBytes);
    EXPECT_EQ(sink.residentBytes(), 0u);
}

TEST(BamRecordSink, RejectsNullRecordsWithoutPoisoningTheSink)
{
    std::vector<int> observed;
    BamRecordSink sink([&](const bam1_t* record) { observed.push_back(record->core.pos); });
    auto invalid = records({ 1 });
    invalid.emplace_back(nullptr);
    EXPECT_THROW(sink.submit(std::move(invalid)), std::invalid_argument);
    sink.submit(records({ 2 }));
    sink.finish();
    EXPECT_EQ(observed, (std::vector<int> { 2 }));
}

struct BatchOutputFixture
{
    BatchOutputFixture()
        : graph(graphtools::makeDeletionGraph("ACGTAC", "C", "TGTA"))
        , contigs({ { "chr1", 1000 } })
    {
        graph.setGraphId("batch-locus");
        catalog.push_back(ehunter::LocusSpecification(
            "batch-locus", ehunter::ChromType::kAutosome, { ehunter::GenomicRegion(0, 0, 100) }, graph,
            { { 0, ehunter::GenomicRegion(0, 10, 16) } }, ehunter::GenotyperParameters()));
    }
    void write(BamletWriter& writer, int index)
    {
        const auto alignment = graphtools::decodeGraphAlignment(0, "0[5M]", &graph);
        writer.write("batch-locus", "repeat-" + std::to_string(index % 7), "aCGTa", index % 2,
            index % 3, index % 5, alignment);
    }
    graphtools::Graph graph;
    ehunter::ReferenceContigInfo contigs;
    ehunter::RegionCatalog catalog;
};

struct BatchOutputFile
{
    BatchOutputFile()
    {
        char name[] = "/tmp/eh-batch-sink-XXXXXX";
        const int fd = mkstemp(name);
        if (fd < 0)
            throw std::runtime_error("Cannot create batch-output test file");
        close(fd);
        path = name;
    }
    ~BatchOutputFile() { std::remove(path.c_str()); }
    std::string path;
};

std::pair<std::string, std::vector<std::string>> readBatchOutput(const std::string& path)
{
    std::unique_ptr<htsFile, decltype(&hts_close)> file(hts_open(path.c_str(), "rb"), hts_close);
    if (!file)
        throw std::runtime_error("Cannot open batch-output test file");
    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header(sam_hdr_read(file.get()), bam_hdr_destroy);
    if (!header)
        throw std::runtime_error("Cannot read batch-output test header");
    BamRecordSink::Record record(bam_init1());
    std::vector<std::string> result;
    int status = 0;
    while ((status = sam_read1(file.get(), header.get(), record.get())) >= 0)
    {
        // Complete uncompressed record bytes retain fields, tags and duplicates.
        result.emplace_back(reinterpret_cast<const char*>(&record->core), sizeof(record->core));
        result.back().append(reinterpret_cast<const char*>(record->data), record->l_data);
    }
    if (status != -1)
        throw std::runtime_error("Cannot read batch-output test records");
    std::sort(result.begin(), result.end());
    return { std::string(header->text, header->l_text), std::move(result) };
}

TEST(BamletWriterBatch, PreservesCompleteOutputAcrossAutomaticBatchBoundaries)
{
    BatchOutputFixture fixture;
    BatchOutputFile direct, batched;
    {
        BamletWriter writer(direct.path, fixture.contigs, fixture.catalog, 0, 1);
        for (int record = 0; record < 300; ++record)
            fixture.write(writer, record);
        writer.finish();
    }
    {
        BamletWriter writer(batched.path, fixture.contigs, fixture.catalog, 0, 1);
        BamletWriter::Batch batch(&writer);
        for (int record = 0; record < 300; ++record)
            fixture.write(writer, record);
        batch.finish();
        batch.finish();
        writer.finish();
    }
    const auto expected = readBatchOutput(direct.path);
    ASSERT_EQ(expected.second.size(), 300u);
    EXPECT_EQ(readBatchOutput(batched.path), expected);
}

TEST(BamletWriterBatch, UnfinishedScopeDiscardsPendingRecordsAndRestoresWriter)
{
    BatchOutputFixture fixture;
    BatchOutputFile actual, expected;
    {
        BamletWriter writer(actual.path, fixture.contigs, fixture.catalog, 0);
        {
            BamletWriter::Batch batch(&writer);
            fixture.write(writer, 1);
        }
        fixture.write(writer, 2);
        writer.finish();
    }
    {
        BamletWriter writer(expected.path, fixture.contigs, fixture.catalog, 0);
        fixture.write(writer, 2);
        writer.finish();
    }
    EXPECT_EQ(readBatchOutput(actual.path), readBatchOutput(expected.path));
}

TEST(BamletWriterBatch, NullScopeIsNoOpAndNestedScopesRestoreTheirParent)
{
    BatchOutputFixture fixture;
    BatchOutputFile actual, expected;
    {
        BamletWriter writer(actual.path, fixture.contigs, fixture.catalog, 0);
        BamletWriter::Batch outer(&writer);
        fixture.write(writer, 1);
        {
            BamletWriter::Batch nullScope(nullptr);
            fixture.write(writer, 2);
            nullScope.finish();
        }
        {
            BamletWriter::Batch inner(&writer);
            fixture.write(writer, 3);
            inner.finish();
        }
        fixture.write(writer, 4);
        outer.finish();
        writer.finish();
    }
    {
        BamletWriter writer(expected.path, fixture.contigs, fixture.catalog, 0);
        for (int record = 1; record <= 4; ++record)
            fixture.write(writer, record);
        writer.finish();
    }
    EXPECT_EQ(readBatchOutput(actual.path), readBatchOutput(expected.path));
}

TEST(BamletWriterBatch, ThreadLocalScopeDoesNotCaptureAnotherWritersRecords)
{
    BatchOutputFixture fixture;
    BatchOutputFile first, second;
    {
        BamletWriter writerA(first.path, fixture.contigs, fixture.catalog, 0);
        BamletWriter writerB(second.path, fixture.contigs, fixture.catalog, 0);
        BamletWriter::Batch batch(&writerA);
        fixture.write(writerA, 1);
        fixture.write(writerB, 2);
        batch.finish();
        writerA.finish();
        writerB.finish();
    }
    EXPECT_EQ(readBatchOutput(first.path).second.size(), 1u);
    EXPECT_EQ(readBatchOutput(second.path).second.size(), 1u);
    EXPECT_NE(readBatchOutput(first.path).second, readBatchOutput(second.path).second);
}

}
