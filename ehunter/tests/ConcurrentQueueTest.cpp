//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
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

#include "core/ConcurrentQueue.hh"
#include "core/ThreadPool.hh"
#include "sample/HtsStreamingReadPairQueue.hh"

#include <chrono>
#include <future>
#include <vector>

#include "gtest/gtest.h"

using namespace ehunter;

TEST(ConcurrentQueueTest, SerialOps)
{
    // Sanity-check non-concurrent operation
    ConcurrentQueue<int> cq;
    cq.push(1);
    cq.push(2);
    int r;
    cq.pop(r);
    EXPECT_FALSE(cq.empty());
    EXPECT_EQ(r, 1);
    cq.pop(r);
    EXPECT_TRUE(cq.empty());
    EXPECT_EQ(r, 2);
}

TEST(ConcurrentQueueTest, WeightedBudgetBlocksUntilSpaceIsReleased)
{
    ConcurrentQueue<int> queue(10);
    queue.push(1, 8);
    std::promise<void> started;
    auto producer = std::async(std::launch::async, [&] {
        started.set_value();
        queue.push(2, 3);
    });
    started.get_future().wait();
    EXPECT_EQ(producer.wait_for(std::chrono::milliseconds(100)), std::future_status::timeout);

    int value = 0;
    queue.pop(value);
    EXPECT_EQ(value, 1);
    EXPECT_EQ(producer.wait_for(std::chrono::seconds(5)), std::future_status::ready);
    producer.get();
    queue.pop(value);
    EXPECT_EQ(value, 2);
    EXPECT_TRUE(queue.empty());
}

TEST(ConcurrentQueueTest, OversizedRecordCanProgressWithoutDiscardingIt)
{
    ConcurrentQueue<int> queue(10);
    queue.push(1, 20);
    std::promise<void> started;
    auto producer = std::async(std::launch::async, [&] {
        started.set_value();
        queue.push(2, 1);
    });
    started.get_future().wait();
    EXPECT_EQ(producer.wait_for(std::chrono::milliseconds(100)), std::future_status::timeout);
    int value = 0;
    queue.pop(value);
    EXPECT_EQ(value, 1);
    producer.get();
    queue.pop(value);
    EXPECT_EQ(value, 2);
}

TEST(ConcurrentQueueTest, ZeroWeightSentinelFollowsRecordsAtCapacity)
{
    ConcurrentQueue<int> queue(10);
    queue.push(42, 10);
    queue.push(-1, 0);
    int value = 0;
    queue.pop(value);
    EXPECT_EQ(value, 42);
    queue.pop(value);
    EXPECT_EQ(value, -1);
    EXPECT_TRUE(queue.empty());
}

TEST(ConcurrentQueueTest, MultipleProducersPreserveEveryRecordAndTheirOrder)
{
    ConcurrentQueue<int> queue(17);
    const int producerCount = 4;
    const int recordsPerProducer = 1000;
    std::vector<std::future<void>> producers;
    for (int producer = 0; producer < producerCount; ++producer)
    {
        producers.emplace_back(std::async(std::launch::async, [&, producer] {
            for (int record = 0; record < recordsPerProducer; ++record)
                queue.push(producer * recordsPerProducer + record, 1 + record % 23);
        }));
    }
    std::vector<int> next(producerCount, 0);
    for (int index = 0; index < producerCount * recordsPerProducer; ++index)
    {
        int value = 0;
        queue.pop(value);
        const int producer = value / recordsPerProducer;
        ASSERT_GE(producer, 0);
        ASSERT_LT(producer, producerCount);
        EXPECT_EQ(value % recordsPerProducer, next[producer]++);
    }
    for (auto& producer : producers)
        producer.get();
    for (int count : next)
        EXPECT_EQ(count, recordsPerProducer);
    EXPECT_TRUE(queue.empty());
}

namespace
{
HtsStreamingReadPairQueue::ReadPair streamingPair(const std::string& name)
{
    return { locus::RegionType::kTarget, AnalyzerInputType::kBothReads,
        Read(ReadId(name, MateNumber::kFirstMate), "ACGT", false),
        Read(ReadId(name, MateNumber::kSecondMate), "TGCA", true) };
}

HtsStreamingReadPairQueue::ReadPair streamingPair150(const std::string& name)
{
    return { locus::RegionType::kTarget, AnalyzerInputType::kBothReads,
        Read(ReadId(name, MateNumber::kFirstMate), std::string(150, 'A'), false),
        Read(ReadId(name, MateNumber::kSecondMate), std::string(150, 'C'), true) };
}
}

TEST(StreamingReadPairQueueTest, BoundedHelpRetainsActiveOwnershipAndAppendedPairOrder)
{
    using Result = HtsStreamingReadPairQueue::InsertResult;
    HtsStreamingReadPairQueue queue(1, 2);
    for (unsigned i = 0; i < 75; ++i)
        ASSERT_EQ(i ? Result::Queued : Result::Activated,
            queue.insertReadPair(0, streamingPair150(std::to_string(i))));
    boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
    std::vector<unsigned> observed;
    const auto process = [&](HtsStreamingReadPairQueue::ReadPair& current) {
        EXPECT_EQ(150u, current.read.sequence().size());
        observed.push_back(std::stoul(current.read.fragmentId()));
    };
    EXPECT_TRUE(queue.processBatch(0, 32, pair, process));
    EXPECT_FALSE(bool(pair));
    ASSERT_EQ(32u, observed.size());
    // The yielded locus remains active. Appending must neither schedule a
    // second owner nor call the reader helper while its continuation is pending.
    EXPECT_EQ(Result::Queued, queue.insertReadPair(0, streamingPair150("75"), [] {
        ADD_FAILURE() << "An active locus must not request another helper";
        return false;
    }));
    EXPECT_FALSE(queue.processBatch(0, 0, pair, process));
    ASSERT_EQ(76u, observed.size());
    for (unsigned i = 0; i < observed.size(); ++i)
        EXPECT_EQ(i, observed[i]);
    EXPECT_EQ(Result::Activated, queue.insertReadPair(1, streamingPair150("capacity-released")));
}

TEST(StreamingReadPairQueueTest, ExactHelpBudgetPublishesOneEmptyContinuationThatReleasesCapacity)
{
    using Result = HtsStreamingReadPairQueue::InsertResult;
    HtsStreamingReadPairQueue queue(1, 2);
    for (unsigned i = 0; i < 32; ++i)
        queue.insertReadPair(0, streamingPair150(std::to_string(i)));
    ctpl::thread_pool pool(0);
    unsigned processed = 0, helped = 0, resumed = 0;
    auto first = pool.push([&](int index) {
        EXPECT_EQ(2, index); // Explicit reader slot, not the last worker slot.
        boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
        const bool continuation = queue.processBatch(0, 32, pair,
            [&](HtsStreamingReadPairQueue::ReadPair&) { ++processed; });
        ASSERT_TRUE(continuation);
        ASSERT_FALSE(bool(pair));
        pool.push([&](int worker) {
            EXPECT_EQ(0, worker);
            ++resumed;
            boost::optional<HtsStreamingReadPairQueue::ReadPair> next;
            EXPECT_FALSE(queue.processBatch(0, 0, next,
                [&](HtsStreamingReadPairQueue::ReadPair&) { ADD_FAILURE() << "No pair should remain"; }));
        });
    });
    auto reader = pool.pop();
    ASSERT_TRUE(bool(reader));
    reader(2);
    first.get();
    EXPECT_EQ(32u, processed);
    EXPECT_EQ(Result::Activated, queue.insertReadPair(1, streamingPair150("next-locus"), [&] {
        ++helped;
        auto continuation = pool.pop();
        EXPECT_TRUE(bool(continuation));
        EXPECT_FALSE(bool(pool.pop())); // Exactly one continuation was published.
        continuation(0);
        return true;
    }));
    EXPECT_EQ(1u, helped);
    EXPECT_EQ(1u, resumed);
}

TEST(StreamingReadPairQueueTest, BoundedHelpFailurePreservesCurrentPairAndDoesNotRequestContinuation)
{
    HtsStreamingReadPairQueue queue(1, 2);
    for (unsigned i = 0; i < 40; ++i)
        queue.insertReadPair(0, streamingPair150(std::to_string(i)));
    boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
    unsigned processed = 0;
    bool continuation = false;
    EXPECT_THROW(continuation = queue.processBatch(0, 32, pair,
        [&](HtsStreamingReadPairQueue::ReadPair&) {
            if (++processed == 7)
                throw std::runtime_error("alignment failure");
        }), std::runtime_error);
    EXPECT_FALSE(continuation);
    ASSERT_TRUE(bool(pair));
    EXPECT_EQ("6", pair->read.fragmentId());
    queue.cancel();
    EXPECT_EQ(HtsStreamingReadPairQueue::InsertResult::Cancelled,
        queue.insertReadPair(1, streamingPair150("cancelled")));
    // Cancellation currently rejects activation but retains existing queues;
    // preserve that established contract while checking the unconsumed prefix.
    queue.getNextReadPair(0, pair);
    ASSERT_TRUE(bool(pair));
    EXPECT_EQ("7", pair->read.fragmentId());
}

TEST(StreamingReadPairQueueTest, ReaderExecutesOneQueuedTaskAndPreservesPerLocusOrder)
{
    using Result = HtsStreamingReadPairQueue::InsertResult;
    HtsStreamingReadPairQueue queue(1, 2);
    EXPECT_EQ(Result::Activated, queue.insertReadPair(0, streamingPair("first")));
    EXPECT_EQ(Result::Queued, queue.insertReadPair(0, streamingPair("second")));
    // No pool threads are needed to verify the exact pop() path used by the
    // reader: only its invocation of the queued task can release capacity.
    ctpl::thread_pool pool(0);
    std::vector<std::string> names;
    int executionIndex = -1;
    auto taskResult = pool.push([&](int index) {
        executionIndex = index;
        boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
        while (true)
        {
            queue.getNextReadPair(0, pair);
            if (!pair)
                break;
            names.push_back(pair->read.fragmentId());
        }
    });
    unsigned helpCount = 0;
    const auto helper = [&] {
        ++helpCount;
        auto task = pool.pop();
        if (!task)
            return false;
        task(2);
        return true;
    };
    EXPECT_EQ(Result::Activated, queue.insertReadPair(1, streamingPair("third"), helper));
    taskResult.get();
    EXPECT_EQ(1u, helpCount);
    EXPECT_EQ(2, executionIndex);
    EXPECT_EQ((std::vector<std::string> { "first", "second" }), names);
    EXPECT_FALSE(bool(pool.pop()));
    boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
    queue.getNextReadPair(1, pair);
    ASSERT_TRUE(bool(pair));
    EXPECT_EQ("third", pair->read.fragmentId());
    queue.getNextReadPair(1, pair);
    EXPECT_FALSE(bool(pair));
}

TEST(StreamingReadPairQueueTest, ActiveQueueAppendDoesNotRunHelper)
{
    using Result = HtsStreamingReadPairQueue::InsertResult;
    HtsStreamingReadPairQueue queue(1, 1);
    EXPECT_EQ(Result::Activated, queue.insertReadPair(0, streamingPair("first")));
    unsigned calls = 0;
    EXPECT_EQ(Result::Queued, queue.insertReadPair(0, streamingPair("second"), [&] {
        ++calls;
        return false;
    }));
    EXPECT_EQ(0u, calls);
}

TEST(StreamingReadPairQueueTest, RechecksCapacityAfterHelperFindsNoPendingTask)
{
    using Result = HtsStreamingReadPairQueue::InsertResult;
    HtsStreamingReadPairQueue queue(1, 2);
    EXPECT_EQ(Result::Activated, queue.insertReadPair(0, streamingPair("first")));
    unsigned calls = 0;
    EXPECT_EQ(Result::Activated, queue.insertReadPair(1, streamingPair("second"), [&] {
        ++calls;
        // Model a worker completing between an unsuccessful pool.pop() and
        // the producer reacquiring the capacity mutex. Its notification has
        // already happened, so waiting without a fresh check would deadlock.
        boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
        queue.getNextReadPair(0, pair);
        EXPECT_TRUE(bool(pair));
        queue.getNextReadPair(0, pair);
        EXPECT_FALSE(bool(pair));
        return false;
    }));
    EXPECT_EQ(1u, calls);
}

TEST(StreamingReadPairQueueTest, CancellationWakesCapacityWaitWithoutLosingQueuedPairs)
{
    using Result = HtsStreamingReadPairQueue::InsertResult;
    HtsStreamingReadPairQueue queue(1, 3);
    EXPECT_EQ(Result::Activated, queue.insertReadPair(0, streamingPair("keep")));
    std::promise<void> helperStarted;
    bool announced = false;
    auto producer = std::async(std::launch::async, [&] {
        return queue.insertReadPair(1, streamingPair("cancelled"), [&] {
            if (!announced)
            {
                announced = true;
                helperStarted.set_value();
            }
            return false;
        });
    });
    helperStarted.get_future().wait();
    queue.cancel();
    EXPECT_EQ(std::future_status::ready, producer.wait_for(std::chrono::seconds(5)));
    EXPECT_EQ(Result::Cancelled, producer.get());
    EXPECT_EQ(Result::Cancelled, queue.insertReadPair(2, streamingPair("also-cancelled")));
    boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
    queue.getNextReadPair(0, pair);
    ASSERT_TRUE(bool(pair));
    EXPECT_EQ("keep", pair->read.fragmentId());
    queue.getNextReadPair(0, pair);
    EXPECT_FALSE(bool(pair));
}

TEST(StreamingReadPairQueueTest, CapturedHelperFailureCancelsWaitAndRemainsObservable)
{
    using Result = HtsStreamingReadPairQueue::InsertResult;
    HtsStreamingReadPairQueue queue(1, 2);
    EXPECT_EQ(Result::Activated, queue.insertReadPair(0, streamingPair("first")));
    ctpl::thread_pool pool(0);
    std::exception_ptr recordedError;
    auto taskResult = pool.push([&](int) {
        try
        {
            throw std::runtime_error("locus failure");
        }
        catch (...)
        {
            recordedError = std::current_exception();
            queue.cancel();
            throw;
        }
    });
    EXPECT_EQ(Result::Cancelled, queue.insertReadPair(1, streamingPair("not-inserted"), [&] {
        auto task = pool.pop();
        if (!task)
            return false;
        task(2); // packaged_task captures the exception just as in production.
        return true;
    }));
    EXPECT_THROW(taskResult.get(), std::runtime_error);
    ASSERT_TRUE(bool(recordedError));
    EXPECT_THROW(std::rethrow_exception(recordedError), std::runtime_error);
}
