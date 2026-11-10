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
