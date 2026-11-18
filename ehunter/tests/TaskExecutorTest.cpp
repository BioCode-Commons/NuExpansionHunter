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
#include "core/TaskExecutor.hh"
#include "core/OrderedBlocks.hh"

#include <atomic>
#include <chrono>
#include <future>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "gtest/gtest.h"

using namespace ehunter;

namespace
{
// Ensure a failed expectation cannot leave a fixture's pool waiting at teardown.
struct ReleaseGate
{
    std::promise<void> promise;
    std::shared_future<void> future { promise.get_future() };
    bool opened = false;
    void open()
    {
        if (!opened)
        {
            promise.set_value();
            opened = true;
        }
    }
    ~ReleaseGate() { open(); }
};
}

TEST(OrderedBlocks, EmptySharedPoolFallsBackToLocalWorkers)
{
    ctpl::thread_pool empty;
    std::vector<size_t> result;
    orderedBlocks(4, 1, 2, [](int, size_t first, size_t) { return first; },
                  [&](size_t value) { result.push_back(value); }, &empty);
    EXPECT_EQ((std::vector<size_t>{0, 1, 2, 3}), result);
}

TEST(OrderedBlocks, SharedPoolDrainsBorrowedWorkBeforePropagatingWriterFailure)
{
    ctpl::thread_pool pool(2);
    std::promise<void> enteredConsumer, release;
    auto gate = release.get_future().share();
    auto entered = enteredConsumer.get_future();
    std::atomic<bool> finished(false);
    auto writer = std::async(std::launch::async, [&] {
        orderedBlocks(2, 1, 2, [&](int, size_t begin, size_t) {
            if (begin == 1)
            {
                gate.wait();
                finished = true;
            }
            return begin;
        }, [&](size_t) {
            enteredConsumer.set_value();
            throw std::runtime_error("output failed");
        }, &pool);
    });
    EXPECT_EQ(std::future_status::ready, entered.wait_for(std::chrono::seconds(3)));
    EXPECT_EQ(std::future_status::timeout, writer.wait_for(std::chrono::milliseconds(10)));
    release.set_value();
    EXPECT_THROW(writer.get(), std::runtime_error);
    EXPECT_TRUE(finished);
    EXPECT_EQ(42, pool.push([](int) { return 42; }).get());
}

TEST(TaskExecutor, RejectsInvalidConfigurationWithoutTakingOwnership)
{
    EXPECT_THROW(TaskExecutor(0), std::invalid_argument);
    TaskExecutor executor(2);
    std::atomic<int> completed(0);
    EXPECT_THROW(executor.submit(TaskExecutor::Task(), 0, [&] { ++completed; }), std::invalid_argument);
    EXPECT_THROW(executor.submit([](size_t) {}, 2, [&] { ++completed; }), std::out_of_range);
    executor.wait();
    EXPECT_EQ(0, completed.load());
}

TEST(TaskExecutor, StableWorkerIdsOwnDistinctThreadLocalState)
{
    TaskExecutor executor(4);
    std::vector<std::thread::id> owners(4);
    std::vector<size_t> counts(4, 0);
    std::mutex mutex;
    std::atomic<bool> consistent(true);
    for (size_t task = 0; task != 2048; ++task)
    {
        executor.submit([&](size_t worker) {
            std::lock_guard<std::mutex> lock(mutex);
            if (owners[worker] == std::thread::id())
            {
                owners[worker] = std::this_thread::get_id();
            }
            if (owners[worker] != std::this_thread::get_id())
            {
                consistent.store(false);
            }
            ++counts[worker];
        });
    }
    executor.wait();
    size_t total = 0;
    for (size_t worker = 0; worker != counts.size(); ++worker)
    {
        total += counts[worker];
        for (size_t other = worker + 1; other != owners.size(); ++other)
        {
            if (counts[worker] && counts[other])
            {
                EXPECT_NE(owners[worker], owners[other]);
            }
        }
    }
    EXPECT_TRUE(consistent.load());
    EXPECT_EQ(2048u, total);
}

TEST(TaskExecutor, StealsSkewedWorkWhileAnotherTaskIsBlocked)
{
    TaskExecutor executor(4);
    std::promise<void> started;
    std::promise<void> release;
    std::shared_future<void> gate(release.get_future());
    std::atomic<size_t> blockedWorker(0);
    executor.submit([&](size_t worker) {
        blockedWorker.store(worker);
        started.set_value();
        gate.wait();
    }, 0);
    started.get_future().wait();
    std::promise<void> finished;
    std::future<void> finishedFuture = finished.get_future();
    std::atomic<int> remaining(512);
    for (size_t task = 0; task != 512; ++task)
    {
        executor.submit([&](size_t) {
            if (remaining.fetch_sub(1) == 1)
            {
                finished.set_value();
            }
        }, blockedWorker.load());
    }
    const bool progressed = finishedFuture.wait_for(std::chrono::seconds(5)) == std::future_status::ready;
    release.set_value();
    executor.wait();
    EXPECT_TRUE(progressed);
    EXPECT_EQ(0, remaining.load());
}

TEST(TaskExecutor, WaitIncludesDependentContinuationsAndCompletionSubmissions)
{
    TaskExecutor executor(3);
    std::atomic<int> bodies(0);
    std::atomic<int> completions(0);
    std::function<void(size_t)> enqueue;
    enqueue = [&](size_t depth) {
        executor.submit([&, depth](size_t) {
            ++bodies;
            if (depth)
            {
                enqueue(depth - 1);
            }
        }, TaskExecutor::anyWorker(), [&] { ++completions; });
    };
    executor.submit([](size_t) {}, TaskExecutor::anyWorker(), [&] { enqueue(499); });
    executor.wait();
    EXPECT_EQ(500, bodies.load());
    EXPECT_EQ(500, completions.load());
}

TEST(TaskExecutor, RepeatedIdlePhasesDoNotLoseWakeups)
{
    TaskExecutor executor(3);
    std::atomic<int> total(0);
    for (size_t phase = 0; phase != 200; ++phase)
    {
        executor.submit([&](size_t) { ++total; });
        executor.wait();
        EXPECT_EQ(static_cast<int>(phase + 1), total.load());
    }
}

TEST(TaskExecutor, ConcurrentProducersDrainAllAcceptedTasks)
{
    TaskExecutor executor(4);
    std::atomic<int> total(0);
    std::vector<std::thread> producers;
    for (size_t producer = 0; producer != 4; ++producer)
    {
        producers.emplace_back([&] {
            for (size_t task = 0; task != 500; ++task)
            {
                executor.submit([&](size_t) { ++total; });
            }
        });
    }
    for (std::thread& producer : producers)
    {
        producer.join();
    }
    executor.wait();
    EXPECT_EQ(2000, total.load());
}

TEST(TaskExecutor, FailureCompletesCancelledTasksAndPreservesFirstException)
{
    TaskExecutor executor(1);
    std::promise<void> started;
    std::promise<void> release;
    std::shared_future<void> gate(release.get_future());
    std::atomic<int> bodies(0);
    std::atomic<int> completions(0);
    executor.submit([&](size_t) {
        started.set_value();
        gate.wait();
        throw std::runtime_error("first task failure");
    }, 0, [&] { ++completions; });
    started.get_future().wait();
    for (size_t task = 0; task != 64; ++task)
    {
        executor.submit([&](size_t) { ++bodies; }, 0, [&] { ++completions; });
    }
    release.set_value();
    try
    {
        executor.wait();
        FAIL() << "Expected the task failure";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_EQ("first task failure", std::string(error.what()));
    }
    EXPECT_TRUE(executor.isCancelled());
    EXPECT_EQ(0, bodies.load());
    EXPECT_EQ(65, completions.load());
    EXPECT_THROW(executor.rethrowIfFailed(), std::runtime_error);
    EXPECT_THROW(executor.submit([](size_t) {}, 0, [&] { ++completions; }), std::runtime_error);
    executor.cancel(std::make_exception_ptr(std::logic_error("later error")));
    EXPECT_THROW(executor.wait(), std::runtime_error);
    EXPECT_EQ(65, completions.load());
}

TEST(TaskExecutor, CancellationReturnsCreditsForEveryAcceptedTask)
{
    TaskExecutor executor(1);
    std::promise<void> started;
    std::promise<void> release;
    std::shared_future<void> gate(release.get_future());
    std::atomic<int> bodies(0);
    std::atomic<int> completions(0);
    executor.submit([&](size_t) {
        started.set_value();
        gate.wait();
        ++bodies;
    }, 0, [&] { ++completions; });
    started.get_future().wait();
    for (size_t task = 0; task != 64; ++task)
    {
        executor.submit([&](size_t) { ++bodies; }, 0, [&] { ++completions; });
    }
    executor.cancel();
    release.set_value();
    EXPECT_THROW(executor.wait(), std::runtime_error);
    EXPECT_EQ(1, bodies.load());
    EXPECT_EQ(65, completions.load());
}

TEST(TaskExecutor, CompletionFailureAlsoCancelsAndDrains)
{
    TaskExecutor executor(1);
    std::atomic<int> completed(0);
    executor.submit([](size_t) {}, 0, [&] {
        ++completed;
        throw std::domain_error("completion failure");
    });
    EXPECT_THROW(executor.wait(), std::domain_error);
    EXPECT_EQ(1, completed.load());
    EXPECT_TRUE(executor.isCancelled());
}

TEST(TaskExecutor, RejectsWaitFromWorkerInsteadOfDeadlocking)
{
    TaskExecutor executor(1);
    std::atomic<int> completed(0);
    executor.submit([&](size_t) { executor.wait(); }, 0, [&] { ++completed; });
    EXPECT_THROW(executor.wait(), std::logic_error);
    EXPECT_EQ(1, completed.load());
}

TEST(TaskExecutor, DestroysBodyCapturesBeforeReturningAdmissionCredit)
{
    TaskExecutor executor(1);
    std::shared_ptr<int> payload(new int(7));
    std::weak_ptr<int> borrowed(payload);
    std::promise<void> release;
    std::shared_future<void> gate(release.get_future());
    std::atomic<bool> released(false);
    executor.submit([payload, gate](size_t) { gate.wait(); }, 0, [&] { released.store(borrowed.expired()); });
    payload.reset();
    release.set_value();
    executor.wait();
    EXPECT_TRUE(released.load());
}

TEST(TaskExecutor, DestructorDrainsContinuationsWithoutExplicitWait)
{
    std::atomic<int> completed(0);
    {
        TaskExecutor executor(2);
        executor.submit([&](size_t) {
            executor.submit([&](size_t) { ++completed; }, TaskExecutor::anyWorker(), [&] { ++completed; });
        }, TaskExecutor::anyWorker(), [&] { ++completed; });
    }
    EXPECT_EQ(3, completed.load());
}

TEST(TaskExecutor, DestructorSuppressesFailureAfterRunningTerminalCallbacks)
{
    std::atomic<int> completed(0);
    {
        TaskExecutor executor(2);
        executor.submit([](size_t) { throw std::runtime_error("failure during teardown"); },
            TaskExecutor::anyWorker(), [&] { ++completed; });
    }
    EXPECT_EQ(1, completed.load());
}
