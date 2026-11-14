#include "core/TaskExecutor.hh"

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

TEST(TaskExecutor, PublisherUsesDistinctIdAndStealsAcrossAllQueues)
{
    TaskExecutor executor(2);
    ReleaseGate release;
    std::promise<void> started;
    std::atomic<int> remaining(2);
    std::vector<std::thread::id> workerThreads(2);
    for (size_t preferred = 0; preferred < 2; ++preferred)
    {
        executor.submit([&](size_t worker) {
            workerThreads[worker] = std::this_thread::get_id();
            if (remaining.fetch_sub(1) == 1)
                started.set_value();
            release.future.wait();
        }, preferred);
    }
    started.get_future().wait();
    const auto publisher = std::this_thread::get_id();
    std::vector<size_t> ids;
    int completions = 0;
    for (size_t preferred = 0; preferred < 2; ++preferred)
    {
        executor.submit([&](size_t worker) {
            EXPECT_EQ(publisher, std::this_thread::get_id());
            ids.push_back(worker);
        }, preferred, [&] {
            EXPECT_EQ(publisher, std::this_thread::get_id());
            ++completions;
        });
    }
    EXPECT_TRUE(executor.tryRunOne());
    EXPECT_TRUE(executor.tryRunOne());
    EXPECT_FALSE(executor.tryRunOne()); // Running pool tasks are not stolen.
    EXPECT_EQ((std::vector<size_t> { 2, 2 }), ids);
    EXPECT_EQ(2, completions);
    EXPECT_NE(workerThreads[0], publisher);
    EXPECT_NE(workerThreads[1], publisher);
    EXPECT_NE(workerThreads[0], workerThreads[1]);
    release.open();
    executor.wait();
    EXPECT_FALSE(executor.tryRunOne());
}

TEST(TaskExecutor, PublisherContinuationsAndCompletionSubmissionsUseRealQueues)
{
    TaskExecutor executor(1);
    ReleaseGate release;
    std::promise<void> started;
    executor.submit([&](size_t) {
        started.set_value();
        release.future.wait();
    });
    started.get_future().wait();
    int bodies = 0;
    int completions = 0;
    auto child = [&](size_t worker) {
        EXPECT_EQ(executor.workerCount(), worker);
        ++bodies;
    };
    executor.submit([&](size_t worker) {
        child(worker);
        executor.submit(child, TaskExecutor::anyWorker(), [&] { ++completions; });
    }, TaskExecutor::anyWorker(), [&] {
        ++completions;
        executor.submit(child, TaskExecutor::anyWorker(), [&] { ++completions; });
    });
    EXPECT_TRUE(executor.tryRunOne());
    EXPECT_TRUE(executor.tryRunOne());
    EXPECT_TRUE(executor.tryRunOne());
    EXPECT_FALSE(executor.tryRunOne());
    EXPECT_EQ(3, bodies);
    EXPECT_EQ(3, completions);
    release.open();
    executor.wait();
}

TEST(TaskExecutor, PublisherFailureReleasesCapturesAndCompletesCancelledTasks)
{
    TaskExecutor executor(1);
    ReleaseGate release;
    std::promise<void> started;
    executor.submit([&](size_t) {
        started.set_value();
        release.future.wait();
    });
    started.get_future().wait();
    std::shared_ptr<int> payload(new int(7));
    std::weak_ptr<int> borrowed(payload);
    bool releasedBeforeCompletion = false;
    int bodies = 0;
    int completions = 0;
    executor.submit([payload](size_t) { throw std::runtime_error("publisher task failure"); },
        TaskExecutor::anyWorker(), [&] {
            releasedBeforeCompletion = borrowed.expired();
            ++completions;
        });
    payload.reset();
    for (size_t task = 0; task < 8; ++task)
        executor.submit([&](size_t) { ++bodies; }, TaskExecutor::anyWorker(), [&] { ++completions; });

    EXPECT_TRUE(executor.tryRunOne());
    EXPECT_TRUE(executor.isCancelled());
    for (size_t task = 0; task < 8; ++task)
        EXPECT_TRUE(executor.tryRunOne());
    EXPECT_FALSE(executor.tryRunOne());
    EXPECT_TRUE(releasedBeforeCompletion);
    EXPECT_EQ(0, bodies);
    EXPECT_EQ(9, completions);
    release.open();
    try
    {
        executor.wait();
        FAIL() << "Expected the publisher task failure";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_EQ("publisher task failure", std::string(error.what()));
    }
}

TEST(TaskExecutor, PublisherCompletionFailureCancelsAndDrains)
{
    TaskExecutor executor(1);
    ReleaseGate release;
    std::promise<void> started;
    executor.submit([&](size_t) {
        started.set_value();
        release.future.wait();
    });
    started.get_future().wait();
    int bodies = 0;
    int completions = 0;
    executor.submit([](size_t) {}, TaskExecutor::anyWorker(), [&] {
        ++completions;
        throw std::domain_error("publisher completion failure");
    });
    executor.submit([&](size_t) { ++bodies; }, TaskExecutor::anyWorker(), [&] { ++completions; });
    EXPECT_TRUE(executor.tryRunOne());
    EXPECT_TRUE(executor.tryRunOne());
    EXPECT_FALSE(executor.tryRunOne());
    EXPECT_EQ(0, bodies);
    EXPECT_EQ(2, completions);
    release.open();
    EXPECT_THROW(executor.wait(), std::domain_error);
}

TEST(TaskExecutor, RejectsWorkerReentrancyAndConcurrentPublishers)
{
    TaskExecutor executor(1);
    std::atomic<bool> workerRejected(false);
    executor.submit([&](size_t) {
        try { executor.tryRunOne(); }
        catch (const std::logic_error&) { workerRejected.store(true); }
    });
    executor.wait();
    EXPECT_TRUE(workerRejected.load());

    ReleaseGate poolRelease;
    ReleaseGate publisherRelease;
    std::promise<void> poolStarted;
    std::promise<void> publisherStarted;
    executor.submit([&](size_t) {
        poolStarted.set_value();
        poolRelease.future.wait();
    });
    poolStarted.get_future().wait();
    std::atomic<bool> reentrantRejected(false);
    std::atomic<bool> assisted(false);
    executor.submit([&](size_t worker) {
        EXPECT_EQ(executor.workerCount(), worker);
        try { executor.tryRunOne(); }
        catch (const std::logic_error&) { reentrantRejected.store(true); }
        publisherStarted.set_value();
        publisherRelease.future.wait();
    });
    std::thread publisher([&] { assisted.store(executor.tryRunOne()); });
    publisherStarted.get_future().wait();
    EXPECT_THROW(executor.tryRunOne(), std::logic_error);
    publisherRelease.open();
    publisher.join();
    poolRelease.open();
    executor.wait();
    EXPECT_TRUE(assisted.load());
    EXPECT_TRUE(reentrantRejected.load());
    EXPECT_FALSE(executor.tryRunOne());
}

TEST(TaskExecutor, WaitIncludesAnAssistedCompletion)
{
    TaskExecutor executor(1);
    ReleaseGate poolRelease;
    ReleaseGate completionRelease;
    std::promise<void> poolStarted;
    std::promise<void> completionStarted;
    executor.submit([&](size_t) {
        poolStarted.set_value();
        poolRelease.future.wait();
    });
    poolStarted.get_future().wait();
    executor.submit([](size_t) {}, TaskExecutor::anyWorker(), [&] {
        completionStarted.set_value();
        completionRelease.future.wait();
    });
    std::thread publisher([&] { EXPECT_TRUE(executor.tryRunOne()); });
    completionStarted.get_future().wait();
    poolRelease.open();
    auto waiter = std::async(std::launch::async, [&] { executor.wait(); });
    EXPECT_EQ(std::future_status::timeout, waiter.wait_for(std::chrono::milliseconds(20)));
    completionRelease.open();
    publisher.join();
    EXPECT_NO_THROW(waiter.get());
}
