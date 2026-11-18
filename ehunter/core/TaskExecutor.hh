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
#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace ehunter
{

// A fixed set of compute workers. The caller separately budgets input/codec
// threads and bounds admitted batches; submit never waits for capacity.
//
// Tasks must not wait for other tasks. Submit a continuation instead. Preferred
// workers improve locality but do not pin tasks: idle workers may steal them.
// Pool worker IDs belong to distinct threads, including during completions.
// Contract [K-36]
// Ensures:  every accepted task's completion runs exactly once; wait() returns only after quiescence.
class TaskExecutor
{
public:
    using Task = std::function<void(size_t)>;
    using Completion = std::function<void()>;

    static size_t anyWorker() { return std::numeric_limits<size_t>::max(); }

    explicit TaskExecutor(size_t workerCount);
    ~TaskExecutor();
    TaskExecutor(const TaskExecutor&) = delete;
    TaskExecutor& operator=(const TaskExecutor&) = delete;

    // On success the executor owns the task and invokes completion exactly once,
    // after its body finishes or is skipped due to cancellation. If submit throws,
    // ownership was not transferred and completion is not invoked. A completion
    // can release admission credit, even when the task body throws. The body and
    // its captures are destroyed before completion is called. Completion must
    // not block on other work; its exception also cancels the executor.
    void submit(Task task, size_t preferredWorker = anyWorker(), Completion completion = Completion());

    // Wait for quiescence, including continuations submitted by accepted tasks or
    // completions, then rethrow the first error. Workers remain available for the
    // next phase. External producers must have stopped for a final drain. Calling
    // wait from one of this executor's workers is an error.
    void wait();
    void cancel(std::exception_ptr reason = std::exception_ptr());
    bool isCancelled() const noexcept { return cancelled_.load(std::memory_order_acquire); }
    void rethrowIfFailed() const;
    size_t workerCount() const noexcept { return queues_.size(); }
    size_t queuedTasks() const noexcept { return ready_.load(std::memory_order_relaxed); }

private:
    struct Work
    {
        Work(Task body, Completion terminal)
            : task(std::move(body))
            , completion(std::move(terminal))
        {
        }
        Work() = default;
        Task task;
        Completion completion;
    };
    struct Queue
    {
        std::mutex mutex;
        std::deque<Work> work;
    };

    bool take(size_t workerId, size_t& nextVictim, Work& work);
    void execute(Work& work, size_t workerId);
    void run(size_t workerId);
    void recordFailure(std::exception_ptr error);
    void drain();
    void stopAndJoin() noexcept;

    std::vector<std::unique_ptr<Queue>> queues_;
    std::vector<std::thread> workers_;
    std::atomic<size_t> nextQueue_ { 0 };
    std::atomic<size_t> ready_ { 0 };
    std::atomic<bool> cancelled_ { false };
    mutable std::mutex stateMutex_;
    std::condition_variable workAvailable_;
    std::condition_variable drained_;
    std::atomic<size_t> outstanding_ { 0 };
    std::atomic<size_t> sleepers_ { 0 };
    std::atomic<bool> stopping_ { false };
    std::exception_ptr firstError_;

    static thread_local TaskExecutor* currentExecutor_;
    static thread_local size_t currentWorker_;
};

}
