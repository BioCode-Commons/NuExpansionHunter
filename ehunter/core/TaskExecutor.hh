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
// Optional publisher assistance uses the separate ID workerCount().
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

    // The single external publisher may execute one queued task synchronously,
    // using worker ID workerCount() without creating another OS thread. Returns
    // false when no queued task can be taken immediately (including contention),
    // true after an accepted task is processed, even if cancellation skips its
    // body. Task/completion errors are recorded for wait()/rethrowIfFailed().
    // The publisher must provide separate mutable workspace for its extra ID.
    // Preferred queue IDs still range from zero to workerCount()-1; assisted
    // continuations can use anyWorker() to select a real queue.
    // Do not call from a task/completion or concurrently from two publishers;
    // these uses throw logic_error. Stop assistance before final destruction.
    bool tryRunOne();

    // Wait for quiescence, including continuations submitted by accepted tasks or
    // completions, then rethrow the first error. Workers remain available for the
    // next phase. External producers must have stopped for a final drain. Calling
    // wait from one of this executor's workers is an error.
    void wait();
    void cancel(std::exception_ptr reason = std::exception_ptr());
    bool isCancelled() const noexcept { return cancelled_.load(std::memory_order_acquire); }
    void rethrowIfFailed() const;
    size_t workerCount() const noexcept { return queues_.size(); }

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
    std::atomic<bool> publisherActive_ { false };
    size_t publisherNextVictim_ = 0; // Protected by publisherActive_.
    mutable std::mutex stateMutex_;
    std::condition_variable workAvailable_;
    std::condition_variable drained_;
    size_t outstanding_ = 0;
    bool stopping_ = false;
    std::exception_ptr firstError_;

    static thread_local TaskExecutor* currentExecutor_;
    static thread_local size_t currentWorker_;
};

}
