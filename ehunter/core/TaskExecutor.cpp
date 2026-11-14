#include "core/TaskExecutor.hh"

#include <stdexcept>
#include <utility>

namespace ehunter
{

thread_local TaskExecutor* TaskExecutor::currentExecutor_ = nullptr;
thread_local size_t TaskExecutor::currentWorker_ = TaskExecutor::anyWorker();

TaskExecutor::TaskExecutor(size_t workerCount)
{
    if (workerCount == 0)
    {
        throw std::invalid_argument("TaskExecutor requires at least one worker");
    }
    queues_.reserve(workerCount);
    workers_.reserve(workerCount);
    for (size_t worker = 0; worker != workerCount; ++worker)
    {
        queues_.emplace_back(new Queue());
    }
    try
    {
        for (size_t worker = 0; worker != workerCount; ++worker)
        {
            workers_.emplace_back(&TaskExecutor::run, this, worker);
        }
    }
    catch (...)
    {
        // A partly constructed object has no destructor. Wake and join every
        // thread already started before propagating the startup failure.
        stopAndJoin();
        throw;
    }
}

TaskExecutor::~TaskExecutor()
{
    // The owner must outlive its tasks and stop external producers before
    // destruction. Waiting from a worker could otherwise join that same thread.
    if (currentExecutor_ == this)
    {
        std::terminate();
    }
    drain();
    stopAndJoin();
}

void TaskExecutor::submit(Task task, size_t preferredWorker, Completion completion)
{
    if (!task)
    {
        throw std::invalid_argument("Cannot submit an empty task");
    }
    if (preferredWorker != anyWorker() && preferredWorker >= queues_.size())
    {
        throw std::out_of_range("TaskExecutor worker preference is out of range");
    }
    const size_t target = preferredWorker != anyWorker() ? preferredWorker
        : currentExecutor_ == this && currentWorker_ < queues_.size() ? currentWorker_
            : nextQueue_.fetch_add(1, std::memory_order_relaxed) % queues_.size();
    {
        // This mutex protects only admission/drain bookkeeping and the sleep
        // predicate; task storage and stealing use separate per-worker locks.
        std::lock_guard<std::mutex> stateLock(stateMutex_);
        if (firstError_)
        {
            std::rethrow_exception(firstError_);
        }
        if (stopping_)
        {
            throw std::logic_error("TaskExecutor is stopping");
        }
        Queue& queue = *queues_[target];
        std::lock_guard<std::mutex> queueLock(queue.mutex);
        queue.work.emplace_back(std::move(task), std::move(completion));
        // Publish the count before unlocking the queue so a thief cannot take a
        // newly accepted task before it has been included in both counts.
        ++outstanding_;
        ready_.fetch_add(1, std::memory_order_release);
    }
    workAvailable_.notify_one();
}

bool TaskExecutor::take(size_t workerId, size_t& nextVictim, Work& work)
{
    // The external publisher has no own queue. It only attempts nonblocking
    // steals; its extra workspace ID must never index queues_.
    if (workerId < queues_.size())
    {
        Queue& own = *queues_[workerId];
        std::lock_guard<std::mutex> lock(own.mutex);
        if (!own.work.empty())
        {
            // Local continuations retain their recent data. Thieves take older
            // work from the opposite end to spread skewed batches.
            work = std::move(own.work.back());
            own.work.pop_back();
            ready_.fetch_sub(1, std::memory_order_acq_rel);
            return true;
        }
    }
    for (size_t scanned = 0; scanned != queues_.size(); ++scanned)
    {
        const size_t victim = nextVictim;
        nextVictim = (nextVictim + 1) % queues_.size();
        if (victim == workerId)
        {
            continue;
        }
        Queue& other = *queues_[victim];
        std::unique_lock<std::mutex> lock(other.mutex, std::try_to_lock);
        if (lock.owns_lock() && !other.work.empty())
        {
            work = std::move(other.work.front());
            other.work.pop_front();
            ready_.fetch_sub(1, std::memory_order_acq_rel);
            return true;
        }
    }
    return false;
}

bool TaskExecutor::tryRunOne()
{
    if (currentExecutor_ != nullptr)
    {
        throw std::logic_error("TaskExecutor publisher assistance cannot run inside a task or completion");
    }
    bool inactive = false;
    if (!publisherActive_.compare_exchange_strong(inactive, true, std::memory_order_acq_rel))
    {
        throw std::logic_error("TaskExecutor permits only one assisting publisher at a time");
    }
    struct PublisherScope
    {
        std::atomic<bool>& active;
        TaskExecutor*& executor;
        size_t& worker;
        TaskExecutor* previousExecutor;
        size_t previousWorker;
        ~PublisherScope()
        {
            executor = previousExecutor;
            worker = previousWorker;
            active.store(false, std::memory_order_release);
        }
    } scope { publisherActive_, currentExecutor_, currentWorker_, currentExecutor_, currentWorker_ };

    Work work;
    if (!take(workerCount(), publisherNextVictim_, work))
    {
        return false;
    }
    currentExecutor_ = this;
    currentWorker_ = workerCount();
    execute(work, currentWorker_);
    return true;
}

void TaskExecutor::execute(Work& work, size_t workerId)
{
    if (!isCancelled())
    {
        try
        {
            work.task(workerId);
        }
        catch (...)
        {
            recordFailure(std::current_exception());
        }
    }
    work.task = Task();
    try
    {
        if (work.completion)
        {
            work.completion();
        }
    }
    catch (...)
    {
        recordFailure(std::current_exception());
    }
    // Release completion captures before reporting the task drained.
    work = Work();
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        --outstanding_;
        if (outstanding_ == 0)
        {
            drained_.notify_all();
        }
    }
}

void TaskExecutor::run(size_t workerId)
{
    currentExecutor_ = this;
    currentWorker_ = workerId;
    size_t nextVictim = (workerId + 1) % queues_.size();
    while (true)
    {
        Work work;
        if (!take(workerId, nextVictim, work))
        {
            std::unique_lock<std::mutex> lock(stateMutex_);
            // Submission changes this predicate while holding stateMutex_, so
            // notification cannot be lost between checking and falling asleep.
            workAvailable_.wait(lock, [this] { return stopping_ || ready_.load(std::memory_order_acquire) != 0; });
            if (stopping_)
            {
                break;
            }
            continue;
        }
        execute(work, workerId);
    }
    currentWorker_ = anyWorker();
    currentExecutor_ = nullptr;
}

void TaskExecutor::recordFailure(std::exception_ptr error)
{
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        if (!firstError_)
        {
            firstError_ = error;
            cancelled_.store(true, std::memory_order_release);
        }
    }
    workAvailable_.notify_all();
}

void TaskExecutor::cancel(std::exception_ptr reason)
{
    if (!reason)
    {
        reason = std::make_exception_ptr(std::runtime_error("TaskExecutor cancelled"));
    }
    recordFailure(reason);
}

void TaskExecutor::rethrowIfFailed() const
{
    std::exception_ptr error;
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        error = firstError_;
    }
    if (error)
    {
        std::rethrow_exception(error);
    }
}

void TaskExecutor::drain()
{
    std::unique_lock<std::mutex> lock(stateMutex_);
    drained_.wait(lock, [this] { return outstanding_ == 0; });
}

void TaskExecutor::wait()
{
    if (currentExecutor_ == this)
    {
        throw std::logic_error("A TaskExecutor worker must submit continuations instead of waiting");
    }
    drain();
    rethrowIfFailed();
}

void TaskExecutor::stopAndJoin() noexcept
{
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        stopping_ = true;
    }
    workAvailable_.notify_all();
    for (std::thread& worker : workers_)
    {
        if (worker.joinable())
        {
            worker.join();
        }
    }
}

}
