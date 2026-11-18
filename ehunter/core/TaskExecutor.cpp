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

#include <stdexcept>
#if defined(__linux__)
#include <fstream>
#include <map>
#include <sched.h>
#include <string>
#endif
#include <utility>

namespace ehunter
{

#if defined(__linux__)
static thread_local std::vector<size_t> cachePeers;
static thread_local size_t nextCachePeer;

static std::vector<cpu_set_t> workerCacheMasks(size_t workers)
{
    cpu_set_t permitted;
    if (sched_getaffinity(0, sizeof(permitted), &permitted) != 0
        || workers * 2 < static_cast<size_t>(CPU_COUNT(&permitted))) return {};
    struct Domain
    {
        cpu_set_t mask {};
        size_t cpuCount = 0;
    };
    std::map<std::string, Domain> groups;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
    {
        if (!CPU_ISSET(cpu, &permitted)) continue;
        const auto path = "/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/cache/index3/shared_cpu_list";
        std::ifstream input(path);
        std::string key;
        if (!std::getline(input, key)) return {};
        auto& domain = groups[key];
        CPU_SET(cpu, &domain.mask);
        ++domain.cpuCount;
    }
    std::vector<cpu_set_t> masks;
    if (groups.size() > 1)
    {
        // Interleave the permitted domains so extra workers are distributed
        // in proportion to their available CPUs, including partial domains.
        for (size_t pass = 0; masks.size() < static_cast<size_t>(CPU_COUNT(&permitted)); ++pass)
            for (const auto& group : groups)
                if (pass < group.second.cpuCount)
                    masks.push_back(group.second.mask);
    }
    return masks;
}
#endif

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
#if defined(__linux__)
        const auto masks = workerCacheMasks(workerCount);
#endif
        for (size_t worker = 0; worker != workerCount; ++worker)
        {
#if defined(__linux__)
            if (!masks.empty())
            {
                const auto mask = masks[worker % masks.size()];
                std::vector<size_t> peers;
                for (size_t offset = 1; offset < workerCount; ++offset)
                {
                    const auto other = (worker + offset) % workerCount;
                    if (CPU_EQUAL(&mask, &masks[other % masks.size()])) peers.push_back(other);
                }
                workers_.emplace_back([this, worker, mask, peers]() mutable {
                    if (sched_setaffinity(0, sizeof(mask), &mask) == 0)
                        cachePeers = std::move(peers);
                    run(worker);
                });
                continue;
            }
#endif
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
        Queue& queue = *queues_[target];
        std::lock_guard<std::mutex> queueLock(queue.mutex);
        rethrowIfFailed();
        if (stopping_.load(std::memory_order_acquire))
            throw std::logic_error("TaskExecutor is stopping");
        queue.work.emplace_back(std::move(task), std::move(completion));
        // Publish the count before unlocking the queue so a thief cannot take a
        // newly accepted task before it has been included in both counts.
        outstanding_.fetch_add(1, std::memory_order_relaxed);
        ready_.fetch_add(1, std::memory_order_seq_cst);
    }
    // Sequentially consistent registration and publication ensure that a
    // worker either sees the task on its final check or receives a wakeup.
    if (sleepers_.load(std::memory_order_seq_cst))
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        workAvailable_.notify_one();
    }
}

bool TaskExecutor::take(size_t workerId, size_t& nextVictim, Work& work)
{
    if (!ready_.load(std::memory_order_acquire))
        return false;
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
    auto steal = [&](size_t victim) {
        Queue& other = *queues_[victim];
        std::unique_lock<std::mutex> lock(other.mutex, std::try_to_lock);
        if (!lock.owns_lock() || other.work.empty())
            return false;
        work = std::move(other.work.front());
        other.work.pop_front();
        ready_.fetch_sub(1, std::memory_order_acq_rel);
        return true;
    };
#if defined(__linux__)
    // Locally queued chunks share catalog and index data. Always fall back
    // to global stealing when a skewed region exhausts the local queues.
    for (size_t scanned = 0; scanned < cachePeers.size(); ++scanned)
    {
        const size_t peer = cachePeers[nextCachePeer];
        if (++nextCachePeer == cachePeers.size()) nextCachePeer = 0;
        if (steal(peer)) return true;
    }
#endif
    for (size_t scanned = 0; scanned != queues_.size(); ++scanned)
    {
        const size_t victim = nextVictim;
        nextVictim = (nextVictim + 1) % queues_.size();
        if (victim != workerId && steal(victim)) return true;
    }
    return false;
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
    if (outstanding_.fetch_sub(1, std::memory_order_acq_rel) == 1)
    {
        // Synchronize with drain's transition into wait before notifying it.
        std::lock_guard<std::mutex> lock(stateMutex_);
        drained_.notify_all();
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
            sleepers_.fetch_add(1, std::memory_order_seq_cst);
            workAvailable_.wait(lock, [this] {
                return stopping_.load(std::memory_order_acquire) || ready_.load(std::memory_order_seq_cst) != 0;
            });
            sleepers_.fetch_sub(1, std::memory_order_seq_cst);
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
    if (!isCancelled())
        return;
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
    drained_.wait(lock, [this] { return outstanding_.load(std::memory_order_acquire) == 0; });
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
