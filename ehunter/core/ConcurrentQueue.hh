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

#pragma once

#include <condition_variable>
#include <cstddef>
#include <memory>
#include <mutex>
#include <queue>

namespace ehunter
{

template <typename T> class ConcurrentQueue
{
public:
    // A zero budget preserves the unbounded queue. Weights can represent bytes
    // rather than item counts. One oversized item is allowed in an empty queue
    // so a large record can always make progress.
    explicit ConcurrentQueue(std::size_t maxWeight = 0) : maxWeight_(maxWeight)
    {
        if (maxWeight_)
            weightedQueue_.reset(new std::queue<Entry>());
    }

    void push(T const& data, std::size_t weight = 1)
    {
        if (!maxWeight_)
        {
            // Keep the legacy unbounded path free of per-record weight storage,
            // budget bookkeeping and notifications to blocked producers.
            {
                std::lock_guard<std::mutex> lock(mutex_);
                queue_.push(data);
            }
            cv_.notify_one();
            return;
        }
        {
            std::unique_lock<std::mutex> lock(mutex_);
            spaceAvailable_.wait(lock, [&] {
                return weightedQueue_->empty()
                    || (weight <= maxWeight_ && queuedWeight_ <= maxWeight_ - weight);
            });
            weightedQueue_->emplace(data, weight);
            queuedWeight_ += weight;
        }
        cv_.notify_one();
    }

    bool empty() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return maxWeight_ ? weightedQueue_->empty() : queue_.empty();
    }

    void pop(T& value)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        if (!maxWeight_)
        {
            while (queue_.empty())
                cv_.wait(lock);
            value = queue_.front();
            queue_.pop();
            return;
        }
        while (weightedQueue_->empty())
            cv_.wait(lock);
        value = weightedQueue_->front().data;
        queuedWeight_ -= weightedQueue_->front().weight;
        weightedQueue_->pop();
        lock.unlock();
        // Producers can have different weights; wake all to let one that fits
        // proceed. The consumer must keep draining, including after I/O errors.
        spaceAvailable_.notify_all();
    }

private:
    struct Entry
    {
        Entry(const T& value, std::size_t itemWeight) : data(value), weight(itemWeight) {}
        T data;
        std::size_t weight;
    };

    const std::size_t maxWeight_;
    std::size_t queuedWeight_ = 0;
    std::queue<T> queue_;
    std::unique_ptr<std::queue<Entry>> weightedQueue_;
    mutable std::mutex mutex_;
    std::condition_variable cv_;
    std::condition_variable spaceAvailable_;
};

}
