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

#include "sample/HtsStreamingReadPairQueue.hh"

namespace ehunter
{

HtsStreamingReadPairQueue::InsertResult HtsStreamingReadPairQueue::insertReadPair(
    const unsigned locusIndex, ReadPair readPair, const WaitHelper& waitHelper)
{
    auto& locusAnalyzerQueue(queues_[locusIndex]);
    std::unique_lock<std::mutex> locusAnalyzerQueueLock(queueMutex(locusIndex));
    const bool wasInActive(not locusAnalyzerQueue);
    if (wasInActive)
    {
        std::unique_lock<std::mutex> globalLock(mutex_);
        while (activeLocusAnalyzerQueues_ >= maxActiveLocusAnalyzerQueues_ && !cancelled_)
        {
            // An active locus can share this stripe. It must be able to drain
            // while this inactive locus waits for capacity. The sole producer
            // is the only thread that can activate the inactive locus.
            locusAnalyzerQueueLock.unlock();
            bool helped = false;
            if (waitHelper)
            {
                globalLock.unlock();
                helped = waitHelper();
                globalLock.lock();
            }
            // Recheck under the mutex: a worker may have freed capacity (or
            // cancelled processing) while the helper was looking for work.
            if (!helped && activeLocusAnalyzerQueues_ >= maxActiveLocusAnalyzerQueues_ && !cancelled_)
                cv_.wait(globalLock);
            // Restore stripe -> global acquisition order after waiting/helping.
            globalLock.unlock();
            locusAnalyzerQueueLock.lock();
            globalLock.lock();
        }
        if (cancelled_)
            return InsertResult::Cancelled;
        if (idleQueues_.empty())
            locusAnalyzerQueue.reset(new ReadPairQueue);
        else
        {
            locusAnalyzerQueue = std::move(idleQueues_.back());
            idleQueues_.pop_back();
        }
        activeLocusAnalyzerQueues_++;
        globalLock.unlock();
    }
    locusAnalyzerQueue->emplace(std::move(readPair));
    return wasInActive ? InsertResult::Activated : InsertResult::Queued;
}

void HtsStreamingReadPairQueue::cancel()
{
    {
        std::lock_guard<std::mutex> lock(mutex_);
        cancelled_ = true;
    }
    cv_.notify_all();
}

void HtsStreamingReadPairQueue::getNextReadPair(const unsigned locusIndex, boost::optional<ReadPair>& readPair)
{
    auto& locusAnalyzerQueue(queues_[locusIndex]);
    std::unique_lock<std::mutex> locusAnalyzerQueueLock(queueMutex(locusIndex));

    if (locusAnalyzerQueue->empty())
    {
        std::unique_lock<std::mutex> globalLock(mutex_);
        idleQueues_.push_back(std::move(locusAnalyzerQueue));
        activeLocusAnalyzerQueues_--;
        globalLock.unlock();
        locusAnalyzerQueueLock.unlock();

        cv_.notify_one();
        readPair = boost::none;
    }
    else
    {
        readPair = std::move(locusAnalyzerQueue->front());
        locusAnalyzerQueue->pop();
    }
}

}
