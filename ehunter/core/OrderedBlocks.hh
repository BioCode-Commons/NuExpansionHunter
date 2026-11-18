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

#include <algorithm>
#include <cstddef>
#include <deque>
#include <future>
#include <memory>

#include "core/ThreadPool.hh"

namespace ehunter
{
// Compute disjoint blocks with worker-private state and consume their results
// in source order. Limit live blocks even when the first one takes longer.
// Contract [K-36]
// Ensures:  consume runs on the calling thread for every block in index order, for every schedule of compute.
template <typename Compute, typename Consume>
void orderedBlocks(size_t size, size_t grain, size_t workers, Compute compute, Consume consume,
                   ctpl::thread_pool* shared = nullptr)
{
    if (!size)
        return;
    workers = std::min(workers, 1 + (size - 1) / grain);
    if (workers <= 1)
    {
        for (size_t begin = 0; begin < size; begin += grain)
            consume(compute(0, begin, std::min(size, begin + grain)));
        return;
    }
    using Result = decltype(compute(0, size, size));
    std::deque<std::future<Result>> pending;
    std::unique_ptr<ctpl::thread_pool> owned;
    if (!shared || !shared->size())
    {
        owned.reset(new ctpl::thread_pool(static_cast<int>(workers)));
        shared = owned.get();
    }
    auto consumeFirst = [&] {
        consume(pending.front().get());
        pending.pop_front();
    };
    try
    {
        for (size_t begin = 0; begin < size; begin += grain)
        {
            // Allocate the future slot before handing borrowed state to a worker.
            pending.emplace_back();
            pending.back() = shared->push(compute, begin, std::min(size, begin + grain));
            if (pending.size() >= 2 * workers) consumeFirst();
        }
        while (!pending.empty()) consumeFirst();
    }
    catch (...)
    {
        // A shared executor can outlive this formatter and its worker scratch.
        for (auto& task : pending)
            if (task.valid()) task.wait();
        throw;
    }
}
}
