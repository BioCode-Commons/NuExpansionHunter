//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Chris Saunders <csaunders@illumina.com>
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
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: streaming computes batches independently and reduces locus evidence in input order.
// Modifications are licensed under the Apache License, Version 2.0.
//

// Streaming execution: independently computed batches, ordered locus reduction.
#include "sample/HtsStreamingSampleAnalysis.hh"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <locale>
#include <mutex>
#include <unordered_map>

#include "core/TaskExecutor.hh"
#include "core/ThreadPool.hh"
#include "io/BamletWriter.hh"
#include "sample/GenomeQueryCollection.hh"
#include "sample/HtsFileStreamer.hh"
#include "sample/StreamingPairSource.hh"
#include "sample/ParallelPairSource.hh"
#include "spdlog/spdlog.h"

namespace ehunter
{
namespace
{
using Analyzer = locus::LocusEvidence;

// Construct a locus on the worker that first needs its indexes. The catalog
// drives read discovery independently, including loci with no admitted pairs.
class StreamingAnalyzers
{
public:
    StreamingAnalyzers(const RegionCatalog& catalog, const HeuristicParameters& params, locus::AlignWriterPtr writer)
        : catalog_(catalog), params_(params), writer_(std::move(writer)), entries_(catalog.size()), initialized_(catalog.size()) {}
    size_t size() const { return catalog_.size(); }

    struct Index
    {
        explicit Index(size_t locus) : locus(locus) {}
        size_t locus;
        std::once_flag ready;
        std::unique_ptr<locus::LocusAligner> aligner;
    };

    Analyzer& get(size_t index)
    {
        std::call_once(initialized_[index], [&, index] {
            entries_[index].reset(new Analyzer(catalog_[index]));
        });
        return *entries_[index];
    }
    const locus::LocusAligner& aligner(Index& index)
    {
        std::call_once(index.ready, [&] { index.aligner = buildAligner(index.locus); });
        return *index.aligner;
    }
    // Evidence rows retained for a locus, or none if no pair reached it.
    boost::optional<uint64_t> evidenceRows(size_t index) const
    {
        if (!entries_[index])
            return boost::none;
        return entries_[index]->evidenceRows();
    }
    LocusFindings finish(size_t index, Sex sex, ctpl::thread_pool* gridWorkers)
    {
        // Construct and validate indexes for loci that received no read batch.
        if (!entries_[index]) buildAligner(index);
        auto findings = get(index).analyze(sex, boost::none, gridWorkers);
        entries_[index].reset();
        return findings;
    }
private:
    std::unique_ptr<locus::LocusAligner> buildAligner(size_t index)
    {
        const auto& spec = get(index).locusSpec();
        return std::unique_ptr<locus::LocusAligner>(
            new locus::LocusAligner(spec.locusId(), &spec.regionGraph(), params_, writer_, uint64_t(index) * 2 + 1));
    }

    const RegionCatalog& catalog_;
    const HeuristicParameters& params_;
    locus::AlignWriterPtr writer_;
    std::vector<std::unique_ptr<Analyzer>> entries_;
    std::vector<std::once_flag> initialized_;
};

template <typename Function>
void forLoci(TaskExecutor& executor, size_t count, Function function, size_t grain = 0);

// Only active chunks have scheduling state. There is no mutex, sequence counter,
// or queue allocated for every catalog entry. Two prerequisites gate reduction:
// this chunk's computation and its preceding chunk's reduction. Workers never
// wait for another worker, even if an early alignment takes a very long time.
// Contract [K-35]
// Ensures:  1. each (pair, locus) item is aligned exactly once; 2. for every locus, items are reduced in the serial
//           order; 3. findings are stored by catalog index; 4. the BAM receives every record exactly once, in
//           unspecified order; 5. no deadlock under the admission limits, and termination at end of input.
//           Both front ends (fewer than 4 threads, and 4 or more) refine the same specification.
// Reference: the same code at --threads 1; model verification/models/pipeline.tla.
class OrderedLocusPipeline
{
    struct Input
    {
        Input(PackedRead& first, PackedRead& second) : read(first.releaseRead()), mate(second.releaseRead()) {}
        Read read, mate;
    };
    struct Item
    {
        std::shared_ptr<const Input> input;
        AnalyzerInputType inputType;
        locus::RegionType regionType;
        boost::optional<Read> read, mate;
        Analyzer::AlignedPair alignments;

        const Read& firstRead() const
        {
            if (read) return *read;
            return inputType == AnalyzerInputType::kMateOnly ? input->mate : input->read;
        }
        const Read* secondRead() const
        {
            if (mate) return &*mate;
            return input && inputType == AnalyzerInputType::kBothReads ? &input->mate : nullptr;
        }
    };
    struct Chunk
    {
        explicit Chunk(size_t locus) : locus(locus) {}
        size_t locus;
        uint64_t number = 0;
        std::vector<Item> items;
        std::shared_ptr<StreamingAnalyzers::Index> index;
        size_t charge = 0;
        std::atomic<unsigned> prerequisites { 2 };
        std::mutex linkMutex;
        std::shared_ptr<Chunk> successor;
        bool committed = false;
    };
    struct Completed
    {
        size_t locus;
        uint64_t number;
    };

public:
    OrderedLocusPipeline(TaskExecutor& executor, StreamingAnalyzers& analyzers,
                         graphtools::AlignerType type, locus::AlignWriterPtr writer, size_t stagingSlots)
        : executor_(executor), analyzers_(analyzers), writer_(dynamic_cast<BamletWriter*>(writer.get())),
          indexes_(analyzers.size()),
          recency_(analyzers.size() + 1), indexLimit_(std::min(analyzers.size(),
              executor.workerCount() > analyzers.size() / 8192 ? analyzers.size() : 8192 * executor.workerCount())),
          maxChunks_(32 * executor.workerCount()),
          cacheBudget_(40000000 * executor.workerCount()), publisherBudget_(cacheBudget_ / 8),
          stagingItems_(maxChunks_ * kChunkItems), stagingBytes_(cacheBudget_ / std::max(size_t(1), stagingSlots))
    {
        recency_.back() = { analyzers.size(), analyzers.size() };
        // AlignerSelector contains shared_ptr scratch. Construct each separately;
        // copying one prototype would introduce races in the dynamic programs.
        outputBuffers_.reserve(writer_ ? executor.workerCount() : 0);
        const std::locale analysisLocale;
        for (size_t worker = 0; worker < executor.workerCount(); ++worker)
        {
            if (writer_)
                outputBuffers_.emplace_back(writer_);
            selectors_.emplace_back(new graphtools::AlignerSelector(type));
            selectors_.back()->setSeedSearchLocale(analysisLocale);
        }
    }

    ~OrderedLocusPipeline()
    {
        if (!finished_)
        {
            executor_.cancel();
            try { executor_.wait(); } catch (...) {}
        }
    }

    // Memory a pair occupies once staged: its items, plus one decoded copy of
    // the reads (shared by all of the pair's loci, or moved into its one item).
    static size_t stagedBytes(const StreamingInputPair& pair)
    {
        return pair.bundles.size() * sizeof(Item) + (pair.bundles.size() > 1 ? sizeof(Input) : 0)
            + pair.read.retainedBytes() + pair.mate.retainedBytes();
    }

    static size_t itemCharge(const StreamingInputPair& pair)
    {
        return 1024 + 16 * (pair.read.length() + pair.mate.length())
            + 4 * (pair.read.fragmentId().size() + pair.mate.fragmentId().size());
    }

    template <class Append>
    static void distribute(StreamingInputPair& pair, Append append)
    {
        const size_t charge = itemCharge(pair);
        std::shared_ptr<const Input> shared;
        if (pair.bundles.size() > 1)
            shared = std::make_shared<Input>(pair.read, pair.mate);
        for (const auto& bundle : pair.bundles)
        {
            Item item {shared, bundle.inputType, bundle.regionType, boost::none, boost::none, {}};
            if (!shared)
            {
                auto& first = bundle.inputType == AnalyzerInputType::kMateOnly ? pair.mate : pair.read;
                item.read.emplace(first.releaseRead());
                if (bundle.inputType == AnalyzerInputType::kBothReads)
                    item.mate.emplace(pair.mate.releaseRead());
            }
            append(bundle.locusIndex, charge, std::move(item));
        }
    }

    void push(StreamingInputPair pair)
    {
        distribute(pair, [&](size_t locus, size_t charge, Item item) {
            auto& chunk = filling_[locus];
            if (!chunk)
                chunk = std::make_shared<Chunk>(locus);
            chunk->items.push_back(std::move(item));
            chunk->charge += charge;
            ++bufferedItems_;
            bufferedCharge_ += charge;
            if (chunk->items.size() == kChunkItems)
            {
                auto full = std::move(chunk);
                filling_.erase(locus);
                submit(std::move(full));
            }
            if (bufferedItems_ >= kStreamingBatchRecords || bufferedCharge_ >= publisherBudget_)
                flush();
        });
    }

    std::function<void()> prepareBatch(const std::vector<StreamingInputPair*>& pairs)
    {
        struct StagedLocus
        {
            size_t remaining = 0;
            std::shared_ptr<Chunk> chunk;
        };
        std::unordered_map<size_t, StagedLocus> partial;
        // Stage the longest prefix of the batch within both limits. The
        // publisher distributes any remainder itself, after the prefix, so
        // each locus still receives its pairs in input order.
        size_t items = 0, bytes = 0, staged = pairs.size();
        for (size_t index = 0; index < pairs.size(); ++index)
        {
            const auto* pair = pairs[index];
            if (!pair)
                continue;
            const size_t pairBytes = stagedBytes(*pair);
            if (pair->bundles.size() > stagingItems_ - items || pairBytes > stagingBytes_ - bytes)
            {
                staged = index;
                break;
            }
            items += pair->bundles.size();
            bytes += pairBytes;
            for (const auto& bundle : pair->bundles)
                ++partial[bundle.locusIndex].remaining;
        }
        using Chunks = std::vector<std::shared_ptr<Chunk>>;
        auto ready = std::make_shared<Chunks>();
        for (size_t index = 0; index < staged; ++index)
        {
            auto* pair = pairs[index];
            if (!pair)
                continue;
            distribute(*pair, [&](size_t locus, size_t charge, Item item) {
                auto& stage = partial.find(locus)->second;
                auto& chunk = stage.chunk;
                if (!chunk)
                {
                    chunk = std::make_shared<Chunk>(locus);
                    const auto capacity = std::min(size_t(kChunkItems), stage.remaining);
                    chunk->items.reserve(capacity);
                    stage.remaining -= capacity;
                }
                chunk->items.push_back(std::move(item));
                chunk->charge += charge;
                if (chunk->items.size() == kChunkItems)
                    ready->push_back(std::move(chunk));
            });
        }
        for (auto& entry : partial)
            if (entry.second.chunk)
                ready->push_back(std::move(entry.second.chunk));
        // Only the publishing thread links per-locus predecessors. The slot
        // stays alive until this closure has transferred every staged chunk
        // and distributed the remainder.
        return [this, ready, staged, &pairs] {
            // Reserve credits once for a bounded group. Workers can complete
            // published chunks while this thread links the new ones.
            const size_t budget = cacheBudget_ - publisherBudget_;
            for (size_t begin = 0; begin < ready->size();)
            {
                size_t end = begin + 1;
                size_t charge = (*ready)[begin]->charge;
                while (end < ready->size() && end - begin < 16 && charge < budget
                       && (*ready)[end]->charge <= budget - charge)
                    charge += (*ready)[end++]->charge;
                admit(end - begin, charge);
                for (; begin < end; ++begin)
                    schedule(std::move((*ready)[begin]));
            }
            if (staged == pairs.size())
                return;
            for (size_t index = staged; index < pairs.size(); ++index)
                if (pairs[index])
                    push(std::move(*pairs[index]));
            flush();
        };
    }

    void finish()
    {
        flush();
        executor_.wait();
        for (auto& buffer : outputBuffers_)
            buffer.flush();
        collectCompleted();
        if (inFlight_ != 0)
            throw std::logic_error("Streaming reduction did not drain");
        finished_ = true;
        selectors_.clear();
        // Active chunks have all drained. Retire index allocations in parallel
        // before graph owners are released during per-locus genotyping.
        forLoci(executor_, (indexes_.size() + 1023) / 1024, [this](size_t block) {
            const size_t end = std::min(indexes_.size(), 1024 * (block + 1));
            for (size_t i = 1024 * block; i < end; ++i) indexes_[i].reset();
        });
    }

private:
    void flush()
    {
        // No worker accesses this map. Preserve its small allocation across
        // flushes; a catalog with millions of cold loci creates no entries here.
        for (auto& entry : filling_)
            submit(std::move(entry.second));
        filling_.clear();
    }

    void collectCompleted()
    {
        collected_.clear();
        {
            std::lock_guard<std::mutex> lock(admissionMutex_);
            collected_.swap(completed_);
        }
        for (const auto& entry : collected_)
        {
            const auto tail = tails_.find(entry.locus);
            if (tail != tails_.end() && tail->second->number == entry.number)
                tails_.erase(tail);
        }
    }

    void submit(std::shared_ptr<Chunk> chunk)
    {
        bufferedItems_ -= chunk->items.size();
        bufferedCharge_ -= chunk->charge;
        admit(1, chunk->charge);
        schedule(std::move(chunk));
    }

    void admit(size_t count, size_t charge)
    {
        const size_t budget = cacheBudget_ - publisherBudget_;
        std::unique_lock<std::mutex> lock(admissionMutex_);
        while (inFlight_ && (inFlight_ > maxChunks_ - count
            || charge > budget - std::min(budget, inFlightCharge_)))
        {
            executor_.rethrowIfFailed();
            // On failure, a computation can terminate before releasing its
            // ordering token. Observe executor failure independently of tokens.
            const size_t burst = inFlight_ > 2 * executor_.workerCount()
                ? std::min(size_t(32), executor_.workerCount() / 2) : 1;
            wakeAtChunks_ = std::min(maxChunks_ - count, inFlight_ - std::min(inFlight_, std::max(size_t(1), burst)));
            wakeAtBytes_ = budget - std::min(budget, charge);
            waitingForSpace_ = true;
            admissionChanged_.wait_for(lock, std::chrono::milliseconds(10), [this] {
                return executor_.isCancelled() || !waitingForSpace_;
            });
            waitingForSpace_ = false;
        }
        executor_.rethrowIfFailed();
        inFlight_ += count;
        inFlightCharge_ += charge;
    }

    void schedule(std::shared_ptr<Chunk> chunk)
    {
        chunk->number = nextNumber_++;
        auto& index = indexes_[chunk->locus];
        // Only the ordered publisher touches the links. A chunk holds its own
        // shared index until reduction, including after eviction from this cache.
        const size_t root = indexes_.size();
        const auto unlink = [this](size_t locus) {
            const auto links = recency_[locus];
            recency_[links.previous].next = links.next;
            recency_[links.next].previous = links.previous;
        };
        if (!index)
        {
            if (cachedIndexes_ == indexLimit_)
            {
                const size_t oldest = recency_[root].previous;
                unlink(oldest);
                indexes_[oldest].reset();
            }
            else ++cachedIndexes_;
            index = std::make_shared<StreamingAnalyzers::Index>(chunk->locus);
        }
        else unlink(chunk->locus);
        const size_t newest = recency_[root].next;
        recency_[chunk->locus] = { root, newest };
        recency_[newest].previous = chunk->locus;
        recency_[root].next = chunk->locus;
        chunk->index = index;
        // A completed tail has already released all item data. Reclaim these
        // small ordering entries in batches to amortize producer/worker locking.
        if (chunk->number % 32 == 0)
            collectCompleted();
        // When workers lack ready work, split a large locus batch without
        // splitting its ordered evidence reduction.
        const size_t workers = executor_.workerCount();
        const size_t spare = workers - std::min(workers - 1, executor_.queuedTasks());
        const size_t segments = std::min(spare, (chunk->items.size() + 31) / 32);
        chunk->prerequisites.store(segments + 1, std::memory_order_relaxed);
        auto previous = tails_.find(chunk->locus);
        if (previous == tails_.end())
            chunk->prerequisites.store(segments, std::memory_order_relaxed);
        else
        {
            std::lock_guard<std::mutex> lock(previous->second->linkMutex);
            if (previous->second->committed)
                chunk->prerequisites.store(segments, std::memory_order_relaxed);
            else
                previous->second->successor = chunk;
        }
        tails_[chunk->locus] = chunk;
        const size_t stride = (chunk->items.size() + segments - 1) / segments;
        const size_t itemCount = chunk->items.size();
        for (size_t segment = 0; segment < segments; ++segment)
        {
            const size_t first = segment * stride, last = std::min(first + stride, itemCount);
            executor_.submit([this, chunk, first, last](size_t worker) {
                const auto& aligner = analyzers_.aligner(*chunk->index);
                BamletWriter::Batch output(writer_, writer_ ? &outputBuffers_[worker] : nullptr);
                for (size_t i = first; i < last; ++i)
                {
                    auto& item = chunk->items[i];
                    if (item.regionType == locus::RegionType::kTarget)
                    {
                        if (item.input)
                            item.alignments = aligner.alignSharedMates(item.firstRead(), item.secondRead(),
                                item.read, item.mate, *selectors_[worker]);
                        else
                            item.alignments = aligner.alignMates(
                                *item.read, item.mate ? &*item.mate : nullptr, *selectors_[worker]);
                        if (writer_)
                            aligner.writeAlignedPair(item.firstRead(), item.secondRead(), item.alignments);
                    }
                }
                output.finish();
                if (chunk->prerequisites.fetch_sub(1, std::memory_order_acq_rel) == 1)
                    reduce(chunk);
            }, (chunk->locus + segment) % workers);
        }
    }

    void reduce(std::shared_ptr<Chunk> chunk)
    {
        // Iterative continuation avoids recursion when many later chunks have
        // already finished. The acquire/release prerequisite counter publishes
        // every computed result to whichever worker performs the reduction.
        while (chunk)
        {
            auto& analyzer = analyzers_.get(chunk->locus);
            for (auto& item : chunk->items)
            {
                if (!writer_ && item.regionType == locus::RegionType::kTarget)
                    chunk->index->aligner->writeAlignedPair(item.firstRead(), item.secondRead(), item.alignments);
                analyzer.consumeEvidence(item.firstRead(), item.secondRead(), item.regionType, item.alignments);
                // Output records and retained evidence own their data. Release
                // transient storage before accumulating the next pair's evidence.
                item.alignments.first = boost::none;
                item.alignments.second = boost::none;
                item.read = boost::none;
                item.mate = boost::none;
                item.input.reset();
            }
            std::vector<Item>().swap(chunk->items);
            std::shared_ptr<Chunk> successor;
            {
                std::lock_guard<std::mutex> lock(chunk->linkMutex);
                chunk->committed = true;
                successor = std::move(chunk->successor);
            }
            bool wakePublisher = false;
            {
                std::lock_guard<std::mutex> lock(admissionMutex_);
                --inFlight_;
                inFlightCharge_ -= chunk->charge;
                completed_.push_back(Completed { chunk->locus, chunk->number });
                // Wake the publisher once for several free slots when workers
                // already have a backlog. Near drain a single completion wakes it.
                if (waitingForSpace_ && inFlight_ <= wakeAtChunks_ && inFlightCharge_ <= wakeAtBytes_)
                {
                    waitingForSpace_ = false;
                    wakePublisher = true;
                }
            }
            if (wakePublisher) admissionChanged_.notify_one();
            chunk.reset();
            if (!successor || successor->prerequisites.fetch_sub(1, std::memory_order_acq_rel) != 1)
                break;
            chunk = std::move(successor);
        }
    }

    TaskExecutor& executor_;
    StreamingAnalyzers& analyzers_;
    BamletWriter* writer_;
    std::vector<BamletWriter::BatchBuffer> outputBuffers_;
    std::vector<std::unique_ptr<graphtools::AlignerSelector>> selectors_;
    std::unordered_map<size_t, std::shared_ptr<Chunk>> filling_, tails_;
    uint64_t nextNumber_ = 0;
    size_t bufferedItems_ = 0, bufferedCharge_ = 0;
    std::vector<std::shared_ptr<StreamingAnalyzers::Index>> indexes_;
    struct Links { size_t previous, next; };
    std::vector<Links> recency_;
    const size_t indexLimit_;
    size_t cachedIndexes_ = 0;
    const size_t maxChunks_, cacheBudget_, publisherBudget_;
    // Per input batch: at most one pipeline capacity of items, and a share of
    // the cache budget such that all concurrently staged batches fit in it.
    const size_t stagingItems_, stagingBytes_;
    std::mutex admissionMutex_;
    std::condition_variable admissionChanged_;
    size_t inFlight_ = 0, inFlightCharge_ = 0;
    size_t wakeAtChunks_ = 0, wakeAtBytes_ = 0;
    bool waitingForSpace_ = false;
    std::vector<Completed> completed_, collected_;
    bool finished_ = false;
    static constexpr size_t kChunkItems = 128;
};

// One task per worker consumes blocks from an atomic cursor. There is no
// unbounded task allocation proportional to the number of catalog loci.
template <typename Function>
void forLoci(TaskExecutor& executor, size_t count, Function function, size_t grain)
{
    std::atomic<size_t> next { 0 };
    if (!grain)
        grain = std::max(size_t(1), std::min(size_t(32), count / (8 * executor.workerCount())));
    try
    {
        for (size_t worker = 0; worker < executor.workerCount(); ++worker)
            executor.submit([&](size_t) {
                while (!executor.isCancelled())
                {
                    const size_t begin = next.fetch_add(grain, std::memory_order_relaxed);
                    if (begin >= count)
                        return;
                    for (size_t index = begin; index < std::min(begin + grain, count); ++index)
                        function(index);
                }
            }, worker);
        executor.wait();
    }
    catch (...)
    {
        const auto error = std::current_exception();
        executor.cancel(error);
        try { executor.wait(); } catch (...) {}
        std::rethrow_exception(error);
    }
}

}

SampleFindings htsStreamingSampleAnalysis(
    const InputPaths& inputPaths, Sex sampleSex, const HeuristicParameters& heuristicParams, int threadCount,
    const RegionCatalog& regionCatalog, locus::AlignWriterPtr alignmentWriter, unsigned decodeThreadCount)
{
    if (threadCount < 1)
        throw std::invalid_argument("Streaming requires at least one compute worker");
    spdlog::info("Initializing all loci");
    TaskExecutor executor(static_cast<size_t>(threadCount));
    StreamingAnalyzers analyzers(regionCatalog, heuristicParams, alignmentWriter);
    GenomeQueryCollection queries(regionCatalog, static_cast<unsigned>(threadCount));

    spdlog::info("Streaming reads");
    {
        htshelpers::HtsFileStreamer input(inputPaths.htsFile(), inputPaths.reference(), decodeThreadCount);
        const size_t laneCount = std::max(2, threadCount / 8);
        OrderedLocusPipeline pipeline(
            executor, analyzers, heuristicParams.alignerType(), alignmentWriter,
            threadCount >= 4 ? streamingWindowSlots(laneCount) : 1);
        auto check = [&executor] { executor.rethrowIfFailed(); };
        if (threadCount >= 4)
        {
            streamParallelPairs(input, queries.targetRegionMask, queries.analyzerFinder, laneCount, check,
                decodeThreadCount,
                [&](const std::vector<StreamingInputPair*>& pairs) { return pipeline.prepareBatch(pairs); });
        }
        else
        {
            StreamingPairSource source(input, queries.targetRegionMask, queries.analyzerFinder,
                StreamingPairSource::defaultPendingRawBytes, check);
            while (auto pair = source.next())
                pipeline.push(std::move(*pair));
        }
        pipeline.finish();
    }

    spdlog::info("Analyzing read evidence");
    SampleFindings findings(regionCatalog.size());
    std::unique_ptr<ctpl::thread_pool> gridWorkers;
    if (threadCount > 1)
        gridWorkers.reset(new ctpl::thread_pool(threadCount));
    // Heaviest loci first, one locus per grab: a few loci with hundreds of
    // thousands of evidence rows take most of the genotyping time. Loci are
    // ordered by the binary magnitude of their rows (stable within a class);
    // loci that no pair reached come last. Findings are stored by index, so
    // the order does not affect output.
    std::vector<uint32_t> order(regionCatalog.size());
    {
        const size_t classes = 66; // log2 classes of a 64-bit count, plus loci without evidence
        std::vector<uint8_t> weight(regionCatalog.size());
        forLoci(executor, regionCatalog.size(), [&](size_t index) {
            const auto rows = analyzers.evidenceRows(index);
            uint8_t magnitude = 0;
            if (rows)
                for (uint64_t value = *rows + 1; value; value >>= 1)
                    ++magnitude;
            weight[index] = magnitude;
        });
        std::vector<size_t> starts(classes + 1, 0);
        for (const auto magnitude : weight)
            ++starts[classes - magnitude];
        size_t position = 0;
        for (auto& start : starts)
        {
            const size_t count = start;
            start = position;
            position += count;
        }
        for (size_t index = 0; index < weight.size(); ++index)
            order[starts[classes - weight[index]]++] = uint32_t(index);
    }
    forLoci(executor, order.size(), [&](size_t position) {
        const size_t index = order[position];
        findings[index] = analyzers.finish(index, sampleSex, gridWorkers.get());
    }, 1);
    return findings;
}
}
