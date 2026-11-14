// Streaming execution: independently computed batches, ordered locus reduction.
#include "sample/HtsStreamingSampleAnalysis.hh"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <limits>
#include <locale>
#include <mutex>
#include <unordered_map>

#include "core/TaskExecutor.hh"
#include "io/BamletWriter.hh"
#include "sample/GenomeQueryCollection.hh"
#include "sample/HtsFileStreamer.hh"
#include "sample/StreamingPairSource.hh"
#include "spdlog/spdlog.h"

namespace ehunter
{
namespace
{
using Clock = std::chrono::steady_clock;
using Analyzer = locus::LocusAnalyzer;

double elapsed(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

// Only active chunks have scheduling state. There is no mutex, sequence counter,
// or queue allocated for every catalog entry. Two prerequisites gate reduction:
// this chunk's computation and its preceding chunk's reduction. Workers never
// wait for another worker, even if an early alignment takes a very long time.
class OrderedLocusPipeline
{
    struct Input
    {
        explicit Input(StreamingInputPair pair) : packed(std::move(pair)) {}
        StreamingInputPair packed;
        std::once_flag decoded;
        std::exception_ptr decodeError;
        boost::optional<Read> read, mate;
        void decode()
        {
            // A heavily overlapping catalog still decodes a pair only once.
            // Every locus gets private oriented Read objects afterward.
            std::call_once(decoded, [this] {
                try
                {
                    read.emplace(packed.read.releaseRead());
                    mate.emplace(packed.mate.releaseRead());
                }
                catch (...)
                {
                    // releaseRead consumes its payload. Publish a terminal
                    // failure instead of letting call_once retry partial input.
                    decodeError = std::current_exception();
                }
            });
            if (decodeError)
                std::rethrow_exception(decodeError);
        }
    };
    struct Item
    {
        std::shared_ptr<Input> input;
        AnalyzerInputType inputType;
        locus::RegionType regionType;
        bool soleUse;
        boost::optional<Read> read, mate;
        Analyzer::AlignedPair alignments;
    };
    struct Chunk
    {
        explicit Chunk(size_t locus, uint64_t number)
            : locus(locus), number(number)
        {
        }
        size_t locus;
        uint64_t number;
        std::vector<Item> items;
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
    OrderedLocusPipeline(TaskExecutor& executor, std::vector<std::unique_ptr<Analyzer>>& analyzers,
                         graphtools::AlignerType type, locus::AlignWriterPtr writer)
        : executor_(executor), analyzers_(analyzers), writer_(dynamic_cast<BamletWriter*>(writer.get())),
          maxChunks_(4 * executor.workerCount()),
          cacheBudget_(20000000 * executor.workerCount()), publisherBudget_(std::min(size_t(8000000), cacheBudget_ / 2))
    {
        // AlignerSelector contains shared_ptr scratch. Construct each separately;
        // copying one prototype would introduce races in the dynamic programs.
        const std::locale analysisLocale;
        for (size_t worker = 0; worker <= executor.workerCount(); ++worker)
        {
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

    void push(StreamingInputPair pair)
    {
        auto input = std::make_shared<Input>(std::move(pair));
        const bool soleUse = input->packed.bundles.size() == 1;
        // Conservative admission charge for retained inputs, oriented reads and
        // ordinary alignment results. Actual process RSS is separately measured:
        // allocator slack, DP scratch and arbitrarily large graph results are not
        // represented by this charge. Never truncate an oversized read/result.
        const size_t bases = input->packed.read.length() + input->packed.mate.length();
        const size_t names = input->packed.read.fragmentId().size() + input->packed.mate.fragmentId().size();
        const size_t charge = 1024 + 16 * bases + 4 * names;
        for (const auto& bundle : input->packed.bundles)
        {
            auto& chunk = filling_[bundle.locusIndex];
            if (!chunk)
                chunk = std::make_shared<Chunk>(bundle.locusIndex, nextNumber_++);
            chunk->items.push_back(Item { input, bundle.inputType, bundle.regionType, soleUse,
                                         boost::none, boost::none, {} });
            chunk->charge += charge;
            ++bufferedItems_;
            bufferedCharge_ += charge;
            ++assignments_;
            if (chunk->items.size() == 64)
            {
                auto full = std::move(chunk);
                filling_.erase(bundle.locusIndex);
                submit(std::move(full));
            }
            // Bound the publisher too, including a pair with very large fanout.
            // Flushing keeps each locus's encounter order; inter-locus order has
            // no observable meaning except the permitted BAM record order.
            if (bufferedItems_ >= 4096 || bufferedCharge_ >= publisherBudget_)
                flush();
        }
        ++pairs_;
    }

    void finish()
    {
        flush();
        executor_.wait();
        collectCompleted();
        if (inFlight_ != 0)
            throw std::logic_error("Streaming reduction did not drain");
        finished_ = true;
        selectors_.clear();
        spdlog::info(
            "Streaming batches: pairs={} assignments={} chunks={} peak_active_chunks={} "
            "peak_charged_cache_bytes={} cache_budget_bytes={} admission_wait_seconds={} "
            "alignment_worker_seconds={} reduction_worker_seconds={} publisher_assisted_chunks={}",
            pairs_, assignments_, chunks_, peakChunks_, peakCharge_, cacheBudget_, admissionWait_,
            double(alignmentNanos_.load()) / 1e9, double(reductionNanos_.load()) / 1e9, assistedChunks_);
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
        std::vector<Completed> completed;
        {
            std::lock_guard<std::mutex> lock(admissionMutex_);
            completed.swap(completed_);
        }
        for (const auto& entry : completed)
        {
            const auto tail = tails_.find(entry.locus);
            if (tail != tails_.end() && tail->second->number == entry.number)
                tails_.erase(tail);
        }
    }

    void submit(std::shared_ptr<Chunk> chunk)
    {
        collectCompleted();
        {
            std::unique_lock<std::mutex> lock(admissionMutex_);
            while (inFlight_ && (inFlight_ >= maxChunks_
                || inFlightCharge_ + chunk->charge > cacheBudget_ - publisherBudget_))
            {
                executor_.rethrowIfFailed();
                lock.unlock();
                const bool assisted = executor_.tryRunOne();
                lock.lock();
                if (assisted)
                {
                    ++assistedChunks_;
                    continue;
                }
                // On failure, a computation can terminate before releasing its
                // ordering token. Observe executor failure independently of tokens.
                const auto started = Clock::now();
                admissionChanged_.wait_for(lock, std::chrono::milliseconds(10), [this, &chunk] {
                    return executor_.isCancelled() || inFlight_ == 0
                        || (inFlight_ < maxChunks_
                            && inFlightCharge_ + chunk->charge <= cacheBudget_ - publisherBudget_);
                });
                admissionWait_ += elapsed(started);
            }
            executor_.rethrowIfFailed();
            ++inFlight_;
            inFlightCharge_ += chunk->charge;
            peakChunks_ = std::max(peakChunks_, inFlight_);
            peakCharge_ = std::max(peakCharge_, inFlightCharge_ + bufferedCharge_ - chunk->charge);
        }
        bufferedItems_ -= chunk->items.size();
        bufferedCharge_ -= chunk->charge;
        auto previous = tails_.find(chunk->locus);
        if (previous == tails_.end())
            chunk->prerequisites.store(1, std::memory_order_relaxed);
        else
        {
            std::lock_guard<std::mutex> lock(previous->second->linkMutex);
            if (previous->second->committed)
                chunk->prerequisites.store(1, std::memory_order_relaxed);
            else
                previous->second->successor = chunk;
        }
        tails_[chunk->locus] = chunk;
        ++chunks_;
        executor_.submit([this, chunk](size_t worker) {
            const auto start = Clock::now();
            auto& analyzer = *analyzers_[chunk->locus];
            for (auto& item : chunk->items)
            {
                const bool mateOnly = item.inputType == AnalyzerInputType::kMateOnly;
                item.input->decode();
                Read& source = mateOnly ? *item.input->mate : *item.input->read;
                if (item.soleUse)
                    item.read.emplace(std::move(source));
                else
                    item.read.emplace(source);
                if (item.inputType == AnalyzerInputType::kBothReads)
                {
                    if (item.soleUse)
                        item.mate.emplace(std::move(*item.input->mate));
                    else
                        item.mate.emplace(*item.input->mate);
                }
                item.input.reset();
                item.alignments = analyzer.computeMates(
                    *item.read, item.mate ? &*item.mate : nullptr, item.regionType, *selectors_[worker]);
            }
            alignmentNanos_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count(),
                                     std::memory_order_relaxed);
            if (chunk->prerequisites.fetch_sub(1, std::memory_order_acq_rel) == 1)
                reduce(chunk);
        }, TaskExecutor::anyWorker(), [this] { admissionChanged_.notify_one(); });
    }

    void reduce(std::shared_ptr<Chunk> chunk)
    {
        // Iterative continuation avoids recursion when many later chunks have
        // already finished. The acquire/release prerequisite counter publishes
        // every computed result to whichever worker performs the reduction.
        const auto start = Clock::now();
        while (chunk)
        {
            auto& analyzer = *analyzers_[chunk->locus];
            BamletWriter::Batch output(writer_);
            for (auto& item : chunk->items)
            {
                analyzer.consumeMates(*item.read, item.mate ? &*item.mate : nullptr,
                                      item.regionType, item.alignments);
                // Output records and retained evidence own their data. Release
                // transient storage before accumulating the next pair's evidence.
                item.alignments.first = boost::none;
                item.alignments.second = boost::none;
                item.read = boost::none;
                item.mate = boost::none;
                item.input.reset();
            }
            output.finish();
            std::vector<Item>().swap(chunk->items);
            std::shared_ptr<Chunk> successor;
            {
                std::lock_guard<std::mutex> lock(chunk->linkMutex);
                chunk->committed = true;
                successor = std::move(chunk->successor);
            }
            {
                std::lock_guard<std::mutex> lock(admissionMutex_);
                --inFlight_;
                inFlightCharge_ -= chunk->charge;
                completed_.push_back(Completed { chunk->locus, chunk->number });
            }
            admissionChanged_.notify_one();
            chunk.reset();
            if (!successor || successor->prerequisites.fetch_sub(1, std::memory_order_acq_rel) != 1)
                break;
            chunk = std::move(successor);
        }
        reductionNanos_.fetch_add(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count(),
                                 std::memory_order_relaxed);
    }

    TaskExecutor& executor_;
    std::vector<std::unique_ptr<Analyzer>>& analyzers_;
    BamletWriter* writer_;
    std::vector<std::unique_ptr<graphtools::AlignerSelector>> selectors_;
    std::unordered_map<size_t, std::shared_ptr<Chunk>> filling_, tails_;
    uint64_t nextNumber_ = 0, pairs_ = 0, assignments_ = 0, chunks_ = 0, assistedChunks_ = 0;
    size_t bufferedItems_ = 0, bufferedCharge_ = 0;
    const size_t maxChunks_, cacheBudget_, publisherBudget_;
    std::mutex admissionMutex_;
    std::condition_variable admissionChanged_;
    size_t inFlight_ = 0, inFlightCharge_ = 0, peakChunks_ = 0, peakCharge_ = 0;
    std::vector<Completed> completed_;
    double admissionWait_ = 0;
    std::atomic<uint64_t> alignmentNanos_ { 0 }, reductionNanos_ { 0 };
    bool finished_ = false;
};

// One task per worker consumes blocks from an atomic cursor. There is no
// unbounded task allocation proportional to the number of catalog loci.
template <typename Function>
void forLoci(TaskExecutor& executor, size_t count, Function function)
{
    std::atomic<size_t> next { 0 };
    const size_t grain = std::max(size_t(1), std::min(size_t(32), count / (8 * executor.workerCount())));
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
    const RegionCatalog& regionCatalog, locus::AlignWriterPtr alignmentWriter, unsigned decodeThreadCount,
    unsigned readDecodeThreadCount)
{
    if (threadCount < 1)
        throw std::invalid_argument("Streaming requires at least one compute worker");
    // The compatibility hint no longer creates an independent preprocessing
    // pool: pure decoding now runs on the shared executor for every input.
    if (readDecodeThreadCount)
        spdlog::info("Streaming preprocessing hint={} (decoding uses shared compute workers)", readDecodeThreadCount);
    const auto indexStart = Clock::now();
    TaskExecutor executor(static_cast<size_t>(threadCount));
    std::vector<std::unique_ptr<Analyzer>> analyzers(regionCatalog.size());
    forLoci(executor, regionCatalog.size(), [&](size_t index) {
        analyzers[index].reset(new Analyzer(regionCatalog[index], heuristicParams, alignmentWriter,
                                             uint64_t(index) * 2 + 1));
    });
    GenomeQueryCollection queries(analyzers, &regionCatalog);
    spdlog::info("Phase catalog_index_seconds={}", elapsed(indexStart));

    spdlog::info("Streaming reads with ordered batch reduction");
    const auto readStart = Clock::now();
    {
        htshelpers::HtsFileStreamer input(inputPaths.htsFile(), inputPaths.reference(), decodeThreadCount);
        StreamingPairSource source(input, queries.targetRegionMask, queries.analyzerFinder,
            StreamingPairSource::defaultPendingRawBytes, [&executor] {
                if (executor.isCancelled())
                    executor.rethrowIfFailed();
            });
        OrderedLocusPipeline pipeline(executor, analyzers, heuristicParams.alignerType(), alignmentWriter);
        while (auto pair = source.next())
            pipeline.push(std::move(*pair));
        pipeline.finish();
        spdlog::info("Streaming pending mates: peak_reads={} peak_raw_bytes={} peak_cold_reads={} cold_records={}",
                     source.peakPendingReadCount(), source.peakPendingRawBytes(),
                     source.peakColdPendingReadCount(), source.coldReadCount());
    }
    spdlog::info("Phase read_processing_seconds={}", elapsed(readStart));

    const auto genotypeStart = Clock::now();
    SampleFindings findings(regionCatalog.size());
    forLoci(executor, regionCatalog.size(), [&](size_t index) {
        findings[index] = analyzers[index]->analyze(sampleSex, boost::none);
        analyzers[index].reset();
    });
    spdlog::info("Phase genotyping_seconds={}", elapsed(genotypeStart));
    return findings;
}
}
