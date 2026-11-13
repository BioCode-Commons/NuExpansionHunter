// Modified for NuExpansionHunter (2026): transfer cached mates without copying sequence buffers.
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

#include "sample/HtsStreamingSampleAnalysis.hh"

#include <memory>

#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "absl/strings/string_view.h"
#include "spdlog/spdlog.h"
#include <boost/optional.hpp>

#include "core/HtsHelpers.hh"
#include "core/ThreadPool.hh"
#include "locus/LocusAnalyzer.hh"
#include "locus/LocusAnalyzerUtil.hh"
#include "sample/GenomeQueryCollection.hh"
#include "sample/HtsFileStreamer.hh"
#include "sample/HtsStreamingReadPairQueue.hh"
#include "sample/DecodedReadStream.hh"
#include <chrono>

using ehunter::locus::initializeLocusAnalyzers;
using ehunter::locus::LocusAnalyzer;
using graphtools::AlignmentWriter;
using std::string;
using std::vector;

namespace ehunter
{

namespace
{

// Looking up a borrowed BAM qname avoids constructing a Read for a second
// record whose pair has no exact target/off-target routing. Both overloads use
// the identical hash; the existing fragment-only pairing rule is unchanged.
struct FragmentHash
{
    using is_transparent = void;
    size_t operator()(absl::string_view name) const { return absl::Hash<absl::string_view>()(name); }
    size_t operator()(const Read& read) const { return (*this)(absl::string_view(read.fragmentId())); }
};

struct FragmentEqual
{
    using is_transparent = void;
    bool operator()(const Read& a, const Read& b) const { return a.fragmentId() == b.fragmentId(); }
    bool operator()(const Read& a, absl::string_view b) const { return a.fragmentId() == b; }
    bool operator()(absl::string_view a, const Read& b) const { return a == b.fragmentId(); }
};

/// \brief Mutable data shared by all LocusAnalyzer-processing threads
///
class LocusAnalyzerThreadSharedData
{
public:
    LocusAnalyzerThreadSharedData(const unsigned maxActiveLocusAnalyzerQueues, const unsigned locusAnalyzerCount)
        : isWorkerThreadException(false)
        , readPairQueue(maxActiveLocusAnalyzerQueues, locusAnalyzerCount)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    HtsStreamingReadPairQueue readPairQueue;
    vector<std::unique_ptr<LocusAnalyzer>> locusAnalyzers;
};

/// \brief Data isolated to each LocusAnalyzer-processing thread
///
struct LocusAnalyzerThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
    std::shared_ptr<graphtools::AlignerSelector> alignerSelectorPtr;
};

/// Process queue for a single LocusAnalyzer on one thread
void processLocusAnalyzerQueue(
    const int threadIndex, LocusAnalyzerThreadSharedData& locusAnalyzerThreadSharedData,
    std::vector<LocusAnalyzerThreadLocalData>& locusAnalyzerThreadLocalDataPool, const unsigned locusIndex,
    ctpl::thread_pool& pool, const int readerThreadIndex)
{
    if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
    {
        return;
    }

    LocusAnalyzerThreadLocalData& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
    auto& locusAnalyzer(*locusAnalyzerThreadSharedData.locusAnalyzers[locusIndex]);

    boost::optional<HtsStreamingReadPairQueue::ReadPair> readPair;

    try
    {
        // Workers retain the original drain-to-empty behavior. The reader may
        // take a heavy pending locus while capacity is full; bound its help so
        // that workers finishing other loci can promptly receive further input.
        // readerThreadIndex is the compute-worker count, not the last storage
        // slot (which is a worker when EH_READER_HELPER_DISABLE is defined).
        const std::size_t limit = threadIndex == readerThreadIndex ? 32 : 0;
        const bool needsContinuation = locusAnalyzerThreadSharedData.readPairQueue.processBatch(
            locusIndex, limit, readPair, [&](HtsStreamingReadPairQueue::ReadPair& pair) {
                processAnalyzerBundleReadPair(locusAnalyzer, pair.regionType, pair.inputType, pair.read, pair.mate,
                    *locusAnalyzerThreadData.alignerSelectorPtr);
            });
        if (needsContinuation && !locusAnalyzerThreadSharedData.isWorkerThreadException.load())
        {
            // All locus access above has ended. Queue publication transfers its
            // sole ownership to one continuation, while its active flag/count
            // remain unchanged. Only the synchronous reader reaches this path,
            // so every continuation is submitted before pool.stop(true).
            pool.push(processLocusAnalyzerQueue, std::ref(locusAnalyzerThreadSharedData),
                std::ref(locusAnalyzerThreadLocalDataPool), locusIndex, std::ref(pool), readerThreadIndex);
        }
    }
    catch (const std::exception& e)
    {
        locusAnalyzerThreadData.threadExceptionPtr = std::current_exception();
        locusAnalyzerThreadSharedData.isWorkerThreadException.store(true);
        locusAnalyzerThreadSharedData.readPairQueue.cancel();

        std::ostringstream oss;
        oss << "Exception caught in thread " << threadIndex << " while processing read pair queue for locus: `"
            << locusAnalyzer.locusId() << "`";
        if (readPair)
        {
            oss << " current readPair: `" << readPair->read.fragmentId() << "`";
        }
        oss << ": " << e.what();
        spdlog::error(oss.str());
        throw;
    }
    catch (...)
    {
        locusAnalyzerThreadData.threadExceptionPtr = std::current_exception();
        locusAnalyzerThreadSharedData.isWorkerThreadException.store(true);
        locusAnalyzerThreadSharedData.readPairQueue.cancel();

        std::ostringstream oss;
        oss << "Exception caught in thread " << threadIndex << " while processing read pair queue for locus: `"
            << locusAnalyzer.locusId() << "`";
        if (readPair)
        {
            oss << " current readPair: `" << readPair->read.fragmentId() << "`";
        }
        spdlog::error(oss.str());
        throw;
    }
}

/// \brief Mutable data shared by all SampleFindings-processing threads
///
class SampleFindingsThreadSharedData
{
public:
    SampleFindingsThreadSharedData()
        : isWorkerThreadException(false)
        , locusIndex(0)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    std::atomic<unsigned> locusIndex;
};

/// \brief Data isolated to each SampleFindings-processing thread
///
struct SampleFindingsThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
};

/// \brief Analyze a series of loci on one thread
///
void analyzeLocus(
    const int threadIndex, const Sex sampleSex, vector<std::unique_ptr<LocusAnalyzer>>& locusAnalyzers,
    SampleFindings& sampleFindings, SampleFindingsThreadSharedData& sampleFindingsThreadSharedData,
    std::vector<SampleFindingsThreadLocalData>& sampleFindingsThreadLocalData)
{
    SampleFindingsThreadLocalData& sampleFindingsThreadData(sampleFindingsThreadLocalData[threadIndex]);
    std::string locusId = "Unknown";

    try
    {

        const unsigned size(locusAnalyzers.size());
        while (true)
        {
            if (sampleFindingsThreadSharedData.isWorkerThreadException.load())
            {
                return;
            }
            const auto locusIndex(sampleFindingsThreadSharedData.locusIndex.fetch_add(1));
            if (locusIndex >= size)
            {
                return;
            }

            auto& locusAnalyzer(*locusAnalyzers[locusIndex]);
            locusId = locusAnalyzer.locusId();
            sampleFindings[locusIndex] = locusAnalyzer.analyze(sampleSex, boost::none);
        }
    }
    catch (const std::exception& e)
    {
        sampleFindingsThreadSharedData.isWorkerThreadException = true;
        sampleFindingsThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Exception caught in thread {} while analyzing locus: {} : {}", threadIndex, locusId, e.what());
        throw;
    }
    catch (...)
    {
        sampleFindingsThreadSharedData.isWorkerThreadException = true;
        sampleFindingsThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Unknown exception caught in thread {} while analyzing locus: {}", threadIndex, locusId);
        throw;
    }
}

}

SampleFindings htsStreamingSampleAnalysis(
    const InputPaths& inputPaths, Sex sampleSex, const HeuristicParameters& heuristicParams, const int threadCount,
    const RegionCatalog& regionCatalog, locus::AlignWriterPtr bamletWriter, unsigned decodeThreadCount,
    unsigned readDecodeThreadCount)
{
    // Setup thread-specific data structures and thread pool
    const unsigned maxActiveLocusAnalyzerQueues(threadCount + 5);
    const unsigned locusAnalyzerCount(regionCatalog.size());
    LocusAnalyzerThreadSharedData locusAnalyzerThreadSharedData(maxActiveLocusAnalyzerQueues, locusAnalyzerCount);
#ifndef EH_READER_HELPER_DISABLE
    const int processingSlotCount = threadCount + 1;
#else
    const int processingSlotCount = threadCount;
#endif
    // Fix the storage before workers start. The reader's additional slot is
    // never shared with a worker and its aligner is allocated only if needed.
    std::vector<LocusAnalyzerThreadLocalData> locusAnalyzerThreadLocalDataPool(processingSlotCount);
    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        auto& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
        locusAnalyzerThreadData.alignerSelectorPtr.reset(
            new graphtools::AlignerSelector(heuristicParams.alignerType()));
    }
    ctpl::thread_pool pool(threadCount);
    std::size_t assistedLocusJobs = 0; // Accessed only by the input reader.
    HtsStreamingReadPairQueue::WaitHelper waitHelper;
#ifndef EH_READER_HELPER_DISABLE
    waitHelper = [&]() {
        if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
            return false;
        // pop() removes a single pending task atomically, using the same FIFO
        // queue as the workers. A locus still has exactly one task owner.
        auto pending = pool.pop();
        if (!pending)
            return false;
        auto& readerData = locusAnalyzerThreadLocalDataPool[threadCount];
        if (!readerData.alignerSelectorPtr)
            readerData.alignerSelectorPtr.reset(new graphtools::AlignerSelector(heuristicParams.alignerType()));
        pending(threadCount);
        ++assistedLocusJobs;
        // The packaged task captures exceptions; processLocusAnalyzerQueue
        // records the exception and cancels capacity waits before rethrowing.
        return true;
    };
#endif

    spdlog::info("Initializing all loci");
    const auto indexStart = std::chrono::steady_clock::now();
    locusAnalyzerThreadSharedData.locusAnalyzers
        = initializeLocusAnalyzers(regionCatalog, heuristicParams, bamletWriter, threadCount);
    GenomeQueryCollection genomeQuery(locusAnalyzerThreadSharedData.locusAnalyzers, &regionCatalog);
    spdlog::info(
        "Phase catalog_index_seconds={}",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - indexStart).count());

    spdlog::info("Streaming reads");
    const auto readsStart = std::chrono::steady_clock::now();

    using ReadCatalog = absl::flat_hash_set<Read, FragmentHash, FragmentEqual>;
    ReadCatalog unpairedReads(1000);

    htshelpers::HtsFileStreamer readStreamer(inputPaths.htsFile(), inputPaths.reference(), decodeThreadCount);
    {
        DecodedReadStream decodedStream(readStreamer, genomeQuery.targetRegionMask, readDecodeThreadCount);
        for (;;)
        {
            const auto* raw = readDecodeThreadCount == 0 ? decodedStream.nextUndecoded() : nullptr;
            auto* current = readDecodeThreadCount != 0 ? decodedStream.next() : nullptr;
            if (!raw && !current)
                break;
            // Stop processing reads if an exception is thrown in the worker pool:
            if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
            {
                break;
            }

            boost::optional<Read> decodedRead;
            if (current)
                decodedRead.emplace(std::move(current->read));
            const auto fragment = raw ? absl::string_view(bam_get_qname(raw))
                                      : absl::string_view(decodedRead->fragmentId());
            const auto mateIterator = unpairedReads.find(fragment);
            if (mateIterator == unpairedReads.end())
            {
                // The first record must remain available even when its own
                // coordinates look irrelevant: the eventual second record
                // supplies the routing coordinates in the original algorithm.
                if (raw)
                    decodedRead.emplace(htshelpers::decodeReadParts(bam_get_qname(raw), raw->core.flag,
                        raw->core.l_qseq, bam_get_seq(raw), bam_get_qual(raw)));
                unpairedReads.emplace(std::move(*decodedRead));
                continue;
            }

            const int32_t contig = raw ? raw->core.tid : current->contig;
            const int32_t mateContig = raw ? raw->core.mtid : current->mateContig;
            const int64_t position = raw ? raw->core.pos : current->position;
            const int64_t matePosition = raw ? raw->core.mpos : current->matePosition;
            const int64_t readEnd = position + (raw ? raw->core.l_qseq : decodedRead->sequence().length());
            const int64_t mateEnd = matePosition + mateIterator->sequence().length();

            vector<AnalyzerBundle> analyzerBundles = genomeQuery.analyzerFinder.queryStreaming(
                contig, position, readEnd, mateContig, matePosition, mateEnd);
            if (analyzerBundles.empty())
            {
                // Consume both records exactly as before; never leave a marker
                // that could alter subsequent repeated-name pairing or grow
                // the existing unpaired cache.
                unpairedReads.erase(mateIterator);
                continue;
            }
            if (raw)
                decodedRead.emplace(htshelpers::decodeReadParts(bam_get_qname(raw), raw->core.flag,
                    raw->core.l_qseq, bam_get_seq(raw), bam_get_qual(raw)));
            Read& read = *decodedRead;
            auto mateNode = unpairedReads.extract(mateIterator);
            Read mate = std::move(mateNode.value());

            const unsigned bundleCount(analyzerBundles.size());
            for (unsigned bundleIndex(0); bundleIndex < bundleCount; ++bundleIndex)
            {
                auto& bundle(analyzerBundles[bundleIndex]);
                auto sendReadPair = [&](HtsStreamingReadPairQueue::ReadPair readPair)
                {
                    if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
                        return;
                    const auto insertion = locusAnalyzerThreadSharedData.readPairQueue.insertReadPair(
                        bundle.locusIndex, std::move(readPair), waitHelper);
                    if (insertion == HtsStreamingReadPairQueue::InsertResult::Cancelled)
                        return;
                    if (insertion == HtsStreamingReadPairQueue::InsertResult::Activated)
                    {
                        pool.push(
                            processLocusAnalyzerQueue, std::ref(locusAnalyzerThreadSharedData),
                            std::ref(locusAnalyzerThreadLocalDataPool), bundle.locusIndex,
                            std::ref(pool), threadCount);
                    }
                };

                if ((bundleIndex + 1) < bundleCount)
                {
                    sendReadPair({ bundle.regionType, bundle.inputType, read, mate });
                }
                else
                {
                    sendReadPair({ bundle.regionType, bundle.inputType, std::move(read), std::move(mate) });
                }
            }
        }

        spdlog::info("Read preprocessing workers={} peak_worker_cache_bytes={} synchronous_large_records={}",
            readDecodeThreadCount, decodedStream.peakWorkerCacheBytes(), decodedStream.synchronousLargeRecords());
    }
    pool.stop(true);
    spdlog::info(
        "Phase read_processing_seconds={}",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - readsStart).count());
    spdlog::info("Streaming reader assisted_locus_jobs={}", assistedLocusJobs);

    // Rethrow exceptions from the pool in thread order:
    if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < processingSlotCount; ++threadIndex)
        {
            const LocusAnalyzerThreadLocalData& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
            if (locusAnalyzerThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(locusAnalyzerThreadData.threadExceptionPtr);
            }
        }
    }

    spdlog::info("Analyzing read evidence");
    const auto genotypeStart = std::chrono::steady_clock::now();

    SampleFindingsThreadSharedData sampleFindingsThreadSharedData;
    std::vector<SampleFindingsThreadLocalData> sampleFindingsThreadLocalDataPool(threadCount);

    const unsigned locusCount(locusAnalyzerThreadSharedData.locusAnalyzers.size());
    SampleFindings sampleFindings(locusCount);

    // Start all sampleFindings worker threads
    std::vector<std::thread> sampleFindingsThreads;
    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        sampleFindingsThreads.emplace_back(
            analyzeLocus, threadIndex, sampleSex, std::ref(locusAnalyzerThreadSharedData.locusAnalyzers),
            std::ref(sampleFindings), std::ref(sampleFindingsThreadSharedData),
            std::ref(sampleFindingsThreadLocalDataPool));
    }

    // Rethrow exceptions from worker pool in thread order:
    if (sampleFindingsThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
        {
            const auto& sampleFindingsThreadData(sampleFindingsThreadLocalDataPool[threadIndex]);
            if (sampleFindingsThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(sampleFindingsThreadData.threadExceptionPtr);
            }
        }
    }

    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        sampleFindingsThreads[threadIndex].join();
    }

    spdlog::info(
        "Phase genotyping_seconds={}",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - genotypeStart).count());
    return sampleFindings;
}

}
