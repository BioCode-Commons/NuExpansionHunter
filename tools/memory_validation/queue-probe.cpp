// Compile against each checkout's queue source and core/Read.cpp with -O3 -pthread.
// This isolates queue metadata and scheduling; it is not an end-to-end WGS benchmark.
#include "sample/HtsStreamingReadPairQueue.hh"
#include "core/ThreadPool.hh"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace ehunter;

static HtsStreamingReadPairQueue::ReadPair makePair(unsigned index)
{
    const auto name = std::to_string(index);
    return { locus::RegionType::kTarget, AnalyzerInputType::kBothReads,
        Read(ReadId(name, MateNumber::kFirstMate), std::string(150, 'A'), false),
        Read(ReadId(name, MateNumber::kSecondMate), std::string(150, 'T'), true) };
}

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "usage: queue-probe idle|serial|parallel LOCI PAIRS WORKERS\n";
        return 2;
    }
    const std::string mode = argv[1];
    const unsigned loci = std::stoul(argv[2]);
    const unsigned pairs = std::stoul(argv[3]);
    const unsigned workers = std::stoul(argv[4]);
    const auto start = std::chrono::steady_clock::now();
    HtsStreamingReadPairQueue queue(workers + 5, loci);
    const auto initialized = std::chrono::steady_clock::now();
    uint64_t checksum = 0;
    if (mode == "serial")
    {
        boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
        for (unsigned i = 0; i < pairs; ++i)
        {
            const unsigned locus = i % loci;
            if (queue.insertReadPair(locus, makePair(i)) != HtsStreamingReadPairQueue::InsertResult::Activated)
                return 3;
            queue.getNextReadPair(locus, pair);
            checksum += std::stoul(pair->read.fragmentId());
            queue.getNextReadPair(locus, pair);
            if (pair)
                return 4;
        }
    }
    else if (mode == "parallel")
    {
        ctpl::thread_pool pool(workers);
        std::vector<uint64_t> sums(workers, 0);
        for (unsigned i = 0; i < pairs; ++i)
        {
            const unsigned locus = i % loci;
            if (queue.insertReadPair(locus, makePair(i)) == HtsStreamingReadPairQueue::InsertResult::Activated)
                pool.push([&, locus](int worker) {
                    boost::optional<HtsStreamingReadPairQueue::ReadPair> pair;
                    queue.processBatch(locus, 0, pair, [&](HtsStreamingReadPairQueue::ReadPair& current) {
                        sums[worker] += std::stoul(current.read.fragmentId());
                    });
                });
        }
        pool.stop(true);
        for (auto sum : sums)
            checksum += sum;
    }
    else if (mode != "idle")
        return 2;
    if (mode != "idle" && checksum != uint64_t(pairs) * (pairs - 1) / 2)
        return 5;
    const auto end = std::chrono::steady_clock::now();
    std::cout << "{\"construction_seconds\":" << std::chrono::duration<double>(initialized - start).count()
              << ",\"processing_seconds\":" << std::chrono::duration<double>(end - initialized).count()
              << ",\"checksum\":" << checksum << "}\n";
}
