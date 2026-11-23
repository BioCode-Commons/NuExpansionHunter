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

#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>

#include "graphalign/PinnedDagAligner.hh"
#include "graphutils/ExtensionTraceCache.hh"

namespace graphtools
{
namespace extension_cache
{
struct TableTestAccess
{
    static size_t bucketIndex(const Table::Key& key, size_t length)
    {
        return Table::bucketIndex(key, length);
    }

    static std::mutex& mutex(Table& table, const Table::Key& key, size_t length)
    {
        return table.buckets_[bucketIndex(key, length)].mutex;
    }

    static Table::Key& lastPublishedKey(Table& table, const Table::Key& key, size_t length)
    {
        auto& bucket = table.buckets_[bucketIndex(key, length)];
        return bucket.donors[(bucket.next + bucket.donors.size() - 1) % bucket.donors.size()].key;
    }
};
}
}

namespace
{
using namespace graphtools;
using namespace graphtools::extension_cache;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void tableCases(Table& table)
{
    const auto token = nextPairToken();
    const Table::Key donor{token, 1, true, false, false};
    const Table::Key recipient{token, 2, true, false, false};
    const std::string query = "ACT", target = "AGTCACGAG", matchingPrefix = "AGTCACGAT";
    // Q=3, S=7 => C=8 and L=8. A tie cannot be excluded at t=7.
    constexpr int optimum = 7;
    std::vector<Cigar> trace(2), observed;
    trace[0].push_back({Cigar::NODE_START, 0});
    trace[0].push_back({Cigar::MISMATCH, 1});
    trace[0].push_back({Cigar::MATCH, 2});
    trace[0].push_back({Cigar::NODE_END, 0});
    trace[1].push_back({Cigar::NODE_START, 0});
    trace[1].push_back({Cigar::MATCH, 2});
    trace[1].push_back({Cigar::MISMATCH, 1});
    trace[1].push_back({Cigar::NODE_END, 0});
    int score = -111;
    const auto find = [&](const Table::Key& key, const std::string& q, const std::string& t) {
        return table.find(key, q, t, observed, score);
    };
    require(!find(recipient, query, matchingPrefix) && score == -111, "hit before publication");
    require(table.publish(donor, query, target, optimum, trace), "donor did not publish");
    require(!find(donor, query, matchingPrefix), "donor used at the same locus");
    const Table::Key otherPair{7038, 2, true, false, false};
    require(otherPair.pair != donor.pair
            && TableTestAccess::bucketIndex(otherPair, query.size()) == TableTestAccess::bucketIndex(donor, query.size()),
        "pair tokens must collide for this fixture");
    require(!find(otherPair, query, matchingPrefix), "pair identity ignored within a bucket");

    // At this query length, flipping one flag changes the hash for a fixed pair.
    // Change the stored flags in place to exercise each comparison in isolation.
    auto& stored = TableTestAccess::lastPublishedKey(table, donor, query.size());
    stored.mateSlot = !donor.mateSlot;
    require(!find(recipient, query, matchingPrefix), "physical mate slot ignored within a bucket");
    stored = donor;
    stored.reversedRead = !donor.reversedRead;
    require(!find(recipient, query, matchingPrefix), "read orientation ignored within a bucket");
    stored = donor;
    stored.reversedExtension = !donor.reversedExtension;
    require(!find(recipient, query, matchingPrefix), "extension direction ignored within a bucket");
    stored = donor;
    require(!find(recipient, "ATC", matchingPrefix), "query bytes ignored");
    require(!find(recipient, query, "AGTCACGCT"), "certified prefix byte ignored");
    require(!find(recipient, query, "AGTCACG"), "short recipient matched");
    require(!find(recipient, query, matchingPrefix + std::string(3277, 'A')), "wide recipient matched");
    require(find(recipient, query, matchingPrefix) && score == optimum && observed == trace,
        "matching prefix lost a score, CIGAR or operation order");
    observed[0].begin()->value_ = 22;
    require(find(recipient, query, "AGTCACGA") && observed == trace, "CIGAR storage leaked a recipient node ID");

    const auto wideToken = nextPairToken();
    table.publish({wideToken, 1, true, false, false}, query, target + std::string(3277, 'A'), optimum, trace);
    require(!find({wideToken, 2, true, false, false}, query, matchingPrefix), "wide donor published");

    std::unique_lock<std::mutex> held(TableTestAccess::mutex(table, donor, query.size()));
    bool hit = true, published = true;
    std::thread contending([&] {
        hit = find(recipient, query, matchingPrefix);
        published = table.publish(donor, query, target, optimum, observed);
    });
    contending.join();
    require(!hit && !published,
        "contending lookup or publish did not fall back");
}

void packedTraceCases(Table& table)
{
    struct Case
    {
        size_t firstLength, secondLength;
        Cigar::Operation lastOperation;
        bool fits;
    };
    for (const auto& test : {
             Case{64, 64, {Cigar::MATCH, 4095}, true},
             Case{64, 65, {Cigar::MATCH, 4095}, true},
             Case{254, 1, {Cigar::MATCH, 4095}, true},
             Case{255, 1, {Cigar::MATCH, 4095}, true},
             Case{255, 2, {Cigar::MATCH, 4095}, false},
             Case{256, 0, {Cigar::MATCH, 4095}, false},
             Case{64, 65, {Cigar::MATCH, 4096}, false},
             Case{64, 65, {static_cast<Cigar::OpCode>(15), 4095}, false},
         })
    {
        // Two distinct ordered CIGARs exercise the codec with every legal
        // opcode and values across the packed bit positions.
        std::vector<Cigar> trace(test.secondLength ? 2 : 1);
        for (size_t i = 0; i < trace.size(); ++i)
        {
            const size_t length = i ? test.secondLength : test.firstLength;
            for (size_t n = 1; n < length; ++n)
            {
                const auto code = static_cast<Cigar::OpCode>((n + i) % (Cigar::NODE_END + 1));
                trace[i].push_back({code, n * 16 + i});
            }
            trace[i].push_back(test.lastOperation);
        }
        const auto token = nextPairToken();
        const std::string query = "ACT", target = "AGTCACGAG";
        const std::vector<Cigar> sentinel(1);
        std::vector<Cigar> observed = sentinel;
        constexpr int optimum = 7;
        int score = -111;
        const bool published = table.publish({token, 1, true, false, false}, query, target, optimum, trace);
        require(published == test.fits, "publication disagrees with the packed trace bounds");
        const bool hit = table.find({token, 2, true, false, false}, query, target, observed, score);
        require(hit == test.fits, "packed trace capacity or field validation differs");
        if (test.fits)
        {
            require(score == optimum && observed == trace, "packed trace changed score or ordered CIGARs");
        }
        else
        {
            require(score == -111 && observed == sentinel,
                "rejected trace affected the cache or lookup result");
        }
    }
}

void wrapperCase(Table& table, bool suffix, const std::string& query, const std::string& target,
    size_t extensionLength)
{
    // The shifted homopolymer has several optimal deletion placements. The
    // donor and recipients use different original node IDs and seed lengths.
    const auto makeGraph = [suffix](NodeId node, const std::string& seed, std::string extension) {
        Graph graph(node + 1);
        if (suffix)
            std::reverse(extension.begin(), extension.end());
        graph.setNodeSeq(node, suffix ? extension + seed : seed + extension);
        return graph;
    };
    const auto makeSeed = [suffix](const Graph& graph, size_t seedLength) {
        const auto node = static_cast<NodeId>(graph.numNodes() - 1);
        const auto length = graph.nodeSeqLength(node);
        return suffix ? Path(&graph, length - seedLength, {node}, length) : Path(&graph, 0, {node}, seedLength);
    };
    const auto extend = [&](PinnedDagAligner& aligner, const Path& seed, int& score) {
        return suffix ? aligner.suffixAlign(seed, std::string(query.rbegin(), query.rend()), extensionLength, score)
                      : aligner.prefixAlign(seed, query, extensionLength, score);
    };
    Graph donorGraph = makeGraph(0, "GT", target + "ACTT");
    Graph receiverGraph = makeGraph(2, "AGC", target + "CCGGAATA");
    Graph differentGraph = makeGraph(3, "AGCT", "T" + target.substr(1) + "CCGGAATA");
    const auto donor = makeSeed(donorGraph, 2), receiver = makeSeed(receiverGraph, 3);
    const auto different = makeSeed(differentGraph, 4);
    int referenceScore = -1, score = -1, differentScore = -1;
    PinnedDagAligner referenceAligner(5, -4, -8, -2), aligner(5, -4, -8, -2);
    const auto expected = extend(referenceAligner, receiver, referenceScore);
    bool cacheHit = false;
    const auto expectedDifferent = extend(referenceAligner, different, differentScore);
    require(expected.size() > 1 && expected != expectedDifferent, "fixture does not distinguish ordered ties");
    const auto print = [](int score, const std::list<PathAndAlignment>& alignments) {
        std::cout << score << ' ' << alignments.size() << '\n';
        for (const auto& alignment : alignments)
            std::cout << alignment.first.encode() << ' ' << alignment.second.referenceStart() << ' '
                      << alignment.second.generateCigar() << '\n';
    };
    const auto token = nextPairToken();
    {
        Scope scope(table, token, 0, true, &cacheHit);
        const auto actual = extend(aligner, donor, score);
        require(!cacheHit, "donor reused an unpublished trace");
        require(actual.size() == expected.size() && score == referenceScore,
            "donor differs before caching");
        print(score, actual);
    }
    {
        Scope scope(table, token, 2, true, &cacheHit);
        const auto actual = extend(aligner, receiver, score);
        require(cacheHit == verify::kExtensionCache, "receiver cache use differs from the build policy");
        cacheHit = false;
        require(actual == expected && score == referenceScore,
            "cached wrapper output differs after recipient node remapping");
        print(score, actual);
    }
    {
        Scope scope(table, token, 3, true, &cacheHit);
        const auto actual = extend(aligner, different, score);
        require(!cacheHit, "changed target used a cached trace");
        require(actual == expectedDifferent && score == differentScore,
            "wrapper ignored changed target content");
        print(score, actual);
    }
    {
        Scope scope(table, token, 2, true, &cacheHit);
        PinnedDagAligner custom(6, -4, -8, -2);
        const auto actual = extend(custom, receiver, score);
        require(!cacheHit, "custom scores used a cached trace");
        require(!actual.empty(), "custom-score fallback produced no paths");
        print(score, actual);
    }
    if (query.size() < 16)
    {
        Graph branching(2);
        const auto length = donorGraph.nodeSeqLength(0);
        const NodeId seedNode = suffix ? 1 : 0;
        branching.setNodeSeq(seedNode, std::string(donorGraph.nodeSeq(0)));
        branching.setNodeSeq(1 - seedNode, "AAAACCCC");
        branching.addEdge(0, 1);
        const Path seed = suffix ? Path(&branching, length - 2, {seedNode}, length)
                                 : Path(&branching, 0, {seedNode}, 2);
        const auto multiNode = [&](PinnedDagAligner& pinned, int& result) {
            return suffix ? pinned.suffixAlign(seed, std::string(query.rbegin(), query.rend()), length + 3, result)
                          : pinned.prefixAlign(seed, query, length + 3, result);
        };
        const auto reference = multiNode(referenceAligner, referenceScore);
        Scope scope(table, token, 4, true, &cacheHit);
        const auto actual = multiNode(aligner, score);
        require(!cacheHit, "multi-node preparation used a cached trace");
        require(actual == reference && score == referenceScore, "multi-node fallback differs");
        print(score, actual);
    }
    require(currentScope() == nullptr, "wrapper did not restore its cache scope");
}
}

int main()
{
    try
    {
        Table table;
        tableCases(table);
        packedTraceCases(table);
        // The homopolymer creates ties, and isolated substitutions make their
        // combined raw CIGARs long enough to exercise the upper half of the packed capacity.
        std::string longQuery = "CAAGT", longTarget = "CAAAGT";
        std::mt19937 bases(3198);
        for (size_t n = 0; n < 250; ++n)
        {
            const auto base = bases() % 4;
            longQuery += "ACGT"[base];
            longTarget += "ACGT"[(base + (n % 10 == 4)) % 4];
        }
        longTarget += std::string(160, 'T');
        for (const bool suffix : {false, true})
        {
            std::cout << (suffix ? "suffix" : "prefix") << '\n';
            wrapperCase(table, suffix, "CAAAAAGT", "CAAAAAAGTGGCTGTACC", 14);
            wrapperCase(table, suffix, "CAAAAAGTCCGAATCTGGTCA", "CAAAAAAGTCCGAATCTGGTCAGGCTGTACC", 25);
            wrapperCase(table, suffix, longQuery, longTarget, longTarget.size() + 3);
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
