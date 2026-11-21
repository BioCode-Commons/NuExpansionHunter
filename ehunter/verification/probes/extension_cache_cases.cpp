// Copyright 2026 NuExpansionHunter contributors
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>

#include "graphalign/PinnedDagAligner.hh"

namespace graphtools
{
namespace extension_cache
{
struct TableTestAccess
{
    static std::mutex& mutex(Table& table, const Table::Key& key, size_t length)
    {
        return table.buckets_[Table::bucketIndex(key, length)].mutex;
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
    Counters counts;
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
        return table.find(key, q, t, counts, observed, score);
    };
    require(!find(recipient, query, matchingPrefix) && score == -111, "hit before publication");
    table.publish(donor, query, target, optimum, trace, counts);
    require(!find(donor, query, matchingPrefix), "donor used at the same locus");
    require(!find({nextPairToken(), 2, true, false, false}, query, matchingPrefix), "pair identity ignored");
    require(!find({token, 2, false, false, false}, query, matchingPrefix), "physical mate slot ignored");
    require(!find({token, 2, true, true, false}, query, matchingPrefix), "read orientation ignored");
    require(!find({token, 2, true, false, true}, query, matchingPrefix), "extension direction ignored");
    require(!find(recipient, "ATC", matchingPrefix), "query bytes ignored");
    require(!find(recipient, query, "AGTCACGCT"), "certified prefix byte ignored");
    require(!find(recipient, query, "AGTCACG"), "short recipient matched");
    require(!find(recipient, query, matchingPrefix + std::string(3277, 'A')), "wide recipient matched");
    require(find(recipient, query, matchingPrefix) && score == optimum && observed == trace,
        "matching prefix lost a score, CIGAR or operation order");
    observed[0].begin()->value_ = 22;
    require(find(recipient, query, "AGTCACGA") && observed == trace, "CIGAR storage leaked a recipient node ID");

    const auto wideToken = nextPairToken();
    table.publish({wideToken, 1, true, false, false}, query, target + std::string(3277, 'A'), optimum, trace, counts);
    require(!find({wideToken, 2, true, false, false}, query, matchingPrefix), "wide donor published");

    std::unique_lock<std::mutex> held(TableTestAccess::mutex(table, donor, query.size()));
    bool hit = true;
    std::thread contending([&] {
        hit = find(recipient, query, matchingPrefix);
        table.publish(donor, query, target, optimum, observed, counts);
    });
    contending.join();
    require(!hit && counts.busyLookups == 1 && counts.busyPublishes == 1 && counts.published == 1,
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
        Counters counts;
        const auto token = nextPairToken();
        const std::string query = "ACT", target = "AGTCACGAG";
        const std::vector<Cigar> sentinel(1);
        std::vector<Cigar> observed = sentinel;
        constexpr int optimum = 7;
        int score = -111;
        table.publish({token, 1, true, false, false}, query, target, optimum, trace, counts);
        const bool hit = table.find({token, 2, true, false, false}, query, target, counts, observed, score);
        require(hit == test.fits, "packed trace capacity or field validation differs");
        if (test.fits)
        {
            require(counts.published == 1 && counts.hits == 1 && counts.traceTooLarge == 0,
                "representable trace was not published and read once");
            require(score == optimum && observed == trace, "packed trace changed score or ordered CIGARs");
            require(counts.expandedHits == (test.firstLength + test.secondLength > 128 ? 1u : 0u),
                "cache hit was counted in the wrong capacity range");
        }
        else
        {
            require(counts.published == 0 && counts.hits == 0 && counts.traceTooLarge == 1
                    && counts.expandedHits == 0 && score == -111 && observed == sentinel,
                "rejected trace affected the cache or lookup result");
        }
    }
}

void wrapperCase(Table& table, bool suffix, const std::string& query, const std::string& target,
    size_t extensionLength, bool expanded = false)
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
    const auto expectedDifferent = extend(referenceAligner, different, differentScore);
    require(expected.size() > 1 && expected != expectedDifferent, "fixture does not distinguish ordered ties");
    const auto token = nextPairToken();
    Counters counts;
    {
        Scope scope(table, counts, token, 0, true);
        require(extend(aligner, donor, score).size() == expected.size() && score == referenceScore,
            "donor differs before caching");
    }
    {
        Scope scope(table, counts, token, 2, true);
        require(extend(aligner, receiver, score) == expected && score == referenceScore,
            "cached wrapper output differs after recipient node remapping");
    }
    {
        Scope scope(table, counts, token, 3, true);
        require(extend(aligner, different, score) == expectedDifferent && score == differentScore,
            "wrapper ignored changed target content");
    }
    const auto guarded = counts.eligible;
    {
        Scope scope(table, counts, token, 2, true);
        PinnedDagAligner custom(6, -4, -8, -2);
        require(!extend(custom, receiver, score).empty(), "custom-score fallback produced no paths");
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
        Scope scope(table, counts, token, 4, true);
        require(multiNode(aligner, score) == reference && score == referenceScore, "multi-node fallback differs");
    }
    require(counts.published >= 1 && counts.hits == 1 && counts.eligible == guarded
            && counts.expandedHits == (expanded ? 1u : 0u)
            && counts.shadowChecks == (verify::kShadow ? 1u : 0u) && currentScope() == nullptr,
        "wrapper cache eligibility or hit count differs");
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
        // combined raw CIGARs large enough to require the expanded capacity.
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
            wrapperCase(table, suffix, "CAAAAAGT", "CAAAAAAGTGGCTGTACC", 14);
            wrapperCase(table, suffix, "CAAAAAGTCCGAATCTGGTCA", "CAAAAAAGTCCGAATCTGGTCAGGCTGTACC", 25);
            wrapperCase(table, suffix, longQuery, longTarget, longTarget.size() + 3, true);
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
