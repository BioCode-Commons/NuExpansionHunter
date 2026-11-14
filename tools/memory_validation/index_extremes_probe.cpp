// Deterministic boundary workloads for exact compact-index and fallback checks.
// No fixtures or public data downloads are required.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "graphalign/KmerIndex.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/KmerEncoding.hh"
#include "graphutils/SequenceOperations.hh"

using namespace graphtools;

namespace
{
std::string randomSequence(size_t length, uint32_t seed)
{
    std::string sequence(length, 'A');
    for (char& base : sequence)
    {
        seed ^= seed << 13;
        seed ^= seed >> 17;
        seed ^= seed << 5;
        base = "ACGT"[seed & 3];
    }
    return sequence;
}

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::string deBruijn8()
{
    std::string sequence;
    std::vector<int> digits(33, 0);
    std::function<void(int, int)> append = [&](int depth, int period) {
        if (depth > 8)
        {
            if (8 % period == 0)
                for (int i = 1; i <= period; ++i)
                    sequence += "ACGT"[digits[i]];
            return;
        }
        digits[depth] = digits[depth - period];
        append(depth + 1, period);
        for (int value = digits[depth - period] + 1; value < 4; ++value)
        {
            digits[depth] = value;
            append(depth + 1, depth);
        }
    };
    append(1, 1);
    require(sequence.size() == 65536, "Invalid de Bruijn fixture");
    sequence += sequence.substr(0, 7);
    return sequence;
}

struct Profile
{
    Graph graph;
    int k = 14;
    bool membership = false;
    bool expectPrefix = false;
};

const std::vector<std::string> names {
    "start2047", "start2048", "node3", "node4", "dense14", "dense10", "skew14", "skew10",
    "entries65535", "entries65536", "complex14", "complex10", "ambiguous14", "ambiguous10",
    "tail14", "tail10", "empty14", "maximum16"
};

Profile makeProfile(const std::string& name)
{
    Profile profile;
    profile.membership = name.size() > 2 && name.substr(name.size() - 2) == "10";
    profile.k = profile.membership ? 10 : 14;
    if (name == "start2047" || name == "start2048")
    {
        profile.graph = Graph(1);
        profile.graph.setNodeSeq(0, randomSequence(name == "start2047" ? 2061 : 2062, 17891));
        profile.expectPrefix = name == "start2047";
    }
    else if (name == "node3" || name == "node4")
    {
        const NodeId node = name == "node3" ? 3 : 4;
        profile.graph = Graph(node + 1);
        for (NodeId i = 0; i < node; ++i)
            profile.graph.setNodeSeq(i, "AC");
        profile.graph.setNodeSeq(node, randomSequence(1000, 39171));
        profile.expectPrefix = node == 3;
    }
    else if (name == "dense14" || name == "dense10")
        profile.graph = makeStrGraph(std::string(1000, 'A'), "A", std::string(1000, 'A'));
    else if (name == "skew14" || name == "skew10")
    {
        auto sequence = randomSequence(1800, 71673);
        for (size_t i = 0; i + 16 < sequence.size(); i += 16)
            sequence.replace(i, 6, "AAAAAA");
        profile.graph = Graph(1);
        profile.graph.setNodeSeq(0, sequence);
    }
    else if (name == "entries65535" || name == "entries65536")
    {
        profile.k = 8;
        const size_t entries = name == "entries65535" ? 65535 : 65536;
        profile.graph = Graph(1);
        profile.graph.setNodeSeq(0, deBruijn8().substr(0, entries + 7));
    }
    else if (name == "complex14" || name == "complex10")
    {
        profile.graph = Graph(4);
        auto left = randomSequence(1000, 1713);
        left[270] = 'N';
        profile.graph.setNodeSeq(0, left);
        profile.graph.setNodeSeq(1, "A");
        profile.graph.setNodeSeq(2, "GTC");
        profile.graph.setNodeSeq(3, randomSequence(1000, 44291));
        for (const auto& edge : std::vector<std::pair<NodeId, NodeId>> {
                 { 0, 1 }, { 0, 2 }, { 1, 1 }, { 1, 3 }, { 2, 2 }, { 2, 3 } })
            profile.graph.addEdge(edge.first, edge.second);
        profile.expectPrefix = true;
    }
    else if (name == "ambiguous14" || name == "ambiguous10")
    {
        auto left = randomSequence(1000, 11613), right = randomSequence(1000, 86171);
        left[440] = 'N';
        right[700] = 'Y';
        profile.graph = makeStrGraph(left, "CAR", right);
        profile.expectPrefix = true;
    }
    else if (name == "tail14" || name == "tail10")
    {
        auto left = randomSequence(1000, 77213), right = randomSequence(1000, 61163);
        // The last prefix bucket contains low suffix 1. Suffix 0 is absent,
        // although padding contains zero words immediately after the table.
        if (profile.membership)
            right.replace(850, 10, "TTTTGAAAAC");
        else
            right.replace(850, 14, "TTTTTAAAAAAAAC");
        profile.graph = makeStrGraph(left, "CAG", right);
        profile.expectPrefix = true;
    }
    else if (name == "empty14")
    {
        profile.graph = Graph(1);
        profile.graph.setNodeSeq(0, "ACG");
    }
    else if (name == "maximum16")
    {
        profile.k = 16;
        profile.graph = Graph(3);
        profile.graph.setNodeSeq(0, std::string(16, 'A'));
        profile.graph.setNodeSeq(1, std::string(16, 'T'));
        profile.graph.setNodeSeq(2, std::string(16, 'T'));
    }
    else
        throw std::invalid_argument("Unknown profile " + name);
    return profile;
}

std::map<std::string, std::vector<Path>> enumerate(const Graph& graph, int k)
{
    std::map<std::string, std::vector<Path>> expected;
    for (NodeId node = 0; node < graph.numNodes(); ++node)
        for (size_t position = 0; position < graph.nodeSeq(node).size(); ++position)
            for (const auto& path : extendPathEnd(Path(&graph, position, { node }, position), k))
            {
                std::vector<std::string> expansions;
                expandReferenceSequence(path.seq(), expansions);
                for (const auto& key : expansions)
                    expected[key].push_back(path);
            }
    return expected;
}

struct ArchiveFile
{
    std::string path;
    ArchiveFile()
    {
        char name[] = "/tmp/eh-index-extremes-XXXXXX";
        const int fd = mkstemp(name);
        require(fd >= 0, "Cannot reserve archive file");
        close(fd);
        path = name;
        unlink(path.c_str());
    }
    ~ArchiveFile() { unlink(path.c_str()); }
};

void audit(const Profile& profile, const KmerIndex& index)
{
    const auto expected = enumerate(profile.graph, profile.k);
    require(index.kmers().size() == expected.size(), "Key count differs from path enumeration");
    std::unique_ptr<KmerIndex> retained;
    {
        ArchiveFile file;
        KmerIndex::prepareArchive(file.path, 1, [&](std::ostream& out, uint64_t) { index.writeRecord(out); });
        auto archive = std::make_shared<KmerIndexArchive>(file.path);
        KmerIndex mapped(profile.graph, profile.k, profile.membership, archive, 0), copy(mapped);
        retained.reset(new KmerIndex(std::move(copy)));
    }
    require(index.encode() == retained->encode(), "Mapped diagnostic traversal differs");
    std::cout << "encode " << index.encode() << '\n';
    const TwoBitKmerEncoder encoder(profile.k);
    for (const auto& item : expected)
    {
        const auto key = encoder.encode(item.first);
        require(index.containsEncodedKmer(key) && retained->containsEncodedKmer(key), "Missing present key");
        std::cout << "key " << item.first;
        if (!profile.membership)
        {
            const auto observed = index.getPaths(item.first);
            require(observed == item.second && retained->getPaths(item.first) == item.second,
                    "Ordered paths or duplicate multiplicity differ from enumeration");
            const auto view = retained->findEncodedPaths(key);
            require(view.size() == item.second.size(), "Path view count differs");
            for (size_t i = 0; i < view.size(); ++i)
            {
                require(view[i] == item.second[i], "Path view content differs");
                std::cout << ' ' << observed[i].encode();
            }
            bool threw = false;
            try { view[view.size()]; } catch (const std::out_of_range&) { threw = true; }
            require(threw, "Out-of-range path access was accepted");
        }
        std::cout << '\n';
    }
    const uint32_t mask = profile.k == 16 ? 0xffffffffU : (uint32_t(1) << (2 * profile.k)) - 1;
    const unsigned lowBits = 2 * profile.k - (profile.membership ? 9 : 10);
    std::vector<uint32_t> queries { 0, mask, mask - 1 };
    for (uint32_t prefix = 0; prefix < (1U << (profile.membership ? 9 : 10)); ++prefix)
    {
        queries.push_back(prefix << lowBits);
        queries.push_back((prefix << lowBits) | ((1U << lowBits) - 1));
    }
    for (const auto& item : expected)
    {
        const auto key = encoder.encode(item.first);
        if (key >= (1U << lowBits))
            queries.push_back(key - (1U << lowBits)); // Same suffix in preceding bucket.
    }
    uint32_t random = 912367;
    for (int i = 0; i < 100000; ++i)
    {
        random ^= random << 13;
        random ^= random >> 17;
        random ^= random << 5;
        queries.push_back(random & mask);
    }
    uint64_t hits = 0;
    for (auto key : queries)
    {
        const bool present = expected.count(encoder.decode(key));
        require(index.containsEncodedKmer(key) == present && retained->containsEncodedKmer(key) == present,
                "Membership differs for boundary or random query");
        hits += present;
        if (!present && !profile.membership)
            require(retained->findEncodedPaths(key).size() == 0, "Missing key has paths");
    }
    std::cout << "queries " << queries.size() << ' ' << hits << '\n';
}
}

int main(int argc, char** argv)
{
    try
    {
        if (argc < 2 || std::string(argv[1]) == "list")
        {
            for (const auto& name : names)
                std::cout << name << '\n';
            return 0;
        }
        require(argc >= 3, "Usage: index_extremes_probe audit|bench PROFILE [hit|miss|paths ROUNDS]");
        const std::string mode = argv[1], name = argv[2];
        const auto profile = makeProfile(name);
        const KmerIndex index(profile.graph, profile.k, profile.membership);
        std::ostringstream recordStream;
        index.writeRecord(recordStream);
        const auto record = recordStream.str();
        uint64_t identityBytes;
        std::memcpy(&identityBytes, record.data() + 24, sizeof(identityBytes));
        uint32_t buckets;
        std::memcpy(&buckets, record.data() + 24 + 16 + identityBytes + 4, sizeof(buckets));
        if (mode == "audit")
        {
            std::cerr << "{\"profile\":\"" << name << "\",\"k\":" << profile.k
                      << ",\"membership\":" << profile.membership << ",\"entries\":" << index.kmers().size()
                      << ",\"record_bytes\":" << record.size() << ",\"prefix\":" << (buckets == 0)
                      << ",\"expected_prefix\":" << profile.expectPrefix << "}\n";
            audit(profile, index);
            return 0;
        }
        require(mode == "bench" && argc >= 5, "Missing benchmark operation or rounds");
        const std::string operation = argv[3];
        const int rounds = std::stoi(argv[4]);
        require(!(profile.membership && operation == "paths"), "Membership indexes have no paths");
        const TwoBitKmerEncoder encoder(profile.k);
        const auto set = index.kmers();
        std::vector<std::string> sorted(set.begin(), set.end());
        std::sort(sorted.begin(), sorted.end());
        std::vector<uint32_t> queries;
        for (const auto& key : sorted)
            queries.push_back(encoder.encode(key));
        std::mt19937 random(71391);
        std::shuffle(queries.begin(), queries.end(), random);
        if (queries.size() > 4096)
            queries.resize(4096);
        if (queries.empty())
            queries.resize(1024, 0);
        const uint32_t mask = profile.k == 16 ? 0xffffffffU : (uint32_t(1) << (2 * profile.k)) - 1;
        if (operation == "miss")
        {
            std::vector<uint32_t> missing;
            uint32_t state = 913571;
            for (size_t attempt = 0; attempt < 100000 && missing.size() < queries.size(); ++attempt)
            {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                const auto key = state & mask;
                if (!index.containsEncodedKmer(key))
                    missing.push_back(key);
            }
            // The 65,535-entry eight-base profile has only one absent key.
            if (missing.empty() && profile.k <= 10)
                for (uint32_t key = 0; key <= mask; ++key)
                    if (!index.containsEncodedKmer(key))
                    {
                        missing.push_back(key);
                        break;
                    }
            require(!missing.empty(), "This profile covers the complete key universe; no valid miss exists");
            for (size_t i = 0; i < queries.size(); ++i)
                queries[i] = missing[i % missing.size()];
        }
        uint64_t count = 0, digest = 14695981039346656037ULL;
        const auto begin = std::chrono::steady_clock::now();
        for (int round = 0; round < rounds; ++round)
            for (auto key : queries)
            {
                if (operation == "paths")
                {
                    const auto view = index.findEncodedPaths(key);
                    count += view.size();
                    if (view.size())
                    {
                        const auto path = view[0];
                        digest += path.startPosition() + path.endPosition() + path.numNodes();
                    }
                }
                else
                {
                    count += index.containsEncodedKmer(key);
                    digest += key;
                }
            }
        const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        std::cout << count << ' ' << digest << '\n';
        std::cerr << "{\"seconds\":" << seconds << ",\"record_bytes\":" << record.size()
                  << ",\"prefix\":" << (buckets == 0) << ",\"queries_per_round\":" << queries.size() << "}\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
