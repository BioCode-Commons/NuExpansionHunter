// Differential membership, validation and owner-lifetime coverage for hasEdge.
#include "graphcore/Graph.hh"
#include <atomic>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

using namespace graphtools;

static unsigned checks = 0;
static uint64_t digest = 14695981039346656037ULL;
static void require(bool value, const char* reason)
{
    if (!value) throw std::runtime_error(reason);
    ++checks;
}
static void record(bool value)
{
    digest = (digest ^ unsigned(value)) * 1099511628211ULL;
}
static void expectFailure(const char* label, const std::function<void()>& call, bool logic)
{
    std::string kind, message;
    try { call(); }
    catch (const std::logic_error& e) { kind = "logic"; message = e.what(); }
    catch (const std::runtime_error& e) { kind = "runtime"; message = e.what(); }
    require(kind == (logic ? "logic" : "runtime"), "Wrong or absent exception type");
    std::cout << label << ' ' << kind << ' ' << message << '\n';
}
static MappedGraph graphRecord(const BlobBuilder& blob, Slice record)
{
    MappedGraph graph;
    std::memcpy(&graph, blob.bytes.data() + record.offset, sizeof(graph));
    return graph;
}
static void replaceSuccessors(BlobBuilder& blob, Slice record, Slice successors)
{
    const auto graph = graphRecord(blob, record);
    MappedNode node;
    std::memcpy(&node, blob.bytes.data() + graph.nodes.offset, sizeof(node));
    node.successors = successors;
    blob.replace(graph.nodes.offset, node);
}
static Graph mapWithOwner(const BlobBuilder& blob, Slice record)
{
    auto bytes = std::make_shared<std::string>(blob.bytes);
    return Graph::mapped(MappedBytes(bytes->data(), bytes->size(), bytes), record);
}
static void allMembership(const Graph& graph, const std::vector<std::vector<NodeId>>& expected)
{
    for (unsigned repeat = 0; repeat < 41; ++repeat)
        for (NodeId from = 0; from < graph.numNodes(); ++from)
            for (NodeId to = 0; to < graph.numNodes(); ++to)
            {
                const auto& row = expected[from];
                const bool member = std::find(row.begin(), row.end(), to) != row.end();
                const bool actual = graph.hasEdge(from, to);
                require(actual == member, "Wrong hasEdge membership");
                record(actual);
            }
    expectFailure("invalid-source", [&] { graph.hasEdge(4, 0); }, true);
    expectFailure("invalid-sink", [&] { graph.hasEdge(0, 4); }, true);
    expectFailure("invalid-both", [&] { graph.hasEdge(4, 5); }, true);
}
int main()
{
    try
    {
        Graph owned(4);
        for (NodeId node = 0; node < 4; ++node) owned.setNodeSeq(node, "AC");
        owned.addEdge(0, 0); owned.addEdge(0, 2); owned.addEdge(1, 3); owned.addEdge(2, 3);
        const std::vector<std::vector<NodeId>> ordinary = {{0, 2}, {3}, {3}, {}};
        allMembership(owned, ordinary);
        BlobBuilder blob;
        const Slice recordSlice = owned.freeze(blob);
        const Graph mapped = mapWithOwner(blob, recordSlice);
        allMembership(mapped, ordinary);

        BlobBuilder unsorted = blob;
        const auto repeatedSlice = unsorted.append(std::vector<NodeId>{3, 1, 3, 0});
        replaceSuccessors(unsorted, recordSlice, repeatedSlice);
        const Graph repeated = mapWithOwner(unsorted, recordSlice);
        allMembership(repeated, {{3, 1, 3, 0}, {3}, {3}, {}});
        require(repeated.successors(0).count(3) == 2, "Repeated view element changed");
        require(std::vector<NodeId>(repeated.successors(0)) == std::vector<NodeId>({3, 1, 3, 0}), "View order changed");

        BlobBuilder empty = blob;
        const auto emptySlice = empty.append(std::vector<NodeId>{});
        replaceSuccessors(empty, recordSlice, emptySlice);
        allMembership(mapWithOwner(empty, recordSlice), {{}, {3}, {3}, {}});

        const std::vector<Slice> corrupt = {
            {uint64_t(blob.bytes.size()) + 8, 0},
            {0, std::numeric_limits<uint64_t>::max()},
            {1, 1},
            {1, 0},
            {0, uint64_t(blob.bytes.size() / sizeof(NodeId)) + 1}
        };
        for (const auto slice : corrupt)
        {
            BlobBuilder invalid = blob;
            replaceSuccessors(invalid, recordSlice, slice);
            const Graph graph = mapWithOwner(invalid, recordSlice);
            expectFailure("bad-successor-slice", [&] { graph.hasEdge(0, 1); }, false);
            expectFailure("invalid-sink-before-slice", [&] { graph.hasEdge(0, 4); }, true);
            expectFailure("invalid-source-before-slice", [&] { graph.hasEdge(4, 1); }, true);
            // Lazy slice validation is retained: an unrelated source is valid.
            require(graph.hasEdge(1, 3), "Unrelated valid source became invalid");
        }

        // The Graph alone owns these bytes while hasEdge executes. Existing
        // returned ranges must still retain ownership after the Graph dies.
        CatalogRange<NodeId> escaped;
        std::weak_ptr<std::string> lifetime;
        {
            auto owner = std::make_shared<std::string>(unsorted.bytes);
            lifetime = owner;
            const Graph graph = Graph::mapped(MappedBytes(owner->data(), owner->size(), owner), recordSlice);
            owner.reset();
            require(!lifetime.expired() && graph.hasEdge(0, 3), "Graph did not retain mapped bytes");
            escaped = graph.successors(0);
        }
        require(!lifetime.expired(), "Range lost mapped owner");
        require(escaped.count(3) == 2, "Escaped range contents changed");
        escaped = CatalogRange<NodeId>();
        require(lifetime.expired(), "Mapped owner leaked");

        // Concurrent immutable reads with a shared Graph/owner are semantically
        // identical; this is a correctness smoke check, not a timing benchmark.
        std::atomic<unsigned> incorrect{0};
        std::vector<std::thread> threads;
        for (unsigned worker = 0; worker < 4; ++worker)
            threads.emplace_back([&] {
                for (unsigned n = 0; n < 2000; ++n)
                    if (!repeated.hasEdge(0, 3) || repeated.hasEdge(0, 2) || !repeated.hasEdge(1, 3)) ++incorrect;
            });
        for (auto& thread : threads) thread.join();
        require(incorrect.load() == 0, "Concurrent membership mismatch");
        std::cout << "checks " << checks << " digest " << digest << '\n';
        return 0;
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
