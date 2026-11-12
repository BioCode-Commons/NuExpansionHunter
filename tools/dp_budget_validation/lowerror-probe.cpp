// Exact 0/1/2/8 substitution fixtures; all bytes checked before timing.
#define main inheritedPinnedProbeMain
#include "pinned-probe.cpp"
#undef main

std::string lowErrorQuery(const Graph& graph, const std::string& shape, size_t length,
                          bool suffix, unsigned sample, unsigned errors)
{
    const auto walk = referenceWalk(graph, shape, length, sample);
    if (walk.size() < length || errors > length) throw std::invalid_argument("invalid low-error fixture");
    const auto exact = suffix ? walk.substr(walk.size() - length) : walk.substr(0, length);
    auto query = exact;
    for (unsigned index = 0; index < errors; ++index)
    {
        const size_t position = (sample * 7 + index * (length / errors)) % length;
        query[position] = query[position] == 'A' ? 'C' : 'A';
    }
    unsigned differences = 0;
    for (size_t i = 0; i < length; ++i) differences += query[i] != exact[i];
    if (differences != errors) throw std::runtime_error("fixture error count mismatch");
    return query;
}

void lowErrorBenchmark(const std::string& storage, const std::string& direction, const std::string& shape,
               int queryLength, int extension, int nodeLength, int errors, int repetitions, bool dump)
{
    if (queryLength <= 0 || extension <= 0 || nodeLength < 0 || repetitions <= 0)
        throw std::invalid_argument("positive workload sizes required");
    if (direction != "prefix" && direction != "suffix" && direction != "mixed")
        throw std::invalid_argument("unknown direction");
    if (storage != "owned" && storage != "mapped") throw std::invalid_argument("unknown storage");
    const Graph owned = makeGraph(shape, queryLength, extension, nodeLength);
    BlobBuilder blob;
    Graph mapped;
    if (storage == "mapped")
    {
        const auto record = owned.freeze(blob);
        mapped = Graph::mapped(MappedBytes(blob.bytes.data(), blob.bytes.size()), record);
    }
    const Graph& graph = storage == "owned" ? owned : mapped;
    const auto prefixSeed = makeSeed(graph, false, false), suffixSeed = makeSeed(graph, true, false);
    std::vector<std::string> prefixQueries, suffixQueries;
    for (unsigned sample = 0; sample < 64; ++sample)
    {
        prefixQueries.push_back(lowErrorQuery(graph, shape, queryLength, false, sample, errors));
        suffixQueries.push_back(lowErrorQuery(graph, shape, queryLength, true, sample, errors));
    }
    PinnedDagAligner aligner(5, -4, -8, -2);
    Digest output; output.dump = dump;
    uint64_t checksum = 0, prefixResults = 0, suffixResults = 0, operations = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int repetition = 0; repetition < repetitions; ++repetition)
    {
        const bool suffix = direction == "suffix" || (direction == "mixed" && repetition % 2);
        const unsigned sample = repetition % 64;
        int score = -773;
        const auto results = suffix ? aligner.suffixAlign(suffixSeed, suffixQueries[sample], extension, score)
                                    : aligner.prefixAlign(prefixSeed, prefixQueries[sample], extension, score);
        if (dump)
        {
            std::cout << "CASE lowerror " << shape << ' ' << queryLength << ' ' << errors << ' ' << suffix << ' ' << sample << '\n';
            output.add(score); output.add(results.size());
            for (const auto& result : results)
            {
                output.text(result.first.encode()); output.text(result.first.seq());
                output.add(result.first.startPosition()); output.add(result.first.endPosition());
                for (NodeId node : result.first.nodeIds()) output.add(node);
                inspectAlignment(result.second, output);
            }
        }
        // Require a real completed alignment and consume every path/operation.
        // This measures the full production alignment call, not only a getter.
        if (results.empty()) throw std::runtime_error("benchmark returned no alignment");
        (suffix ? suffixResults : prefixResults) += results.size();
        checksum = checksum * 33 + static_cast<uint64_t>(score);
        for (const auto& result : results)
        {
            const auto& path = result.first;
            const auto& alignment = result.second;
            checksum = checksum * 33 + path.startPosition(); checksum = checksum * 33 + path.endPosition();
            checksum = checksum * 33 + path.numNodes();
            for (NodeId node : path.nodeIds()) checksum = checksum * 33 + node;
            checksum = checksum * 33 + alignment.referenceStart();
            for (const auto& operation : alignment)
            {
                checksum = checksum * 33 + static_cast<int>(operation.type());
                checksum = checksum * 33 + operation.length(); ++operations;
            }
        }
    }
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (!operations) throw std::runtime_error("benchmark returned no operations");
    std::cout << "{\"mode\":\"pinned\",\"storage\":\"" << storage << "\",\"node_length\":" << nodeLength << ",\"direction\":\"" << direction << "\",\"shape\":\"" << shape
              << "\",\"query_length\":" << queryLength << ",\"extension\":" << extension
              << ",\"repetitions\":" << repetitions << ",\"prefix_results\":" << prefixResults
              << ",\"suffix_results\":" << suffixResults << ",\"operations\":" << operations
              << ",\"checksum\":" << checksum << "}\n";
    std::cerr << "{\"kernel_seconds\":" << elapsed << "}\n";
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 10 || (std::string(argv[1]) != "bench" && std::string(argv[1]) != "verify"))
            throw std::invalid_argument("expected bench|verify owned|mapped prefix|suffix|mixed linear|join|repeat QLEN EXTENSION NODE_LENGTH ERRORS REP");
        lowErrorBenchmark(argv[2], argv[3], argv[4], std::stoi(argv[5]), std::stoi(argv[6]),
                          std::stoi(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]), std::string(argv[1]) == "verify");
#ifdef GRAPHTOOLS_BUDGET_DP_INSTRUMENT
        const auto& s = graphalign::dagAligner::PinnedBudgetAlignMatrix::stats();
        std::cerr << "{\"budget_dp\":{\"calls\":" << s.calls << ",\"attempts\":" << s.attempts
                  << ",\"succeeded\":" << s.succeeded << ",\"fallbacks\":" << s.fallbacks
                  << ",\"dense_cells\":" << s.denseCells << ",\"band_cells\":" << s.bandCells
                  << ",\"success_dense_cells\":" << s.successDenseCells
                  << ",\"fallback_dense_cells\":" << s.fallbackDenseCells << "}}\n";
#endif
        return 0;
    }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
