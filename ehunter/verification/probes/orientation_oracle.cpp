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
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "alignment/OrientationPredictor.hh"
#include "graphalign/PinnedDagAligner.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphutils/SequenceOperations.hh"

using namespace graphtools;
using namespace ehunter;

static uint64_t digest = 14695981039346656037ULL;
static void mix(uint64_t x) { digest = (digest ^ x) * 1099511628211ULL; }

static std::string repeat(const std::string& motif, size_t n)
{
    std::string result;
    result.reserve(n);
    for (size_t i = 0; i < n; ++i) result += motif[i % motif.size()];
    return result;
}

static int referenceCount(const std::string& query, const KmerIndex& index)
{
    int count = 0;
    const size_t k = index.kmerLength();
    for (size_t pos = 0; pos + k <= query.size();)
    {
        std::string word = query.substr(pos, k);
        std::transform(word.begin(), word.end(), word.begin(), ::toupper);
        if (index.contains(word)) { ++count; pos += k; }
        else ++pos;
    }
    return count;
}

static OrientationPrediction referencePredict(const std::string& q, const KmerIndex& index, int threshold)
{
    const int f = referenceCount(q, index);
    const int r = referenceCount(reverseComplement(q), index);
    if (std::max(f, r) < threshold) return OrientationPrediction::kDoesNotAlign;
    return f >= r ? OrientationPrediction::kAlignsInOriginalOrientation
                  : OrientationPrediction::kAlignsInReverseComplementOrientation;
}

static std::vector<std::string> reads(const std::string& workload)
{
    std::vector<size_t> lengths = { 35, 75, 100, 150, 250, 300 };
    if (workload == "short") lengths = { 1, 9, 10, 15, 16, 17, 25 };
    if (workload == "long") lengths = { 500, 1000, 5000 };
    std::mt19937 rng(78191);
    std::vector<std::string> result;
    const std::string alphabet = "ACGTacgtRYSWKMBDHVNryswkmbdhvnXx?";
    for (size_t n : lengths)
        for (int i = 0; i < 32; ++i)
        {
            const int kind = workload == "forward" ? 0 : workload == "reverse" ? 1
                : workload == "invalid" ? 3 : i % 4;
            std::string q = repeat(kind == 1 ? "CTG" : "CAG", n);
            if (kind == 2)
            {
                for (size_t j = 0; j < n; j += 13) q[j] = "ACGTN"[rng() % 5];
                for (size_t j = 1; j < n; j += 5) q[j] = static_cast<char>(::tolower(q[j]));
            }
            if (kind == 3) for (char& c : q) c = alphabet[rng() % alphabet.size()];
            result.push_back(q);
        }
    return result;
}

static size_t checkOrientation()
{
    size_t cases = 0;
    const Graph graph = makeStrGraph("ACGTTGCA", "CAG", "GATCCGTA");
    auto queries = reads("short");
    const auto regular = reads("mixed");
    queries.insert(queries.end(), regular.begin(), regular.end());
    queries.push_back("");
    for (int byte = 0; byte < 128; ++byte)
    {
        std::string q = repeat("ctgCTG", 51);
        q[17] = static_cast<char>(byte);
        queries.push_back(q);
    }
    for (int k : { 1, 2, 3, 10, 16 })
    {
        const KmerIndex index(graph, k);
        for (int threshold : { 0, 1, 3, 50 })
        {
            OrientationPredictor predictor(&graph, k, threshold);
            for (const auto& q : queries)
            {
                const auto expected = referencePredict(q, index, threshold);
                const auto actual = predictor.predict(q);
                if (expected != actual) throw std::runtime_error("orientation mismatch k=" + std::to_string(k));
                mix(static_cast<uint64_t>(actual));
                ++cases;
            }
        }
    }
    return cases;
}


int main()
{
    try { const auto cases = checkOrientation(); std::cout << cases << " " << digest << "\n"; }
    catch (const std::exception& error) { std::cerr << error.what() << "\n"; return 1; }
}
