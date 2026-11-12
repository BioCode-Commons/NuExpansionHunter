#include "alignment/OrientationPredictor.hh"

#include <algorithm>
#include <random>

#include "gtest/gtest.h"
#include "graphcore/GraphBuilders.hh"
#include "graphutils/SequenceOperations.hh"

using namespace ehunter;

namespace
{
int countMatches(const std::string& query, const graphtools::KmerIndex& index)
{
    int count = 0;
    const size_t k = index.kmerLength();
    for (size_t pos = 0; pos + k <= query.size();)
    {
        std::string kmer = query.substr(pos, k);
        std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
        if (index.contains(kmer))
        {
            ++count;
            pos += k;
        }
        else
            ++pos;
    }
    return count;
}

OrientationPrediction exhaustivePrediction(const std::string& query, const graphtools::KmerIndex& index, int minimum)
{
    const int forward = countMatches(query, index);
    const int reverse = countMatches(graphtools::reverseComplement(query), index);
    if (std::max(forward, reverse) < minimum)
        return OrientationPrediction::kDoesNotAlign;
    return forward >= reverse ? OrientationPrediction::kAlignsInOriginalOrientation
                              : OrientationPrediction::kAlignsInReverseComplementOrientation;
}
}

TEST(OrientationPrediction, UpperBoundPreservesThresholdsAndForwardTies)
{
    const auto graph = graphtools::makeStrGraph("ACGTAC", "CAG", "GTACGT");
    for (int k : { 1, 2, 3, 10, 16 })
    {
        const graphtools::KmerIndex index(graph, k, true);
        for (int minimum : { 0, 1, 3, 10 })
        {
            const OrientationPredictor predictor(&graph, k, minimum);
            for (const std::string query : { "", "A", "NNN", "ACGTACGT", "CAGCAGCAG", "cagCAGcagN", "CTGCTGCTG" })
                EXPECT_EQ(exhaustivePrediction(query, index, minimum), predictor.predict(query)) << query;
        }
    }
}

TEST(OrientationPrediction, MatchesFullTwoStrandScan)
{
    std::mt19937 random(71423);
    const auto graph = graphtools::makeStrGraph("ACGTAC", "CAG", "GTACGT");
    const std::string alphabet = "ACGTacgtN";
    for (int k : { 1, 2, 3, 10, 16 })
    {
        const graphtools::KmerIndex index(graph, k, true);
        for (int minimum : { 0, 1, 3, 10 })
        {
            const OrientationPredictor predictor(&graph, k, minimum);
            for (int trial = 0; trial < 200; ++trial)
            {
                std::string query(random() % 80, 'A');
                for (char& base : query)
                    base = alphabet[random() % alphabet.size()];
                EXPECT_EQ(exhaustivePrediction(query, index, minimum), predictor.predict(query)) << query;
            }
        }
    }
}

TEST(OrientationPrediction, ReverseAdapterPreservesAmbiguityCaseAndInvalidWindows)
{
    const auto graph = graphtools::makeStrGraph("ACGTTGCA", "CAG", "GATCCGTA");
    const std::string alphabet = "ACGTacgtRYSWKMBDHVNryswkmbdhvnXx? -";
    std::mt19937 random(8675309);
    for (int k : { 1, 2, 3, 10, 16 })
    {
        const graphtools::KmerIndex index(graph, k, true);
        for (int minimum : { 0, 1, 3, 50 })
        {
            const OrientationPredictor predictor(&graph, k, minimum);
            for (size_t length : { 0u, 1u, 9u, 10u, 15u, 16u, 17u, 75u, 150u, 300u })
            {
                for (int trial = 0; trial < 10; ++trial)
                {
                    std::string query(length, 'A');
                    for (size_t pos = 0; pos < length; ++pos)
                        query[pos] = trial < 4 ? "CTGctg"[pos % 6] : alphabet[random() % alphabet.size()];
                    if (length && trial < 4)
                        query[random() % length] = alphabet[random() % alphabet.size()];
                    EXPECT_EQ(exhaustivePrediction(query, index, minimum), predictor.predict(query));
                }
            }
            // Embedded control bytes and punctuation must invalidate exactly the
            // same windows as the materialized reverse-complement reference.
            for (int byte = 0; byte < 128; ++byte)
            {
                std::string query = "ctgCTGctgCTGctgCTGctg";
                query[8] = static_cast<char>(byte);
                EXPECT_EQ(exhaustivePrediction(query, index, minimum), predictor.predict(query));
            }
        }
    }
}
