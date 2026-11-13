#include "alignment/OrientationPredictor.hh"

#include <algorithm>
#include <limits>
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

TEST(OrientationPrediction, BoundedScansPreserveLateHitsAndStrandTies)
{
    const auto graph = graphtools::makeStrGraph("AAAAAA", "CAG", "AAAAAA");
    const graphtools::KmerIndex index(graph, 3, true);
    for (int minimum : { -1, 0, 1, 2, 3, 4, 50, std::numeric_limits<int>::max() })
    {
        const OrientationPredictor predictor(&graph, 3, minimum);
        // Mix both strand-specific motifs with invalid windows. In particular,
        // forward may stop below threshold while reverse still passes, or both
        // strands may pass with equal counts and must retain the forward tie.
        for (int forwardRepeats = 0; forwardRepeats <= 5; ++forwardRepeats)
        {
            for (int reverseRepeats = 0; reverseRepeats <= 5; ++reverseRepeats)
            {
                for (size_t padding : { 0u, 1u, 2u, 3u, 17u, 150u })
                {
                    std::string query(padding, 'N');
                    for (int i = 0; i < forwardRepeats; ++i)
                        query += "cag";
                    query += 'N';
                    for (int i = 0; i < reverseRepeats; ++i)
                        query += "CTG";
                    query += std::string(padding, 'N');
                    EXPECT_EQ(exhaustivePrediction(query, index, minimum), predictor.predict(query))
                        << "minimum=" << minimum << " query=" << query;
                }
            }
        }
    }
}

TEST(OrientationPrediction, BoundedScansMatchExhaustiveRandomizedCatalogs)
{
    std::mt19937 random(921632);
    const std::string alphabet = "ACGTacgtNRYry";
    for (int graphNumber = 0; graphNumber < 5; ++graphNumber)
    {
        std::string left(40, 'A'), motif(1 + graphNumber, 'A'), right(40, 'A');
        for (auto sequence : { &left, &motif, &right })
            for (char& base : *sequence)
                base = "ACGT"[random() % 4];
        const auto graph = graphtools::makeStrGraph(left, motif, right);
        for (int k = 1; k <= 16; ++k)
        {
            const graphtools::KmerIndex index(graph, k, true);
            for (int minimum : { -1, 0, 1, 2, 3, 5, 20 })
            {
                const OrientationPredictor predictor(&graph, k, minimum);
                for (int trial = 0; trial < 30; ++trial)
                {
                    std::string query(random() % 200, 'A');
                    for (char& base : query)
                        base = alphabet[random() % alphabet.size()];
                    // Half the reads contain genuine matching substrings, in
                    // either orientation and potentially near the bound.
                    if (trial % 2 == 0)
                    {
                        std::string matching = left;
                        for (int repeat = 0; repeat < trial; ++repeat)
                            matching += motif;
                        matching += right;
                        if (trial % 4 == 0)
                            matching = graphtools::reverseComplement(matching);
                        query.insert(query.empty() ? 0 : random() % query.size(), matching);
                    }
                    EXPECT_EQ(exhaustivePrediction(query, index, minimum), predictor.predict(query))
                        << "graph=" << graphNumber << " k=" << k << " minimum=" << minimum << " query=" << query;
                }
            }
        }
    }
}
