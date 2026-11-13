//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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

#include "graphalign/DagAlignerAffine.hh"
#include "graphalign/dagAligner/BaseMatchingPenaltyMatrix.hh"

#include "gtest/gtest.h"

using std::string;

using namespace graphalign;
using namespace graphalign::dagAligner;

template <typename T> std::string toString(const T& obj)
{
    std::stringstream ss;
    ss << obj;
    return ss.str();
}

TEST(SimpleAlignment, Short_to_short)
{
    DagAligner<true> aligner({ 5, -4 }, 0, -8);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 151, 151 } }), std::vector<int>({ 0 }));

    const string query = "tgCccgcCCcCCCCcccC";

    const string reference = "TGCAGTCCCGCCCCGTCCC";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score secondBestScore = 0;
    std::vector<Cigar> cigars;
    const Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);

    EXPECT_EQ(std::size_t(3), cigars.size());
    EXPECT_EQ("0[3=3X3=1X4=1X1D3=]", toString(cigars.at(2)));
    EXPECT_EQ("0[3=3X3=1X4=1D1X3=]", toString(cigars.at(1)));
    EXPECT_EQ("0[3=3X3=1D4=2X3=]", toString(cigars.at(0)));
    EXPECT_EQ(37, bestScore);
}

TEST(SimpleAlignment, Long_to_Long)
{
    DagAligner<true> aligner({ 5, -4 }, 0, -8);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 151, 151 } }), std::vector<int>({ 0 }));

    const string query = "TCTCGCCCCGCCCCTCAGGCGGCCTCCCTGCtgtgCCCCGCCCCGGCCcCGCCCCgCCCCcCCCCCcCCaCgCCCCCCcCccCcCCCCgCCCC"
                         "CCCCctCcCCCCccctCCCCtccCCtgCccgcCCcCCCCcccC";
    const string reference = "CCGCCCCGCCCCCGTCTCGCCCCGCCCCTCAGGCGGCCTCCCTGCTGTGCCCCGCCCCGGCCTCGCCACGCCCCTACCTCACCACGCCC"
                             "CCCGCATCGCCACGCCCCCCGCATCGCCACGCCTCCCTTACCATGCAGTCCCGCCCCGTCCC";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score secondBestScore = 0;
    std::vector<Cigar> cigars;
    const Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(3), cigars.size());

    EXPECT_EQ(
        "0[14D48=1X4=1X6=2X2=1X1=1X11=1X1=2X1=1X2=1X8=1X1=1X2=1X2=1X1=1X6=1X1=1X2=1X3=3X3=1X4=1X1D3=]",
        toString(cigars.at(2)));
    EXPECT_EQ(
        "0[14D48=1X4=1X6=2X2=1X1=1X11=1X1=2X1=1X2=1X8=1X1=1X2=1X2=1X1=1X6=1X1=1X2=1X3=3X3=1X4=1D1X3=]",
        toString(cigars.at(1)));
    EXPECT_EQ(
        "0[14D48=1X4=1X6=2X2=1X1=1X11=1X1=2X1=1X2=1X8=1X1=1X2=1X2=1X1=1X6=1X1=1X2=1X3=3X3=1D4=2X3=]",
        toString(cigars.at(0)));
    EXPECT_EQ(344, bestScore);
}

TEST(SimpleAlignment, AAAC_to_AGC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 3, 3 } }), std::vector<int>({ 0 }));

    const string query = "AAAC";
    const string reference = "AGC";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score secondBestScore = 0;

    std::vector<Cigar> cigars;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(3), cigars.size());
    EXPECT_EQ("0[1S1=1X1=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=1I1X1=]", toString(cigars.at(1)));
    EXPECT_EQ("0[1=1X1I1=]", toString(cigars.at(2)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1S1=1X1=]", toString(cigar));
}

TEST(SimpleAlignment, ATGC_to_AGC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "ATGC";
    const string reference = "AGC";

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 3, 3 } }), std::vector<int>({ 0 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=1I2=]", toString(cigar));
}

TEST(SimpleAlignment, TAACTTTTGGG_to_TGCTTTTAA)
{
    // query:     TAACTTTTGGG
    //            |x |||||xx
    // reference: TG-CTTTTAA-

    const string query = "TAACTTTTGGG";
    const string reference = "TGCTTTTAA";

    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 9, 9 } }), std::vector<int>({ 0 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=1I1X5=3S]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=1X1I5=3S]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);

    EXPECT_EQ("0[1=1I1X5=3S]", toString(cigar));
}

TEST(SimpleAlignment, TCACGGAGA_to_TACGAGAG)
{
    // TCACGGAGA
    // | ||| |||
    // T-ACG-AGAG

    const string query = "TCACGGAGA";
    const string reference = "TACGAGAG";

    DagAligner<true> aligner({ 5, -4 }, 0, -8);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 8, 8 } }), std::vector<int>({ 0 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=1I2=1I4=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=1I3=1I3=]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=1I2=1I4=]", toString(cigar));
}

TEST(ForkAlignment, AAAC_to_AAC_fork_AAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAAC";
    const string reference = "AAC"
                             "AAC";

    EdgeMap edges(std::vector<std::pair<int, int>>({ { -1, 3 }, { 6, 6 } }), std::vector<int>({ 0, 1 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(6), cigars.size());
    EXPECT_EQ("0[1S3=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=1I2=]", toString(cigars.at(1)));
    EXPECT_EQ("0[2=1I1=]", toString(cigars.at(2)));
    EXPECT_EQ("1[1S3=]", toString(cigars.at(3)));
    EXPECT_EQ("1[1=1I2=]", toString(cigars.at(4)));
    EXPECT_EQ("1[2=1I1=]", toString(cigars.at(5)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1S3=]", toString(cigar));
}

TEST(ForkAlignment, AAAC_to_AGC_fork_AAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAAC";
    const string reference = "AGC"
                             "AAC";
    EdgeMap edges(std::vector<std::pair<int, int>>({ { -1, 3 }, { 6, 6 } }), std::vector<int>({ 0, 1 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(3), cigars.size());
    EXPECT_EQ("1[1S3=]", toString(cigars.at(0)));
    EXPECT_EQ("1[1=1I2=]", toString(cigars.at(1)));
    EXPECT_EQ("1[2=1I1=]", toString(cigars.at(2)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("1[1S3=]", toString(cigar));
}

TEST(Fork, AAC_to_AGC_AAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAC";
    const string reference = "AGC"
                             "AAC";

    EdgeMap edges(std::vector<std::pair<int, int>>({ { -1, 3 }, { 6, 6 } }), std::vector<int>({ 0, 1 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("1[3=]", toString(cigar));
}

/*
 *    A
 *   / \
 *  A   C
 *   \ /
 *    G
 */
TEST(ForkJoin1Base, AAC_to_AGC_AAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAC";
    const string reference = "A"
                             "A"
                             "G"
                             "C";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 }, { 0, 2 }, { 1, 3 }, { 2, 3 }, { 4, 4 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=]1[1=]3[1=]", toString(cigar));
}

/*
 *    AA
 *   /  \
 *  A    C
 *   \  /
 *    AG
 */
TEST(ForkJoin2Base, AAGC_to_AAGC_AAAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAGC";
    const string reference = "A"
                             "AA"
                             "AG"
                             "C";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 }, { 0, 3 }, { 2, 5 }, { 4, 5 }, { 6, 6 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=]2[2=]3[1=]", toString(cigar));
}

/*
 *    AA
 *   /  \
 *  A    C
 *   \  /
 *    AC
 */
TEST(ForkJoin2Base, AAGC_to_AAAC_AACC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAGC";
    const string reference = "A"
                             "AA"
                             "AC"
                             "C";
    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 }, { 0, 3 }, { 2, 5 }, { 4, 5 }, { 6, 6 } }),
        std::vector<int>({ 0, 1, 2, 3 }));
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=]1[1=1X]3[1=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=]2[1=1X]3[1=]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=]2[1=1X]3[1=]", toString(cigar));
}

/*
 *    AA
 *   /  \
 *  A    C
 *   \  /
 *    AC
 */
TEST(ForkJoin2Base, AAC_to_AAAC_AACC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAC";
    const string reference = "A"
                             "AA"
                             "AC"
                             "C";
    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 }, { 0, 3 }, { 2, 5 }, { 4, 5 }, { 6, 6 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=]2[2=]", toString(cigar));
}

/*
 *  AAAA
 *      \
 *       T
 *      /
 *  GACC
 */
TEST(JoinStartAtOffset, ACCT_to_AAAAT_GACCT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "ACCT";
    const string reference = "AAAA"
                             "GACC"
                             "T";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { -1, 4 }, { 3, 8 }, { 7, 8 }, { 9, 9 } }), std::vector<int>({ 0, 1, 2 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("1[1D3=]2[1=]", toString(cigar));
}

/*
 *  AA
 *    \
 *     T
 *    /
 *   C
 */
TEST(JoinStartAtOffset, AAT_to_AAT_CT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAT";
    const string reference = "AA"
                             "C"
                             "T";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { -1, 2 }, { 1, 3 }, { 2, 3 }, { 4, 4 } }), std::vector<int>({ 0, 1, 2 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    EXPECT_EQ(
        "Aligner(AffineAlignMatrix(AlignMatrix(\n"
        "[0\t-2\t-4\t-6]\n"
        "[-2\t1\t-1\t-3]\n"
        "[-4\t-1\t2\t0]\n"
        "[-2\t-1\t-3\t-5]\n"
        "[-4\t-3\t0\t3]\n"
        ")))",
        toString(aligner));

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[2=]2[1=]", toString(cigar));
}

/*
 *   C
 *    \
 *     T
 *    /
 *  AA
 */
TEST(JoinStartAtOffset, AAT_to_CT_AAT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAT";
    const string reference = "C"
                             "AA"
                             "T";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { -1, 1 }, { 0, 3 }, { 2, 3 }, { 4, 4 } }), std::vector<int>({ 0, 1, 2 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    EXPECT_EQ(
        "Aligner(AffineAlignMatrix(AlignMatrix(\n"
        "[0\t-2\t-4\t-6]\n"
        "[-2\t-1\t-3\t-5]\n"
        "[-2\t1\t-1\t-3]\n"
        "[-4\t-1\t2\t0]\n"
        "[-4\t-3\t0\t3]\n"
        ")))",
        toString(aligner));

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("1[2=]2[1=]", toString(cigar));
}

/*
 *           AA
 *          /  \
 *  TCGTGTAA    CCCCCCCCTTTTT
 *          \  /
 *           GC
 */
TEST(ForkJoinLong, AAGCCCCCCCCCTTTTT_to_TCGTGTAACCCCCCCCTTTTT_TCGTGTGCCCCCCCCCTTTTT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAGCCCCCCCCCTTTTT";
    const string reference = "TCGTGTAA"
                             "AA"
                             "GC"
                             "CCCCCCCCTTTTT";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 7, 8 }, { 7, 10 }, { 9, 12 }, { 11, 12 }, { 25, 25 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[6D2=]2[2=]3[13=]", toString(cigar));
}

/*
 *           AAA
 *          /   \
 *  TCGTGTAA     CCCCCCCCTTTTT
 *          \   /
 *           GC
 */
TEST(ForkJoinLong, AAGCCCCCCCCCTTTTT_to_TCGTGTAAAAACCCCCCCCTTTTT_TCGTGTAAGCCCCCCCCCTTTTT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAGCCCCCCCCCTTTTT";
    const string reference = "TCGTGTAA"
                             "AAA"
                             "GC"
                             "CCCCCCCCTTTTT";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 7, 8 }, { 7, 11 }, { 10, 13 }, { 12, 13 }, { 26, 26 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[6D2=]2[2=]3[13=]", toString(cigar));
}

/*
 *  A---   TT
 *      \ /
 *       T
 *      / \
 *  GACC   C
 */
class SimpleGraphTest : public testing::Test
{
protected:
    string reference;
    EdgeMap edges;
    DagAligner<true> aligner;
    Score bestScore = 0;
    Score secondBestScore = 0;

    SimpleGraphTest()
        : reference("A"
                    "GACC"
                    "T"
                    "TT"
                    "C")
        , edges(
              std::vector<std::pair<int, int>>({ { -1, 1 }, { 0, 5 }, { 4, 5 }, { 5, 6 }, { 5, 8 }, { 9, 9 } }),
              std::vector<int>({ 0, 1, 2, 3, 4 }))
        , aligner({ 1, -1 }, 0, -2)
    {
    }
};

TEST_F(SimpleGraphTest, OffEnd)
{
    string query = "ATCTG";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=]2[1=1I]3[1=1S]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=]2[1=]3[1X1=1S]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(0, bestScore);
    EXPECT_EQ("0[1=]2[1=1I]3[1=1S]", toString(cigar));
}
//
TEST_F(SimpleGraphTest, QueryAllN)
{
    string query = "N";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    bestScore = aligner.backtrackAllPaths<true>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=]", toString(cigars.at(0)));
    EXPECT_EQ("1[1=]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<true>(edges, bestScore, secondBestScore);
    EXPECT_EQ(1, bestScore);
    EXPECT_EQ("0[1=]", toString(cigar));

    cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(1, bestScore);
    EXPECT_EQ("0[1=]", toString(cigar));

    query = "NNNNNNNNNN";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    cigars.clear();
    bestScore = aligner.backtrackAllPaths<true>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(1), cigars.size());
    EXPECT_EQ("1[4=]2[1=]3[2=3S]", toString(cigars.at(0)));

    cigar = aligner.backtrackBestPath<true>(edges, bestScore, secondBestScore);
    EXPECT_EQ(7, bestScore);
    EXPECT_EQ("1[4=]2[1=]3[2=3S]", toString(cigar));

    cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(1, bestScore);
    EXPECT_EQ("1[3S4=]2[1=]3[2=]", toString(cigar));
}

TEST_F(SimpleGraphTest, QuerySomeN)
{
    string query = "GANCNC";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(6, bestScore);
    EXPECT_EQ("1[4=]2[1=]4[1=]", toString(cigar));
}

TEST_F(SimpleGraphTest, EmptyQuery)
{
    string query = "";
    EXPECT_ANY_THROW(aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges));
}

TEST_F(SimpleGraphTest, BadQualities)
{
    string query = "gACc";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(4, bestScore);
    EXPECT_EQ("1[4=]", toString(cigar));
}

/*
 *     _
 *    / \
 *  G-TCC-AAAAA
 */
TEST(RepeatExpansion, SimpleRepeat)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string reference = "G"
                             "TCC"
                             "TCC"
                             "TCC"
                             "AAAAA";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 },
                                           { 3, 4 },
                                           { 6, 7 },
                                           { 3, 10 },
                                           { 6, 10 },
                                           { 9, 10 },
                                           { reference.length(), reference.length() } }),
        std::vector<int>({ 0, 1, 2, 3, 4 }));

    string query = "TCCTCCAA";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);
    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(6, bestScore);
    EXPECT_EQ("0[1D]1[3=]2[3=]4[2=]", toString(cigar));

    query = "GTCTCCCCAA";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=]1[1=1D1=]2[3=]3[1D2=]4[2=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=]1[2=1D]2[3=]3[1D2=]4[2=]", toString(cigars.at(1)));

    cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(6, bestScore);
    EXPECT_EQ("0[1=]1[1=1D1=]2[3=]3[1D2=]4[2=]", toString(cigar));
}

/*
 *  |\     |\
 *  G--TCC-C
 */
TEST(RepeatExpansion, HomoPolymers)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "GGTCCGC";
    const string reference = "G"
                             "G"
                             "G"
                             "TCC"
                             "C"
                             "C";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 },
                                           { 1, 2 },
                                           { 0, 3 },
                                           { 1, 3 },
                                           { 2, 3 },
                                           { 5, 6 },
                                           { 6, 7 },
                                           { reference.length(), reference.length() } }),
        std::vector<int>({ 0, 1, 2, 3, 4, 5 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(5, bestScore);
    EXPECT_EQ("0[1=]1[1=]3[3=]4[1X]5[1=]", toString(cigar));
}

/*
 *              ---
 *             /   \
 * C-GCC-GACAAC-GAC-CTTCCTGAACT
 *   \ /        \ /
 *    -          -
 */
TEST(RepeatExpansion, two_repeats)
{
    DagAligner<true> aligner({ 5, -4 }, 0, -8);

    const string query = "CgCCGCCA";
    const string reference = "C"
                             "GCC"
                             "GCC"
                             "GCC"
                             "GCC"
                             "GCC"
                             "GCC"
                             "GACAAC"
                             "GAC"
                             "GAC"
                             "GAC"
                             "GAC"
                             "CTTCCTGAACT";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 },
                                           { 3, 4 },
                                           { 6, 7 },
                                           { 9, 10 },
                                           { 12, 13 },
                                           { 15, 16 },
                                           { 0, 19 },
                                           { 3, 19 },
                                           { 6, 19 },
                                           { 9, 19 },
                                           { 12, 19 },
                                           { 15, 19 },
                                           { 18, 19 },
                                           { 24, 25 },
                                           { 27, 28 },
                                           { 30, 31 },
                                           { 33, 34 },
                                           { 27, 37 },
                                           { 30, 37 },
                                           { 33, 37 },
                                           { 36, 37 },
                                           { 24, 37 },
                                           { reference.length(), reference.length() } }),
        std::vector<int>({ 6, 7, 8, 9, 10, 11, 12, 5, 1, 2, 3, 4, 0 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score secondBestScore = 0;
    std::vector<Cigar> cigars;
    const Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);

    EXPECT_EQ(32, bestScore);
    EXPECT_EQ("6[1=]7[3=]8[3=]5[1D1=]", toString(cigars.at(0)));
}

namespace
{
template <int pad>
typename PaddedAlignMatrix<pad>::const_iterator legacyNextBestAlign(
    const PaddedAlignMatrix<pad>& matrix, typename PaddedAlignMatrix<pad>::const_iterator start, Score& bestScore)
{
    auto it = start;
    if (matrix.end() == it)
    {
        return matrix.end();
    }
    if (!(std::distance(matrix.begin(), it) % matrix.paddedRowLen()))
    {
        ++it;
    }
    if (matrix.rowLen() == std::distance(matrix.begin(), it) % matrix.paddedRowLen())
    {
        it += matrix.paddedRowLen() - matrix.rowLen() + 1;
    }
    auto best = it;
    bestScore = *(it++);
    while (matrix.end() != it)
    {
        if (matrix.rowLen() == std::distance(matrix.begin(), it) % matrix.paddedRowLen())
        {
            it += matrix.paddedRowLen() - matrix.rowLen() + 1;
        }
        if (bestScore < *it)
        {
            bestScore = *it;
            best = it;
        }
        ++it;
    }
    return best;
}

template <int pad> void checkBestAlignmentScan()
{
    PaddedAlignMatrix<pad> matrix;
    for (int qLen : { 1, 2, 15, 16, 17, 31, 32, 33, 63, 64, 65 })
    {
        matrix.reset(qLen, 7);
        for (int t = 0; t < 7; ++t)
        {
            // Include tied maxima, boundary maxima and negative interior scores.
            for (int q = -1; q < qLen; ++q)
            {
                matrix.at(q, t) = q < 0 ? 3 : (q * 7 + t * 3) % 5 - 3;
            }
        }
        for (int t = 0; t < 7; ++t)
        {
            for (int q = -1; q < matrix.paddedRowLen() - 1; ++q)
            {
                const auto start = matrix.cellOneOne() + t * matrix.paddedRowLen() + q;
                if (start > matrix.end())
                {
                    continue;
                }
                Score expectedScore = 123;
                Score actualScore = 123;
                const auto expected = legacyNextBestAlign(matrix, start, expectedScore);
                const auto actual = matrix.nextBestAlign(start, actualScore);
                EXPECT_EQ(expected, actual) << "qLen=" << qLen << ", q=" << q << ", t=" << t;
                EXPECT_EQ(expectedScore, actualScore);
            }
        }
    }
}

template <bool penalizeMove> void checkSinglePredecessorRecurrence()
{
    const FixedPenaltyMatrix<> penalties(5, -4);
    AffineAlignMatrix<FixedPenaltyMatrix<>, penalizeMove> reference(penalties, -3, -1);
    AffineAlignMatrixVectorized<FixedPenaltyMatrix<>, penalizeMove> optimized(penalties, -3, -1);
    const string target = "ACGTNACGT";
    // Single predecessor rows before and after a join, including a second root.
    const EdgeMap edges(
        { 0, 0, 1, 1, 2, 2, 3, 4, 4 }, { 0, 1, 2, 3, 4, 6, 7, 8, 10, 11 },
        { -1, 0, -1, 2, 1, 3, 4, 4, 5, 6, 7 });
    for (int qLen : { 1, 15, 16, 17, 63, 64, 65, 129, 2, 32 })
    {
        string query;
        for (int q = 0; q < qLen; ++q)
        {
            query += "ANTGXC"[(q * 7 + qLen) % 6];
        }
        reference.init(query.begin(), query.end(), target.begin(), target.end(), edges);
        optimized.init(query.begin(), query.end(), target.begin(), target.end(), edges);
        // The public dump excludes padding and compares every V score.
        EXPECT_EQ(toString(reference), toString(optimized));
        for (int t = -1; t < int(target.size()); ++t)
        {
            for (int q = -1; q < qLen; ++q)
            {
                EXPECT_EQ(reference.isInsertion(q, t), optimized.isInsertion(q, t));
                if (t < 0)
                {
                    continue;
                }
                for (auto p = edges.prevNodesBegin(t); p != edges.prevNodesEnd(t); ++p)
                {
                    EXPECT_EQ(reference.isDeletion(q, t, *p), optimized.isDeletion(q, t, *p));
                    EXPECT_EQ(reference.isMatch(q, t, *p), optimized.isMatch(q, t, *p));
                    EXPECT_EQ(reference.isMismatch(q, t, *p), optimized.isMismatch(q, t, *p));
                }
            }
        }
    }
}
template <bool penalizeMove> void checkScoreLimitRecurrence()
{
    const string target = "ACGTNACGT";
    const std::vector<int> nodes = { 0, 0, 1, 1, 2, 2, 3, 4, 4 };
    // Duplicate each edge to force the general multiple-predecessor recurrence.
    // Max is idempotent, so this must retain every score and traceback predicate.
    // The disconnected row also exercises the zero-predecessor fallback.
    const std::vector<std::size_t> index = { 0, 1, 2, 3, 4, 6, 6, 7, 9, 10 };
    const std::vector<int> predecessors = { -1, 0, -1, 2, 1, 3, 4, 5, 6, 7 };
    std::vector<std::size_t> duplicateIndex;
    std::vector<int> duplicatePredecessors;
    for (auto i : index)
    {
        duplicateIndex.push_back(i * 2);
    }
    for (auto predecessor : predecessors)
    {
        duplicatePredecessors.push_back(predecessor);
        duplicatePredecessors.push_back(predecessor);
    }
    const EdgeMap edges(nodes, index, predecessors);
    const EdgeMap duplicateEdges(nodes, duplicateIndex, duplicatePredecessors);
    const Score settings[][4] = { { 5, -4, -3, -1 }, { 0, 0, 0, 0 }, { 2, -2, 2, 1 },
                                 { 32767, -32768, 32767, 32767 }, { -32768, 32767, -32768, -32768 },
                                 { 30000, -30000, 30000, -30000 }, { 1, -1, -32768, -1 }, { -1, 1, 32767, 1 } };
    for (const auto& scores : settings)
    {
        const FixedPenaltyMatrix<> penalties(scores[0], scores[1]);
        AffineAlignMatrixVectorized<FixedPenaltyMatrix<>, penalizeMove> optimized(penalties, scores[2], scores[3]);
        AffineAlignMatrixVectorized<FixedPenaltyMatrix<>, penalizeMove> reference(penalties, scores[2], scores[3]);
        for (int qLen : { 1, 15, 16, 17, 63, 64, 65, 129, 2, 32 })
        {
            string query;
            for (int q = 0; q < qLen; ++q)
            {
                query += "ANTGXC"[(q * 7 + qLen) % 6];
            }
            reference.init(query.begin(), query.end(), target.begin(), target.end(), duplicateEdges);
            optimized.init(query.begin(), query.end(), target.begin(), target.end(), edges);
            EXPECT_EQ(toString(reference), toString(optimized));
            for (int t = -1; t < int(target.size()); ++t)
            {
                for (int q = -1; q < qLen; ++q)
                {
                    EXPECT_EQ(reference.isInsertion(q, t), optimized.isInsertion(q, t));
                    if (t < 0)
                    {
                        continue;
                    }
                    for (auto p = edges.prevNodesBegin(t); p != edges.prevNodesEnd(t); ++p)
                    {
                        EXPECT_EQ(reference.isDeletion(q, t, *p), optimized.isDeletion(q, t, *p));
                        EXPECT_EQ(reference.isMatch(q, t, *p), optimized.isMatch(q, t, *p));
                        EXPECT_EQ(reference.isMismatch(q, t, *p), optimized.isMismatch(q, t, *p));
                    }
                }
            }
        }
    }
}
} // namespace

TEST(PaddedAlignmentMatrix, RowScanPreservesLegacyScoresAndTieOrder)
{
    checkBestAlignmentScan<1>();
    checkBestAlignmentScan<16>();
}

TEST(AffineAlignmentMatrix, SinglePredecessorRowsMatchScalarRecurrence)
{
    checkSinglePredecessorRecurrence<true>();
    checkSinglePredecessorRecurrence<false>();
}

TEST(AffineAlignmentMatrix, SinglePredecessorMatchesGeneralRecurrenceAtScoreLimits)
{
    checkScoreLimitRecurrence<true>();
    checkScoreLimitRecurrence<false>();
}

namespace
{
template <bool penalizeMove> void checkReusedTargetPenaltyProfiles()
{
    const BaseMatchingPenaltyMatrix penalties(5, -4);
    AffineAlignMatrix<BaseMatchingPenaltyMatrix, penalizeMove> reference(penalties, -3, -1);
    AffineAlignMatrixVectorized<BaseMatchingPenaltyMatrix, penalizeMove> optimized(penalties, -3, -1);
    const std::vector<string> targets = { "ACGTXBDHKMNRSVWY", "A", "CGTXBDHKMNRSVWY", "CGTXBDHKMNRSVWY",
                                         "NNNNN", "acgt", "Y", "ACGTXBDHKMNRSVWYACGT", "X", "ACACAC" };
    const int queryLengths[] = { 16, 48, 32, 48, 17, 16, 31, 129, 2, 7 };
    for (std::size_t iteration = 0; iteration < targets.size(); ++iteration)
    {
        const auto& target = targets[iteration];
        const int qLen = queryLengths[iteration];
        string query(qLen, 'A');
        for (int q = 0; q < qLen; ++q) query[q] = "ANTGXacgt"[(q + iteration) % 9];
        const EdgeMap edges({ { int(target.size()), int(target.size()) } }, { 0 });
        // The scalar recurrence computes penalties directly, without profiles.
        reference.init(query.begin(), query.end(), target.begin(), target.end(), edges);
        optimized.init(query.begin(), query.end(), target.begin(), target.end(), edges);
        EXPECT_EQ(toString(reference), toString(optimized));
        Score expectedScore = 0, actualScore = 0;
        const auto expected = reference.template nextBestAlign<false>(reference.alignBegin(), expectedScore);
        const auto actual = optimized.template nextBestAlign<false>(optimized.alignBegin(), actualScore);
        EXPECT_EQ(expectedScore, actualScore);
        EXPECT_EQ(reference.targetOffset(expected), optimized.targetOffset(actual));
        EXPECT_EQ(reference.queryOffset(expected), optimized.queryOffset(actual));
        for (int t = -1; t < int(target.size()); ++t)
        {
            for (int q = -1; q < qLen; ++q)
            {
                EXPECT_EQ(reference.isInsertion(q, t), optimized.isInsertion(q, t));
                if (t < 0) continue;
                for (auto p = edges.prevNodesBegin(t); p != edges.prevNodesEnd(t); ++p)
                {
                    EXPECT_EQ(reference.isDeletion(q, t, *p), optimized.isDeletion(q, t, *p));
                    EXPECT_EQ(reference.isMatch(q, t, *p), optimized.isMatch(q, t, *p));
                    EXPECT_EQ(reference.isMismatch(q, t, *p), optimized.isMismatch(q, t, *p));
                }
            }
        }
    }
}
} // namespace

TEST(AffineAlignmentMatrix, ReusedTargetAlphabetAndLengthsMatchScalarPenalties)
{
    checkReusedTargetPenaltyProfiles<true>();
    checkReusedTargetPenaltyProfiles<false>();
}

TEST(AffineAlignmentMatrix, PaddedQueryLifetimePreservesLogicalLengthAndEmptyInputState)
{
    AffineAlignMatrixVectorized<BaseMatchingPenaltyMatrix, true> matrix({ 6, -4 }, -8, -2);
    const string target(180, 'A');
    const EdgeMap edges({ { 180, 180 } }, { 0 });
    for (int queryLength : { 150, 17, 1, 149, 16, 150 })
    {
        const string query(queryLength, 'A');
        matrix.init(query.begin(), query.end(), target.begin(), target.end(), edges);
        for (int stage = 0; stage != 3; ++stage)
        {
            if (stage == 1)
                EXPECT_THROW(matrix.init(query.begin(), query.begin(), target.begin(), target.end(), edges),
                             std::logic_error);
            if (stage == 2)
                EXPECT_THROW(matrix.init(query.begin(), query.end(), target.begin(), target.begin(), edges),
                             std::logic_error);
            EXPECT_EQ(queryLength, matrix.queryLen());
            Score score = 0;
            const auto best = matrix.nextBestAlign<false>(matrix.alignBegin(), score);
            EXPECT_EQ(6 * queryLength, score);
            EXPECT_EQ(queryLength - 1, matrix.queryOffset(best));
            EXPECT_EQ(queryLength - 1, matrix.targetOffset(best));
            EXPECT_TRUE(matrix.isMatch(queryLength - 1, queryLength - 1, queryLength - 2));
        }
    }
}
