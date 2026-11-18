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

#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include "genotyping/AlignMatrix.hh"
#include "genotyping/AlignMatrixFiltering.hh"
#include "reference/genotyping/AlignMatrix_v5.hh"
#include "reference/genotyping/AlignMatrixFiltering_v5.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphBuilders.hh"

namespace product = ehunter::strgt;
namespace reference = ehunter::strgt::reference;

// Contract [K-25]: each public operation preserves every abstract row and best cell.
template <typename A, typename B> bool equalCell(const A& a, const B& b)
{
    return int(a.type()) == int(b.type()) && a.numMotifs() == b.numMotifs()
        && a.score() == b.score() && a.numIndels() == b.numIndels();
}
void compare(const product::AlignMatrix& actual, const reference::AlignMatrix& expected, unsigned step)
{
    const auto fail = [step]() { throw std::runtime_error("K-25 difference at step " + std::to_string(step)); };
    if (actual.numReads() != expected.numReads() || actual.getMaxMotifCount() != expected.getMaxMotifCount()) fail();
    for (int row = 0; row < actual.numReads(); ++row)
    {
        const auto& cells = expected.matrix()[row];
        if (size_t(actual.rowSize(row)) != cells.size() || !equalCell(actual.getBestAlign(row), expected.getBestAlign(row))) fail();
        for (size_t col = 0; col < cells.size(); ++col)
            if (!equalCell(actual.rowData(row)[col], cells[col]) || !equalCell(actual.getAlign(row, col), cells[col])) fail();
        for (int col : {int(cells.size()), int(cells.size()) + 7})
            if (!equalCell(actual.getAlign(row, col), expected.getAlign(row, col))) fail();
    }
}

graphtools::GraphAlignment alignment(const graphtools::Graph& graph, std::mt19937& random)
{
    const unsigned kind = random() % 5, motifs = 1 + random() % 6;
    if (kind == 0) return graphtools::decodeGraphAlignment(0, "0[6M]", &graph);
    std::string cigar = kind == 1 || kind == 2 ? "0[6M]" : "";
    for (unsigned i = 0; i < motifs; ++i)
    {
        const auto op = random() % 12;
        cigar += op == 0 ? "1[1M1I2M]" : op == 1 ? "1[1M1D1M]" : op == 2 ? "1[1M1X1M]" : "1[3M]";
    }
    if (kind == 1 || kind == 3) cigar += "2[6M]";
    return graphtools::decodeGraphAlignment(0, cigar, &graph);
}

int main()
{
    try
    {
        std::mt19937 random(0x62451u);
        const auto graph = graphtools::makeStrGraph("ATTCGA", "CAG", "ATGTCG");
        unsigned step = 0;
        for (unsigned sequence = 0; sequence < 100; ++sequence)
        {
            product::AlignMatrix actual(1);
            reference::AlignMatrix expected(1);
            compare(actual, expected, step++);
            for (unsigned operation = 0; operation < 200; ++operation)
            {
                const unsigned choice = random() % 5;
                if (choice <= 1 && actual.numReads() < 64)
                {
                    const auto read = alignment(graph, random), mate = alignment(graph, random);
                    actual.add(read, mate); expected.add(read, mate);
                }
                else if (choice == 2 && actual.numReads() >= 2)
                {
                    const int row = 2 * (random() % (actual.numReads() / 2));
                    actual.remove(row + 1); expected.remove(row + 1);
                    compare(actual, expected, step++);
                    actual.remove(row); expected.remove(row);
                }
                else if (choice == 3)
                {
                    product::filter(actual); reference::filter(expected);
                }
                else if (actual.numReads() < 64)
                {
                    const int maxMotifs = 1 + random() % 6, pairs = random() % 4;
                    product::addIrrPairsIfPossibleExpansion(maxMotifs, actual, pairs);
                    reference::addIrrPairsIfPossibleExpansion(maxMotifs, expected, pairs);
                }
                compare(actual, expected, step++);
            }
        }
        std::cout << "K-25 states " << step << '\n';
    }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
