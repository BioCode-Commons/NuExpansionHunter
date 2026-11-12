//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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
//

#include "core/CountTable.hh"

#include "gtest/gtest.h"

using std::map;
using std::vector;

using namespace ehunter;

TEST(InitializationOfCountTable, TypicalCountTable_Initialized)
{
    const map<int32_t, int32_t> elementsAndCounts = { { 1, 2 }, { 3, 5 } };
    CountTable countTable(elementsAndCounts);
    EXPECT_EQ(2, countTable.countOf(1));
    EXPECT_EQ(0, countTable.countOf(2));
    EXPECT_EQ(5, countTable.countOf(3));
}

TEST(ManipulatingCountTable, TypicalOperations_TableUpdated)
{
    CountTable countTable;
    countTable.incrementCountOf(4);
    EXPECT_EQ(1, countTable.countOf(4));

    countTable.setCountOf(4, 3);
    EXPECT_EQ(3, countTable.countOf(4));
}

TEST(ObtainingElementsWithNonzeroCounts, TypicalCountTable_ElementsObtained)
{
    const map<int32_t, int32_t> elementsAndCounts = { { 1, 2 }, { 3, 5 }, { 7, 15 } };
    CountTable countTable(elementsAndCounts);

    countTable.setCountOf(3, 0);

    vector<int32_t> expectedElements = { 1, 7 };
    EXPECT_EQ(expectedElements, countTable.getElementsWithNonzeroCounts());
}

TEST(TruncatingCounts, TypicalCountTables_CountsTruncated)
{
    const map<int32_t, int32_t> elementsAndCounts = { { 1, 2 }, { 3, 5 }, { 7, 15 }, { 10, 2 } };
    CountTable countTable(elementsAndCounts);

    {
        const map<int32_t, int32_t> expectedElementsAndCounts = { { 1, 2 }, { 3, 5 }, { 5, 17 } };
        const CountTable expectedCountTable(expectedElementsAndCounts);

        EXPECT_EQ(expectedCountTable, collapseTopElements(countTable, 5));
    }

    {
        const map<int32_t, int32_t> expectedElementsAndCounts = { { 1, 2 }, { 3, 22 } };
        const CountTable expectedCountTable(expectedElementsAndCounts);

        EXPECT_EQ(expectedCountTable, collapseTopElements(countTable, 3));
    }
}

TEST(MovingCountTable, PreservesTargetAndAllowsSourceReuse)
{
    const CountTable expected({ { -2, 3 }, { 7, 9 } });
    CountTable source = expected;
    CountTable target(std::move(source));
    EXPECT_EQ(expected, target);
    source.clear();
    source.incrementCountOf(11);
    EXPECT_EQ(1, source.countOf(11));
    CountTable assigned({ { 100, 2 } });
    assigned = std::move(target);
    EXPECT_EQ(expected, assigned);
    target.clear();
    target.setCountOf(17, 4);
    EXPECT_EQ(4, target.countOf(17));
}

TEST(MovingCountTable, SelfMovePreservesExistingPublicBehavior)
{
    CountTable table({ { 1, 3 }, { 8, 7 } });
    const CountTable expected = table;
    CountTable& alias = table;
    table = std::move(alias);
    EXPECT_EQ(expected, table);
}
