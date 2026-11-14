//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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

#include "sample/GenomeMask.hh"

#include "gtest/gtest.h"
#include <limits>

using namespace ehunter;

TEST(GenomeMaskTest, covered)
{
    const int binSize = 1 << 16;
    GenomeMask mask;
    mask.addRegion(1, 10 * binSize + 1000, 11 * binSize + 10);
    mask.addRegion(1, 10 * binSize, 10 * binSize + 100);

    ASSERT_TRUE(mask.query(1, 10 * binSize));
    ASSERT_TRUE(mask.query(1, 10 * binSize + 50));
    ASSERT_TRUE(mask.query(1, 11 * binSize + 10));
    ASSERT_FALSE(mask.query(1, 12 * binSize));
    ASSERT_FALSE(mask.query(1, 10 * binSize - 1));
}

TEST(GenomeMaskTest, notCovered)
{
    const int binSize = 1 << 16;
    GenomeMask mask;
    mask.addRegion(0, 10 * binSize, 10 * binSize + 100);
    mask.addRegion(2, 10 * binSize + 1000, 11 * binSize + 10);

    ASSERT_FALSE(mask.query(1, 10 * binSize));
    ASSERT_FALSE(mask.query(3, 10 * binSize));
    ASSERT_FALSE(mask.query(2, 12 * binSize));
    ASSERT_FALSE(mask.query(2, 100 * binSize));
    ASSERT_TRUE(mask.query(0, 10 * binSize + 50));
}

TEST(GenomeMaskTest, outOfBounds)
{
    const int binSize = 1 << 16;
    GenomeMask mask;
    mask.addRegion(0, 10 * binSize, 10 * binSize + 100);
    mask.addRegion(2, 10 * binSize + 1000, 11 * binSize + 10);

    ASSERT_FALSE(mask.query(100, 10));
    ASSERT_FALSE(mask.query(1, 0));
}

TEST(GenomeMaskTest, BinBoundariesAndNegativeCoordinatesRemainExact)
{
    GenomeMask mask;
    mask.addRegion(0, 1024, 1024);
    mask.addRegion(0, 65536, 65536);
    EXPECT_FALSE(mask.query(0, 1023));
    EXPECT_TRUE(mask.query(0, 1024));
    EXPECT_TRUE(mask.query(0, 2047));
    EXPECT_FALSE(mask.query(0, 2048));
    EXPECT_FALSE(mask.query(0, 65535));
    EXPECT_TRUE(mask.query(0, 65536));
    EXPECT_TRUE(mask.query(0, 66559));
    EXPECT_FALSE(mask.query(0, 66560));
    EXPECT_FALSE(mask.query(-1, 1024));
    EXPECT_FALSE(mask.query(0, -1));
    EXPECT_FALSE(mask.query(0, std::numeric_limits<int64_t>::min()));
    EXPECT_FALSE(mask.query(0, std::numeric_limits<int64_t>::max()));
}

TEST(GenomeMaskTest, PreparedBoundsValidationKeepsOriginalShortCircuitOrder)
{
    graphtools::BlobBuilder builder;
    const auto bits = builder.append(uint64_t(2));
    PreparedContigIntervals root {};
    root.mask = bits;
    const auto roots = builder.append(root);
    GenomeMask mask;
    mask.usePrepared({ builder.bytes.data(), builder.bytes.size() }, roots);
    EXPECT_FALSE(mask.query(0, 1023));
    EXPECT_TRUE(mask.query(0, 1024));
    EXPECT_TRUE(mask.query(0, 2047));
    EXPECT_FALSE(mask.query(0, 2048));
    EXPECT_FALSE(mask.query(0, 65536));

    // Prepared bytes are read on each query: a changed backing bit is visible.
    builder.replace(bits.offset, uint64_t(4));
    EXPECT_FALSE(mask.query(0, 1024));
    EXPECT_TRUE(mask.query(0, 2048));

    // Setting an invalid root still does not throw until an in-range query.
    EXPECT_NO_THROW(mask.usePrepared(
        { builder.bytes.data(), builder.bytes.size() }, { builder.bytes.size() + 8, 1 }));
    EXPECT_FALSE(mask.query(-1, 0));
    EXPECT_FALSE(mask.query(0, -1));
    EXPECT_FALSE(mask.query(1, 0));
    EXPECT_THROW(mask.query(0, 0), std::runtime_error);

    root.mask = { builder.bytes.size() + 8, 1 };
    builder.replace(roots.offset, root);
    mask.usePrepared({ builder.bytes.data(), builder.bytes.size() }, roots);
    // The word-count check precedes validation of this invalid mask slice.
    EXPECT_FALSE(mask.query(0, 65536));
    EXPECT_THROW(mask.query(0, 0), std::runtime_error);
}
