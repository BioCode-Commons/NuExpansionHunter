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
#include "gtest/gtest.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <string>
#include <vector>

#include "core/HtsHelpers.hh"

namespace
{
TEST(HtsReadParts, MatchesScalarDecoderForAllPackedBasesAndQualityBoundaries)
{
    // Exercise all 16 BAM symbols, both nibbles, odd lengths, and the inclusive
    // quality cutoff (including absent BAM qualities, represented as 255).
    const std::vector<uint8_t> qualities = { 0, 19, 20, 21, 40, 255 };
    for (unsigned packed = 0; packed != 256; ++packed)
    {
        for (const auto firstQuality : qualities)
        {
            for (const auto secondQuality : qualities)
            {
                const uint8_t sequence[] = { static_cast<uint8_t>(packed), static_cast<uint8_t>(255 - packed) };
                const uint8_t quality[] = { firstQuality, secondQuality, firstQuality };
                for (int length = 1; length <= 3; ++length)
                {
                    std::string expected;
                    for (int base = 0; base < length; ++base)
                    {
                        const char value = seq_nt16_str[bam_seqi(sequence, base)];
                        expected += quality[base] <= 20 ? static_cast<char>(std::tolower(value)) : value;
                    }
                    const auto read = ehunter::htshelpers::decodeReadParts(
                        "packed-bases", BAM_FPAIRED | BAM_FREAD1 | BAM_FREVERSE, length, sequence, quality);
                    ASSERT_EQ(read.sequence(), expected);
                    ASSERT_EQ(read.fragmentId(), "packed-bases");
                    ASSERT_TRUE(read.isFirstMate());
                    ASSERT_TRUE(read.isReversed());
                }
            }
        }
    }
}

std::string scalarReadBases(int32_t length, const uint8_t* sequence, const uint8_t* qualities)
{
    std::string expected;
    for (int32_t index = 0; index < length; ++index)
    {
        const char base = seq_nt16_str[bam_seqi(sequence, index)];
        expected += qualities[index] <= 20 ? static_cast<char>(std::tolower(base)) : base;
    }
    return expected;
}

TEST(HtsReadParts, VectorChunksMatchAllPackedBytesAndAllUnsignedQualities)
{
    // Two vector chunks plus an odd scalar tail. Every packed byte occurs
    // against every unsigned quality, including all high-bit values.
    const int32_t length = 65;
    std::vector<uint8_t> sequence((length + 1) / 2);
    std::vector<uint8_t> qualities(length);
    for (unsigned packed = 0; packed < 256; ++packed)
    {
        std::fill(sequence.begin(), sequence.end(), static_cast<uint8_t>(packed));
        for (unsigned quality = 0; quality < 256; ++quality)
        {
            for (int32_t index = 0; index < length; ++index)
                qualities[index] = static_cast<uint8_t>(quality + 37 * index);
            const uint32_t flags = BAM_FPAIRED | (packed % 2 ? BAM_FREAD1 : BAM_FREAD2)
                | (quality % 2 ? BAM_FREVERSE : 0);
            const auto read = ehunter::htshelpers::decodeReadParts(
                "vector-bases", flags, length, sequence.data(), qualities.data());
            ASSERT_EQ(read.sequence(), scalarReadBases(length, sequence.data(), qualities.data()));
            ASSERT_EQ(read.fragmentId(), "vector-bases");
            ASSERT_EQ(read.isFirstMate(), bool(flags & BAM_FREAD1));
            ASSERT_EQ(read.isReversed(), bool(flags & BAM_FREVERSE));
        }
    }
}

TEST(HtsReadParts, VectorBoundariesAndMisalignedInputsMatchScalarOracle)
{
    for (const int32_t length : { 1, 2, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 149, 150, 151, 255, 256, 257 })
    {
        for (std::size_t offset = 0; offset < 32; ++offset)
        {
            // Allocate exactly through the last logical input byte. Sanitized
            // builds can therefore detect an accidental full-vector tail load.
            std::vector<uint8_t> sequence(offset + (length + 1) / 2);
            const std::size_t qualityOffset = 31 - offset;
            std::vector<uint8_t> qualities(qualityOffset + length);
            for (std::size_t i = offset; i < sequence.size(); ++i)
                sequence[i] = static_cast<uint8_t>(i * 73 + offset * 17);
            for (std::size_t i = qualityOffset; i < qualities.size(); ++i)
                qualities[i] = static_cast<uint8_t>(i * 31 + offset * 19);
            const auto* packed = sequence.data() + offset;
            const auto* quality = qualities.data() + qualityOffset;
            const auto read = ehunter::htshelpers::decodeReadParts(
                "unaligned-bases", BAM_FPAIRED | BAM_FREAD1, length, packed, quality);
            ASSERT_EQ(read.sequence(), scalarReadBases(length, packed, quality));
        }
    }
}
}
