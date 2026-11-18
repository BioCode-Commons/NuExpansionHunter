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
#pragma once
#include "graphutils/MappedData.hh"

namespace ehunter
{
struct PreparedInterval
{
    uint64_t start, stop, locus, rank;
    uint32_t regionType, reserved;
};
struct PreparedIntervalNode
{
    uint64_t center;
    graphtools::Slice intervals, left, right;
};
struct PreparedContigIntervals
{
    int32_t contig;
    uint32_t reserved;
    graphtools::Slice tree, sorted, mask;
};
}
