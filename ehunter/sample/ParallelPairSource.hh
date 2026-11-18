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

#include <functional>
#include "sample/StreamingPairSource.hh"

namespace ehunter
{
// Admitted records per input batch, and buffered assignments per publisher flush.
constexpr size_t kStreamingBatchRecords = 16384;

// Input batches in flight between the reader and the publisher: each can be
// staged at the same time.
size_t streamingWindowSlots(size_t laneCount);

// Partition complete names across independent pending-mate tables. Publication
// follows the input ordinal of the arriving mate regardless of lane completion.
using PrepareStreamingBatch = std::function<std::function<void()>(const std::vector<StreamingInputPair*>&)>;
// Each batch is prepared on a staging thread; the returned commit runs on the calling thread in input order.
void streamParallelPairs(htshelpers::HtsFileStreamer& input, const GenomeMask& mask, const AnalyzerFinder& finder,
    size_t laneCount, const std::function<void()>& check, unsigned decodeWorkers, const PrepareStreamingBatch& prepare);
}
