//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#pragma once

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
}

#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"

namespace ehunter
{

namespace htshelpers
{

LinearAlignmentStats decodeAlignmentStats(bam1_t* htsAlignPtr);
bool isPrimaryAlignment(bam1_t* htsAlignPtr);
Read decodeRead(bam1_t* htsAlignPtr);
Read decodeReadParts(const char* qname, uint32_t samFlag, int32_t readLength,
    const uint8_t* packedSequence, const uint8_t* qualities);
#ifdef EH_READ_DECODE_TEST_INSTRUMENT_SIMD
namespace detail
{
uint64_t readDecodeSimdBlocksForTesting();
bool readDecodeSimdSupportedForTesting();
}
#endif
ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr);

} // namespace htshelpers

}
