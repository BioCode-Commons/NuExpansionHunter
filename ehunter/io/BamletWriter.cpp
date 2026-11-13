//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "io/BamletWriter.hh"

#include <cstring>

#if defined(__aarch64__) && !defined(EH_DISABLE_BAM_ENCODING_SIMD)
#include <arm_neon.h>
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(EH_DISABLE_BAM_ENCODING_SIMD)
#include <immintrin.h>
#endif

// cppcheck-suppress missingInclude
#include "htslib/bgzf.h"
// cppcheck-suppress missingInclude
#include "htslib/kseq.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"

using graphtools::GraphAlignment;
using graphtools::GraphReferenceMapping;
using graphtools::ReferenceInterval;
using std::string;
using std::to_string;

namespace ehunter
{

// Establish the owned HTS compression pool before starting the existing writer
// thread. A failed setup then unwinds the file without a joinable std::thread.
static htsFile* openBamletFile(const string& path, int compressionThreads)
{
    std::unique_ptr<htsFile, decltype(&hts_close)> file(hts_open(path.c_str(), "wb"), hts_close);
    if (!file)
    {
        throw std::runtime_error("Failed to open BAMlet output " + path);
    }
    if (compressionThreads > 0 && hts_set_threads(file.get(), compressionThreads) != 0)
    {
        throw std::runtime_error("Failed to initialize BAMlet compression threads");
    }
    return file.release();
}

static void validateMapping(const ReferenceContigInfo& contigInfo, const LocusSpecification& locusSpec)
{
    for (const auto& nodeAndRegion : locusSpec.referenceProjectionOfNodes())
    {
        auto nodeId = nodeAndRegion.first;
        const auto& region = nodeAndRegion.second;
        const auto& regionEncoding = encode(contigInfo, region);
        ReferenceInterval referenceInterval = ReferenceInterval::parseRegion(regionEncoding);
        graphtools::NodeReferenceMapping mapping(locusSpec.regionGraph(), nodeId, referenceInterval);
        (void)mapping;
    }
}

BamletWriter::BamletWriter(
    const string& bamletPath, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
    int compressionThreads, std::size_t maxQueueBytes)
    : filePtr_(openBamletFile(bamletPath, compressionThreads), hts_close)
    , bamHeader_(bam_hdr_init(), bam_hdr_destroy)
    , contigInfo_(contigInfo)
    , catalog_(regionCatalog)
    , writeQueue_(maxQueueBytes)
{
    if (!regionCatalog.mappedBytes().data)
        for (const auto& locusSpec : regionCatalog)
        {
            validateMapping(contigInfo, locusSpec);
        }

    writeHeader();
    bamContigIds_.reserve(contigInfo_.numContigs());
    for (int index = 0; index != contigInfo_.numContigs(); ++index)
    {
        const auto& name = contigInfo_.getContigName(index);
        bamContigIds_.push_back(name.empty() ? -1 : bam_name2id(bamHeader_.get(), name.c_str()));
    }
    writeThread_ = std::thread(&BamletWriter::writeHtsAlignments, this);
}

BamletWriter::~BamletWriter()
{
    // Explicit finish() reports errors on the main path. Unwinding still owns cleanup.
    try
    {
        finish();
    }
    catch (...)
    {
    }
}

void BamletWriter::finish()
{
    if (writeThread_.joinable())
    {
        writeQueue_.push(nullptr, 0);
        writeThread_.join();
    }
    const int closeResult = filePtr_ ? hts_close(filePtr_.release()) : 0;
    if (writeError_)
        std::rethrow_exception(writeError_);
    if (closeResult != 0)
        throw std::runtime_error("Cannot close alignment output");
}

void BamletWriter::writeHeader()
{
    const string initHeader = "@HD\tVN:1.4\tSO:unknown\n";
    bamHeader_->l_text = strlen(initHeader.c_str());
    bamHeader_->text = strdup(initHeader.c_str());
    bamHeader_->n_targets = contigInfo_.numContigs();

    // All this memory gets freed by the header (bam_hdr_destroy) using free
    bamHeader_->target_len = (uint32_t*)calloc(contigInfo_.numContigs(), sizeof(uint32_t));
    bamHeader_->target_name = (char**)calloc(contigInfo_.numContigs(), sizeof(char*));

    for (int index = 0; index != contigInfo_.numContigs(); ++index)
    {
        const string& contigName = contigInfo_.getContigName(index);
        bamHeader_->target_name[index] = (char*)malloc(contigName.length() + 1);
        memcpy(bamHeader_->target_name[index], contigName.c_str(), contigName.length() + 1);
        bamHeader_->target_len[index] = contigInfo_.getContigSize(index);
    }

    if (bam_hdr_write(filePtr_->fp.bgzf, bamHeader_.get()) != 0)
    {
        throw std::logic_error("Failed to write header");
    }
}

// generateSubpathForEachNode used to allocate a list and a node vector for every
// node, solely to obtain its start position. Check the same per-node bounds
// without constructing those Paths. Invalid inputs retain the original Path
// constructor's exception behavior through the slow fallback below.
static bool nodeSubpathsAreValid(const graphtools::Path& path)
{
    try
    {
        for (std::size_t index = 0; index < path.numNodes(); ++index)
        {
            const int32_t start = path.getStartPositionOnNodeByIndex(index);
            const int32_t end = path.getEndPositionOnNodeByIndex(index);
            const auto length = path.graphRawPtr()->nodeSeqLength(path.getNodeIdByIndex(index));
            if (start < 0 || end < 0 || static_cast<unsigned>(start) > length
                || static_cast<unsigned>(end) > length || start > end)
                return false;
        }
        return true;
    }
    catch (const std::exception&)
    {
        return false;
    }
}

void BamletWriter::write(
    graphtools::TextView locusId, const string& fragmentName, const string& query, bool isFirstMate, bool isReversed,
    bool isMateReversed, const GraphAlignment& alignment)
{
    const auto locusIndex = catalog_.findLocusIndex(locusId);
    const auto& path = alignment.path();
    if (!nodeSubpathsAreValid(path))
        (void)graphtools::generateSubpathForEachNode(path);
    int32_t contigId = -1;
    int64_t position = -1;
    for (std::size_t index = 0; index < path.numNodes(); ++index)
    {
        const auto node = path.getNodeIdByIndex(index);
        if (node >= catalog_.nodeCount(locusIndex))
            throw std::logic_error("Invalid node " + std::to_string(node));
        const auto region = catalog_.referenceRegion(locusIndex, node);
        if (region)
        {
            const auto start = path.getStartPositionOnNodeByIndex(index);
            if (start >= region->length())
                throw std::logic_error("Cannot map position outside node");
            const auto& contigName = contigInfo_.getContigName(region->contigIndex());
            contigId = bamContigIds_[region->contigIndex()];
            if (contigId == -1 && !contigName.empty())
                throw std::logic_error("Unknown contig name " + contigName);
            position = region->start() + start;
            break;
        }
    }
    write(contigId, position, fragmentName, query, isFirstMate, isReversed, isMateReversed, alignment);
}

// TODO: Consider moving to graph tools
static string summarizeAlignment(const GraphAlignment& alignment)
{
    string encoding = alignment.path().graphRawPtr()->graphId();
    encoding += ',';
    encoding += to_string(alignment.path().startPosition());
    encoding += ',';
    // Operation::generateCigar constructs a stream just to append its type.
    // Keep the same order and decimal formatting while writing directly into
    // the XG string instead of allocating per-operation/per-node strings.
    for (int32_t index = 0; index != static_cast<int32_t>(alignment.size()); ++index)
    {
        encoding += to_string(static_cast<int32_t>(alignment.getNodeIdByIndex(index)));
        encoding += '[';
        for (const auto& operation : alignment[index].operations())
        {
            encoding += to_string(operation.length());
            switch (operation.type())
            {
            case graphtools::OperationType::kMatch: encoding += 'M'; break;
            case graphtools::OperationType::kMismatch: encoding += 'X'; break;
            case graphtools::OperationType::kInsertionToRef: encoding += 'I'; break;
            case graphtools::OperationType::kDeletionFromRef: encoding += 'D'; break;
            case graphtools::OperationType::kSoftclip: encoding += 'S'; break;
            case graphtools::OperationType::kMissingBases: encoding += 'N'; break;
            }
        }
        encoding += ']';
    }
    return encoding;
}

static const string graphAlignmentBamTag = "XG";

// ACGTN, including lower-case low-quality bases, covers the usual BAM input.
// Other characters use HTSlib's complete mapping and the original quality rule.
#if defined(__aarch64__) && !defined(EH_DISABLE_BAM_ENCODING_SIMD)
static std::size_t encodeCommonBases(const string& query, uint8_t* sequence, uint8_t* quality)
{
    std::size_t index = 0;
    for (; index + 16 <= query.size(); index += 16)
    {
        const auto bases = vld1q_u8(reinterpret_cast<const uint8_t*>(query.data() + index));
        const auto upper = vandq_u8(bases, vdupq_n_u8(0xdf));
        const auto a = vceqq_u8(upper, vdupq_n_u8('A'));
        const auto c = vceqq_u8(upper, vdupq_n_u8('C'));
        const auto g = vceqq_u8(upper, vdupq_n_u8('G'));
        const auto t = vceqq_u8(upper, vdupq_n_u8('T'));
        const auto n = vceqq_u8(upper, vdupq_n_u8('N'));
        const auto valid = vorrq_u8(vorrq_u8(a, c), vorrq_u8(vorrq_u8(g, t), n));
        if (vminvq_u8(valid) != 255)
            break;
        auto codes = vorrq_u8(vandq_u8(a, vdupq_n_u8(1)), vandq_u8(c, vdupq_n_u8(2)));
        codes = vorrq_u8(codes, vorrq_u8(vandq_u8(g, vdupq_n_u8(4)), vandq_u8(t, vdupq_n_u8(8))));
        codes = vorrq_u8(codes, vandq_u8(n, vdupq_n_u8(15)));
        const auto even = vget_low_u8(vuzp1q_u8(codes, codes));
        const auto odd = vget_low_u8(vuzp2q_u8(codes, codes));
        vst1_u8(sequence + index / 2, vorr_u8(vshl_n_u8(even, 4), odd));
        vst1q_u8(quality + index, vandq_u8(vceqq_u8(bases, upper), vdupq_n_u8(40)));
    }
    return index;
}
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(EH_DISABLE_BAM_ENCODING_SIMD)
__attribute__((target("avx2"))) static std::size_t encodeCommonBasesAvx2(
    const string& query, uint8_t* sequence, uint8_t* quality)
{
    std::size_t index = 0;
    for (; index + 32 <= query.size(); index += 32)
    {
        const auto bases = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(query.data() + index));
        const auto upper = _mm256_and_si256(bases, _mm256_set1_epi8(char(0xdf)));
        const auto a = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8('A'));
        const auto c = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8('C'));
        const auto g = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8('G'));
        const auto t = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8('T'));
        const auto n = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8('N'));
        const auto valid = _mm256_or_si256(_mm256_or_si256(a, c), _mm256_or_si256(_mm256_or_si256(g, t), n));
        if (_mm256_movemask_epi8(valid) != -1)
            break;
        auto codes = _mm256_or_si256(_mm256_and_si256(a, _mm256_set1_epi8(1)), _mm256_and_si256(c, _mm256_set1_epi8(2)));
        codes = _mm256_or_si256(codes, _mm256_or_si256(_mm256_and_si256(g, _mm256_set1_epi8(4)), _mm256_and_si256(t, _mm256_set1_epi8(8))));
        codes = _mm256_or_si256(codes, _mm256_and_si256(n, _mm256_set1_epi8(15)));
        const auto packedWords = _mm256_maddubs_epi16(codes, _mm256_set1_epi16(0x0110));
        const auto packed = _mm256_packus_epi16(packedWords, _mm256_setzero_si256());
        const auto bytes = _mm_unpacklo_epi64(_mm256_castsi256_si128(packed), _mm256_extracti128_si256(packed, 1));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(sequence + index / 2), bytes);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(quality + index),
            _mm256_and_si256(_mm256_cmpeq_epi8(bases, upper), _mm256_set1_epi8(40)));
    }
    return index;
}

static std::size_t encodeCommonBases(const string& query, uint8_t* sequence, uint8_t* quality)
{
    static const bool hasAvx2 = __builtin_cpu_supports("avx2");
    return hasAvx2 ? encodeCommonBasesAvx2(query, sequence, quality) : 0;
}
#else
static std::size_t encodeCommonBases(const string&, uint8_t*, uint8_t*) { return 0; }
#endif

void BamletWriter::write(
    int32_t contigId, int64_t position, const string& fragmentName, const string& query, bool isFirstMate,
    bool isReversed, bool isMateReversed, const GraphAlignment& alignment)
{

    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> ownedAlignment(bam_init1(), bam_destroy1);
    if (!ownedAlignment)
        throw std::bad_alloc();
    bam1_t* htsAlignmentPtr = ownedAlignment.get();

    htsAlignmentPtr->core.tid = contigId;
    htsAlignmentPtr->core.pos = position;
    htsAlignmentPtr->core.mtid = -1;
    htsAlignmentPtr->core.mpos = -1;
    htsAlignmentPtr->core.flag = BAM_FUNMAP;

    htsAlignmentPtr->core.flag += BAM_FPAIRED + BAM_FMUNMAP;

    if (isReversed)
        htsAlignmentPtr->core.flag += BAM_FREVERSE;
    if (isMateReversed)
        htsAlignmentPtr->core.flag += BAM_FMREVERSE;

    htsAlignmentPtr->core.flag += isFirstMate ? BAM_FREAD1 : BAM_FREAD2;

    const string alignmentEncoding = summarizeAlignment(alignment);

    htsAlignmentPtr->core.l_qname = fragmentName.length() + 1; // +1 includes the tailing '\0'
    htsAlignmentPtr->core.l_qseq = query.length();
    htsAlignmentPtr->core.n_cigar = 0; // we have no cigar sequence

    //`q->data` structure: qname-cigar-seq-qual-aux
    int seqQualLength = (int)(1.5 * query.length() + (query.length() % 2 != 0));
    htsAlignmentPtr->l_data = htsAlignmentPtr->core.l_qname + seqQualLength;
    // Include XG's two-byte name, type byte and terminating NUL in the first
    // allocation. bam_aux_append still owns tag serialization and l_data.
    htsAlignmentPtr->m_data = htsAlignmentPtr->l_data + 3 + alignmentEncoding.length() + 1;
    kroundup32(htsAlignmentPtr->m_data);
    auto* data = static_cast<uint8_t*>(realloc(htsAlignmentPtr->data, htsAlignmentPtr->m_data));
    if (!data)
        throw std::bad_alloc();
    htsAlignmentPtr->data = data;
    memcpy(htsAlignmentPtr->data, fragmentName.c_str(), htsAlignmentPtr->core.l_qname); // first set qname

    uint8_t* const htsSequencePtr = bam_get_seq(htsAlignmentPtr);
    uint8_t* const htsQualityPtr = bam_get_qual(htsAlignmentPtr);
    std::size_t index = encodeCommonBases(query, htsSequencePtr, htsQualityPtr);
    // Write each packed byte once and derive qualities directly, without an
    // intermediate vector<int>. Retain HTSlib's complete nucleotide mapping
    // and the original isupper-based quality rule (including non-ACGT bases).
    for (; index + 1 < query.length(); index += 2)
    {
        const auto first = static_cast<unsigned char>(query[index]);
        const auto second = static_cast<unsigned char>(query[index + 1]);
        htsSequencePtr[index / 2] = (seq_nt16_table[first] << 4) | seq_nt16_table[second];
        htsQualityPtr[index] = isupper(query[index]) ? 40 : 0;
        htsQualityPtr[index + 1] = isupper(query[index + 1]) ? 40 : 0;
    }
    if (index < query.length())
    {
        htsSequencePtr[index / 2] = seq_nt16_table[static_cast<unsigned char>(query[index])] << 4;
        htsQualityPtr[index] = isupper(query[index]) ? 40 : 0;
    }

    const int appendResult = bam_aux_append(
        htsAlignmentPtr, graphAlignmentBamTag.c_str(), 'Z', alignmentEncoding.length() + 1,
        reinterpret_cast<const uint8_t*>(alignmentEncoding.c_str()));
    if (appendResult != 0)
        throw std::runtime_error("Cannot append graph alignment BAM tag");

    writeQueue_.push(htsAlignmentPtr, sizeof(bam1_t) + std::size_t(htsAlignmentPtr->m_data));
    ownedAlignment.release();
}

void BamletWriter::writeHtsAlignments()
{
    bam1_t* htsAlignmentPtr(nullptr);
    while (true)
    {
        writeQueue_.pop(htsAlignmentPtr);
        if (not htsAlignmentPtr)
            return;
        if (!writeError_)
        {
            try
            {
                if (bam_write1(filePtr_->fp.bgzf, htsAlignmentPtr) <= 0)
                    throw std::logic_error("Cannot write alignment");
            }
            catch (...)
            {
                writeError_ = std::current_exception();
            }
        }
        // Continue draining after failure so locus producers cannot deadlock.
        bam_destroy1(htsAlignmentPtr);
    }
}

}
