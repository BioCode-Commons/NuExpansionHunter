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
//
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from ExpansionHunter v5.0.0: BAM serialization uses producer batches, asynchronous sinks, optional codecs, and explicit completion.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "io/BamletWriter.hh"
#include "core/OrderedBlocks.hh"
#include "core/VerifyPolicy.hh"
#include "core/ThreadPool.hh"

#include <algorithm>
#include <cerrno>
#ifdef __linux__
#include <fcntl.h>
#include <atomic>
#include <boost/align/aligned_allocator.hpp>
#include <sys/stat.h>
#include <unistd.h>
#include "htslib/hfile.h"
#endif
#include <cstring>
#include <deque>
#include <limits>

#ifdef EH_HAVE_ISAL
#include <isa-l/igzip_lib.h>
#endif
#include <utility>

#if defined(__aarch64__)
#include <arm_neon.h>
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#endif

// cppcheck-suppress missingInclude
#include "htslib/bgzf.h"
#include "htslib/hts_endian.h"
// cppcheck-suppress missingInclude
#include "htslib/kseq.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"
#include "io/BamEncodingKernels.hh"

using graphtools::GraphAlignment;
using graphtools::ReferenceInterval;
using std::string;
using std::to_string;

namespace ehunter
{

// Open ordinary local files once. HTSlib owns the descriptor after hdopen.
static htsFile* openBamOutput(const string& path, bool fast, int& writebackFd)
{
#ifdef __linux__
    if (fast && path != "-" && path.find(':') == string::npos && path.find(HTS_IDX_DELIM) == string::npos)
    {
        const int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0666);
        if (fd < 0) return nullptr;
        auto* stream = hdopen(fd, "w");
        if (!stream)
        {
            close(fd);
            return nullptr;
        }
        auto* file = hts_hopen(stream, path.c_str(), "wb1");
        if (!file) hclose_abruptly(stream);
        else
        {
            struct stat status;
            if (fstat(fd, &status) == 0 && S_ISREG(status.st_mode)) writebackFd = fd;
        }
        return file;
    }
#else
    (void)writebackFd;
#endif
    return hts_open(path.c_str(), fast ? "wb1" : "wb");
}

// Establish compression and the header before starting the output consumer.
static htsFile* openBamletFile(const string& path, int compressionThreads, bool fastCompression, int& writebackFd)
{
    std::unique_ptr<htsFile, decltype(&hts_close)> file(openBamOutput(path, fastCompression, writebackFd), hts_close);
    if (!file)
    {
        throw std::runtime_error("Failed to open BAMlet output " + path);
    }
    if (writebackFd >= 0 && hts_set_opt(file.get(), HTS_OPT_BLOCK_SIZE, 2 * 1024 * 1024) != 0)
        throw std::runtime_error("Cannot configure the local BAM output buffer");
    if (compressionThreads > 0 && hts_set_threads(file.get(), compressionThreads) != 0)
    {
        throw std::runtime_error("Failed to initialize BAMlet compression threads");
    }
    return file.release();
}

static size_t bgzfWorkspace(size_t bytes)
{
    const size_t blocks = bytes / BGZF_BLOCK_SIZE + (bytes % BGZF_BLOCK_SIZE != 0);
    if (blocks > std::numeric_limits<size_t>::max() / BGZF_MAX_BLOCK_SIZE)
        throw std::length_error("BAM batch is too large to compress");
    return blocks * BGZF_MAX_BLOCK_SIZE;
}

static int compressFastBlock(uint8_t* out, size_t* capacity, const uint8_t* input, size_t length)
{
#ifdef EH_HAVE_ISAL
    // ISA-L writes the raw DEFLATE payload and standard gzip trailer. The
    // BGZF extra field stores the complete block size minus one, little endian.
    static constexpr uint8_t header[] = { 31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0 };
    static thread_local std::unique_ptr<isal_zstream> codec;
    if (!codec)
        codec.reset(new isal_zstream);
    auto& stream = *codec;
    isal_deflate_stateless_init(&stream);
    // Each encoder sees many BGZF blocks made of the same BAM record layout.
    // Estimate one Huffman table from the first substantial block. ISA-L assigns
    // codes to every symbol, including symbols absent from that first block.
    static thread_local std::unique_ptr<isal_hufftables> tables;
    static thread_local bool trained = false;
    if (!trained && length >= BGZF_BLOCK_SIZE / 2)
    {
        trained = true;
        std::unique_ptr<isal_huff_histogram> histogram(new isal_huff_histogram {});
        tables.reset(new isal_hufftables);
        isal_update_histogram(const_cast<uint8_t*>(input), static_cast<int>(length), histogram.get());
        if (isal_create_hufftables(tables.get(), histogram.get()) != 0) tables.reset();
    }
    if (!tables || isal_deflate_set_hufftables(&stream, tables.get(), IGZIP_HUFFTABLE_CUSTOM) != COMP_OK)
    {
        static thread_local std::vector<uint8_t> levelBuffer(ISAL_DEF_LVL1_DEFAULT);
        stream.level = 1;
        stream.level_buf = levelBuffer.data();
        stream.level_buf_size = static_cast<uint32_t>(levelBuffer.size());
    }
    stream.gzip_flag = IGZIP_GZIP_NO_HDR;
    stream.next_in = const_cast<uint8_t*>(input);
    stream.avail_in = static_cast<uint32_t>(length);
    stream.next_out = out + sizeof(header);
    stream.avail_out = static_cast<uint32_t>(*capacity - sizeof(header));
    if (isal_deflate_stateless(&stream) == COMP_OK && stream.avail_in == 0)
    {
        const size_t bytes = sizeof(header) + stream.total_out;
        std::memcpy(out, header, sizeof(header));
        out[16] = static_cast<uint8_t>(bytes - 1);
        out[17] = static_cast<uint8_t>((bytes - 1) >> 8);
        *capacity = bytes;
        return 0;
    }
#endif
    return bgzf_compress(out, capacity, input, length, 1);
}

// Contract [K-32]
// Ensures:  the concatenated blocks decode to the submitted batches in submission order; file offsets are
//           contiguous. The codecs (zlib, libdeflate, ISA-L) are trusted.
static BamRecordSink::Bytes compressBlocks(const BamRecordSink::Bytes& input)
{
    BamRecordSink::Bytes output(bgzfWorkspace(input.size()), boost::container::default_init);
    size_t written = 0;
    for (size_t start = 0; start < input.size();)
    {
        const size_t blockSize = std::min(size_t(BGZF_BLOCK_SIZE), input.size() - start);
        size_t size = BGZF_MAX_BLOCK_SIZE;
        if (compressFastBlock(output.data() + written, &size, input.data() + start, blockSize) != 0)
            throw std::runtime_error("Cannot compress BAM block");
        written += size;
        start += blockSize;
    }
    output.resize(written);
    return output;
}

static void validateMapping(const ReferenceContigInfo& contigInfo, const LocusSpecification& locusSpec)
{
    for (const auto& nodeAndRegion : locusSpec.referenceProjectionOfNodes())
    {
        auto nodeId = nodeAndRegion.first;
        const auto& region = nodeAndRegion.second;
        const auto& contig = contigInfo.getContigName(region.contigIndex());
        // parseRegion consumes nonnegative decimal int coordinates and a
        // single colon. Keep its validation for other encodings.
        const bool direct = region.start() >= 0 && region.end() >= 0
            && region.start() <= std::numeric_limits<int>::max() && region.end() <= std::numeric_limits<int>::max()
            && contig.find(':') == string::npos;
        const auto referenceInterval = direct ? ReferenceInterval(contig, region.start(), region.end())
            : ReferenceInterval::parseRegion(encode(contigInfo, region));
        graphtools::NodeReferenceMapping mapping(locusSpec.regionGraph(), nodeId, referenceInterval);
        (void)mapping;
    }
}

#ifdef __linux__
// Assign compressed byte offsets in FIFO order. Independent positioned writes
// can complete in any order. The next ordinary HTSlib write drains the jobs
// and resumes at their exact end, including the final BGZF EOF block.
class PositionalBamOutput
{
public:
    static constexpr size_t blockSize = 1024 * 1024;

    PositionalBamOutput(BGZF* stream, int fd, unsigned writers, size_t slots)
        : stream_(stream), fd_(fd), pool_(writers), slots_(slots)
    {
        // Open the same inode through the original descriptor, without changing
        // HTSlib's buffered file description or reopening a replaceable pathname.
        const auto name = "/proc/self/fd/" + std::to_string(fd);
        direct_.fd = open(name.c_str(), O_WRONLY | O_DIRECT | O_CLOEXEC);
    }

    void append(const uint8_t* data, size_t size)
    {
        if (!active_)
        {
            if (bgzf_flush(stream_) != 0 || (next_ = hseek(stream_->fp, 0, SEEK_CUR)) < 0)
                throw std::runtime_error("Cannot position BAM output after the header");
            active_ = true;
        }
        while (size)
        {
            if (staged_.empty()) staged_.reserve(blockSize);
            // The prefix after the BGZF header and every partial suffix use
            // buffered pwrite. All interior ranges begin at MiB boundaries.
            const size_t limit = blockSize - static_cast<uint64_t>(next_) % blockSize;
            const size_t count = std::min(size, limit - staged_.size());
            staged_.insert(staged_.end(), data, data + count);
            data += count;
            size -= count;
            if (staged_.size() == limit) submit();
        }
    }

    void synchronize()
    {
        if (!active_) return;
        std::exception_ptr error;
        try { submit(); } catch (...) { error = std::current_exception(); }
        while (!pending_.empty())
        {
            try { retire(); } catch (...) { if (!error) error = std::current_exception(); }
        }
        active_ = false;
        if (error) std::rethrow_exception(error);
        if (hseek(stream_->fp, next_, SEEK_SET) != next_)
            throw std::runtime_error("Cannot resume buffered BAM output");
    }

private:
    void retire()
    {
        auto first = std::move(pending_.front());
        pending_.pop_front();
        first.get();
    }

    void submit()
    {
        if (staged_.empty()) return;
        if (pending_.size() == slots_) retire();
        const size_t size = staged_.size();
        if (next_ > INT64_MAX - static_cast<int64_t>(size)) throw std::length_error("BAM file too large");
        auto data = std::make_shared<AlignedBytes>(std::move(staged_));
        pending_.emplace_back();
        const auto offset = next_;
        try
        {
            pending_.back() = pool_.push([this, data, offset](int) {
                bool direct = direct_.fd >= 0 && direct_.supported.load(std::memory_order_relaxed)
                    && offset % blockSize == 0 && data->size() == blockSize;
                size_t written = 0;
                while (written < data->size())
                {
                    const int fd = direct ? direct_.fd : fd_;
                    const auto n = pwrite(fd, data->data() + written, data->size() - written, offset + written);
                    if (n < 0 && direct && (errno == EINVAL || errno == EOPNOTSUPP))
                    {
                        direct_.supported.store(false, std::memory_order_relaxed);
                        direct = false;
                        continue;
                    }
                    if (n < 0 && errno == EINTR) continue;
                    if (n <= 0) throw std::runtime_error("Cannot write compressed BAM range");
                    written += static_cast<size_t>(n);
                    if (written % 4096) direct = false;
                }
                return written;
            });
        }
        catch (...) { pending_.pop_back(); throw; }
        next_ += size;
    }

    using AlignedBytes = std::vector<uint8_t, boost::alignment::aligned_allocator<uint8_t, 4096>>;
    struct DirectFile
    {
        ~DirectFile() { if (fd >= 0) close(fd); }
        int fd = -1;
        std::atomic<bool> supported { true };
    };
    BGZF* stream_;
    int fd_;
    DirectFile direct_;
    ctpl::thread_pool pool_;
    size_t slots_;
    int64_t next_ = 0;
    bool active_ = false;
    AlignedBytes staged_;
    std::deque<std::future<size_t>> pending_;
};
constexpr size_t PositionalBamOutput::blockSize;
#endif

BamletWriter::BamletWriter(
    const string& bamletPath, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
    int compressionThreads, std::size_t maxQueueBytes, bool fastCompression, unsigned validationThreads)
    : filePtr_(openBamletFile(bamletPath, fastCompression ? 0 : compressionThreads, fastCompression, fileDescriptor_), hts_close)
    , bamHeader_(bam_hdr_init(), bam_hdr_destroy)
    , contigInfo_(contigInfo)
    , catalog_(regionCatalog)
{
    if (!regionCatalog.mappedBytes().data)
    {
        auto validate = [&](int, size_t begin, size_t end) {
            for (size_t index = begin; index < end; ++index)
                validateMapping(contigInfo, regionCatalog[index]);
            return true;
        };
        orderedBlocks(regionCatalog.size(), 512, validationThreads, validate, [](bool) {});
    }

    writeHeader();
    bamContigIds_.reserve(contigInfo_.numContigs());
    for (int index = 0; index != contigInfo_.numContigs(); ++index)
    {
        const auto& name = contigInfo_.getContigName(index);
        bamContigIds_.push_back(name.empty() ? -1 : bam_name2id(bamHeader_.get(), name.c_str()));
    }
    BamRecordSink::PrepareAccounted prepare;
    if (fastCompression && compressionThreads > 0)
    {
        encoderPool_.reset(new ctpl::thread_pool(compressionThreads));
        prepare = [this](BamRecordSink::Bytes bytes, BamRecordSink::SetPreparedCapacity release) {
            auto input = std::make_shared<BamRecordSink::Bytes>(std::move(bytes));
            return encoderPool_->push([input, release](int) {
                auto encoded = compressBlocks(*input);
                const size_t inputBytes = input->capacity();
                BamRecordSink::Bytes().swap(*input);
                // The released raw allocation pays for the temporary copy while
                // both compressed buffers exist. Retain the existing slab for
                // incompressible batches or when compaction would save little.
                if (encoded.size() <= inputBytes && encoded.size() < encoded.capacity() / 2)
                {
                    BamRecordSink::Bytes compact(encoded.begin(), encoded.end());
                    encoded.swap(compact);
                }
                release(encoded.capacity());
                return encoded;
            });
        };
    }
#ifdef __linux__
    const size_t slots = maxQueueBytes ? std::min(size_t(32), maxQueueBytes / (4 * PositionalBamOutput::blockSize)) : 16;
    if (encoderPool_ && fileDescriptor_ >= 0 && slots >= 2)
    {
        const unsigned writers = std::max(1, std::min(4, compressionThreads / 2));
        positional_.reset(new PositionalBamOutput(filePtr_->fp.bgzf, fileDescriptor_, writers, slots));
        if (maxQueueBytes) maxQueueBytes -= (slots + 1) * PositionalBamOutput::blockSize;
    }
#endif
    sink_.reset(new BamRecordSink([this](const bam1_t* record) {
#ifdef __linux__
        if (positional_) positional_->synchronize();
#endif
        if (bam_write1(filePtr_->fp.bgzf, record) <= 0)
            throw std::logic_error("Cannot write alignment");
    }, maxQueueBytes, [this](const uint8_t* bytes, size_t size) {
        auto* file = filePtr_->fp.bgzf;
        // Preserve the relative order of the header, individual bam_write1
        // calls and pre-compressed batches on the one file-owning thread.
#ifdef __linux__
        if (positional_)
        {
            positional_->append(bytes, size);
            return;
        }
#endif
        ssize_t written;
        if (encoderPool_)
        {
            if (bgzf_flush(file) != 0) throw std::logic_error("Cannot flush BAM block");
            written = bgzf_raw_write(file, bytes, size);
        }
        else written = bgzf_write(file, bytes, size);
        if (written != static_cast<ssize_t>(size))
            throw std::logic_error("Cannot write alignment batch");
    }, std::move(prepare)));
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
    if (sink_)
    {
        try
        {
            sink_->finish();
        }
        catch (...)
        {
            if (!finishError_)
                finishError_ = std::current_exception();
        }
    }
    encoderPool_.reset();
#ifdef __linux__
    if (positional_)
    {
        try { positional_->synchronize(); }
        catch (...) { if (!finishError_) finishError_ = std::current_exception(); }
        positional_.reset();
    }
#endif
    const int closeResult = filePtr_ ? hts_close(filePtr_.release()) : 0;
    if (closeResult != 0 && !finishError_)
        finishError_ = std::make_exception_ptr(std::runtime_error("Cannot close alignment output"));
    if (finishError_)
        std::rethrow_exception(finishError_);
}

thread_local BamletWriter::Batch* BamletWriter::activeBatch_ = nullptr;

BamletWriter::Batch::Batch(BamletWriter* writer, BatchBuffer* buffer)
    : writer_(writer), owned_(writer), buffer_(buffer ? buffer : &owned_)
{
    if (writer_)
    {
        if (buffer_->owner_ != writer_)
            throw std::logic_error("BAM batch buffer belongs to a different writer");
        previous_ = activeBatch_;
        activeBatch_ = this;
    }
}

BamletWriter::Batch::~Batch() noexcept
{
    if (writer_ && !finished_ && activeBatch_ == this)
        activeBatch_ = previous_;
}

void BamletWriter::BatchBuffer::flush()
{
    if (!count)
        return;
    const auto workspace = owner_->encoderPool_ ? bgzfWorkspace(data.size()) : 0;
    owner_->sink_->submitEncoded(std::move(data), count, workspace);
    count = 0;
}

uint8_t* BamletWriter::Batch::append(size_t recordBytes)
{
    const size_t maxBytes = 256 * 1024;
    auto& data = buffer_->data;
    if (buffer_->count && (recordBytes > maxBytes || data.size() > maxBytes - recordBytes))
        buffer_->flush();
    const size_t offset = data.size();
    // A long-lived worker fills this buffer to the batch limit. Reserve it
    // once per flush, avoiding repeated copies as records are appended.
    if (buffer_ != &owned_ && data.capacity() == 0)
        data.reserve(std::max(maxBytes, recordBytes));
    // The serializer below assigns every byte, including odd-sequence padding.
    data.resize(offset + recordBytes, boost::container::default_init);
    ++buffer_->count;
    return data.data() + offset;
}

size_t BamletWriter::Batch::locusIndex(graphtools::TextView id)
{
    if (!buffer_->hasIndex || graphtools::TextView(buffer_->lastId) != id)
    {
        const auto& catalog = writer_->catalog_;
        buffer_->lastIndex = catalog.findLocusIndex(id);
        buffer_->nodeCount = catalog.nodeCount(buffer_->lastIndex);
        buffer_->cachedNodes = 0;
        buffer_->lastId.assign(id.data(), id.size());
        buffer_->hasIndex = true;
    }
    return buffer_->lastIndex;
}

boost::optional<GenomicRegion> BamletWriter::Batch::projection(graphtools::NodeId node) const
{
    if (node >= buffer_->projections.size())
        return writer_->catalog_.referenceRegion(buffer_->lastIndex, node);
    const uint8_t bit = static_cast<uint8_t>(1u << node);
    if (!(buffer_->cachedNodes & bit))
    {
        buffer_->projections[node] = writer_->catalog_.referenceRegion(buffer_->lastIndex, node);
        buffer_->cachedNodes |= bit;
    }
    return buffer_->projections[node];
}

void BamletWriter::Batch::finish()
{
    if (!writer_ || finished_)
        return;
    if (activeBatch_ != this)
        throw std::logic_error("BAM output batches must finish in stack order");
    activeBatch_ = previous_;
    finished_ = true;
    if (buffer_ == &owned_)
        owned_.flush();
}

void BamletWriter::submitRecord(BamRecordSink::Record record)
{
    BamRecordSink::Records single;
    single.push_back(std::move(record));
    sink_->submit(std::move(single), BamRecordSink::WriteErrorPolicy::kDeferUntilFinish);
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

// Check per-node bounds without constructing subpaths. Invalid inputs use
// generateSubpathForEachNode below to preserve v5.0.0's Path constructor
// exceptions.
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
    auto* batch = activeBatch_ && activeBatch_->writer_ == this ? activeBatch_ : nullptr;
    const auto locusIndex = batch ? batch->locusIndex(locusId) : catalog_.findLocusIndex(locusId);
    const auto nodeCount = batch ? batch->nodeCount() : catalog_.nodeCount(locusIndex);
    const auto& path = alignment.path();
    if (!nodeSubpathsAreValid(path))
        (void)graphtools::generateSubpathForEachNode(path);
    int32_t contigId = -1;
    int64_t position = -1;
    for (std::size_t index = 0; index < path.numNodes(); ++index)
    {
        const auto node = path.getNodeIdByIndex(index);
        if (node >= nodeCount)
            throw std::logic_error("Invalid node " + std::to_string(node));
        const auto region = batch ? batch->projection(node) : catalog_.referenceRegion(locusIndex, node);
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
static void summarizeAlignment(const GraphAlignment& alignment, string& encoding)
{
    const auto id = alignment.path().graphRawPtr()->graphId();
    encoding.assign(id.data(), id.size());
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
}

static const string graphAlignmentBamTag = "XG";

// Letters of either case use the vector packers; other characters use HTSlib's complete mapping and v5.0.0's
// quality rule in the tail.
#if defined(__aarch64__) && defined(__ARM_NEON)
static std::size_t encodeCommonBases(const string& query, uint8_t* sequence, uint8_t* quality)
{
    return be_encode_neon_blocks(query.data(), query.size(), sequence, quality, seq_nt16_table);
}
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
static std::size_t encodeCommonBases(const string& query, uint8_t* sequence, uint8_t* quality)
{
    static const bool hasAvx2 = !verify::kNoAvx2 && __builtin_cpu_supports("avx2");
    return hasAvx2 ? be_encode_avx2_blocks(query.data(), query.size(), sequence, quality, seq_nt16_table) : 0;
}
#else
static std::size_t encodeCommonBases(const string&, uint8_t*, uint8_t*) { return 0; }
#endif

// Contract [K-31]: the vector block (be_encode_avx2_blocks, or the NEON block) followed by be_encode_tail
// equals be_encode_reference (BamEncodingKernels.hh). For an odd length the unused low nibble of the last
// sequence byte is 0, where v5.0.0 left it unset; decoded records are equal.
static void encodeReadPayload(const string& query, uint8_t* htsSequencePtr, uint8_t* htsQualityPtr)
{
    const std::size_t index = verify::kScalarOnly ? 0 : encodeCommonBases(query, htsSequencePtr, htsQualityPtr);
    be_encode_tail(index, query.data(), query.size(), htsSequencePtr, htsQualityPtr, seq_nt16_table);
    if (verify::kShadow)
    {
        std::vector<uint8_t> sequence((query.size() + 1) / 2), quality(query.size());
        be_encode_reference(query.data(), query.size(), sequence.data(), quality.data(), seq_nt16_table);
        if (!std::equal(sequence.begin(), sequence.end(), htsSequencePtr)
            || !std::equal(quality.begin(), quality.end(), htsQualityPtr))
            throw std::logic_error("BAM payload encoding differs from the scalar reference");
    }
}

void BamletWriter::write(
    int32_t contigId, int64_t position, const string& fragmentName, const string& query, bool isFirstMate,
    bool isReversed, bool isMateReversed, const GraphAlignment& alignment)
{

    auto* batch = activeBatch_ && activeBatch_->writer_ == this ? activeBatch_ : nullptr;
    string temporary;
    auto& alignmentEncoding = batch ? batch->buffer_->alignmentEncoding : temporary;
    summarizeAlignment(alignment, alignmentEncoding);
    const size_t nameSize = fragmentName.size() + 1;
    const size_t seqSize = (query.size() + 1) / 2;
    const size_t dataSize = nameSize + seqSize + query.size() + 4 + alignmentEncoding.size();
    const unsigned flags = BAM_FUNMAP | BAM_FPAIRED | BAM_FMUNMAP | (isReversed ? BAM_FREVERSE : 0)
        | (isMateReversed ? BAM_FMREVERSE : 0) | (isFirstMate ? BAM_FREAD1 : BAM_FREAD2);

    // The streaming writer emits an unmapped BAM core, QNAME, SEQ, QUAL and a
    // single XG:Z tag. Encode this fixed layout once in producer-owned storage.
    // Keep bam_write1's size checks and explicit little-endian representation.
    if (batch)
    {
        if (nameSize > 255 || position > INT32_MAX || dataSize > INT32_MAX - 32)
            throw std::logic_error("Cannot write alignment");
        uint8_t* record = batch->append(36 + dataSize);
        u32_to_le(32 + dataSize, record);
        u32_to_le(contigId, record + 4);
        u32_to_le(position, record + 8);
        u32_to_le(nameSize, record + 12);
        u32_to_le(flags << 16, record + 16);
        u32_to_le(query.size(), record + 20);
        u32_to_le(uint32_t(-1), record + 24);
        u32_to_le(uint32_t(-1), record + 28);
        u32_to_le(0, record + 32);
        uint8_t* data = record + 36;
        std::memcpy(data, fragmentName.c_str(), nameSize);
        auto* sequence = data + nameSize;
        auto* qualities = sequence + seqSize;
        encodeReadPayload(query, sequence, qualities);
        auto* tag = qualities + query.size();
        std::memcpy(tag, "XGZ", 3);
        std::memcpy(tag + 3, alignmentEncoding.c_str(), alignmentEncoding.size() + 1);
        return;
    }

    BamRecordSink::Record ownedAlignment(bam_init1());
    if (!ownedAlignment)
        throw std::bad_alloc();
    auto* htsAlignmentPtr = ownedAlignment.get();
    htsAlignmentPtr->core.tid = contigId;
    htsAlignmentPtr->core.pos = position;
    htsAlignmentPtr->core.mtid = -1;
    htsAlignmentPtr->core.mpos = -1;
    htsAlignmentPtr->core.flag = flags;
    htsAlignmentPtr->core.l_qname = nameSize;
    htsAlignmentPtr->core.l_qseq = query.size();
    htsAlignmentPtr->l_data = nameSize + seqSize + query.size();
    htsAlignmentPtr->m_data = dataSize;
    kroundup32(htsAlignmentPtr->m_data);
    htsAlignmentPtr->data = static_cast<uint8_t*>(malloc(htsAlignmentPtr->m_data));
    if (!htsAlignmentPtr->data)
        throw std::bad_alloc();
    std::memcpy(htsAlignmentPtr->data, fragmentName.c_str(), nameSize);
    encodeReadPayload(query, bam_get_seq(htsAlignmentPtr), bam_get_qual(htsAlignmentPtr));

    const int appendResult = bam_aux_append(
        htsAlignmentPtr, graphAlignmentBamTag.c_str(), 'Z', alignmentEncoding.length() + 1,
        reinterpret_cast<const uint8_t*>(alignmentEncoding.c_str()));
    if (appendResult != 0)
        throw std::runtime_error("Cannot append graph alignment BAM tag");

    submitRecord(std::move(ownedAlignment));
}

}
