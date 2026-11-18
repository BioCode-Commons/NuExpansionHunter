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

// Compare complete BAM records at WGS scale with compressed disk partitions.
// Build: c++ -O3 -std=c++17 bam_multiset_compare.cpp -pthread -lz -o bam_multiset_compare
// partition DIR [WORKERS] reads uncompressed BAM bytes on stdin (e.g. bgzip -dc -@4).
// compare REFERENCE_DIR CANDIDATE_DIR requires identical headers and records.
// --ignore-unused-nibble [WORKERS] permits just the unused low nibble of an odd sequence.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string_view>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <memory>
#include <memory_resource>
#include <mutex>
#include <thread>
#include <future>
#include <zlib.h>

namespace fs = std::filesystem;
using std::string;
constexpr size_t buckets = 256;

static uint32_t u32(const char* p)
{
    const auto* b = reinterpret_cast<const uint8_t*>(p);
    return uint32_t(b[0]) | uint32_t(b[1]) << 8 | uint32_t(b[2]) << 16 | uint32_t(b[3]) << 24;
}

static void readExact(std::istream& in, char* p, size_t size)
{
    if (!in.read(p, size))
        throw std::runtime_error("truncated BAM");
}

static void copyHeader(std::istream& in, std::ostream& out)
{
    auto number = [&] {
        char b[4];
        readExact(in, b, 4);
        out.write(b, 4);
        if (u32(b) > INT32_MAX)
            throw std::runtime_error("negative BAM header size");
        return u32(b);
    };
    auto copy = [&](size_t size) {
        char b[65536];
        while (size)
        {
            auto n = std::min(sizeof(b), size);
            readExact(in, b, n);
            out.write(b, n);
            size -= n;
        }
    };
    char magic[4];
    readExact(in, magic, 4);
    if (std::memcmp(magic, "BAM\1", 4) != 0)
        throw std::runtime_error("invalid BAM magic");
    out.write(magic, 4);
    copy(number());
    const auto refs = number();
    for (uint32_t i = 0; i < refs; ++i)
    {
        copy(number());
        number();
    }
}

// Each worker is the sole owner of its gzip streams. Queue whole blocks, so
// input ordering and per-bucket append order are independent of compression.
class Compressor
{
    struct Block { size_t bucket; string bytes; };
    fs::path directory_;
    std::mutex mutex_;
    std::condition_variable changed_;
    std::deque<Block> queue_;
    bool closed_ = false;
    std::exception_ptr error_;
    std::thread thread_;

    void run() noexcept
    {
        std::array<gzFile, buckets> writers{};
        try
        {
            for (;;)
            {
                std::unique_lock<std::mutex> lock(mutex_);
                changed_.wait(lock, [&] { return closed_ || !queue_.empty(); });
                if (queue_.empty()) break;
                auto block = std::move(queue_.front());
                queue_.pop_front();
                lock.unlock();
                changed_.notify_all();
                auto& writer = writers[block.bucket];
                if (!writer)
                {
                    writer = gzopen((directory_ / std::to_string(block.bucket)).c_str(), "wb1");
                    if (!writer) throw std::runtime_error("cannot open partition");
                    gzbuffer(writer, 65536);
                }
                size_t offset = 0;
                while (offset < block.bytes.size())
                {
                    const auto length = std::min(size_t(1) << 26, block.bytes.size() - offset);
                    if (gzwrite(writer, block.bytes.data() + offset, length) != int(length))
                        throw std::runtime_error("cannot write partition");
                    offset += length;
                }
            }
        }
        catch (...)
        {
            std::lock_guard<std::mutex> lock(mutex_);
            error_ = std::current_exception();
            changed_.notify_all();
        }
        for (auto writer : writers)
            if (writer && gzclose(writer) != Z_OK)
            {
                std::lock_guard<std::mutex> lock(mutex_);
                error_ = std::make_exception_ptr(std::runtime_error("partition close failed"));
            }
    }
public:
    explicit Compressor(const fs::path& dir) : directory_(dir), thread_(&Compressor::run, this) {}
    void put(size_t bucket, string bytes)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        changed_.wait(lock, [&] { return queue_.size() < 16 || error_; });
        if (error_) std::rethrow_exception(error_);
        queue_.push_back({bucket, std::move(bytes)});
        lock.unlock();
        changed_.notify_all();
    }
    void close()
    {
        { std::lock_guard<std::mutex> lock(mutex_); closed_ = true; }
        changed_.notify_all();
        if (thread_.joinable()) thread_.join();
    }
    void finish() { close(); if (error_) std::rethrow_exception(error_); }
    ~Compressor() { close(); }
};

static void partition(const fs::path& dir, unsigned jobs)
{
    if (!fs::create_directory(dir))
        throw std::runtime_error("partition destination must be new");
    std::ofstream header(dir / "header", std::ios::binary);
    header.exceptions(std::ios::badbit | std::ios::failbit);
    copyHeader(std::cin, header);
    header.close();
    std::vector<std::unique_ptr<Compressor>> encoders;
    for (unsigned j = 0; j < jobs; ++j) encoders.emplace_back(new Compressor(dir));
    std::array<string, buckets> buffers;
    std::array<uint64_t, buckets> counts{};
    char length[4], core[8];
    auto flush = [&](size_t bucket) {
        auto& buffer = buffers[bucket];
        if (!buffer.empty()) encoders[bucket % jobs]->put(bucket, std::move(buffer));
        buffer.clear();
    };
    while (std::cin.peek() != EOF)
    {
        readExact(std::cin, length, 4);
        const auto size = u32(length);
        if (size < 32 || size > INT32_MAX)
            throw std::runtime_error("invalid BAM record size");
        readExact(std::cin, core, 8);
        const auto bucket = (u32(core) * 97u + (u32(core + 4) >> 18)) % buckets;
        auto& buffer = buffers[bucket];
        const auto start = buffer.size();
        buffer.resize(start + 4 + size);
        std::memcpy(&buffer[start], length, 4);
        std::memcpy(&buffer[start + 4], core, 8);
        readExact(std::cin, &buffer[start + 12], size - 8);
        if (buffer.size() >= 256 * 1024) flush(bucket);
        ++counts[bucket];
    }
    if (std::cin.bad()) throw std::runtime_error("input read failed");
    for (size_t i = 0; i < buckets; ++i) flush(i);
    for (auto& encoder : encoders) encoder->finish();
    std::ofstream manifest(dir / "counts");
    manifest.exceptions(std::ios::badbit | std::ios::failbit);
    for (auto count : counts) manifest << count << '\n';
    manifest.close();
}

template <class Visit>
static void visitPartition(const fs::path& dir, size_t bucket, uint64_t count, Visit visit)
{
    if (!count) return;
    auto stream = gzopen((dir / std::to_string(bucket)).c_str(), "rb");
    if (!stream) throw std::runtime_error("cannot read partition");
    std::unique_ptr<gzFile_s, decltype(&gzclose)> owner(stream, gzclose);
    gzbuffer(stream, 256 * 1024);
    char length[4];
    string record;
    for (uint64_t i = 0; i < count; ++i)
    {
        if (gzread(stream, length, 4) != 4) throw std::runtime_error("truncated partition");
        const auto size = u32(length);
        if (size < 32 || size > INT32_MAX) throw std::runtime_error("invalid partition record length");
        record.resize(size);
        if (gzread(stream, &record[0], size) != int(size)) throw std::runtime_error("truncated record");
        visit(record);
    }
    if (gzread(stream, length, 1) != 0 || gzclose(owner.release()) != Z_OK)
        throw std::runtime_error("invalid partition end");
}

static size_t paddingByte(std::string_view r)
{
    const auto length = u32(r.data() + 16);
    if (!(length & 1))
        return string::npos;
    const size_t cigarCount = u32(r.data() + 12) & 65535;
    return 32 + uint8_t(r[8]) + 4 * cigarCount + length / 2;
}

static uint8_t takePadding(string& record, bool allowPadding)
{
    const auto offset = allowPadding ? paddingByte(record) : string::npos;
    if (offset == string::npos) return 0;
    auto& byte = record.at(offset);
    const uint8_t nibble = uint8_t(byte) & 15;
    byte = char(uint8_t(byte) & 0xf0);
    return nibble;
}

struct Multiplicity
{
    uint64_t remaining = 0, commonRemaining = 0;
    uint8_t common = 0;
    std::unique_ptr<std::array<uint64_t, 16>> mixed;
    void add(uint8_t nibble)
    {
        if (!remaining) common = nibble;
        ++remaining;
        if (nibble == common)
            ++commonRemaining;
        else
        {
            if (!mixed) mixed.reset(new std::array<uint64_t, 16> {});
            ++(*mixed)[nibble];
        }
    }
    bool takeRaw(uint8_t nibble)
    {
        if (nibble != common && !mixed) return false;
        auto& available = nibble == common ? commonRemaining : (*mixed)[nibble];
        if (!available) return false;
        --available;
        return true;
    }
};

static uint64_t compareBucket(const fs::path& ref, const fs::path& cand, size_t i, uint64_t count, bool allowPadding, bool zeroCandidate)
{
    // Hashing only selects a lookup bucket. Its key is the entire record and
    // std::string equality resolves collisions. Equal partition sizes plus one
    // successful removal for every candidate consume the reference exactly.
    std::unordered_map<string, Multiplicity> records;
    records.reserve(count);
    visitPartition(ref, i, count, [&](string& record) {
        const auto nibble = takePadding(record, allowPadding);
        auto& entry = records.try_emplace(std::move(record)).first->second;
        entry.add(nibble);
    });
    uint64_t unmatched = 0;
    visitPartition(cand, i, count, [&](string& record) {
        const auto nibble = takePadding(record, allowPadding);
        if (zeroCandidate && nibble) throw std::runtime_error("candidate contains nonzero padding");
        auto found = records.find(record);
        if (found == records.end() || !found->second.remaining)
            throw std::runtime_error("BAM records differ in partition " + std::to_string(i));
        --found->second.remaining;
        if (!found->second.takeRaw(nibble)) ++unmatched;
    });
    return unmatched;
}

static void compare(const fs::path& ref, const fs::path& cand, bool allowPadding, bool zeroCandidate, unsigned jobs)
{
    auto header = [](const fs::path& dir) {
        std::ifstream in(dir / "header", std::ios::binary);
        if (!in) throw std::runtime_error("missing header");
        return string(std::istreambuf_iterator<char>(in), {});
    };
    if (header(ref) != header(cand)) throw std::runtime_error("BAM headers differ");
    std::ifstream ca(ref / "counts"), cb(cand / "counts");
    uint64_t total = 0;
    std::array<uint64_t, buckets> counts;
    for (size_t i = 0; i < buckets; ++i)
    {
        uint64_t nb;
        if (!(ca >> counts[i]) || !(cb >> nb)) throw std::runtime_error("missing or incomplete partition manifest");
        if (counts[i] != nb) throw std::runtime_error("record counts differ in partition " + std::to_string(i));
        total += counts[i];
    }
    std::atomic<size_t> next {0};
    std::vector<std::future<uint64_t>> workers;
    for (unsigned j = 0; j < jobs; ++j)
        workers.push_back(std::async(std::launch::async, [&] {
            uint64_t different = 0;
            for (size_t i = next++; i < buckets; i = next++)
                different += compareBucket(ref, cand, i, counts[i], allowPadding, zeroCandidate);
            return different;
        }));
    uint64_t different = 0;
    for (auto& worker : workers) different += worker.get();
    std::cout << "identical_header=true identical_defined_record_multiset=true raw_record_multiset_equal="
              << std::boolalpha << (different == 0) << " records=" << total
              << " raw_unmatched_records_per_file=" << different << '\n';
}

// Near-coordinate streams can be compared while they are decoded. Each lane
// owns an exact counted difference of full bytes, removing matching copies as
// they arrive. General reorderings can use the disk-partitioned verifier.
struct RecordSpan { size_t offset, size; };
struct Wave
{
    std::vector<char> data;
    std::vector<std::vector<RecordSpan>> lanes;
    explicit Wave(size_t workers) : lanes(workers) {}
    uint64_t nonzeroPadding = 0;
    size_t read(std::istream& input, bool normalizePadding)
    {
        data.clear();
        for (auto& lane : lanes) lane.clear();
        size_t count = 0;
        while (count < 131072 && data.size() < 64 * 1024 * 1024)
        {
            char length[4];
            input.read(length, 4);
            if (!input && input.eof() && input.gcount() == 0) break;
            if (!input) throw std::runtime_error("truncated record length");
            const uint32_t size = u32(length);
            if (size < 32 || size > INT32_MAX) throw std::runtime_error("invalid record size");
            size_t start = data.size();
            data.resize(start + size);
            char* record = data.data() + start;
            readExact(input, record, size);
            if (normalizePadding)
            {
                const auto offset = paddingByte(std::string_view(record, size));
                if (offset != string::npos)
                {
                    if (offset >= size) throw std::runtime_error("invalid BAM sequence length");
                    nonzeroPadding += (uint8_t(record[offset]) & 15) != 0;
                    record[offset] = char(uint8_t(record[offset]) & 0xf0);
                }
            }
            const size_t lane = (u32(record + 4) + 31u * u32(record)) % lanes.size();
            lanes[lane].push_back({start, size});
            ++count;
        }
        return count;
    }
};
struct Difference
{
    // Each lane owns the storage behind its full-byte keys. Look up the
    // incoming span directly, and copy only records that lack a matching key.
    std::pmr::unsynchronized_pool_resource storage;
    std::pmr::unordered_map<std::string_view, int64_t> counts { &storage };
    size_t bytes = 0;
    void add(const Wave& wave, RecordSpan span, int direction)
    {
        const std::string_view key(wave.data.data() + span.offset, span.size);
        const auto found = counts.find(key);
        if (found == counts.end())
        {
            char* copy = static_cast<char*>(storage.allocate(key.size(), alignof(char)));
            std::memcpy(copy, key.data(), key.size());
            counts.emplace(std::string_view(copy, key.size()), direction);
            bytes += key.size();
        }
        else if ((found->second += direction) == 0)
        {
            const auto owned = found->first;
            counts.erase(found);
            storage.deallocate(const_cast<char*>(owned.data()), owned.size(), alignof(char));
            bytes -= owned.size();
        }
    }
};
static void compareStreams(const fs::path& first, const fs::path& second, size_t workers, bool zeroPaddedCandidate)
{
    std::vector<char> leftBuffer(1 << 20), rightBuffer(1 << 20);
    std::ifstream left, right;
    left.rdbuf()->pubsetbuf(leftBuffer.data(), leftBuffer.size());
    right.rdbuf()->pubsetbuf(rightBuffer.data(), rightBuffer.size());
    left.open(first, std::ios::binary);
    right.open(second, std::ios::binary);
    std::ostringstream ha, hb;
    copyHeader(left, ha);
    copyHeader(right, hb);
    if (ha.str() != hb.str()) throw std::runtime_error("BAM headers differ");
    Wave a(workers), b(workers), nextA(workers), nextB(workers);
    std::vector<Difference> differences(workers);
    uint64_t leftTotal = 0, rightTotal = 0;
    size_t peak = 0;
    auto readLeft = [&] { return nextA.read(left, zeroPaddedCandidate); };
    auto readRight = [&] { return nextB.read(right, zeroPaddedCandidate); };
    auto leftJob = std::async(std::launch::async, readLeft);
    auto rightJob = std::async(std::launch::async, readRight);
    for (;;)
    {
        const auto na = leftJob.get(), nb = rightJob.get();
        std::swap(a, nextA);
        std::swap(b, nextB);
        if (b.nonzeroPadding)
            throw std::runtime_error("candidate contains nonzero padding; use the partitioned comparison");
        leftTotal += na;
        rightTotal += nb;
        if (!na && !nb) break;
        leftJob = std::async(std::launch::async, readLeft);
        rightJob = std::async(std::launch::async, readRight);
        std::vector<std::future<void>> jobs;
        for (size_t lane = 0; lane < workers; ++lane)
            jobs.push_back(std::async(std::launch::async, [&, lane] {
                const auto& aa = a.lanes[lane];
                const auto& bb = b.lanes[lane];
                for (size_t i = 0; i < std::max(aa.size(), bb.size()); ++i)
                {
                    if (i < aa.size() && i < bb.size() && aa[i].size == bb[i].size
                        && !std::memcmp(a.data.data() + aa[i].offset, b.data.data() + bb[i].offset, aa[i].size))
                        continue;
                    if (i < aa.size()) differences[lane].add(a, aa[i], 1);
                    if (i < bb.size()) differences[lane].add(b, bb[i], -1);
                }
            }));
        for (auto& job : jobs) job.get();
        size_t bytes = 0;
        for (const auto& difference : differences) bytes += difference.bytes;
        peak = std::max(peak, bytes);
        if (bytes > 1024u * 1024u * 1024u)
            throw std::runtime_error("use partitioned comparison for this degree of reordering or difference");
    }
    if (leftTotal != rightTotal) throw std::runtime_error("record counts differ");
    for (const auto& difference : differences)
        if (!difference.counts.empty()) throw std::runtime_error("BAM record multisets differ");
    // When every candidate pad is zero, each nonzero reference pad leaves
    // exactly one unmatched raw record on each side, including duplicate keys.
    std::cout << "identical_header=true identical_defined_record_multiset=true raw_record_multiset_equal="
              << std::boolalpha << (a.nonzeroPadding + nextA.nonzeroPadding == 0) << " records=" << leftTotal
              << " raw_unmatched_records_per_file=" << a.nonzeroPadding + nextA.nonzeroPadding << " peak_unmatched_bytes=" << peak << '\n';
}

int main(int argc, char** argv)
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    try
    {
        if ((argc == 3 || argc == 4) && string(argv[1]) == "partition")
        {
            int jobs = argc == 4 ? std::stoi(argv[3]) : 1;
            if (jobs < 1 || jobs > 256) throw std::runtime_error("invalid compression worker count");
            partition(argv[2], jobs);
        }
        else if ((argc == 5 || (argc == 6 && string(argv[5]) == "--zero-padded-candidate")) && string(argv[1]) == "stream")
        {
            int jobs = std::stoi(argv[4]);
            if (jobs < 1 || jobs > 256) throw std::runtime_error("invalid stream worker count");
            compareStreams(argv[2], argv[3], jobs, argc == 6);
        }
        else if ((argc == 4 || ((argc == 5 || argc == 6) && (string(argv[4]) == "--ignore-unused-nibble" || string(argv[4]) == "--zero-padded-candidate")))
                 && string(argv[1]) == "compare")
            {
            int jobs = argc == 6 ? std::stoi(argv[5]) : 1;
            if (jobs < 1 || jobs > 256) throw std::runtime_error("invalid compare worker count");
            compare(argv[2], argv[3], argc >= 5, argc >= 5 && string(argv[4]) == "--zero-padded-candidate", jobs);
        }
        else throw std::runtime_error("expected partition DIR [WORKERS], compare REF_DIR CAND_DIR [--ignore-unused-nibble [WORKERS]], or stream REF CAND WORKERS [--zero-padded-candidate]");
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
