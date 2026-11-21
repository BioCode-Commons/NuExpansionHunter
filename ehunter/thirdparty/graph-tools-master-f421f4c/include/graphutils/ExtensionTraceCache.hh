// Copyright 2026 NuExpansionHunter contributors
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "graphalign/dagAligner/Details.hh"

namespace graphtools
{
namespace extension_cache
{
using graphalign::dagAligner::Cigar;

struct Counters
{
    uint64_t eligible = 0, hits = 0, published = 0, evicted = 0, traceTooLarge = 0;
    uint64_t busyLookups = 0, busyPublishes = 0, shadowChecks = 0, expandedHits = 0;
};

inline uint64_t nextPairToken()
{
    static std::atomic<uint64_t> next{1};
    auto token = next.load(std::memory_order_relaxed);
    while (token != std::numeric_limits<uint64_t>::max())
        if (next.compare_exchange_weak(token, token + 1, std::memory_order_relaxed))
            return token;
    // Disable reuse after saturation. A retired token must never name another input.
    return 0;
}

class Table
{
public:
    struct Key
    {
        uint64_t pair;
        size_t locus;
        bool mateSlot, reversedRead, reversedExtension;
    };

    bool find(const Key& key, const std::string& query, const std::string& target, Counters& counters,
        std::vector<Cigar>& cigars, int& score)
    {
        if (!key.pair || !fitsScoreWidth(query.size(), target.size()))
            return false;
        ++counters.eligible;
        PackedTrace trace;
        if (!snapshot(key, query, target, counters, trace, score))
            return false;
        // Only bounded byte copies happen under the bucket lock. The caller
        // owns the decoded CIGARs before it changes their mapped node IDs.
        cigars.clear();
        cigars.resize(trace.cigarCount);
        size_t nextOperation = 0;
        for (size_t i = 0; i < trace.cigarCount; ++i)
        {
            auto& cigar = cigars[i];
            cigar.reserve(trace.cigarLengths[i]);
            for (size_t n = 0; n < trace.cigarLengths[i]; ++n)
            {
                const auto op = trace.ops[nextOperation++];
                cigar.push_back(Cigar::Operation(static_cast<Cigar::OpCode>(op >> kValueBits), op & kValueMask));
            }
        }
        ++counters.hits;
        counters.expandedHits += trace.operationCount > kOriginalOperations;
        return true;
    }

    void publish(const Key& key, const std::string& query, const std::string& target, int score,
        const std::vector<Cigar>& cigars, Counters& counters)
    {
        if (!key.pair || !fitsScoreWidth(query.size(), target.size()) || cigars.empty()
            || score < std::numeric_limits<int16_t>::min() || score > 5 * static_cast<int>(query.size()))
            return;
        // In a single linear node an endpoint consuming t target bases
        // costs at least 2*(t-Q) versus the all-match score. At L and
        // beyond that lower bound is strictly worse than the donor optimum.
        // All optimal endpoints and their ordered traces are inside [0,L).
        const size_t prefix = query.size() + (5 * static_cast<int>(query.size()) - score) / 2 + 1;
        if (prefix > target.size())
            return;
        PackedTrace trace;
        if (!pack(cigars, trace))
        {
            ++counters.traceTooLarge;
            return;
        }
        auto& bucket = buckets_[bucketIndex(key, query.size())];
        const std::unique_lock<std::mutex> lock(bucket.mutex, std::try_to_lock);
        if (!lock.owns_lock())
        {
            ++counters.busyPublishes;
            return;
        }
        auto& slot = bucket.donors[bucket.next];
        if (slot.key.pair)
            ++counters.evicted;
        slot.key = key;
        slot.queryLength = query.size();
        slot.prefixLength = prefix;
        slot.score = score;
        std::memcpy(slot.sequence.data(), query.data(), query.size());
        std::memcpy(slot.sequence.data() + query.size(), target.data(), prefix);
        copyTrace(trace, slot.trace);
        bucket.next = (bucket.next + 1) % kWays;
        ++counters.published;
    }

private:
    friend struct TableTestAccess;
    // The full-width score guard implies Q+L <= (32767-10)/10 when
    // L <= target length. Packed operations store a 4-bit opcode and a
    // 12-bit value. Traces outside the format use ordinary alignment.
    enum : size_t
    {
        kBuckets = 4096,
        kWays = 4,
        kSequenceBytes = 3275,
        kCigars = 16,
        kOriginalOperations = 128,
        kOperations = 256,
        kValueBits = 12,
        kValueMask = (1u << kValueBits) - 1
    };
    static_assert(Cigar::NODE_END <= (std::numeric_limits<uint16_t>::max() >> kValueBits),
        "CIGAR opcode exceeds its packed width");

    struct PackedTrace
    {
        size_t cigarCount, operationCount;
        std::array<uint8_t, kCigars> cigarLengths;
        std::array<uint16_t, kOperations> ops;
    };
    static_assert(sizeof(PackedTrace::ops) == kOriginalOperations * sizeof(uint32_t),
        "Packed trace operation storage changed size");
    struct Donor
    {
        Key key{};
        size_t queryLength, prefixLength;
        int score;
        std::array<char, kSequenceBytes> sequence;
        PackedTrace trace;
    };
    struct Bucket
    {
        std::mutex mutex;
        std::array<Donor, kWays> donors;
        size_t next = 0;
    };
    static_assert(kBuckets * sizeof(Bucket) <= 64 * 1024 * 1024, "Cache exceeds its memory budget");

    static bool fitsScoreWidth(size_t query, size_t target)
    {
        return query && target && query < 32768 && target < 32768 && 15 * query + 10 * target + 10 < 32768;
    }

    static size_t bucketIndex(const Key& key, size_t queryLength)
    {
        const uint64_t part = 8 * queryLength + 4 * key.mateSlot + 2 * key.reversedRead + key.reversedExtension;
        const uint64_t mixed = (key.pair * 0x9e3779b97f4a7c15ULL) ^ (part * 0xbf58476d1ce4e5b9ULL);
        return (mixed >> 32) & (kBuckets - 1);
    }

    static bool pack(const std::vector<Cigar>& cigars, PackedTrace& packed)
    {
        if (cigars.size() > kCigars)
            return false;
        packed.cigarCount = cigars.size();
        packed.operationCount = 0;
        for (size_t i = 0; i < cigars.size(); ++i)
        {
            const auto& cigar = cigars[i];
            const size_t length = cigar.length();
            if (length > std::numeric_limits<uint8_t>::max() || length > kOperations - packed.operationCount)
                return false;
            packed.cigarLengths[i] = static_cast<uint8_t>(length);
            for (const auto& op : cigar)
            {
                // All legal CIGAR opcodes are contiguous from ALIGN to NODE_END.
                const auto code = static_cast<unsigned>(op.code_);
                if (code > Cigar::NODE_END || op.value_ > kValueMask)
                    return false;
                packed.ops[packed.operationCount++] = static_cast<uint16_t>((code << kValueBits) | op.value_);
            }
        }
        return true;
    }

    static void copyTrace(const PackedTrace& from, PackedTrace& to)
    {
        to.cigarCount = from.cigarCount;
        to.operationCount = from.operationCount;
        std::copy_n(from.cigarLengths.begin(), from.cigarCount, to.cigarLengths.begin());
        std::copy_n(from.ops.begin(), from.operationCount, to.ops.begin());
    }

    bool snapshot(const Key& key, const std::string& query, const std::string& target, Counters& counters,
        PackedTrace& trace, int& score)
    {
        auto& bucket = buckets_[bucketIndex(key, query.size())];
        const std::unique_lock<std::mutex> lock(bucket.mutex, std::try_to_lock);
        if (!lock.owns_lock())
        {
            ++counters.busyLookups;
            return false;
        }
        for (const auto& donor : bucket.donors)
        {
            if (donor.key.pair == key.pair && donor.key.locus != key.locus && donor.key.mateSlot == key.mateSlot
                && donor.key.reversedRead == key.reversedRead && donor.key.reversedExtension == key.reversedExtension
                && donor.queryLength == query.size() && donor.prefixLength <= target.size()
                && std::memcmp(donor.sequence.data(), query.data(), donor.queryLength) == 0
                && std::memcmp(donor.sequence.data() + donor.queryLength, target.data(), donor.prefixLength) == 0)
            {
                copyTrace(donor.trace, trace);
                score = donor.score;
                return true;
            }
        }
        return false;
    }

    std::unique_ptr<Bucket[]> buckets_{new Bucket[kBuckets]};
};

class Scope;

inline Scope*& currentScope()
{
    static thread_local Scope* scope = nullptr;
    return scope;
}

// Scoped to one target item backed by a shared input. Non-streaming callers
// and single-locus inputs have no active scope.
class Scope
{
public:
    Scope(Table& table, Counters& counters, uint64_t pair, size_t locus, bool mateSlot)
        : table(table), counters(counters), key{pair, locus, mateSlot, false, false}, previous_(currentScope())
    {
        currentScope() = this;
    }

    ~Scope() { currentScope() = previous_; }
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;

    Table& table;
    Counters& counters;
    Table::Key key;

private:
    Scope* previous_;
};
}
}
