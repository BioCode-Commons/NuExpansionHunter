//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Roman Petrovski <RPetrovski@illumina.com>
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
// Modified by NuExpansionHunter contributors, 2026.
// This file differs from the GraphTools copy shipped with ExpansionHunter v5.0.0: seed searches use rolling encoded k-mers, borrowed path views, and optional locale snapshots.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "graphalign/GappedAligner.hh"

#include <algorithm>
#include <limits>
#include <locale>
#include <stdexcept>
#include <vector>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/VerifyPolicy.hh"

using boost::optional;
using std::list;
using std::logic_error;
using std::make_pair;
using std::string;
using std::to_string;

namespace graphtools
{

/**
 * Trim prefix of a path if it is close to node edge
 *
 * @param requested_trim_len: length of the prefix to attempt to trim by
 * @param min_path_len: minimal length of the trimmed path
 * @param[in|out] path: any path
 * @return: actual trim length
 */
static int32_t trimPrefixNearNodeEdge(int32_t requested_trim_len, int32_t min_path_len, Path& path)
{
    if (path.numNodes() == 1 || static_cast<int32_t>(path.length()) <= min_path_len)
    {
        return 0;
    }

    const int32_t overlap_len_with_first_node = path.getNodeOverlapLengthByIndex(0);
    if (overlap_len_with_first_node > requested_trim_len)
    {
        return 0;
    }

    const bool can_be_fully_trimmed = static_cast<int32_t>(path.length()) >= requested_trim_len + min_path_len;
    const int32_t actual_trim_len = can_be_fully_trimmed ? requested_trim_len : path.length() - min_path_len;

    path.shrinkStartBy(actual_trim_len);
    return actual_trim_len;
}

/**
 * Trim suffix of a path if it is close to node edge
 *
 * @param requested_trim_len: length of the suffix to attempt to trim by
 * @param min_path_len: minimal length of the shrank path
 * @param[in|out] path: any path
 * @return: the length by which the path's end was trimmed
 */
static int32_t trimSuffixNearNodeEdge(int32_t requested_trim_len, int32_t min_path_len, Path& path)
{
    if (path.numNodes() == 1 || static_cast<int32_t>(path.length()) <= min_path_len)
    {
        return 0;
    }

    const int32_t overlap_len_with_last_node = path.getNodeOverlapLengthByIndex(path.numNodes() - 1);
    if (overlap_len_with_last_node > requested_trim_len)
    {
        return 0;
    }

    const bool can_be_fully_trimmed = static_cast<int32_t>(path.length()) >= requested_trim_len + min_path_len;
    const int32_t actual_trim_len = can_be_fully_trimmed ? requested_trim_len : path.length() - min_path_len;

    path.shrinkEndBy(actual_trim_len);
    return actual_trim_len;
}

// Return the exhaustive alignments in comparator order.
list<GraphAlignment> GappedGraphAligner::align(const string& query, AlignerSelector& alignerSelector) const
{
    return alignImpl(query, alignerSelector, nullptr);
}

// Contract [K-39]: return the exhaustive list or a certified singleton with the same canonical result.
std::pair<list<GraphAlignment>, bool>
GappedGraphAligner::alignForCanonicalSelection(const string& query, AlignerSelector& alignerSelector) const
{
    bool selectedProduct = false;
    auto alignments = alignImpl(query, alignerSelector, &selectedProduct);
    return { std::move(alignments), selectedProduct };
}

// Find and trim a seed, then select the requested exhaustive or canonical extension.
list<GraphAlignment>
GappedGraphAligner::alignImpl(const string& query, AlignerSelector& alignerSelector, bool* selectedProduct) const
{
    try
    {
        const auto* seedFacet = alignerSelector.seedSearchFacet();
        optional<AlignmentSeed> optional_seed
            = seedFacet ? searchForAlignmentSeed(query, seedFacet) : searchForAlignmentSeed(query);

        if (optional_seed)
        {
            Path& seed_path = optional_seed->path;
            int seed_start_on_query = optional_seed->start_on_query;

            const int kMinPathLength = 2;
            trimSuffixNearNodeEdge(seed_affix_trim_len_, kMinPathLength, seed_path);
            const int trimmed_prefix_len = trimPrefixNearNodeEdge(seed_affix_trim_len_, kMinPathLength, seed_path);
            return extendSeedToFullAlignmentsImpl(
                std::move(seed_path), query, seed_start_on_query + trimmed_prefix_len, alignerSelector, selectedProduct);
        }
        else
        {
            return {};
        }
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to align " + query + ": " + e.what());
    }
}

// Seed search advances to the right end of an extended path. That end cannot
// precede the original kmer's end: extending left changes its start and length
// by opposite amounts. A cursor can therefore encode visited windows in one
// forward pass, while preserving the search's original skip and tie rules.
class SeedKmerCursor
{
public:
    SeedKmerCursor(const string& query, size_t kmer_length, size_t encoded_length)
        : query_(query)
        , kmer_length_(kmer_length)
        , encoded_length_(encoded_length)
        , mask_(static_cast<uint32_t>((uint64_t(1) << (2 * encoded_length)) - 1))
    {
    }

    bool seek(size_t position)
    {
        // Preserve numPaths' original length validation before examining bases.
        // The index constructor, unlike the size_t search setting, has validated
        // the encoder width. Use that width for the mask even when no window fits.
        if (kmer_length_ != encoded_length_)
        {
            return false;
        }
        // Skipped characters precede the next window and cannot affect its key.
        scanned_ = std::max(scanned_, position);
        while (scanned_ < position + kmer_length_)
        {
            const char base = query_[scanned_++];
            uint32_t bits = 0;
            switch (base)
            {
            case 'A': break;
            case 'C': bits = 1; break;
            case 'G': bits = 2; break;
            case 'T': bits = 3; break;
            default:
                valid_suffix_ = 0;
                key_ = (key_ << 2) & mask_;
                continue;
            }
            key_ = ((key_ << 2) | bits) & mask_;
            valid_suffix_ = std::min(kmer_length_, valid_suffix_ + 1);
        }
        return valid_suffix_ == kmer_length_;
    }

    uint32_t key() const { return key_; }

private:
    const string& query_;
    const size_t kmer_length_;
    const size_t encoded_length_;
    const uint32_t mask_;
    size_t scanned_ = 0;
    size_t valid_suffix_ = 0;
    uint32_t key_ = 0;
};

static void uppercaseSeedQuery(string& query, const std::ctype<char>& facet)
{
    // The classic facet applies the same translation to a range and to each
    // individual byte. Custom facets may override those operations separately.
    static const auto* classic = &std::use_facet<std::ctype<char>>(std::locale::classic());
    if (&facet == classic)
        facet.toupper(&query[0], &query[0] + query.size());
    else
        for (char& base : query)
            base = facet.toupper(base);
}

optional<GappedGraphAligner::AlignmentSeed> GappedGraphAligner::searchForAlignmentSeed(const string& query) const
{
    return searchForAlignmentSeed(query, nullptr);
}

optional<GappedGraphAligner::AlignmentSeed>
// Contract [K-16]
// Requires: the global locale does not change after the selector captured it.
// Ensures:  for all queries and indexes, the chosen seed (path, query offset) equals v5.0.0's searchForAlignmentSeed.
GappedGraphAligner::searchForAlignmentSeed(const string& query, const std::ctype<char>* snapshotFacet) const
{
    string upperQuery = query;
    if (!upperQuery.empty())
    {
        if (snapshotFacet)
            uppercaseSeedQuery(upperQuery, *snapshotFacet);
        else
        {
            const std::locale locale;
            uppercaseSeedQuery(upperQuery, std::use_facet<std::ctype<char>>(locale));
        }
    }

    optional<GappedGraphAligner::AlignmentSeed> optional_seed;
    size_t seed_length = 0;

    bool found_multipath_kmer = false;
    size_t kmer_start_position = 0;
    SeedKmerCursor unique_cursor(upperQuery, kmer_len_, kmer_index_.kmerLength());
    while (kmer_start_position + kmer_len_ <= upperQuery.length())
    {
        const auto paths = unique_cursor.seek(kmer_start_position)
            ? kmer_index_.findEncodedPaths(unique_cursor.key()) : KmerIndex::PathView();

        // Initiate seed construction from a unique kmer
        auto num_kmer_paths = paths.size();
        if (num_kmer_paths > 1)
        {
            found_multipath_kmer = true;
        }

        if (num_kmer_paths == 1)
        {
            // This call updates kmer_start_position to the start of the extended path
            Path extended_path = extendPathMatching(paths[0], upperQuery, kmer_start_position);
            const size_t extended_length = extended_path.length();

            if (!optional_seed || extended_length > seed_length)
            {
                optional_seed = AlignmentSeed(std::move(extended_path), kmer_start_position);
                seed_length = extended_length;
            }

            kmer_start_position += extended_length;
        }
        else
        {
            ++kmer_start_position;
        }
    }

    if (optional_seed || !found_multipath_kmer)
    {
        return optional_seed;
    }

    // If the search for unique kmer failed, consider kmers that correspond to multiple paths
    const int kMaxPathCount = 10;
    kmer_start_position = 0;
    SeedKmerCursor multipath_cursor(upperQuery, kmer_len_, kmer_index_.kmerLength());
    while (kmer_start_position + kmer_len_ <= upperQuery.length())
    {
        const auto paths = multipath_cursor.seek(kmer_start_position)
            ? kmer_index_.findEncodedPaths(multipath_cursor.key()) : KmerIndex::PathView();

        const int numPaths = paths.size();
        if (0 < numPaths && numPaths <= kMaxPathCount)
        {
            size_t longest_kmer_path_extension = 0;
            size_t kmer_start_position_for_longest_extension = 0;
            for (size_t path_index = 0; path_index < paths.size(); ++path_index)
            {
                size_t kmer_start_position_for_kmer_path = kmer_start_position;
                Path extended_path = extendPathMatching(paths[path_index], upperQuery, kmer_start_position_for_kmer_path);
                const size_t extended_length = extended_path.length();

                if (longest_kmer_path_extension < extended_length)
                {
                    longest_kmer_path_extension = extended_length;
                    kmer_start_position_for_longest_extension = kmer_start_position_for_kmer_path;
                }

                if (!optional_seed || extended_length > seed_length)
                {
                    optional_seed = AlignmentSeed(std::move(extended_path), kmer_start_position_for_kmer_path);
                    seed_length = extended_length;
                }
            }

            kmer_start_position = kmer_start_position_for_longest_extension + longest_kmer_path_extension;
        }
        else
        {
            ++kmer_start_position;
        }
    }

    return optional_seed;
}

namespace
{
using Extensions = list<PathAndAlignment>;

// Accept only the shared path and exact consumptions with a valid positive-seed seam.
bool validExtension(const PathAndAlignment& extension, const Path& commonPath, size_t expectedReferenceLength,
    size_t expectedQueryLength, bool isPrefix)
{
    const auto& alignment = extension.second;
    if (!(extension.first == commonPath) || alignment.operations().empty()
        || (isPrefix && alignment.referenceStart() != 0))
        return false;
    const auto& seam = isPrefix ? alignment.back() : alignment.front();
    if (seam.type() == OperationType::kSoftclip && alignment.numOperations() != 1)
        return false;

    // Check the true lengths before using the uint32_t alignment counts and seam sums.
    uint64_t referenceLength = 0;
    uint64_t queryLength = 0;
    for (const auto& op : alignment)
    {
        const auto reference = op.referenceLength();
        const auto query = op.queryLength();
        if (reference == 0 && query == 0)
            return false;
        referenceLength += reference;
        queryLength += query;
        if (referenceLength > expectedReferenceLength || queryLength > expectedQueryLength)
            return false;
    }
    return referenceLength == expectedReferenceLength && queryLength == expectedQueryLength;
}

struct ProductScratch
{
    std::vector<uint32_t> nodeEnds;
    std::vector<uint64_t> key;
    std::vector<uint64_t> bestKey;
};

struct CanonicalProduct
{
    Path path;
    const Alignment& prefix;
    Alignment& suffix;
};

class ProjectedOperationKey
{
public:
    // Start an empty comparison key at a position and node in the common path.
    ProjectedOperationKey(ProductScratch& scratch, uint32_t position, size_t node)
        : key_(scratch.key)
        , nextNodeEnd_(scratch.nodeEnds.data() + node)
        , position_(position)
    {
        key_.clear();
    }

    // Encode the projected pieces in comparator order, leaving zero-reference operations on the active node.
    void append(const Operation& op)
    {
        if (op.referenceLength() == 0)
        {
            appendPiece(op.type(), op.length());
            return;
        }

        uint32_t remaining = op.length();
        while (position_ + remaining > *nextNodeEnd_)
        {
            const uint32_t pieceLength = *nextNodeEnd_ - position_;
            if (pieceLength != 0)
                appendPiece(op.type(), pieceLength);
            key_.push_back(0);
            remaining -= pieceLength;
            position_ = *nextNodeEnd_++;
        }
        appendPiece(op.type(), remaining);
        position_ += remaining;
    }

    // Merge a leading match into the seed match before encoding the suffix operations.
    void appendSuffix(uint32_t matchLength, const Alignment& suffix)
    {
        auto op = suffix.begin();
        if (op->type() == OperationType::kMatch)
            matchLength += (op++)->length();
        append(Operation(OperationType::kMatch, matchLength));
        for (; op != suffix.end(); ++op)
            append(*op);
    }

private:
    // Order node boundaries before operations, and operations by type and length.
    void appendPiece(OperationType type, uint32_t length)
    {
        // Alignment comparisons order by operation type, then by length. A
        // node end sorts before any operation, as in a lexicographic list comparison.
        key_.push_back(((static_cast<uint64_t>(type) + 1) << 32) | length);
    }

    std::vector<uint64_t>& key_;
    const uint32_t* nextNodeEnd_;
    uint32_t position_;
};

// Contract [K-39]: select the first sorted exhaustive product when all complete products share a path.
optional<CanonicalProduct> selectCanonicalProduct(const Extensions& prefixes, const Path& seed, size_t seedLength,
    Extensions& suffixes, size_t queryPrefixLength, size_t queryLength, size_t querySuffixLength)
{
    if (prefixes.empty() || suffixes.empty() || (prefixes.size() == 1 && suffixes.size() == 1)
        || seedLength == 0 || queryLength > std::numeric_limits<int32_t>::max()
        || queryPrefixLength > queryLength || seedLength > queryLength - queryPrefixLength
        || querySuffixLength != queryLength - queryPrefixLength - seedLength)
        return {};

    const Path& prefixPath = prefixes.front().first;
    const Path& suffixPath = suffixes.front().first;
    const size_t prefixLength = prefixPath.length();
    const size_t suffixLength = suffixPath.length();
    constexpr size_t maxLength = std::numeric_limits<int32_t>::max();
    if (prefixLength > maxLength || seedLength > maxLength || suffixLength > maxLength
        || static_cast<uint64_t>(prefixLength) + seedLength + suffixLength > maxLength
        || prefixPath.graphRawPtr() != seed.graphRawPtr() || suffixPath.graphRawPtr() != seed.graphRawPtr()
        || !checkIfBookended(prefixPath, seed) || !checkIfBookended(seed, suffixPath))
        return {};
    for (const auto& prefix : prefixes)
        if (!validExtension(prefix, prefixPath, prefixLength, queryPrefixLength, true))
            return {};
    for (const auto& suffix : suffixes)
        if (!validExtension(suffix, suffixPath, suffixLength, querySuffixLength, false))
            return {};

    Path fullPath = concatenatePaths(concatenatePaths(prefixPath, seed), suffixPath);
    thread_local ProductScratch scratch;
    scratch.nodeEnds.clear();
    uint32_t position = 0;
    for (size_t node = 0; node < fullPath.numNodes(); ++node)
    {
        const size_t length = fullPath.getNodeOverlapLengthByIndex(node);
        if (length == 0 || length > maxLength - position
            || fullPath.graphRawPtr()->nodeSeqLength(fullPath.getNodeIdByIndex(node)) > maxLength)
            return {};
        position += length;
        scratch.nodeEnds.push_back(position);
    }
    if (position != prefixLength + seedLength + suffixLength)
        return {};

    // The positive seed makes the second merge seam a match for every prefix.
    // Start a suffix key at the final seed base. A real prefix can only add a
    // constant to that first match piece, shared by every suffix candidate.
    // Later pieces and the ownership of zero-reference operations are unchanged.
    auto bestSuffix = suffixes.begin();
    if (suffixes.size() > 1)
    {
        const uint32_t finalSeedBase = prefixLength + seedLength - 1;
        const size_t seedNode = std::upper_bound(scratch.nodeEnds.begin(), scratch.nodeEnds.end(), finalSeedBase)
            - scratch.nodeEnds.begin();
        for (auto it = suffixes.begin(); it != suffixes.end(); ++it)
        {
            ProjectedOperationKey key(scratch, finalSeedBase, seedNode);
            key.appendSuffix(1, it->second);
            if (it == suffixes.begin() || scratch.key < scratch.bestKey)
            {
                bestSuffix = it;
                scratch.key.swap(scratch.bestKey);
            }
        }
    }

    auto bestPrefix = prefixes.begin();
    if (prefixes.size() > 1)
    {
        for (auto it = prefixes.begin(); it != prefixes.end(); ++it)
        {
            const Alignment& prefix = it->second;
            auto end = prefix.end();
            uint32_t matchLength = seedLength;
            if (prefix.back().type() == OperationType::kMatch)
                matchLength += (--end)->length();

            ProjectedOperationKey key(scratch, 0, 0);
            for (auto op = prefix.begin(); op != end; ++op)
                key.append(*op);
            key.appendSuffix(matchLength, bestSuffix->second);
            if (it == prefixes.begin() || scratch.key < scratch.bestKey)
            {
                bestPrefix = it;
                scratch.key.swap(scratch.bestKey);
            }
        }
    }

    return CanonicalProduct{ std::move(fullPath), bestPrefix->second, bestSuffix->second };
}

// Append a positive seed match to the prefix using the exhaustive seam rules.
Alignment mergeWithSeed(const Alignment& prefix, size_t seedLength)
{
    // Preserve the CIGAR parser's range error for unusually long paths.
    Alignment seed = seedLength <= static_cast<size_t>(std::numeric_limits<int>::max())
        ? Alignment(prefix.referenceLength(), { Operation(OperationType::kMatch, seedLength) })
        : Alignment(prefix.referenceLength(), to_string(seedLength) + "M");
    return mergeAlignments(prefix, seed);
}

// Merge at the given suffix offset and project the result onto the supplied graph path.
GraphAlignment mergeAndProject(const Alignment& prefixWithSeed, int32_t suffixStart, Alignment& suffix, Path path)
{
    suffix.setReferenceStart(suffixStart);
    Alignment fullAlignment = mergeAlignments(prefixWithSeed, suffix);
    return projectAlignmentOntoGraph(std::move(fullAlignment), std::move(path));
}
}

// Contract [K-17]
// Requires (fast path): the seed covers the whole query on one node.
// Ensures:  the single full-length match equals the general prefix/suffix extension below for the same seed.
// Reference: the general extension (EH_VERIFY_FAST_PATHS=OFF).
list<GraphAlignment> GappedGraphAligner::extendSeedToFullAlignments(
    Path seed_path, const string& query, size_t seed_start_on_query, AlignerSelector& alignerSelector) const
{
    return extendSeedToFullAlignmentsImpl(std::move(seed_path), query, seed_start_on_query, alignerSelector, nullptr);
}

// Build the ordered exhaustive product, or return a certified canonical singleton when requested.
list<GraphAlignment> GappedGraphAligner::extendSeedToFullAlignmentsImpl(
    Path seed_path, const string& query, size_t seed_start_on_query, AlignerSelector& alignerSelector,
    bool* selectedProduct) const
{
    assert(seed_path.length() > 1);

    // A full-query seed inside one node has no unaligned ends. The two
    // artificial one-base matches below would merge back into this same CIGAR.
    if (verify::kFastPaths && seed_start_on_query == 0 && seed_path.numNodes() == 1 && query.size() >= 2
        && query.size() <= size_t(std::numeric_limits<int32_t>::max()) && seed_path.length() == query.size())
    {
        std::vector<Alignment> projected;
        projected.emplace_back(seed_path.startPosition(),
            list<Operation>{Operation(OperationType::kMatch, static_cast<uint32_t>(query.size()))});
        list<GraphAlignment> matches;
        matches.push_back(GraphAlignment::fromOwnedPath(std::move(seed_path), std::move(projected)));
        return matches;
    }

    // Generate prefix extensions
    list<PathAndAlignment> prefix_extensions;
    size_t query_prefix_len = seed_start_on_query;
    if (query_prefix_len != 0)
    {
        const string query_prefix = query.substr(0, query_prefix_len);
        Path prefix_seed_path = seed_path;
        prefix_seed_path.shrinkEndBy(seed_path.length());
        prefix_extensions
            = extendAlignmentPrefix(prefix_seed_path, query_prefix, query_prefix_len + padding_len_, alignerSelector);
    }
    else
    {
        // Because (a) empty alignments are currently disallowed and (b) we don't want to deal with an empty list of
        // prefix_extensions we create a 1bp prefix artificially.
        query_prefix_len = 1;
        Path prefix_path = seed_path;
        prefix_path.shrinkEndBy(prefix_path.length() - 1);
        prefix_extensions.emplace_back(
            std::move(prefix_path), Alignment(0, { Operation(OperationType::kMatch, 1) }));
        seed_path.shrinkStartBy(1);
    }

    // Generate suffix extensions
    list<PathAndAlignment> suffix_extensions;
    size_t query_suffix_len = query.length() - seed_path.length() - query_prefix_len;
    if (query_suffix_len != 0)
    {
        const string query_suffix = query.substr(query_prefix_len + seed_path.length(), query_suffix_len);
        Path suffix_seed_path = seed_path;
        suffix_seed_path.shrinkStartBy(seed_path.length());
        suffix_extensions
            = extendAlignmentSuffix(suffix_seed_path, query_suffix, query_suffix_len + padding_len_, alignerSelector);
    }
    else
    {
        // Because (a) empty alignments are currently disallowed and (b) we don't want to deal with an empty list of
        // suffix_extensions we create a 1bp suffix artificially.
        Path suffix_path = seed_path;
        suffix_path.shrinkStartBy(suffix_path.length() - 1);
        suffix_extensions.emplace_back(
            std::move(suffix_path), Alignment(0, { Operation(OperationType::kMatch, 1) }));
        seed_path.shrinkEndBy(1);
        query_suffix_len = 1;
    }

    // Merge alignments together
    const size_t seed_length = seed_path.length();
    if (verify::kFastPaths && selectedProduct)
    {
        auto product = selectCanonicalProduct(
            prefix_extensions, seed_path, seed_length, suffix_extensions, query_prefix_len, query.size(),
            query_suffix_len);
        if (product)
        {
            Alignment prefixWithSeed = mergeWithSeed(product->prefix, seed_length);
            list<GraphAlignment> selected;
            selected.push_back(mergeAndProject(prefixWithSeed, prefixWithSeed.referenceLength(), product->suffix,
                std::move(product->path)));
            *selectedProduct = true;
            return selected;
        }
    }
    list<GraphAlignment> top_graph_alignments;
    for (PathAndAlignment& prefix_path_and_alignment : prefix_extensions)
    {
        Path& prefix_path = prefix_path_and_alignment.first;
        Path prefix_plus_seed_path = concatenatePaths(prefix_path, seed_path);
        const int32_t prefix_plus_seed_length = prefix_plus_seed_path.length();

        Alignment& prefix_alignment = prefix_path_and_alignment.second;
        Alignment prefix_plus_kmer_alignment = mergeWithSeed(prefix_alignment, seed_length);

        for (PathAndAlignment& suffix_path_and_alignment : suffix_extensions)
        {
            Path& suffix_path = suffix_path_and_alignment.first;
            Alignment& suffix_alignment = suffix_path_and_alignment.second;
            Path full_path = concatenatePaths(prefix_plus_seed_path, suffix_path);

            top_graph_alignments.push_back(mergeAndProject(prefix_plus_kmer_alignment, prefix_plus_seed_length,
                suffix_alignment, std::move(full_path)));
        }
    }

    top_graph_alignments.sort();
    top_graph_alignments.unique();
    return top_graph_alignments;
}

list<PathAndAlignment> GappedGraphAligner::extendAlignmentPrefix(
    const Path& seed_path, const string& query_piece, size_t extension_len, AlignerSelector& alignerSelector) const
{
    assert(seed_path.length() == 0);

    int32_t top_alignment_score = INT32_MIN;
    list<PathAndAlignment> top_paths_and_alignments
        = alignerSelector.suffixAlign(seed_path, query_piece, extension_len, top_alignment_score);

    for (PathAndAlignment& path_and_alignment : top_paths_and_alignments)
    {
        Path& path = path_and_alignment.first;
        Alignment& alignment = path_and_alignment.second;
        alignment.setReferenceStart(0);

        const int32_t overhang = path.length() - alignment.referenceLength();
        path.shrinkStartBy(overhang);

        if (!checkConsistency(path_and_alignment.second, path.seq(), query_piece))
        {
            throw std::logic_error("Inconsistent prefix");
        }
    }

    return top_paths_and_alignments;
}

list<PathAndAlignment> GappedGraphAligner::extendAlignmentSuffix(
    const Path& seed_path, const string& query_piece, size_t extension_len, AlignerSelector& alignerSelector) const
{
    assert(seed_path.length() == 0);

    int32_t top_alignment_score = INT32_MIN;
    list<PathAndAlignment> top_paths_and_alignments
        = alignerSelector.prefixAlign(seed_path, query_piece, extension_len, top_alignment_score);

    for (PathAndAlignment& path_and_alignment : top_paths_and_alignments)
    {
        if (!checkConsistency(path_and_alignment.second, path_and_alignment.first.seq(), query_piece))
        {
            throw std::logic_error("Inconsistent suffix");
        }

        Path& path = path_and_alignment.first;
        Alignment& alignment = path_and_alignment.second;

        const int32_t overhang = path.length() - alignment.referenceLength();
        path.shrinkEndBy(overhang);
    }

    return top_paths_and_alignments;
}
}
