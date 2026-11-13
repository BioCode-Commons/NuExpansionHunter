//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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

#include <array>
#ifdef GRAPHTOOLS_PREPARED_GRAPH_INSTRUMENT
#include <cstdio>
#endif
#include <list>
#include <map>
#include <string>

#include <boost/throw_exception.hpp>

#include "graphalign/DagAlignerAffine.hh"
#include "graphalign/dagAligner/PinnedBudgetAlignMatrix.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphalign/Operation.hh"
#include "graphalign/dagAligner/BaseMatchingPenaltyMatrix.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

#ifndef GRAPHTOOLS_PREPARED_GRAPH_CACHE_SLOTS
#define GRAPHTOOLS_PREPARED_GRAPH_CACHE_SLOTS 4
#endif

namespace graphtools
{

using PathAndAlignment = std::pair<Path, Alignment>;

template <bool penalizeMove, bool clipFront = true>
class BaseMatchingDagAligner : public graphalign::dagAligner::Aligner<
                                   graphalign::dagAligner::AffineAlignMatrixVectorized<
                                       graphalign::dagAligner::BaseMatchingPenaltyMatrix, penalizeMove>,
                                   clipFront>
{
    typedef graphalign::dagAligner::BaseMatchingPenaltyMatrix PenaltyMatrix;
    typedef graphalign::dagAligner::Score Score;

public:
    BaseMatchingDagAligner(const PenaltyMatrix& penaltyMatrix, Score gapOpen, Score gapExt)
        : graphalign::dagAligner::Aligner<
              graphalign::dagAligner::AffineAlignMatrixVectorized<PenaltyMatrix, penalizeMove>, clipFront>(
              penaltyMatrix, gapOpen, gapExt)
    {
    }

    BaseMatchingDagAligner(graphalign::dagAligner::Score match, Score mismatch, Score gapOpen, Score gapExt)
        : graphalign::dagAligner::Aligner<
              graphalign::dagAligner::AffineAlignMatrixVectorized<PenaltyMatrix, penalizeMove>, clipFront>(
              PenaltyMatrix(match, mismatch), gapOpen, gapExt)
    {
    }
};

/**
 * Performs alignment of query pieces that start or end at the seed in the graph
 */
class PinnedDagAligner
{
#ifdef GRAPHTOOLS_PREPARED_GRAPH_INSTRUMENT
public:
    struct PreparationStats
    {
        uint64_t calls = 0, hits = 0, misses = 0, oversized = 0, evictions = 0, maxRetainedBytes = 0;
        ~PreparationStats()
        {
            if (calls)
                std::fprintf(stderr, "PREPARED_GRAPH calls=%llu hits=%llu misses=%llu oversized=%llu evictions=%llu max_retained_bytes=%llu\n",
                    static_cast<unsigned long long>(calls), static_cast<unsigned long long>(hits),
                    static_cast<unsigned long long>(misses), static_cast<unsigned long long>(oversized),
                    static_cast<unsigned long long>(evictions), static_cast<unsigned long long>(maxRetainedBytes));
        }
    };
    static PreparationStats& preparationStats()
    {
        static thread_local PreparationStats value;
        return value;
    }
private:
#endif
#ifdef GRAPHTOOLS_UNIQUE_UNGAPPED_INSTRUMENT
public:
    struct UniqueUngappedStats
    {
        uint64_t attempts = 0, searched = 0, hits = 0, forward = 0, reverse = 0, mismatched = 0;
        uint64_t multiNode = 0, emptyBoundary = 0, ambiguous = 0, exhausted = 0;
        uint64_t provedNoCost9 = 0;
    };
    static UniqueUngappedStats& uniqueUngappedStats()
    {
        static thread_local UniqueUngappedStats value;
        return value;
    }
private:
#endif
    typedef std::pair<int, int> Edge;
    typedef std::vector<Edge> Edges;
    typedef graphalign::dagAligner::Cigar Cigar;
    graphalign::dagAligner::Aligner<graphalign::dagAligner::PinnedBudgetAlignMatrix, false> aligner_;
    const bool defaultScores_;
    unsigned adaptiveHighCostStreak_ = 0;

    int preferredInitialBudget(int provenMinimumCost) const
    {
#ifndef GRAPHTOOLS_ADAPTIVE_BUDGET_HINT_DISABLE
        // Only qualifying DP calls train this hint. The initial search width
        // is a performance preference, not an additional cost lower bound.
        return provenMinimumCost > 9 && adaptiveHighCostStreak_ >= 3 ? 36 : 0;
#else
        return 0;
#endif
    }

    void observeAlignmentCost(size_t queryLength, int bestScore, int provenMinimumCost)
    {
#ifndef GRAPHTOOLS_ADAPTIVE_BUDGET_HINT_DISABLE
        if (provenMinimumCost > 9)
        {
            // The proof's default-score and overflow guards make this exact.
            const int cost = 5 * static_cast<int>(queryLength) - bestScore;
            adaptiveHighCostStreak_ = cost > 18 ? std::min(3u, adaptiveHighCostStreak_ + 1) : 0;
        }
#endif
    }

    struct UngappedSearch
    {
        const Graph& graph;
        const std::string& query;
        const bool forward;
        std::array<NodeId, 257> trail, bestTrail;
        size_t bestDepth = 0, bestPosition = 0, bestMismatchPosition = 0;
        unsigned bestMismatches = 2, ties = 0, visits = 0;
        bool exhausted = false;

        UngappedSearch(const Graph& graph, const std::string& query, bool forward)
            : graph(graph), query(query), forward(forward)
        {
        }

        void visit(NodeId node, size_t position, size_t queryPosition, unsigned mismatches,
                   size_t mismatchPosition, size_t depth)
        {
            if (exhausted || (bestMismatches == 0 && ties > 1) || mismatches > bestMismatches)
                return;
            if (++visits > 512 || depth == trail.size())
            {
                exhausted = true;
                return;
            }
            trail[depth++] = node;
            const auto sequence = graph.nodeSeq(node);
            while (queryPosition < query.size() && (forward ? position < sequence.size() : position != 0))
            {
                const size_t q = forward ? queryPosition : query.size() - queryPosition - 1;
                const size_t t = forward ? position++ : --position;
                const auto qc = codes::kQueryBaseEncodingTable[static_cast<unsigned char>(query[q])];
                const auto tc = codes::kReferenceBaseEncodingTable[static_cast<unsigned char>(sequence[t])];
                if (!codes::kReferenceQueryCodeMatchLookupTable[tc][qc])
                {
                    if (++mismatches > bestMismatches || mismatches > 1)
                        return;
                    mismatchPosition = q;
                }
                ++queryPosition;
            }
            if (queryPosition == query.size())
            {
                if (mismatches < bestMismatches)
                {
                    bestMismatches = mismatches;
                    ties = 0;
                }
                if (++ties == 1)
                {
                    std::copy(trail.begin(), trail.begin() + depth, bestTrail.begin());
                    bestDepth = depth;
                    bestPosition = position;
                    bestMismatchPosition = mismatchPosition;
                }
                return;
            }
            for (NodeId next : forward ? graph.successors(node) : graph.predecessors(node))
                visit(next, forward ? 0 : graph.nodeSeqLength(next), queryPosition, mismatches, mismatchPosition, depth);
        }
    };

    bool tryUniqueUngappedAlignment(
        const Path& seedPath, const std::string& queryPiece, size_t extensionLen, bool extendEnd,
        std::list<PathAndAlignment>& result, int& score, int* minimumCost = nullptr) const
    {
#if defined(GRAPHTOOLS_BUDGET_DP_DISABLE) || defined(GRAPHTOOLS_UNIQUE_UNGAPPED_DISABLE)
        return false;
#else
        const size_t queryLength = queryPiece.size();
        if (!defaultScores_ || queryLength < 16 || queryLength > 256 || extensionLen < queryLength
            || extensionLen > 1024)
            return false;
#ifdef GRAPHTOOLS_UNIQUE_UNGAPPED_INSTRUMENT
        ++uniqueUngappedStats().attempts;
#endif
        const Graph& graph = *seedPath.graphRawPtr();
        if (graph.numNodes() > 64)
            return false;
        const NodeId firstNode = extendEnd ? seedPath.lastNodeId() : seedPath.firstNodeId();
        const size_t boundary = extendEnd ? seedPath.endPosition() : seedPath.startPosition();
        const size_t firstLength = graph.nodeSeqLength(firstNode);
        if (boundary > firstLength || firstLength >= 32768)
            return false;
        const size_t startOffset = extendEnd ? boundary : firstLength - boundary;
        const size_t unrollingLength = startOffset + extensionLen;
        size_t targetUpperBound = 0;
        for (NodeId node = 0; node < graph.numNodes(); ++node)
        {
            const size_t length = graph.nodeSeqLength(node);
            if (!length || length >= 32768)
                return false;
            bool repeats = false;
            for (NodeId next : graph.successors(node))
            {
                // Catalog graphs are ordered DAGs apart from self-loops. Other
                // graph topologies retain the original discovery behavior.
                if (next < node)
                    return false;
                repeats = repeats || next == node;
            }
            targetUpperBound += repeats ? (unrollingLength + length - 1) / length * length : length;
            if (15 * queryLength + 10 * targetUpperBound + 10 >= 32768)
                return false;
        }

        // With default penalties, a mismatch costs 9 and every gap costs at
        // least 10. Enumerating all paths of cost <=9 therefore proves the
        // best alignment when it is unique. The complete target-size bound
        // above also preserves legacy int16 arithmetic outside the shortcut.
        UngappedSearch search(graph, queryPiece, extendEnd);
        search.visit(firstNode, boundary, 0, 0, 0, 0);
#ifndef GRAPHTOOLS_PROVEN_COST9_SKIP_DISABLE
        // An exhaustive unsuccessful search proves that no full-query path
        // costs <=9: every such path would be ungapped with <=1 mismatch.
        // Guard/cap rejection and tied or outer-mismatch paths prove nothing.
        if (minimumCost && !search.exhausted && search.ties == 0)
        {
            *minimumCost = 10;
#ifdef GRAPHTOOLS_UNIQUE_UNGAPPED_INSTRUMENT
            ++uniqueUngappedStats().provedNoCost9;
#endif
        }
#endif
#ifdef GRAPHTOOLS_UNIQUE_UNGAPPED_INSTRUMENT
        ++uniqueUngappedStats().searched;
        uniqueUngappedStats().ambiguous += search.ties > 1;
        uniqueUngappedStats().exhausted += search.exhausted;
#endif
        if (search.exhausted || search.ties != 1)
            return false;
        // Leave outer mismatch/softclip and node-boundary details to traceback.
        const size_t mismatch = search.bestMismatches ? search.bestMismatchPosition : queryLength;
        if (search.bestMismatches && (mismatch == 0 || mismatch + 1 == queryLength))
            return false;

        std::vector<NodeId> nodes;
        nodes.reserve(seedPath.numNodes() + search.bestDepth - 1);
        if (extendEnd)
        {
            nodes.insert(nodes.end(), seedPath.nodeIds().begin(), seedPath.nodeIds().end());
            nodes.insert(nodes.end(), search.bestTrail.begin() + 1, search.bestTrail.begin() + search.bestDepth);
        }
        else
        {
            for (size_t i = search.bestDepth; i != 0; --i)
                nodes.push_back(search.bestTrail[i - 1]);
            nodes.insert(nodes.end(), seedPath.nodeIds().begin() + 1, seedPath.nodeIds().end());
        }
        Path path(&graph, extendEnd ? seedPath.startPosition() : search.bestPosition, std::move(nodes),
                  extendEnd ? search.bestPosition : seedPath.endPosition());
        std::list<Operation> operations;
        if (mismatch)
            operations.emplace_back(OperationType::kMatch, mismatch);
        if (search.bestMismatches)
        {
            operations.emplace_back(OperationType::kMismatch, 1);
            operations.emplace_back(OperationType::kMatch, queryLength - mismatch - 1);
        }
        result.emplace_back(std::move(path), Alignment(extendEnd ? seedPath.length() : 0, std::move(operations)));
        score = 5 * queryLength - 9 * search.bestMismatches;
#ifdef GRAPHTOOLS_UNIQUE_UNGAPPED_INSTRUMENT
        auto& stats = uniqueUngappedStats();
        ++stats.hits;
        stats.forward += extendEnd;
        stats.reverse += !extendEnd;
        stats.mismatched += search.bestMismatches;
        stats.multiNode += search.bestDepth > 1;
        stats.emptyBoundary += extendEnd ? boundary == firstLength : boundary == 0;
#endif
        return true;
#endif
    }

    static bool hasSingleGapCompetitor(
        const TextView& nodeSequence, const std::string& query, size_t boundary, bool extendEnd,
        size_t firstMismatch, size_t secondMismatch)
    {
        const size_t length = query.size();
        const size_t prefix = extendEnd ? firstMismatch : length - secondMismatch - 1;
        const auto matches = [&](size_t queryPosition, size_t targetPosition) {
            const size_t q = extendEnd ? queryPosition : length - queryPosition - 1;
            const size_t t = extendEnd ? boundary + targetPosition : boundary - targetPosition - 1;
            const auto qc = codes::kQueryBaseEncodingTable[static_cast<unsigned char>(query[q])];
            const auto tc = codes::kReferenceBaseEncodingTable[static_cast<unsigned char>(nodeSequence[t])];
            return codes::kReferenceQueryCodeMatchLookupTable[tc][qc];
        };
        // A single deletion of length d costs 8+2d. A compatible gap position
        // exists exactly when the matching prefix and right-aligned suffix
        // cover the query. Leading legacy deletions can cost more, so including
        // them here only causes conservative fallback.
        for (size_t deletion = 1; deletion <= 5; ++deletion)
        {
            size_t suffix = 0;
            while (suffix < length && matches(length - suffix - 1, length + deletion - suffix - 1))
                ++suffix;
            if (prefix + suffix >= length)
                return true;
        }
        // One inserted query base costs 15; two cost at least 22.
        size_t suffix = 0;
        while (suffix + 1 < length && matches(length - suffix - 1, length - suffix - 2))
            ++suffix;
        return prefix + suffix + 1 >= length;
    }

    bool trySingleNodeAlignment(
        const Path& seedPath, const std::string& queryPiece, size_t extensionLen, bool extendEnd,
        std::list<PathAndAlignment>& result, int& score) const
    {
#ifdef GRAPHTOOLS_BUDGET_DP_DISABLE
        return false;
#else
        const size_t queryLength = queryPiece.size();
        if (!defaultScores_ || queryLength < 16 || !extensionLen)
            return false;

        const Graph& graph = *seedPath.graphRawPtr();
        const NodeId nodeId = extendEnd ? seedPath.lastNodeId() : seedPath.firstNodeId();
        const size_t nodeLength = graph.nodeSeqLength(nodeId);
        const size_t boundary = extendEnd ? seedPath.endPosition() : seedPath.startPosition();
        if (boundary > nodeLength)
            return false;
        const size_t targetLength = extendEnd ? nodeLength - boundary : boundary;

        // These are exactly the single-node discovery and compact-matrix
        // guards. In particular, the score bound uses the entire remaining
        // node, not just the part consumed by the query.
        if (queryLength > targetLength || targetLength < 8
            || queryLength >= 32768 || targetLength >= 32768
            || 15 * queryLength + 10 * targetLength + 10 >= 32768)
            return false;
        if (extensionLen > targetLength
            && !(extendEnd ? graph.successors(nodeId) : graph.predecessors(nodeId)).empty())
            return false;

        const auto nodeSequence = graph.nodeSeq(nodeId);
        const size_t targetStart = extendEnd ? boundary : boundary - queryLength;
        size_t mismatch = queryLength;
        size_t secondMismatch = queryLength;
        for (size_t i = 0; i != queryLength; ++i)
        {
            // Use the same asymmetric reference/query encodings as the DP.
            // Comparing both sequences forward also handles suffix alignment:
            // reversing both preserves matches, and the returned operations
            // must ultimately be in this original query order.
            const auto q = codes::kQueryBaseEncodingTable[static_cast<unsigned char>(queryPiece[i])];
            const auto t = codes::kReferenceBaseEncodingTable[static_cast<unsigned char>(nodeSequence[targetStart + i])];
            if (!codes::kReferenceQueryCodeMatchLookupTable[t][q])
            {
                if (secondMismatch != queryLength)
                    return false;
                if (mismatch == queryLength)
                    mismatch = i;
                else
                {
#ifdef GRAPHTOOLS_TWO_MISMATCH_DISABLE
                    return false;
#else
                    secondMismatch = i;
#endif
                }
            }
        }

        if (secondMismatch != queryLength)
        {
            // With two strictly internal mismatches the diagonal costs 18.
            // Every competing cost<=18 path is either a perfect single D1..5
            // or I1, or the legacy boundary I1+D1 / D1+I1 of cost 17. Those
            // boundary pairs skip only the first query/target base and leave
            // both internal mismatches, so cannot compete. A gap+mismatch
            // costs at least 19, longer boundary pairs at least 19, and
            // ordinary multiple gaps at least 20.
            // Five spare target bases keep every possible competitor in this
            // same node, including all cost-18 deletion ties.
            if (mismatch == 0 || secondMismatch + 1 == queryLength || targetLength - queryLength < 5
                || hasSingleGapCompetitor(nodeSequence, queryPiece, boundary, extendEnd, mismatch, secondMismatch))
                return false;
        }

        // With at most one mismatch, cost <=9 is below every gap opening;
        // the guarded two-mismatch case above also proves a unique diagonal.
        // Construct its final representation directly instead of building
        // EdgeMap, Cigar and traceback objects.
        // Legacy traceback represents an outer-end mismatch as a softclip,
        // retaining its mismatch score but not consuming that reference base.
        const bool softclip = secondMismatch == queryLength
            && (extendEnd ? mismatch + 1 == queryLength : mismatch == 0);
        std::list<Operation> operations;
        if (mismatch)
            operations.emplace_back(OperationType::kMatch, mismatch);
        if (mismatch != queryLength)
        {
            if (secondMismatch != queryLength)
            {
                if (secondMismatch == mismatch + 1)
                    operations.emplace_back(OperationType::kMismatch, 2);
                else
                {
                    operations.emplace_back(OperationType::kMismatch, 1);
                    operations.emplace_back(OperationType::kMatch, secondMismatch - mismatch - 1);
                    operations.emplace_back(OperationType::kMismatch, 1);
                }
                operations.emplace_back(OperationType::kMatch, queryLength - secondMismatch - 1);
            }
            else
            {
                operations.emplace_back(softclip ? OperationType::kSoftclip : OperationType::kMismatch, 1);
                if (mismatch + 1 != queryLength)
                    operations.emplace_back(OperationType::kMatch, queryLength - mismatch - 1);
            }
        }
        Path path = seedPath;
        const size_t referenceLength = queryLength - (softclip ? 1 : 0);
        if (extendEnd)
            path.shiftEndAlongNode(referenceLength);
        else
            path.shiftStartAlongNode(referenceLength);
        result.emplace_back(
            std::move(path), Alignment(extendEnd ? seedPath.length() : 0, std::move(operations)));
        score = 5 * queryLength - 9 * ((mismatch != queryLength) + (secondMismatch != queryLength));
        return true;
#endif
    }

    static void appendOperation(OperationType type, uint32_t length, std::list<Operation>& operations)
    {
        if (operations.empty() || operations.back().type() != type)
        {
            operations.push_back(Operation(type, length));
        }
        else
        {
            operations.back() = Operation(type, operations.back().length() + length);
        }
    }

    template <typename GraphT, typename PathT>
    static void parseGraphCigar(
        const GraphT& graph, const graphalign::dagAligner::Cigar& cigar, PathT& path, std::list<Operation>& operations)
    {
        using namespace graphalign::dagAligner;
        for (const Cigar::Operation& op : cigar)
        {
            if (Cigar::NODE_START == op.code_)
            {
                if (path.lastNodeId() != op.value_
                    || std::size_t(path.endPosition()) == graph.nodeSeqLength(op.value_))
                {
                    path.extendEndToNode(op.value_);
                }
            }
            else if (Cigar::NODE_END == op.code_)
            {
            }
            else
                switch (op.code_)
                {
                case Cigar::MATCH:
                {
                    appendOperation(OperationType::kMatch, op.value_, operations);
                    path.shiftEndAlongNode(op.value_);
                    break;
                }
                case Cigar::MISMATCH:
                {
                    appendOperation(OperationType::kMismatch, op.value_, operations);
                    path.shiftEndAlongNode(op.value_);
                    break;
                }
                case Cigar::INSERT:
                {
                    appendOperation(OperationType::kInsertionToRef, op.value_, operations);
                    break;
                }
                case Cigar::SOFT_CLIP:
                {
                    appendOperation(OperationType::kSoftclip, op.value_, operations);
                    break;
                }
                case Cigar::DELETE:
                {
                    appendOperation(OperationType::kDeletionFromRef, op.value_, operations);
                    path.shiftEndAlongNode(op.value_);
                    break;
                }
                default:
                {
                    throw std::logic_error(std::string("Unexpected graph cigar operation:") + std::to_string(op.code_));
                }
                }
        }
    }

    typedef int MappedId;

    struct PreparedGraph
    {
        uint64_t revision = 0;
        uint64_t lastUsed = 0;
        NodeId startNode = 0;
        size_t startOffset = 0, extension = 0;
        bool reverse = false;
        std::vector<MappedId> nodeIds;
        Edges edges;
        std::string target;
        std::vector<NodeId> originalIds;
        graphalign::dagAligner::EdgeMap alignerEdges;

        size_t retainedBytes() const
        {
            return sizeof(*this) + nodeIds.capacity() * sizeof(MappedId) + edges.capacity() * sizeof(Edge)
                + target.capacity() + 1 + originalIds.capacity() * sizeof(NodeId)
                + alignerEdges.retainedCapacityBytes();
        }
    };
    static_assert(GRAPHTOOLS_PREPARED_GRAPH_CACHE_SLOTS >= 1 && GRAPHTOOLS_PREPARED_GRAPH_CACHE_SLOTS <= 4,
                  "Prepared graph cache supports one to four bounded slots");
    std::array<PreparedGraph, GRAPHTOOLS_PREPARED_GRAPH_CACHE_SLOTS> preparedGraphs_;
    uint64_t preparationClock_ = 0;

    const PreparedGraph& prepareGraph(
        const Graph& graph, NodeId startNode, size_t startOffset, size_t extension, bool reverse,
        PreparedGraph& temporary)
    {
#ifdef GRAPHTOOLS_PREPARED_GRAPH_INSTRUMENT
        auto& stats = preparationStats();
        ++stats.calls;
#endif
        const uint64_t revision = graph.alignmentRevision();
#ifndef GRAPHTOOLS_PREPARED_GRAPH_CACHE_DISABLE
        if (++preparationClock_ == 0)
        {
            // Clock rollover affects replacement order only, never identity.
            for (auto& entry : preparedGraphs_)
                entry.lastUsed = 0;
            preparationClock_ = 1;
        }
        for (size_t slot = 0; slot < preparedGraphs_.size(); ++slot)
        {
            auto& entry = preparedGraphs_[slot];
            if (entry.revision == revision && entry.startNode == startNode && entry.startOffset == startOffset
                && entry.extension == extension && entry.reverse == reverse)
            {
                entry.lastUsed = preparationClock_;
#ifdef GRAPHTOOLS_PREPARED_GRAPH_INSTRUMENT
                ++stats.hits;
#endif
                return entry;
            }
        }
        PreparedGraph& entry = *std::min_element(preparedGraphs_.begin(), preparedGraphs_.end(),
            [](const PreparedGraph& a, const PreparedGraph& b) { return a.lastUsed < b.lastUsed; });
#ifdef GRAPHTOOLS_PREPARED_GRAPH_INSTRUMENT
        stats.evictions += entry.revision != 0;
#endif
#else
        PreparedGraph& entry = temporary;
#endif
#ifdef GRAPHTOOLS_PREPARED_GRAPH_INSTRUMENT
        ++stats.misses;
#endif
        entry.revision = 0;
        entry.nodeIds.clear();
        entry.edges.clear();
        entry.target.clear();
        entry.originalIds.clear();
        try
        {
            if (reverse)
                bfsDiscoverEdges(ReverseGraph(graph), startNode, startOffset, extension,
                                 entry.nodeIds, entry.edges, entry.target, entry.originalIds);
            else
                bfsDiscoverEdges(graph, startNode, startOffset, extension,
                                 entry.nodeIds, entry.edges, entry.target, entry.originalIds);
            entry.edges.emplace_back(entry.target.length(), entry.target.length());
            if (!entry.target.empty())
                entry.alignerEdges.reset(entry.edges, entry.nodeIds);
        }
        catch (...)
        {
            // Failed preparation must not keep an oversized allocation alive.
            entry = PreparedGraph();
            throw;
        }
        entry.revision = revision;
        entry.startNode = startNode;
        entry.startOffset = startOffset;
        entry.extension = extension;
        entry.reverse = reverse;
        entry.lastUsed = preparationClock_;
#ifndef GRAPHTOOLS_PREPARED_GRAPH_CACHE_DISABLE
        // Count retained capacity, not size. Oversized results remain available
        // for this alignment, but their buffers die with the local temporary.
        if (entry.retainedBytes() > 64 * 1024)
        {
#ifdef GRAPHTOOLS_PREPARED_GRAPH_INSTRUMENT
            ++stats.oversized;
#endif
            temporary = std::move(entry);
            entry = PreparedGraph();
            return temporary;
        }
#ifdef GRAPHTOOLS_PREPARED_GRAPH_INSTRUMENT
        size_t retained = 0;
        for (const auto& cached : preparedGraphs_)
            retained += cached.retainedBytes();
        stats.maxRetainedBytes = std::max<uint64_t>(stats.maxRetainedBytes, retained);
#endif
#endif
        return entry;
    }

    void unmapNodeIds(const std::vector<NodeId>& originalIds, graphalign::dagAligner::Cigar& cigar)
    {
        using namespace graphalign::dagAligner;
        for (Cigar::Operation& op : cigar)
        {
            if (Cigar::NODE_START == op.code_ || Cigar::NODE_END == op.code_)
            {
                op.value_ = originalIds.at(op.value_);
            }
        }
    }

    /**
     * \brief Tests if the cigar first node is a repeat expansion and corrects
     *        cigar to ensure that fraction of the first expansion is interpreted correctly
     */
    template <typename PathT>
    void fixFirstNodeExpansion(
        const std::vector<MappedId>& nodeIds, const std::vector<NodeId>& originalIds, const PathT& seedPath,
        graphalign::dagAligner::Cigar& cigar)
    {
        using namespace graphalign::dagAligner;

        if (seedPath.lastNodeId() != originalIds.at(cigar.firstNode()) || MappedId(cigar.firstNode()) != nodeIds[0])
        {
            const int skipLen
                = seedPath.graphRawPtr()->nodeSeqLength(seedPath.lastNodeId()) - seedPath.endPosition();
            if (0 > skipLen)
            {
                throw std::logic_error("invalid distance to skip");
            }
            if (skipLen)
            {
                const Cigar::Operation skipFirstNode[]
                    = { Cigar::Operation(
                            Cigar::NODE_START, nodeIds[0]), // endPosition is on the base that belongs to the path...
                        Cigar::Operation(Cigar::DELETE, skipLen), Cigar::Operation(Cigar::NODE_END, nodeIds[0]) };

                cigar.insert(
                    cigar.begin(), skipFirstNode, skipFirstNode + sizeof(skipFirstNode) / sizeof(skipFirstNode[0]));
            }
        }
    }

public:
    explicit PinnedDagAligner(
        const int32_t matchScore, const int32_t mismatchScore, const int32_t gapOpenScore, const int32_t gapExtendScore)
        : aligner_(graphalign::dagAligner::BaseMatchingPenaltyMatrix(matchScore, mismatchScore), gapOpenScore, gapExtendScore)
        , defaultScores_(
              graphalign::dagAligner::Score(matchScore) == 5 && graphalign::dagAligner::Score(mismatchScore) == -4
              && graphalign::dagAligner::Score(gapOpenScore) == -8 && graphalign::dagAligner::Score(gapExtendScore) == -2)
    {
    }

    std::list<PathAndAlignment>
    prefixAlign(const Path& seedPath, const std::string& queryPiece, size_t extensionLen, int& score)
    {
        using namespace graphalign::dagAligner;
        std::list<PathAndAlignment> ret;
        if (trySingleNodeAlignment(seedPath, queryPiece, extensionLen, true, ret, score))
            return ret;
        int minimumCost = 0;
        if (tryUniqueUngappedAlignment(seedPath, queryPiece, extensionLen, true, ret, score, &minimumCost))
            return ret;
        PreparedGraph temporary;
        const auto& prepared = prepareGraph(*seedPath.graphRawPtr(), seedPath.lastNodeId(),
                                            seedPath.endPosition(), extensionLen, false, temporary);
        const auto& target = prepared.target;
        const auto& nodeIds = prepared.nodeIds;
        const auto& originalIds = prepared.originalIds;

        if (!target.empty())
        {
            const auto& alignerEdges = prepared.alignerEdges;

            aligner_.alignWithCostGuidance(queryPiece.begin(), queryPiece.end(), target.begin(), target.end(),
                                           alignerEdges, minimumCost, preferredInitialBudget(minimumCost));

            std::vector<Cigar> cigars;
            Score secondBestScore = 0;
            Score bestScore = aligner_.backtrackAllPaths<false>(alignerEdges, cigars, secondBestScore);

            score = bestScore;
            observeAlignmentCost(queryPiece.size(), bestScore, minimumCost);
            for (Cigar& cigar : cigars)
            {
                fixFirstNodeExpansion(nodeIds, originalIds, seedPath, cigar);

                unmapNodeIds(originalIds, cigar);

                Path path = seedPath;
                std::list<Operation> operations;
                parseGraphCigar(*seedPath.graphRawPtr(), cigar, path, operations);

                ret.emplace_back(std::move(path), Alignment(seedPath.length(), std::move(operations)));
            }
        }

        return ret;
    }

    std::list<PathAndAlignment>
    suffixAlign(const Path& seedPath, std::string queryPiece, size_t extensionLen, int& score)
    {
        using namespace graphalign::dagAligner;
        std::list<PathAndAlignment> ret;
        if (trySingleNodeAlignment(seedPath, queryPiece, extensionLen, false, ret, score))
            return ret;
        int minimumCost = 0;
        if (tryUniqueUngappedAlignment(seedPath, queryPiece, extensionLen, false, ret, score, &minimumCost))
            return ret;
        ReverseGraph rg(*seedPath.graphRawPtr());
        PreparedGraph temporary;
        const auto& prepared = prepareGraph(*seedPath.graphRawPtr(), seedPath.firstNodeId(),
                                            ConstReversePath(seedPath).endPosition(), extensionLen, true, temporary);
        const auto& target = prepared.target;
        const auto& nodeIds = prepared.nodeIds;
        const auto& originalIds = prepared.originalIds;

        if (!target.empty())
        {
            const auto& alignerEdges = prepared.alignerEdges;

            std::reverse(queryPiece.begin(), queryPiece.end());
            aligner_.alignWithCostGuidance(queryPiece.begin(), queryPiece.end(), target.begin(), target.end(),
                                           alignerEdges, minimumCost, preferredInitialBudget(minimumCost));

            std::vector<Cigar> cigars;
            Score secondBestScore = 0;
            Score bestScore = aligner_.backtrackAllPaths<false>(alignerEdges, cigars, secondBestScore);

            score = bestScore;
            observeAlignmentCost(queryPiece.size(), bestScore, minimumCost);
            for (Cigar& cigar : cigars)
            {
                fixFirstNodeExpansion(nodeIds, originalIds, ConstReversePath(seedPath), cigar);

                unmapNodeIds(originalIds, cigar);

                Path path = seedPath;
                ReversePath rp(path);
                std::list<Operation> operations;
                parseGraphCigar(rg, cigar, rp, operations);
                operations.reverse();

                // reversed alignments always start at the beginning of the path because
                // the seed path gets start-extended to incorporate them
                ret.emplace_back(std::move(path), Alignment(0, std::move(operations)));
            }
        }

        return ret;
    }

private:
    template <typename GraphT>
    static std::map<NodeId, int> extractSubgraph(
        const GraphT& graph, const NodeId startNodeId, const std::size_t startNodeOffset, const std::size_t seqLen)
    {
        std::map<NodeId, int> nodeStartSeqOffset;

        if (graph.nodeSeqLength(startNodeId) == startNodeOffset)
        {
            // flag empty start node in a special way
            nodeStartSeqOffset[startNodeId] = -1;
        }
        else
        {
            // first node start is at the start of the sequence
            nodeStartSeqOffset[startNodeId] = 0;
        }

        // nodes to be visited by bfs
        std::deque<NodeId> shouldVisit(1, startNodeId);

        // extract longest subgraph of nodes such that each node begins within
        // seqLen from the start of the start node
        while (!shouldVisit.empty())
        {
            const NodeId currentNodeId = shouldVisit.front();
            shouldVisit.pop_front();

            const auto currentNodeSeqLength = graph.nodeSeqLength(currentNodeId);
            const int currentNodeSeqOffset = nodeStartSeqOffset[currentNodeId];

            // avoid dealing with individual node startNodeOffset (only start node has it)
            // by pretending the sequence starts at the start node start
            if (seqLen + startNodeOffset - std::max(0, currentNodeSeqOffset) > currentNodeSeqLength)
            {
                const int successorSeqOffset
                    = -1 == currentNodeSeqOffset ? 0 : currentNodeSeqOffset + currentNodeSeqLength;
                // sequence does not terminate at this node, enqueue successors
                const auto successors = graph.successors(currentNodeId);
                for (const NodeId successorId : successors)
                {
                    const auto seenSuccessor = nodeStartSeqOffset.find(successorId);
                    if (nodeStartSeqOffset.end() == seenSuccessor || seenSuccessor->second > successorSeqOffset)
                    {
                        // some successors will end up listed in shouldVisit more than once at a time
                        shouldVisit.push_back(successorId);
                        nodeStartSeqOffset[successorId] = successorSeqOffset;
                    }
                }
            }
        }

        return nodeStartSeqOffset;
    }

    typedef std::pair<MappedId, MappedId> IdEdge;
    typedef std::vector<IdEdge> IdEdges;

    struct MappedNodeRange
    {
        NodeId originalId;
        MappedId begin, end;
    };
    typedef std::vector<MappedNodeRange> MappedNodeRanges;

    /**
     * \brief expands repeats up to remainder of sequnce length
     * \return pairs of mapped node ids indicating an edge between them.
     * \postcondition result array is ordered by successor id then by predecessor id
     */
    template <typename GraphT>
    static IdEdges unrollRepeats(
        const GraphT& graph, const std::size_t seqLen, const std::map<NodeId, int>& nodeStartSeqOffset,
        std::vector<NodeId>& originalIds, MappedNodeRanges& mappedRanges)
    {
        IdEdges idEdges;
        mappedRanges.reserve(nodeStartSeqOffset.size());
        for (const auto& nodeIdOffset : nodeStartSeqOffset)
        {
            if (int(seqLen) <= nodeIdOffset.second)
            {
                throw std::logic_error("node should not be in the subgraph");
            }

            const MappedId firstMappedId = originalIds.size();
            const auto successors = graph.successors(nodeIdOffset.first);
            if (successors.end() != successors.find(nodeIdOffset.first))
            {
                const std::size_t nodeSeqLen = graph.nodeSeqLength(nodeIdOffset.first);
                // don't forget empty repeat start node has special sequence offset -1
                for (std::size_t lenLeft = seqLen - std::max(0, nodeIdOffset.second); lenLeft;)
                {
                    // chain unrolled repeat nodes together
                    IdEdge edge(originalIds.size(), originalIds.size() + 1);
                    idEdges.push_back(edge);
                    originalIds.push_back(nodeIdOffset.first);
                    lenLeft -= std::min(lenLeft, nodeSeqLen);
                }

                // since edges point forward, above loop always produces one more edge than we need
                idEdges.pop_back();
            }
            else
            {
                originalIds.push_back(nodeIdOffset.first);
            }
            // All copies of each original node are generated consecutively,
            // and original nodes are visited in ascending NodeId order.
            mappedRanges.push_back({ nodeIdOffset.first, firstMappedId, MappedId(originalIds.size()) });
        }

        linkPredecessors(graph, mappedRanges, idEdges);

        // group edges by successor node
        std::sort(
            idEdges.begin(), idEdges.end(), [](const IdEdge& left, const IdEdge& right)
            { return left.second < right.second || (left.second == right.second && left.first < right.first); });

        return idEdges;
    }

    static void dfsExtractOrderedNodeIds(
        const MappedId currentId, const IdEdges& idEdges, const std::vector<std::size_t>& idEdgesIndex,
        std::vector<char>& seenNodes, std::vector<MappedId>& nodeIds)
    {
        if (!seenNodes.at(currentId))
        {
            seenNodes[currentId] = true;
            for (std::size_t predOffset = idEdgesIndex.at(currentId); idEdgesIndex.at(currentId + 1) != predOffset;
                 ++predOffset)
            {
                if (idEdges.at(predOffset).second != currentId)
                {
                    throw std::logic_error(
                        "dfsExtractOrderedNodeIds: Invalid edge for node " + std::to_string(currentId));
                }

                dfsExtractOrderedNodeIds(idEdges[predOffset].first, idEdges, idEdgesIndex, seenNodes, nodeIds);
            }
            nodeIds.push_back(currentId);
        }
    }
    /**
     * \brief dfs in order to produce the topological ordering with start node on top
     */
    static std::vector<MappedId>
    extractOrderedNodeIds(const IdEdges& idEdges, const std::vector<std::size_t>& idEdgesIndex)
    {
        std::vector<MappedId> nodeIds;
        std::vector<char> seenNodes(idEdgesIndex.size(), false);
        for (std::size_t mappedId = 0; idEdgesIndex.size() - 1 != mappedId; ++mappedId)
        {
            dfsExtractOrderedNodeIds(mappedId, idEdges, idEdgesIndex, seenNodes, nodeIds);
        }

        return nodeIds;
    }

    /**
     * \brief self-repeat edges already in idEdges, add all mapped predecessor edges to each first expansion
     *        and non-repeat nodes
     */
    template <typename GraphT>
    static void linkPredecessors(
        const GraphT& graph, const MappedNodeRanges& mappedRanges, IdEdges& idEdges)
    {
        for (const MappedNodeRange& range : mappedRanges)
        {
            // Only the first expansion receives external predecessor edges;
            // repeat-to-repeat edges were already added by unrollRepeats.
            for (const NodeId& predecessorId : graph.predecessors(range.originalId))
            {
                if (predecessorId != range.originalId)
                {
                    const auto predecessor = std::lower_bound(
                        mappedRanges.begin(), mappedRanges.end(), predecessorId,
                        [](const MappedNodeRange& entry, NodeId id) { return entry.originalId < id; });
                    if (predecessor == mappedRanges.end() || predecessor->originalId != predecessorId)
                        continue;

                    // The integer range has the same order as equal_range in
                    // the former multimap, without allocating a tree node for
                    // every unrolled repeat copy.
                    for (MappedId mappedPredId = predecessor->begin; mappedPredId != predecessor->end; ++mappedPredId)
                    {
                        idEdges.emplace_back(mappedPredId, range.begin);
                    }
                }
            }
        }
    }

    static std::vector<std::size_t> indexEdges(const IdEdges& idEdges, const std::vector<NodeId>& originalIds)
    {
        std::vector<std::size_t> index;
        index.reserve(idEdges.size() + 1);
        index.push_back(0);
        if (!idEdges.size())
        {
            // no edges. Single-node graph, just close the index
            index.push_back(0);
        }
        else
        {
            MappedId lastId = -1;
            for (const IdEdge& edge : idEdges)
            {
                while (lastId != edge.second)
                {
                    index.push_back(index.back());
                    ++lastId;
                }
                ++index.back();
            }
            // add nodes without predecessors to index
            while (lastId != MappedId((originalIds.size() - 1)))
            {
                index.push_back(index.back());
                ++lastId;
            }
        }

        return index;
    }

    template <typename GraphT>
    static std::string buildTargetSequence(
        const GraphT& graph, const std::size_t startNodeOffset, const std::vector<MappedId>& nodeIds,
        const std::vector<NodeId>& originalIds, const std::map<NodeId, int>& nodeStartSeqOffset,
        const std::vector<std::size_t>& idEdgesIndex, const IdEdges& idEdges, std::vector<Edge>& edges)
    {
        std::string target;
        std::vector<int> mappedIdEndOffset(nodeIds.size(), 0);
        // when first node is a repeat expansion fully consumed by seed, just pretend that query start
        // at the beginning of the node rather than after the end...
        std::size_t startOffset = startNodeOffset;
        for (const MappedId mappedId : nodeIds)
        {
            const NodeId originalId = originalIds.at(mappedId);
            const auto nodeSeq = graph.nodeSeq(originalId);
            if (!startOffset || -1 != nodeStartSeqOffset.at(originalId))
            {
                const std::size_t nodePartLen = nodeSeq.length() - startOffset;
                if (!nodePartLen)
                {
                    throw std::logic_error("Empty node in expanded subgraph and it is not the first one");
                }

                for (std::size_t predOffset = idEdgesIndex[mappedId]; idEdgesIndex[mappedId + 1] != predOffset;
                     ++predOffset)
                {
                    if (idEdges[predOffset].second != mappedId)
                    {
                        throw std::logic_error("bfsDiscoverEdges: Invalid edge");
                    }
                    Edge edge(mappedIdEndOffset.at(idEdges[predOffset].first), target.length());
                    edges.push_back(edge);
                }
                target += nodeSeq.substr(startOffset, nodePartLen);
            }
            mappedIdEndOffset[mappedId] = target.length() - 1;
            startOffset = 0;
        }

        return target;
    }

    /**
     * \param seqLen    extensionLen + offset in the first node. This simplifies
     * the process
     *                  by assuming that sequence starts at the node start
     */
    template <typename GraphT>
    static void bfsDiscoverEdges(
        const GraphT& graph, const NodeId startNodeId, const std::size_t startNodeOffset, const std::size_t seqLen,
        std::vector<MappedId>& nodeIds, std::vector<Edge>& edges, std::string& target,
        std::vector<NodeId>& originalIds)
    {
        const auto startNodeLength = graph.nodeSeqLength(startNodeId);
        if (seqLen && startNodeOffset < startNodeLength
            && (seqLen <= startNodeLength - startNodeOffset || graph.successors(startNodeId).empty()))
        {
            // Extraction cannot leave this node, so the general BFS/unrolling
            // produces exactly one mapped node and no edges. This also holds
            // past a terminal node's end. Keep its complete remaining sequence:
            // truncating it to seqLen would change DP.
            nodeIds.assign(1, 0);
            originalIds.assign(1, startNodeId);
            target = graph.nodeSeq(startNodeId).substr(startNodeOffset);
            return;
        }
        // length of shortest path to the first node character
        const std::map<NodeId, int> nodeStartSeqOffset = extractSubgraph(graph, startNodeId, startNodeOffset, seqLen);

        // Repeat expansions need to be unrolled, so we create unique id for each unrolled instance and map
        // them to original ids
        MappedNodeRanges mappedRanges;
        const IdEdges idEdges
            = unrollRepeats(graph, startNodeOffset + seqLen, nodeStartSeqOffset, originalIds, mappedRanges);

        // index the edge array so that for each mappedId there is a pair of entries index[mappedId],
        // index[mappedId+1] which contains offsets of all predecessor edges in the idEdges;
        const std::vector<std::size_t> idEdgesIndex = indexEdges(idEdges, originalIds);

        nodeIds = extractOrderedNodeIds(idEdges, idEdgesIndex);
        if (nodeIds.size() != originalIds.size())
        {
            throw std::logic_error(
                "Invalid number of nodeIds entries " + std::to_string(nodeIds.size()) + " expected "
                + std::to_string(originalIds.size()));
        }

        if (originalIds.at(nodeIds.at(0)) != startNodeId)
        {
            // nodeIds must be topologically sorted. Since we're discovering our subgraph from the
            // start node, its first expansion must sort to the top.
            throw std::logic_error("First expansion of start node must be the first node");
        }

        // extract target sequence in the proper order
        target = buildTargetSequence(
            graph, startNodeOffset, nodeIds, originalIds, nodeStartSeqOffset, idEdgesIndex, idEdges, edges);

        if (-1 == nodeStartSeqOffset.at(startNodeId))
        {
            // if the first node is an empty repeat expansion (consumed by seed), remove it as it has no
            // corresponding sequence in the target.
            nodeIds.erase(nodeIds.begin());
        }
    }
};
}
