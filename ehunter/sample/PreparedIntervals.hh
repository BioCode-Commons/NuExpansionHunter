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
