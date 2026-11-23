# GraphTools modifications

This directory contains the GraphTools copy shipped with Illumina ExpansionHunter
v5.0.0, with changes covered by the notices in the source files and
[NOTICE.txt](../../../NOTICE.txt). The version and license are in
[LICENSE.txt](LICENSE.txt).

## Modified files

- `CMakeLists.txt`
- `external/sparsepp/sparsepp-e40d7a0/sparsepp/spp.h`
- `include/graphalign/DagAlignerAffine.hh`
- `include/graphalign/GappedAligner.hh`
- `include/graphalign/GraphAlignment.hh`
- `include/graphalign/KmerIndex.hh`
- `include/graphalign/OperationOperations.hh`
- `include/graphalign/PinnedDagAligner.hh`
- `include/graphalign/dagAligner/AffineAlignMatrixVectorized.hh`
- `include/graphalign/dagAligner/Details.hh`
- `include/graphcore/Graph.hh`
- `include/graphcore/Path.hh`
- `include/graphcore/PathOperations.hh`
- `include/graphio/AlignmentWriter.hh`
- `include/graphutils/BaseMatching.hh`
- `include/graphutils/SequenceOperations.hh`
- `src/graphalign/GappedAligner.cpp`
- `src/graphalign/GraphAlignmentOperations.cpp`
- `src/graphalign/KmerIndex.cpp`
- `src/graphalign/LinearAlignment.cpp`
- `src/graphalign/LinearAlignmentOperations.cpp`
- `src/graphalign/OperationOperations.cpp`
- `src/graphcore/Graph.cpp`
- `src/graphcore/GraphCoordinates.cpp`
- `src/graphcore/Path.cpp`
- `src/graphcore/PathOperations.cpp`
- `src/graphio/AlignmentWriter.cpp`
- `src/graphio/GraphJson.cpp`
- `src/graphutils/SequenceOperations.cpp`
- `tests/BaseMatchingTest.cpp`
- `tests/DagAlignerTest.cpp`
- `tests/GappedAlignerTest.cpp`
- `tests/GraphAlignmentOperationsTest.cpp`
- `tests/GraphAlignmentTest.cpp`
- `tests/GraphJsonTest.cpp`
- `tests/GraphTest.cpp`
- `tests/KmerIndexTest.cpp`
- `tests/PathOperationsTest.cpp`
- `tests/PathTest.cpp`
- `tests/PinnedDagAlignerTest.cpp`

## Added files

- `include/graphalign/dagAligner/DpKernels.hh`
- `include/graphalign/dagAligner/PinnedBudgetAlignMatrix.hh`
- `include/graphutils/BaseMatchingKernels.hh`
- `include/graphutils/ExtensionTraceCache.hh`
- `include/graphutils/MappedData.hh`
- `include/graphutils/VerifyPolicy.hh`
