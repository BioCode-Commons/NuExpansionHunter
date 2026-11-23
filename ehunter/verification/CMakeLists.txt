#
# NuExpansionHunter
# Copyright 2026 NuExpansionHunter contributors
#
# This file is part of NuExpansionHunter, a modified version of
# Illumina ExpansionHunter v5.0.0 (Copyright 2016-2019 Illumina, Inc.).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# Included only by the opt-in EH_BUILD_VERIFICATION branch of the parent project.
add_library(VerificationReferences STATIC
    reference/genotyping/AlignMatrix_v5.cpp
    reference/genotyping/AlignMatrixFiltering_v5.cpp
    reference/genotyping/StrAlign_v5.cpp
    reference/graphalign/KmerIndex_v5.cpp)
target_include_directories(VerificationReferences PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# The unmodified sparse hash dependency is external to the reference code.
target_include_directories(VerificationReferences SYSTEM PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/graph-tools-master-f421f4c/external/sparsepp/sparsepp-e40d7a0)
target_link_libraries(VerificationReferences PUBLIC ExpansionHunterLib graphtools)

set(VERIFICATION_PROBES
    insertion_scan_oracle pinned_align_cases ungapped_fastpath_cases dense_dp_cases
    budget_matrix_state pinned_wrapper_cases graph_has_edge_oracle kmer_index_oracle
    orientation_oracle pinned_suffix_vs_v5 extension_cache_cases)
foreach(probe IN LISTS VERIFICATION_PROBES)
    add_executable(${probe} probes/${probe}.cpp)
    target_link_libraries(${probe} VerificationReferences ExpansionHunterLib graphtools)
    add_test(NAME verification_${probe} COMMAND ${probe})
endforeach()
if(TARGET extension_cache_cases)
    set_tests_properties(verification_extension_cache_cases PROPERTIES TIMEOUT 30)
endif()
# Not registered with ctest: it needs a BAM, a reference and catalogs (see README).
add_executable(catalog_row_order probes/catalog_row_order.cpp)
target_link_libraries(catalog_row_order ExpansionHunterLib graphtools)
add_executable(K-25_align_matrix harness/K-25_align_matrix.cpp)
target_link_libraries(K-25_align_matrix VerificationReferences ExpansionHunterLib graphtools)
add_test(NAME verification_K-25_align_matrix COMMAND K-25_align_matrix)
