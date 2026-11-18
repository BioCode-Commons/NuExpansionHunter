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

# eh_optional_library(<option> <header> <include var> <library var> <library name> [LINK_ONLY])
#
# <option> takes AUTO (the default), ON or OFF.
#   AUTO  use the library when its development files are installed: the header
#         <header> and the library <library name>. A runtime-only shared object
#         (for example libjemalloc.so.2 without libjemalloc-dev) is not used.
#   ON    the same search, but configuration fails when the library is missing.
#   OFF   do not use the library.
# <include var> and <library var> may be set on the command line to give the
# locations. LINK_ONLY marks a library whose header the build does not use: a
# <library var> set on the command line is then used without the header.
# Sets <option>_FOUND to ON or OFF in the caller's scope.
function(eh_optional_library OPTION HEADER INCLUDE_VAR LIBRARY_VAR LIBRARY_NAME)
    set(${OPTION} AUTO CACHE STRING "Use the optional library: AUTO, ON or OFF")
    set_property(CACHE ${OPTION} PROPERTY STRINGS AUTO ON OFF)
    string(TOUPPER "${${OPTION}}" mode)
    if(NOT mode MATCHES "^(AUTO|ON|OFF)$")
        message(FATAL_ERROR "${OPTION} must be AUTO, ON or OFF, not '${${OPTION}}'")
    endif()
    set(found OFF)
    if(NOT mode STREQUAL "OFF")
        find_path(${INCLUDE_VAR} ${HEADER})
        # The library is searched for only when the header is found, so that a
        # runtime-only package is never picked up.
        if(${INCLUDE_VAR})
            find_library(${LIBRARY_VAR} NAMES ${LIBRARY_NAME})
        endif()
        if(${INCLUDE_VAR} AND ${LIBRARY_VAR})
            set(found ON)
        elseif(${LIBRARY_VAR} AND "LINK_ONLY" IN_LIST ARGN)
            set(found ON)
        elseif(mode STREQUAL "ON")
            message(FATAL_ERROR "${OPTION}=ON, but ${HEADER} or the ${LIBRARY_NAME} library was not found. "
                    "Install its development package, or set ${INCLUDE_VAR} and ${LIBRARY_VAR}.")
        endif()
    endif()
    set(${OPTION}_FOUND ${found} PARENT_SCOPE)
endfunction()
