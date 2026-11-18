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

#pragma once
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <string>
#include "graphalign/PinnedDagAligner.hh"
#include "reference/graphalign/PinnedDagAligner_v5.hh"

// The same graph and scoring values feed both wrappers, including mapped graphs.
class CheckedPinnedAligner
{
    graphtools::PinnedDagAligner actual_;
    graphtools::reference::PinnedDagAligner expected_;

    template <typename Call> std::list<graphtools::PathAndAlignment> check(Call call, int& score)
    {
        int expectedScore = score;
        std::list<graphtools::PathAndAlignment> expected, actual;
        std::string expectedError, actualError;
        bool expectedThrew = false, actualThrew = false;
        try { expected = call(expected_, expectedScore); }
        catch (const std::exception& error) { expectedThrew = true; expectedError = error.what(); }
        try { actual = call(actual_, score); }
        catch (const std::exception& error) { actualThrew = true; actualError = error.what(); }
        if (actualThrew != expectedThrew || actualError != expectedError || score != expectedScore || actual != expected)
        { std::cerr << "Pinned wrapper differs from reference\n"; std::exit(1); }
        if (actualThrew) throw std::runtime_error(actualError);
        return actual;
    }
    // C++11 callable avoids generic lambdas while sharing the comparison body.
    struct Call
    {
        const graphtools::Path& seed;
        const std::string& query;
        size_t extension;
        bool suffix;
        template <typename Aligner> std::list<graphtools::PathAndAlignment> operator()(Aligner& aligner, int& score) const
        {
            return suffix ? aligner.suffixAlign(seed, query, extension, score)
                          : aligner.prefixAlign(seed, query, extension, score);
        }
    };
public:
    CheckedPinnedAligner(int match, int mismatch, int opening, int extension)
        : actual_(match, mismatch, opening, extension), expected_(match, mismatch, opening, extension) {}
    std::list<graphtools::PathAndAlignment> prefixAlign(const graphtools::Path& seed, const std::string& query, size_t extension, int& score)
    { return check(Call{seed, query, extension, false}, score); }
    std::list<graphtools::PathAndAlignment> suffixAlign(const graphtools::Path& seed, const std::string& query, size_t extension, int& score)
    { return check(Call{seed, query, extension, true}, score); }
};
