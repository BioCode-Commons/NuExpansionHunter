//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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
// This file differs from the GraphTools copy shipped with ExpansionHunter v5.0.0: blank alignment output accepts borrowed locus identifiers.
// Modifications are licensed under the Apache License, Version 2.0.
//

#include "graphio/AlignmentWriter.hh"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/algorithm/string/join.hpp>

using std::string;

namespace graphtools
{

void BlankAlignmentWriter::write(
    TextView /*locusId*/, const std::string& /*fragmentName*/, const std::string& /*query*/, bool /*isFirstMate*/,
    bool /*isReversed*/, bool /*isMateReversed*/, const GraphAlignment& /*alignment*/)
{
}
}
