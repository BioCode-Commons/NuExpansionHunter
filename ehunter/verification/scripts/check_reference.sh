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

# Each row is reference path | v5.0.0 path | extraction (whole or writer).
# The JSON header hosts just JsonWriter::write; variant visitors remain in the product.
REFERENCES=(
'reference/graphalign/PinnedDagAligner_v5.hh|ehunter/thirdparty/graph-tools-master-f421f4c/include/graphalign/PinnedDagAligner.hh|whole'
'reference/graphalign/KmerIndex_v5.hh|ehunter/thirdparty/graph-tools-master-f421f4c/include/graphalign/KmerIndex.hh|whole'
'reference/graphalign/KmerIndex_v5.cpp|ehunter/thirdparty/graph-tools-master-f421f4c/src/graphalign/KmerIndex.cpp|whole'
'reference/genotyping/AlignMatrix_v5.hh|ehunter/genotyping/AlignMatrix.hh|whole'
'reference/genotyping/AlignMatrix_v5.cpp|ehunter/genotyping/AlignMatrix.cpp|whole'
'reference/genotyping/AlignMatrixFiltering_v5.hh|ehunter/genotyping/AlignMatrixFiltering.hh|whole'
'reference/genotyping/AlignMatrixFiltering_v5.cpp|ehunter/genotyping/AlignMatrixFiltering.cpp|whole'
'reference/genotyping/StrAlign_v5.hh|ehunter/genotyping/StrAlign.hh|whole'
'reference/genotyping/StrAlign_v5.cpp|ehunter/genotyping/StrAlign.cpp|whole'
'reference/io/JsonWriter_v5.hh|ehunter/io/JsonWriter.cpp|writer'
)
set -eu
verification=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python3 - "$verification" "${REFERENCES[@]}" <<'PYTHON'
import difflib
from pathlib import Path
import subprocess
import sys
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
verification = Path(sys.argv[1])
for entry in sys.argv[2:]:
    target, source, extraction = entry.split("|")
    original = subprocess.check_output(["git", "-C", str(verification.parents[1]), "show", "v5.0.0:" + source], text=True)
    if extraction == "writer":
        # Compare the selected upstream function without renaming body variables.
        header = original[:original.index("#include")]
        body = original[original.index("void JsonWriter::write"):original.index("\nstatic string encodeGenotype")]
        original = header + '#include "io/JsonWriter.hh"\n#include <iomanip>\n#include <sstream>\nnamespace ehunter\n{\nusing Json = nlohmann::json;\nusing std::string;\n' + body + '\n}\n'
    actual = (verification / target).read_text()
    sys.stdout.writelines(difflib.unified_diff(original.splitlines(True), actual.splitlines(True),
        fromfile="v5.0.0:" + source + (" (JsonWriter::write)" if extraction == "writer" else ""), tofile=target))
PYTHON
