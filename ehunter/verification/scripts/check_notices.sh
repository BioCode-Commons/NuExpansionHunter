#!/bin/sh
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

set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec taskset -c 80-83,164-167 python3 - "$script_dir/../../.." <<'PY'
from pathlib import Path
import os
import re
import subprocess
import sys

os.chdir(Path(sys.argv[1]).resolve())
failures = []

def git(*args):
    return subprocess.check_output(["git", *args])

def paths(*args):
    return [os.fsdecode(p) for p in git(*args).split(b"\0") if p]

def require(condition, message):
    if not condition:
        failures.append(message)

def upstream_header(data):
    if data.startswith(b"/*"):
        end = data.index(b"*/") + 2
        if data[end:end + 2] == b"\r\n":
            end += 2
        elif data[end:end + 1] == b"\n":
            end += 1
        return data[:end]
    prefix = b"//" if data.startswith(b"//") else b"#" if data.startswith((b"# ", b"#\n")) else None
    lines = []
    if prefix:
        for line in data.splitlines(keepends=True):
            if not line.startswith(prefix):
                break
            lines.append(line)
    return b"".join(lines)

def source(path):
    return path.name == "CMakeLists.txt" or path.suffix in {
        ".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".sh", ".py", ".cmake"
    }

def prefix(path):
    return b"#" if path.name == "CMakeLists.txt" or path.suffix in {".sh", ".py", ".cmake"} else b"//"

def change_notice(path, text, graph=False, restricted=False, reference=False):
    c = prefix(path)
    base = b"the GraphTools copy shipped with ExpansionHunter v5.0.0" if graph else b"ExpansionHunter v5.0.0"
    start = (b"" if restricted else c + b"\n")
    start += c + b" Modified by NuExpansionHunter contributors, 2026.\n"
    start += c + b" This file differs from " + base + b": "
    if not text.startswith(start):
        return False
    sentence, separator, tail = text[len(start):].partition(b"\n")
    if not separator or not sentence.endswith(b".") or len(sentence) <= 1:
        return False
    if reference and sentence != b"re-hosted in namespace `reference` for equivalence checking.":
        return False
    if re.search(rb"\b(now|no longer|previously|optimized version)\b", sentence, re.I):
        return False
    if not restricted:
        return tail.startswith(c + b" Modifications are licensed under the Apache License, Version 2.0.\n" + c + b"\n")
    return not tail.startswith(c + b" Modifications are licensed")

new_header = b"""//
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
"""
graph_root = "ehunter/thirdparty/graph-tools-master-f421f4c/"
modified = paths("diff", "--name-only", "-z", "--diff-filter=M", "v5.0.0", "--",
                 "CMakeLists.txt", "ehunter", "cmake")
modified_count = illumina_count = 0
for name in modified:
    path = Path(name)
    if not source(path):
        continue
    original = git("show", "v5.0.0:" + name)
    header = upstream_header(original)
    current = path.read_bytes()
    modified_count += 1
    if b"Illumina" in header:
        illumina_count += 1
    require(current.startswith(header), name + ": upstream header differs")
    require(change_notice(path, current[len(header):], name.startswith(graph_root),
                          name.endswith(("/IntervalTree.h", "/spp.h"))),
            name + ": missing or invalid D.3 notice immediately after upstream header")

added = set(paths("diff", "--name-only", "-z", "--diff-filter=A", "v5.0.0", "--", "ehunter", "cmake"))
added.update(paths("ls-files", "--others", "--exclude-standard", "-z", "--", "ehunter", "cmake"))
new_count = reference_count = 0
for name in sorted(added):
    path = Path(name)
    # Markdown and generated result logs are not source; patches use D.6.
    if not source(path):
        continue
    current = path.read_bytes()
    if name.startswith("ehunter/verification/reference/"):
        reference_count += 1
        relative = name.removeprefix("ehunter/verification/reference/").replace("_v5", "")
        candidate = graph_root + "include/" + relative if relative.startswith("graphalign/") and path.suffix == ".hh" else (
            graph_root + "src/" + relative if relative.startswith("graphalign/") else "ehunter/" + relative)
        original = git("show", "v5.0.0:" + candidate)
        header = upstream_header(original)
        require(b"Illumina" in header and current.startswith(header), name + ": reference upstream header differs")
        require(change_notice(path, current[len(header):].lstrip(b"\r\n"), candidate.startswith(graph_root), reference=True),
                name + ": invalid reference D.3 notice")
    else:
        new_count += 1
        expected = re.sub(rb"^//", b"#", new_header, flags=re.M) if prefix(path) == b"#" else new_header
        body = current.split(b"\n", 1)[1] if current.startswith(b"#!") else current
        require(body.startswith(expected), name + ": missing verbatim D.4 header")

patch_count = 0
for path in sorted(Path("cmake/patches").glob("*.patch")):
    patch_count += 1
    data = path.read_bytes()
    lines = data.splitlines()
    first = next((i for i, line in enumerate(lines) if line.startswith(b"--- ")), None)
    require(first is not None and first > 0, str(path) + ": missing patch preamble")
    if first is None:
        continue
    preamble = lines[:first]
    require(all(not line or line.startswith(b"#") for line in preamble), str(path) + ": non-comment preamble")
    require(not any(line.startswith((b"---", b"+++", b"diff ", b"@@")) for line in preamble),
            str(path) + ": patch syntax in preamble")
    prose = b" ".join(line.removeprefix(b"# ").strip() for line in preamble if line.strip())
    require(b"Applies to htslib 1.10.2." in prose, str(path) + ": missing target version")
    require(b"This patch modifies htslib, Copyright (C) Genome Research Ltd., and is provided under htslib's MIT/Expat license." in prose,
            str(path) + ": missing htslib attribution")

require(git("show", "v5.0.0:LICENSE.txt") == Path("LICENSE.txt").read_bytes(), "LICENSE.txt differs from v5.0.0")
original = git("show", "v5.0.0:COPYRIGHT.txt")
copyright_text = Path("COPYRIGHT.txt").read_bytes()
begin = b"----- begin ExpansionHunter v5.0.0 COPYRIGHT.txt -----\n"
end = b"----- end ExpansionHunter v5.0.0 COPYRIGHT.txt -----"
if begin in copyright_text:
    embedded = copyright_text.split(begin, 1)[1]
    require(embedded[:len(original)] == original, "COPYRIGHT.txt upstream bytes differ")
    separator = b"" if original.endswith(b"\n") else b"\n"
    require(embedded[len(original):].startswith(separator + end), "COPYRIGHT.txt end marker misplaced")
else:
    require(False, "COPYRIGHT.txt begin marker missing")

# Split the obsolete wording so this check does not match its own source.
obsolete = "Modified for " + "NuExpansionHunter"
scan = subprocess.run(["grep", "-rn", obsolete, "."], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
require(scan.returncode == 1, "obsolete notice scan failed or found matches")
if scan.stdout:
    failures.append(scan.stdout.decode(errors="replace").rstrip())
if scan.stderr:
    failures.append(scan.stderr.decode(errors="replace").rstrip())

if failures:
    for failure in failures:
        print("FAIL: " + failure)
    sys.exit(1)
print(f"PASS: {modified_count} modified source files carry D.3 notices.")
print(f"PASS: {illumina_count} modified Illumina headers are preserved byte for byte at the start.")
print(f"PASS: {new_count} added source files carry verbatim D.4 headers.")
print(f"PASS: {reference_count} reference files preserve Illumina headers and carry reference D.3 notices.")
print(f"PASS: {patch_count} htslib patches carry valid D.6 preambles.")
print("PASS: LICENSE.txt is byte-identical to v5.0.0.")
print("PASS: COPYRIGHT.txt embeds the upstream file byte for byte.")
print("PASS: obsolete first-line notice scan produced no output.")
PY
