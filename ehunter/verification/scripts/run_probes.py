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

"""Build probes with matching supplied libraries, run, and compare complete stdout.

--compiler/CXX selects the compiler; --cxxflags/CXXFLAGS and --ldflags/LDFLAGS
are shell-split flag strings (no shell is invoked). --include/CPATH supplies
repeatable dependency include directories. --package-build/EH_PACKAGE_BUILD and
--v5-build/EH_V5_BUILD name existing builds of this package and v5.0.0.
--verification-build SETTINGS=DIRECTORY is repeatable: each directory must have
been configured with -DEH_VERIFY="SETTINGS". EH_VERIFICATION_BUILDS can contain a
JSON object mapping settings to directories. Supply matching static libraries;
every supplied policy must match a probe. The script never builds ExpansionHunter.
--library and --v5-library override the
package and v5 library lists, respectively (repeat for each library); environment
alternatives EH_PACKAGE_LIBRARIES and EH_V5_LIBRARIES are path-separated lists.
--dependency-library/EH_DEPENDENCY_LIBRARIES adds libraries shared by all builds.
Policy builds use the same relative library locations as the package build.
Library paths outside the package build must be supplied as dependency libraries.

Default library locations are libExpansionHunterLib.a and
thirdparty/graph-tools-master-f421f4c/libgraphtools.a under each supplied build.
For a static dependency use --dependency-library or --ldflags as needed, e.g.
--ldflags='-lz -ldl -pthread'. Include paths and flags must match the library ABI.
The package-only probes compare an in-process reference in each build as well as
comparing output between default and the documented policy. v5-era probes use a
temporary detached v5.0.0 worktree, removed even on failure. --out keeps commands,
compiler diagnostics, complete stdout/stderr and executables for review.

Run under the shared build/run lock and CPU affinity appropriate for the machine.
"""
import argparse
import filecmp
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

GT = Path("ehunter/thirdparty/graph-tools-master-f421f4c")
PROBES = (
    ("insertion_scan_oracle", "SIMD=SCALAR"),
    ("pinned_align_cases", "BUDGET_SEARCH=NEVER;FAST_PATHS=OFF;EXTENSION_CACHE=OFF"),
    ("ungapped_fastpath_cases", "FAST_PATHS=OFF;EXTENSION_CACHE=OFF;SIMD=SCALAR"),
    ("dense_dp_cases", None),
    ("budget_matrix_state", "BUDGET_SEARCH=ALWAYS;FUSED_ROWS=OFF"),
    ("pinned_wrapper_cases", "GRAPH_CACHE=OFF;FAST_PATHS=OFF;EXTENSION_CACHE=OFF"),
    ("extension_cache_cases", "EXTENSION_CACHE=OFF;BUDGET_SEARCH=ON"),
    ("graph_has_edge_oracle", "LAYOUT=REFERENCE"),
    ("kmer_index_oracle", "LAYOUT=REFERENCE"),
    ("orientation_oracle", None),
    ("pinned_suffix_vs_v5", None),
)


def env_paths(name):
    return [Path(value) for value in os.environ.get(name, "").split(os.pathsep) if value]


def settings(value):
    return dict(item.split("=", 1) for item in value.split(";") if item)


def cache_policy(build):
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if line.startswith("EH_VERIFY:"):
            return settings(line.split("=", 1)[1])
    return {}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--compiler", default=os.environ.get("CXX", "c++"))
    parser.add_argument("--cxxflags", default=os.environ.get("CXXFLAGS", "-O2 -ffp-contract=off"))
    parser.add_argument("--ldflags", default=os.environ.get("LDFLAGS", "-pthread"))
    parser.add_argument("--include", action="append", type=Path, default=env_paths("CPATH"))
    parser.add_argument("--package-build", type=Path, default=os.environ.get("EH_PACKAGE_BUILD"))
    parser.add_argument("--v5-build", type=Path, default=os.environ.get("EH_V5_BUILD"))
    parser.add_argument("--verification-build", action="append", default=[], metavar="SETTINGS=DIRECTORY")
    parser.add_argument("--library", action="append", type=Path, default=env_paths("EH_PACKAGE_LIBRARIES"))
    parser.add_argument("--v5-library", action="append", type=Path, default=env_paths("EH_V5_LIBRARIES"))
    parser.add_argument("--dependency-library", action="append", type=Path, default=env_paths("EH_DEPENDENCY_LIBRARIES"))
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=1800)
    args = parser.parse_args()
    if args.package_build is None or args.v5_build is None:
        parser.error("Supply --package-build and --v5-build (or their environment variables)")
    required_policies = {frozenset(settings(policy).items()): policy for _, policy in PROBES if policy}
    policies = {}
    try:
        supplied = list(json.loads(os.environ.get("EH_VERIFICATION_BUILDS", "{}")).items())
        supplied += [value.rsplit("=", 1) for value in args.verification_build]
        for key, path in supplied:
            policy = required_policies.get(frozenset(settings(key).items()))
            if policy is None:
                parser.error("Unused verification build: " + key)
            if policy in policies:
                parser.error("Duplicate verification build: " + key)
            policies[policy] = Path(path).resolve()
    except (AttributeError, TypeError, ValueError) as error:
        parser.error("Invalid verification build: " + str(error))
    missing = [policy for policy in required_policies.values() if policy not in policies]
    if missing:
        parser.error("Missing verification builds: " + ", ".join(missing))
    source = Path(__file__).resolve().parents[3]
    verification = source / "ehunter/verification"
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    package_build, v5_build = args.package_build.resolve(), args.v5_build.resolve()
    default_names = [Path("libExpansionHunterLib.a"), GT.relative_to("ehunter") / "libgraphtools.a"]
    package_libraries = [p.resolve() for p in args.library] or [package_build / p for p in default_names]
    v5_libraries = [p.resolve() for p in args.v5_library] or [v5_build / p for p in default_names]
    dependencies = [str(p.resolve()) for p in args.dependency_library]
    compiler = shlex.split(args.compiler)
    flags = shlex.split(args.cxxflags)

    def execute(command, prefix, timeout=None):
        prefix.with_suffix(".command.json").write_text(json.dumps(command, indent=2) + "\n")
        with prefix.with_suffix(".stdout").open("wb") as stdout, prefix.with_suffix(".stderr").open("wb") as stderr:
            try:
                return subprocess.run(command, stdout=stdout, stderr=stderr, timeout=timeout).returncode
            except (OSError, subprocess.TimeoutExpired) as error:
                stderr.write((str(error) + "\n").encode())
                return 1

    all_equal = True
    with tempfile.TemporaryDirectory(prefix="eh-reference-") as temporary:
        v5_source = Path(temporary) / "source"
        created = False
        try:
            created = execute(["git", "-C", str(source), "worktree", "add", "--detach", str(v5_source), "v5.0.0"], out / "worktree-add") == 0
            for name, policy in PROBES:
                build_failed = False
                runs = []
                for side in ("package", "reference"):
                    tree, libraries, defines = source, package_libraries, []
                    try:
                        if side == "package":
                            if cache_policy(package_build):
                                raise ValueError("Package build must have empty EH_VERIFY")
                        elif policy is None:
                            if not created:
                                raise ValueError("Cannot create v5.0.0 worktree")
                            tree, libraries = v5_source, v5_libraries
                        else:
                            build = policies[policy]
                            if cache_policy(build) != settings(policy):
                                raise ValueError("Policy build EH_VERIFY does not match " + policy)
                            libraries = [build / p.relative_to(package_build) for p in package_libraries]
                            defines = ["-DEH_VERIFY_" + key + "=EH_VERIFY_" + value for key, value in settings(policy).items()]
                        for library in libraries:
                            if not library.is_file():
                                raise ValueError("Missing library " + str(library))
                    except (ValueError, OSError) as error:
                        (out / (name + "-" + side + "-configuration.stderr")).write_text(str(error) + "\n")
                        build_failed = True
                        continue
                    include_dirs = [tree / "ehunter", tree / "ehunter/thirdparty", tree / GT / "include",
                                    tree / GT / "external/include", verification] + [p.resolve() for p in args.include]
                    binary = out / (name + "-" + side)
                    sources = [verification / "probes" / (name + ".cpp")]
                    if name == "kmer_index_oracle":
                        sources.append(verification / "reference/graphalign/KmerIndex_v5.cpp")
                    command = compiler + ["-std=c++11"] + flags + defines + ["-I" + str(p) for p in include_dirs]
                    command += ["-isystem", str(tree / GT / "external/sparsepp/sparsepp-e40d7a0")]
                    command += [str(p) for p in sources] + [str(p) for p in libraries] + dependencies
                    command += shlex.split(args.ldflags) + ["-o", str(binary)]
                    if execute(command, out / (name + "-" + side + "-build")):
                        build_failed = True
                        continue
                    prefix = out / (name + "-" + side + "-run")
                    code = execute([str(binary)], prefix, args.timeout)
                    runs.append((code, prefix.with_suffix(".stdout")))
                equal = (not build_failed and len(runs) == 2 and all(code == 0 for code, _ in runs)
                         and filecmp.cmp(runs[0][1], runs[1][1], shallow=False))
                print(name + ": " + ("BUILD FAILED" if build_failed else "IDENTICAL" if equal else "DIFFERS"), flush=True)
                all_equal = all_equal and equal
        finally:
            if created:
                if execute(["git", "-C", str(source), "worktree", "remove", "--force", str(v5_source)], out / "worktree-remove"):
                    raise RuntimeError("Could not remove temporary reference worktree; see worktree-remove.stderr")
    return 0 if all_equal else 1


if __name__ == "__main__":
    raise SystemExit(main())
