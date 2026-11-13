#!/usr/bin/env python3
"""Compile identical probes against two existing GraphTools builds and compare stdout.

Example (paths may be absolute or relative):
  python3 tools/round3_validation/verify.py \
    --baseline-source /path/pr-baseline --candidate-source /path/candidate \
    --baseline-build /path/baseline-build --candidate-build /path/candidate-build \
    --extra-include /path/dependencies/include --compiler clang++ --out /path/results

Supply matching Release libraries in each build's
thirdparty/graph-tools-master-f421f4c/libgraphtools.a. --extra-include is repeatable.
Only probes are compiled; the application and supplied libraries are not rebuilt.
All four probe sources come from --candidate-source, including the source included
by ungapped-probe.cpp. Complete stdout is compared byte-for-byte, not by hash.
The existing endpoint probe emits individual endpoint/score and ordered-path data;
its traceback-flag digest remains a diagnostic, not a raw flag-by-flag comparison.
Results, commands, manifests, complete stdout/stderr and binaries remain in --out.
Exit status is zero only when every compile/run and every comparison succeeds.
This Release differential is separate from sanitizer and pipeline-output checks.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


GRAPH_TOOLS = Path("ehunter/thirdparty/graph-tools-master-f421f4c")
LIBRARY = Path("thirdparty/graph-tools-master-f421f4c/libgraphtools.a")
PROBES = [
    ("pinned-probe", "tools/dp_budget_validation/pinned-probe.cpp", ["verify", "--dump"]),
    ("ungapped-probe", "tools/dp_budget_validation/ungapped-probe.cpp", []),
    ("pinned-wrapper-probe", "tools/round3_validation/pinned-wrapper-probe.cpp", []),
    ("endpoint-state-probe", "tools/round3_validation/endpoint-state-probe.cpp", []),
]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare_bytes(first, second):
    """Return exact equality and the first differing byte, reading bounded chunks."""
    offset = 0
    with first.open("rb") as a, second.open("rb") as b:
        while True:
            left, right = a.read(1024 * 1024), b.read(1024 * 1024)
            if left != right:
                index = next((i for i, pair in enumerate(zip(left, right)) if pair[0] != pair[1]),
                             min(len(left), len(right)))
                return False, offset + index
            if not left:
                return True, None
            offset += len(left)


def summary(path, probe):
    lines = path.read_text(encoding="utf-8").splitlines()
    if probe == "endpoint-state-probe":
        states = [line.split()[0] for line in lines if re.fullmatch(r"\S+ \d+ \d+", line)]
        return {
            "states": len(states), "state_names": states,
            "next_endpoint_checks": sum(line.startswith("next ") for line in lines),
            "ordered_path_groups": sum(line.startswith("ordered-paths ") for line in lines),
            "empty_input_exceptions": sum(line.startswith(("empty-query-error ", "empty-target-error "))
                                          for line in lines),
            "traceback_flag_digests": sum(line.startswith("traceback-flags ") for line in lines),
        }
    if not lines:
        raise ValueError("Probe produced empty stdout")
    if lines[-1].startswith("TOTAL "):
        return {"cases": int(lines[-1].split()[1])}
    return json.loads(lines[-1])


def file_record(path, out):
    return {"path": str(path.relative_to(out)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def execute(command, stdout, stderr):
    started = time.monotonic()
    with stdout.open("wb") as output, stderr.open("wb") as errors:
        completed = subprocess.run(command, stdout=output, stderr=errors, check=False)
    return {"command": command, "returncode": completed.returncode,
            "elapsed_seconds": time.monotonic() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("baseline-source", "candidate-source", "baseline-build", "candidate-build", "out"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--extra-include", action="append", default=[], type=Path,
                        help="Additional dependency include directory, repeatable")
    parser.add_argument("--compiler", default="clang++")
    args = parser.parse_args()
    compiler = shutil.which(args.compiler)
    if compiler is None:
        parser.error("Compiler not found: " + args.compiler)
    # Preserve the driver name: clang++ may be a symlink to clang, whose argv[0]
    # changes automatic C++ runtime linkage on some platforms.
    compiler = str(Path(compiler).absolute())
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    trees = {side: getattr(args, side + "_source").resolve() for side in ("baseline", "candidate")}
    builds = {side: getattr(args, side + "_build").resolve() for side in trees}
    extra_includes = [path.resolve() for path in args.extra_include]
    for side in trees:
        for path in (trees[side] / GRAPH_TOOLS / "include", builds[side] / LIBRARY):
            if not path.exists():
                parser.error("Required path does not exist: " + str(path))
    for path in extra_includes:
        if not path.is_dir():
            parser.error("Include directory does not exist: " + str(path))
    for _, relative, _ in PROBES:
        if not (trees["candidate"] / relative).is_file():
            parser.error("Probe source does not exist: " + str(trees["candidate"] / relative))

    report = {
        "complete": False, "passed": False,
        "scope": "Identical probe sources, original/candidate production headers and supplied matching GraphTools libraries; complete stdout byte comparison.",
        "compiler": compiler, "sources": {key: str(value) for key, value in trees.items()},
        "build_directories": {key: str(value) for key, value in builds.items()},
        "extra_include": [str(path) for path in extra_includes],
        "source_sha256": {}, "library_sha256": {}, "builds": [], "comparisons": [],
        "endpoint_traceback_flags": "The unchanged endpoint probe emits a digest for flags; endpoint values and ordered paths are emitted individually. The runner compares all emitted bytes.",
    }
    version = execute([compiler, "--version"], out / "compiler.stdout", out / "compiler.stderr")
    version["stdout"] = file_record(out / "compiler.stdout", out)
    version["stderr"] = file_record(out / "compiler.stderr", out)
    report["compiler_version"] = version
    monitored = {}
    for side, tree in trees.items():
        gt = tree / GRAPH_TOOLS
        paths = sorted(path for directory in ("include", "src") for path in (gt / directory).rglob("*") if path.is_file())
        report["source_sha256"][side] = {str(path.relative_to(tree)): sha256(path) for path in paths}
        library = builds[side] / LIBRARY
        report["library_sha256"][side] = sha256(library)
        monitored.update({path: sha256(path) for path in paths + [library]})
    probe_paths = [trees["candidate"] / relative for _, relative, _ in PROBES]
    report["probe_source_sha256"] = {str(path.relative_to(trees["candidate"])): sha256(path) for path in probe_paths}
    monitored.update({path: sha256(path) for path in probe_paths})
    report["runner_sha256"] = sha256(Path(__file__).resolve())

    def save():
        (out / "results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    save()
    try:
        for name, relative, arguments in PROBES:
            records = {}
            for side, tree in trees.items():
                gt = tree / GRAPH_TOOLS
                binary = out / (side + "-" + name)
                include_dirs = [gt / "include", gt / "external/include"] + extra_includes
                command = [compiler, "-std=c++11", "-O3", "-DNDEBUG", "-Wno-deprecated-declarations",
                           "-Wno-deprecated-literal-operator", "-Wno-deprecated-copy"]
                command += ["-I" + str(path) for path in include_dirs]
                command += [str(trees["candidate"] / relative), str(builds[side] / LIBRARY), "-o", str(binary)]
                build_stdout = out / (side + "-" + name + "-build.stdout")
                build_stderr = out / (side + "-" + name + "-build.stderr")
                record = {"side": side, "probe": name, "build": execute(command, build_stdout, build_stderr)}
                record["build"].update(stdout=file_record(build_stdout, out), stderr=file_record(build_stderr, out))
                report["builds"].append(record)
                records[side] = record
                save()
                if record["build"]["returncode"]:
                    continue
                record["binary"] = file_record(binary, out)
                stdout = out / (side + "-" + name + ".stdout")
                stderr = out / (side + "-" + name + ".stderr")
                record["run"] = execute([str(binary)] + arguments, stdout, stderr)
                record["run"].update(stdout=file_record(stdout, out), stderr=file_record(stderr, out))
                if record["run"]["returncode"] == 0:
                    try:
                        record["summary"] = summary(stdout, name)
                    except (ValueError, IndexError) as error:
                        record["summary_error"] = str(error)
                save()

            comparison = {"probe": name, "full_stdout_identical": False}
            successful = all(record.get("run", {}).get("returncode") == 0 and "summary" in record
                             for record in records.values())
            if successful:
                first, second = [out / (side + "-" + name + ".stdout") for side in ("baseline", "candidate")]
                equal, position = compare_bytes(first, second)
                comparison.update(full_stdout_identical=equal, first_difference_byte=position,
                                  output_bytes={side: records[side]["run"]["stdout"]["bytes"] for side in records},
                                  summaries={side: records[side]["summary"] for side in records})
            else:
                comparison["error"] = "At least one compile/run or summary failed; see retained logs."
            report["comparisons"].append(comparison)
            save()
            print(json.dumps(comparison), flush=True)
        report["changed_inputs_during_verification"] = [str(path) for path, digest in monitored.items()
                                                       if sha256(path) != digest]
        report["complete"] = True
        report["passed"] = all(item["full_stdout_identical"] for item in report["comparisons"]) and not report["changed_inputs_during_verification"]
        save()
        return 0 if report["passed"] else 1
    except (OSError, subprocess.SubprocessError) as error:
        report["error"] = str(error)
        save()
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
