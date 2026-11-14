#!/usr/bin/env python3
"""Native ARM/AVX2 pinned-DP correctness and fair, serialized performance checks."""
import argparse
import concurrent.futures
import difflib
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

from measure import files_identical

GT = Path('ehunter/thirdparty/graph-tools-master-f421f4c')
GRAPH = GT / 'src/graphcore/Graph.cpp'
MATRIX = GT / 'include/graphalign/dagAligner/AffineAlignMatrixVectorized.hh'
BUDGET = GT / 'include/graphalign/dagAligner/PinnedBudgetAlignMatrix.hh'
PINNED = GT / 'include/graphalign/PinnedDagAligner.hh'
CPP = ['graphalign/Operation.cpp', 'graphalign/OperationOperations.cpp',
       'graphalign/LinearAlignment.cpp', 'graphalign/LinearAlignmentOperations.cpp',
       'graphcore/Graph.cpp', 'graphcore/Path.cpp', 'graphutils/SequenceOperations.cpp']
PUBLISHED_SIMD_REF = '7fa386484fa7c6f783dbe98239ad77fd097d057b'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def capture(command):
    try:
        result = subprocess.run([str(v) for v in command], text=True, capture_output=True)
    except OSError as error:
        return {'command': [str(v) for v in command], 'returncode': None, 'stdout': '', 'stderr': str(error)}
    return {'command': [str(v) for v in command], 'returncode': result.returncode,
            'stdout': result.stdout, 'stderr': result.stderr}


def graph_function(text):
    signature = 'bool Graph::hasEdge(NodeId source_id, NodeId sink_id) const\n'
    if text.count(signature) != 1:
        raise RuntimeError('Expected exactly one Graph::hasEdge definition')
    begin = text.index(signature)
    end = text.index('\n}\n', begin) + len('\n}\n')
    return text[:begin], text[begin:end], text[end:]


def relevant_sources(tree):
    files = set(p for p in (tree / GT / 'include').rglob('*') if p.is_file())
    files.update(tree / GT / 'src' / name for name in CPP)
    return {str(p.relative_to(tree)): sha(p) for p in sorted(files)}


def prepare_control(baseline, candidate, output):
    manifests = {'baseline': relevant_sources(baseline), 'candidate': relevant_sources(candidate)}
    changed = sorted(path for path in set(manifests['baseline']) | set(manifests['candidate'])
                     if manifests['baseline'].get(path) != manifests['candidate'].get(path))
    allowed = {str(GRAPH), str(PINNED), str(BUDGET)}
    unexpected = set(changed) - allowed
    if unexpected:
        raise RuntimeError('Unexpected compiled-source differences: ' + ', '.join(sorted(unexpected)))
    manifests['changed_paths'] = changed
    before, after = (baseline / GRAPH).read_text(), (candidate / GRAPH).read_text()
    old, new = graph_function(before), graph_function(after)
    if old[0] != new[0] or old[2] != new[2]:
        raise RuntimeError('Graph.cpp differs outside hasEdge; refusing an unfair control')
    old_body = '    if (mapped_.data)\n        return successors(source_id).count(sink_id) != 0;\n'
    new_body = ('    if (mapped_.data)\n    {\n'
                '        const auto slice = mappedNode(source_id).successors;\n'
                '        const auto* ids = mapped_.array<NodeId>(slice);\n'
                '        return std::find(ids, ids + slice.count, sink_id) != ids + slice.count;\n'
                '    }\n')
    if old[1].count(old_body) != 1 or old[1].replace(old_body, new_body) != new[1]:
        raise RuntimeError('hasEdge change is not the reviewed raw-slice membership patch')
    if sha(baseline / MATRIX) != sha(candidate / MATRIX):
        raise RuntimeError('SIMD kernels differ; this protocol requires the same SIMD fallback')
    if 'PinnedBudgetAlignMatrix' in (baseline / PINNED).read_text():
        raise RuntimeError('Baseline already uses budget DP; supply the published SIMD checkout')
    if not (candidate / BUDGET).is_file():
        raise RuntimeError('Candidate budget DP header is missing')
    control = output / 'simd-edge-control'
    shutil.copytree(baseline, control, ignore=shutil.ignore_patterns(
        '.git', 'build', 'build-*', '__pycache__', '*.o', '*.a', '*.dylib', '*.so'))
    shutil.copyfile(candidate / GRAPH, control / GRAPH)
    patch = ''.join(difflib.unified_diff(before.splitlines(True), after.splitlines(True),
                                        fromfile='published-SIMD/Graph.cpp', tofile='fair-control/Graph.cpp'))
    (output / 'control-patch.log').write_text(patch)
    return control, manifests


def run_logged(command, prefix):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    stdout, stderr = Path(str(prefix) + '.stdout'), Path(str(prefix) + '.stderr')
    with stdout.open('wb') as out, stderr.open('wb') as err:
        result = subprocess.run([str(v) for v in command], stdout=out, stderr=err)
    record = {'command': [str(v) for v in command], 'returncode': result.returncode,
              'stdout_path': str(stdout), 'stderr_path': str(stderr),
              'output_bytes': stdout.stat().st_size, 'output_sha256': sha(stdout)}
    write_json(str(prefix) + '.json', record)
    if result.returncode:
        raise RuntimeError('Command failed: %s; see %s' % (command, stderr))
    return record


def json_lines(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('{'):
            records.append(json.loads(line))
    return records


def assert_equal(records, label):
    first = records[0]['stdout_path']
    if not all(files_identical(first, item['stdout_path']) for item in records[1:]):
        raise RuntimeError('Complete stdout differs: ' + label)
    return {'byte_equal': True, 'output_bytes': records[0]['output_bytes'],
            'output_sha256': records[0]['output_sha256'], 'runs': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True, type=Path, help='Published SIMD checkout')
    parser.add_argument('--candidate', required=True, type=Path, help='Budget DP + hasEdge candidate checkout')
    parser.add_argument('--out', required=True, type=Path, help='New output directory, outside baseline')
    parser.add_argument('--boost-include', type=Path, help='Omit when Boost is in compiler system headers')
    parser.add_argument('--cxx', default=os.environ.get('CXX', 'c++'), help='Compiler executable, not a shell command')
    parser.add_argument('--jobs', type=int, default=2, help='Build concurrency; measurements are always serialized')
    parser.add_argument('--require-avx2', action='store_true', help='Fail unless running natively on an AVX2 x86_64 host')
    parser.add_argument('--skip-benchmarks', action='store_true', help='Correctness only; no performance claim')
    parser.add_argument('--prepare-only', action='store_true', help='Create and validate fair-control source; no compilation')
    parser.add_argument('--pairs', type=int, default=9)
    parser.add_argument('--target-seconds', type=float, default=0.7)
    args = parser.parse_args()
    baseline, candidate, output = args.baseline.resolve(), args.candidate.resolve(), args.out.resolve()
    if output == baseline or baseline in output.parents:
        parser.error('--out must be outside the baseline tree (it is copied into the output)')
    if baseline == candidate or args.jobs < 1 or args.pairs < 1 or args.target_seconds <= 0:
        parser.error('Distinct trees and positive jobs/pairs/target-seconds are required')
    compiler = shutil.which(args.cxx)
    if compiler is None:
        parser.error('Compiler executable not found: ' + args.cxx)
    compiler = str(Path(compiler).absolute())  # Keep the c++/clang++ driver name when it is a symlink.
    output.mkdir(parents=True, exist_ok=False)
    here = Path(__file__).resolve().parent
    report = {'complete': False, 'protocol': 'Native full PinnedDagAligner; fair SIMD+hasEdge control',
              'expected_published_simd_ref': PUBLISHED_SIMD_REF, 'baseline': str(baseline), 'candidate': str(candidate),
              'platform': platform.platform(), 'machine': platform.machine(),
              'compiler': capture([compiler, '--version']), 'compiler_sha256': sha(compiler),
              'baseline_git': capture(['git', '-C', baseline, 'rev-parse', 'HEAD']),
              'candidate_git': capture(['git', '-C', candidate, 'rev-parse', 'HEAD']),
              'source_sha256': {str(p): sha(candidate / p) for p in [GRAPH, MATRIX, BUDGET, PINNED]},
              'probe_sha256': {p.name: sha(p) for p in sorted(here.iterdir()) if p.suffix in ['.py', '.cpp']},
              'builds': {}, 'correctness': {}}
    try:
        if sys.platform == 'darwin':
            report['cpu'] = {key: capture(['sysctl', '-n', key]) for key in
                             ['hw.model', 'hw.ncpu', 'hw.physicalcpu', 'hw.logicalcpu', 'machdep.cpu.brand_string']}
        else:
            report['cpu'] = {'logical_cpu_count': os.cpu_count()}
            if sys.platform.startswith('linux') and Path('/proc/cpuinfo').exists():
                allowed = {'model name', 'vendor_id', 'flags', 'cpu family', 'model', 'stepping'}
                fields = {}
                for line in Path('/proc/cpuinfo').read_text().splitlines():
                    if ':' not in line:
                        continue
                    key, value = (part.strip() for part in line.split(':', 1))
                    if key in allowed and value not in fields.setdefault(key, []):
                        fields[key].append(value)
                report['cpu']['cpuinfo_fields'] = fields
        baseline_head = report['baseline_git']['stdout'].strip() if report['baseline_git']['returncode'] == 0 else None
        report['baseline_reference_status'] = ('matches_expected_commit' if baseline_head == PUBLISHED_SIMD_REF
                                                else 'different_head_recorded' if baseline_head else 'unversioned_snapshot')
        control, source_comparison = prepare_control(baseline, candidate, output)
        report['source_comparison'] = source_comparison
        report['control'] = str(control)
        report['control_graph_sha256'] = sha(control / GRAPH)
        write_json(output / 'results.json', report)
        if args.prepare_only:
            report['prepared_only'] = True
            write_json(output / 'results.json', report)
            return

        flags = ['-std=c++11', '-O3', '-DNDEBUG']
        if args.boost_include:
            flags += ['-isystem', str(args.boost_include.resolve())]
        feature_source = output / 'native-features.cpp'
        feature_source.write_text(r'''#include <iostream>
#if defined(__APPLE__) && defined(__x86_64__)
#include <sys/types.h>
#include <sys/sysctl.h>
#include <cerrno>
#endif
int main() {
// Rosetta translates x86_64. An arm64 helper matching an arm64 host is native
// without consulting a sysctl that may be unavailable inside a sandbox.
#if defined(__APPLE__) && defined(__x86_64__)
  int translated = 0;
  size_t length = sizeof(translated);
  const int translationStatus = sysctlbyname("sysctl.proc_translated", &translated, &length, nullptr, 0);
  if (translationStatus != 0 && errno != ENOENT) {
    std::cout << "{\"arch\":\"unknown\",\"translation_check_errno\":" << errno << "}\n";
    return 0;
  }
  if (translationStatus == 0 && translated == 1) {
    std::cout << "{\"arch\":\"translated\",\"translated\":true}\n";
    return 0;
  }
#endif
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  __builtin_cpu_init();
  std::cout << "{\"arch\":\"x86_64\",\"avx2\":" << (__builtin_cpu_supports("avx2") ? "true" : "false") << "}\n";
#elif defined(__aarch64__) && defined(__ARM_NEON)
  std::cout << "{\"arch\":\"aarch64\",\"neon\":true}\n";
#else
  std::cout << "{\"arch\":\"unsupported\"}\n";
#endif
}
''')
        feature_binary = output / 'native-features'
        feature_command = [compiler] + flags + [str(feature_source), '-o', str(feature_binary)]
        with (output / 'native-features-build.log').open('wb') as log:
            subprocess.run(feature_command, check=True, stdout=log, stderr=subprocess.STDOUT)
        report['native_features'] = json.loads(subprocess.check_output([str(feature_binary)]))
        report['native_feature_command'] = feature_command
        features = report['native_features']
        host_arch = {'arm64': 'aarch64', 'aarch64': 'aarch64', 'amd64': 'x86_64',
                     'x86_64': 'x86_64'}.get(platform.machine().lower(), 'unsupported')
        if features.get('translated') or features.get('arch') != host_arch:
            raise RuntimeError('Compiled helper architecture differs from native host, or executes translated')
        translated = False
        if sys.platform == 'darwin' and host_arch == 'x86_64':
            translation = capture(['sysctl', '-in', 'sysctl.proc_translated'])
            report['process_translation'] = translation
            translated = translation['stdout'].strip() == '1'
        if translated:
            raise RuntimeError('Translated execution is not native performance evidence')
        if args.require_avx2 and not (platform.machine().lower() in ['x86_64', 'amd64']
                                       and features.get('avx2') is True):
            raise RuntimeError('Native AVX2 required; runtime feature check did not pass')
        if not (features.get('avx2') is True or features.get('neon') is True):
            raise RuntimeError('This protocol requires native ARM NEON or native x86_64 AVX2')

        binaries = output / 'binaries'
        binaries.mkdir()
        jobs = [
            ('published-has-edge', baseline, 'has-edge-probe.cpp', []),
            ('candidate-has-edge', candidate, 'has-edge-probe.cpp', []),
            ('control-pinned', control, 'pinned-probe.cpp', []),
            ('candidate-pinned', candidate, 'pinned-probe.cpp', []),
            ('control-lowerror', control, 'lowerror-probe.cpp', []),
            ('candidate-lowerror', candidate, 'lowerror-probe.cpp', []),
            ('instrumented-lowerror', candidate, 'lowerror-probe.cpp', ['GRAPHTOOLS_BUDGET_DP_INSTRUMENT']),
            ('control-activation', control, 'activation-proof.cpp', []),
            ('candidate-activation', candidate, 'activation-proof.cpp', ['GRAPHTOOLS_BUDGET_DP_INSTRUMENT']),
            ('scalar-insertion', candidate, 'insertion-guard.cpp', ['GRAPHTOOLS_DP_TEST_FORCE_SCALAR']),
            ('auto-insertion', candidate, 'insertion-guard.cpp', ['GRAPHTOOLS_DP_TEST_INSTRUMENT_AVX2']),
        ]
        if features.get('avx2'):
            jobs.append(('forced-avx2-insertion', candidate, 'insertion-guard.cpp',
                         ['GRAPHTOOLS_DP_TEST_FORCE_AVX2', 'GRAPHTOOLS_DP_TEST_INSTRUMENT_AVX2']))

        def build(job):
            name, tree, source, defines = job
            gt = tree / GT
            command = [compiler] + flags + ['-I' + str(gt / 'include')]
            command += ['-D' + value for value in defines] + [str(here / source)]
            if source == 'has-edge-probe.cpp':
                command += ['-pthread'] + [str(gt / 'src' / p) for p in
                                           ['graphcore/Graph.cpp', 'graphutils/SequenceOperations.cpp']]
            elif source != 'insertion-guard.cpp':
                command += [str(gt / 'src' / p) for p in CPP]
            binary = binaries / name
            command += ['-o', str(binary)]
            with (output / (name + '-build.log')).open('wb') as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            return name, {'binary': str(binary), 'binary_sha256': sha(binary), 'command': command}

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for name, result in pool.map(build, jobs):
                report['builds'][name] = result
                write_json(output / 'results.json', report)
        def run(name, arguments, label):
            return run_logged([report['builds'][name]['binary']] + list(map(str, arguments)), output / 'correctness' / label)

        records = [run(name, [], name) for name in ['published-has-edge', 'candidate-has-edge']]
        membership = assert_equal(records, 'published SIMD / candidate hasEdge ownership and bounds')
        final_words = Path(records[0]['stdout_path']).read_text().splitlines()[-1].split()
        if len(final_words) != 4 or final_words[0] != 'checks' or final_words[2] != 'digest':
            raise RuntimeError('Unexpected hasEdge correctness summary')
        membership['checks'], membership['digest'] = int(final_words[1]), int(final_words[3])
        if membership['checks'] != 2663:
            raise RuntimeError('Unexpected hasEdge correctness check count')
        report['correctness']['has_edge'] = membership
        write_json(output / 'results.json', report)

        records = [run(name, ['verify', '--dump'], name) for name in ['control-pinned', 'candidate-pinned']]
        pinned = assert_equal(records, '9,380 full pinned cases')
        final_line = Path(records[0]['stdout_path']).read_text().splitlines()[-1]
        pinned['cases'] = json.loads(final_line)['cases']
        if pinned['cases'] != 9380:
            raise RuntimeError('Unexpected pinned correctness case count')
        report['correctness']['pinned'] = pinned
        records = [run(name, [], name) for name in ['control-activation', 'candidate-activation']]
        activation = assert_equal(records, 'explicit sparse activation and fallback witnesses')
        activation['witnesses'] = json_lines(records[1]['stderr_path'])
        if len(activation['witnesses']) != 28:
            raise RuntimeError('Expected 28 individually asserted activation witnesses')
        report['correctness']['activation'] = activation

        lowerror_cases = []
        report['correctness']['lowerror'] = {}
        for shape in ['linear', 'join', 'repeat']:
            for length in [150, 250]:
                for errors in [0, 1, 2, 8]:
                    name = 'lowerror-%s%d-e%d' % (shape, length, errors)
                    common = ['bench', 'mapped' if shape == 'join' else 'owned', 'mixed',
                              shape, length, length + 12, 0, errors, 100]
                    lowerror_cases.append([name, common])
                    verify = common[:]; verify[0] = 'verify'; verify[-1] = 64
                    records = [run(side, verify, name + '-' + side)
                               for side in ['control-lowerror', 'candidate-lowerror', 'instrumented-lowerror']]
                    result = assert_equal(records, name)
                    result['stats'] = next(r['budget_dp'] for r in json_lines(records[-1]['stderr_path']) if 'budget_dp' in r)
                    if result['stats']['calls'] != 64:
                        raise RuntimeError('Unexpected low-error sample count: ' + name)
                    report['correctness']['lowerror'][name] = result
                    write_json(output / 'results.json', report)

        insertion_names = ['scalar-insertion', 'auto-insertion']
        if features.get('avx2'):
            insertion_names.append('forced-avx2-insertion')
        records = [run(name, [], name) for name in insertion_names]
        insertion = assert_equal(records, 'independent scalar insertion oracle / native SIMD')
        insertion['witnesses'] = {name: json_lines(record['stderr_path'])[-1]
                                  for name, record in zip(insertion_names, records)}
        if features.get('avx2'):
            for name in ['auto-insertion', 'forced-avx2-insertion']:
                if insertion['witnesses'][name].get('avx2_blocks', 0) <= 0:
                    raise RuntimeError('AVX2 insertion kernel was not exercised: ' + name)
        report['correctness']['insertion'] = insertion
        write_json(output / 'results.json', report)
        # The complete insertion streams were byte-compared. They are large
        # binary canonical data, so retain only hashes/logs for CI artifacts.
        for record in records:
            Path(record['stdout_path']).unlink()
            record['stdout_removed_after_comparison'] = True

        selected = [case for case in lowerror_cases if case[0] in
                    ['lowerror-linear150-e1', 'lowerror-join150-e1', 'lowerror-repeat150-e1']]
        selected_file = output / 'selected-lowerror-cases.json'
        write_json(selected_file, selected)
        if not args.skip_benchmarks:
            report['benchmark'] = {}
            for group, bname, cname, extra in [
                ('ordinary', 'control-pinned', 'candidate-pinned', []),
                ('lowerror', 'control-lowerror', 'candidate-lowerror', ['--case-file', str(selected_file)]),
            ]:
                bench_output = output / ('benchmark-' + group)
                command = [sys.executable, str(here / 'benchmark.py'), '--baseline', report['builds'][bname]['binary'],
                           '--candidate', report['builds'][cname]['binary'], '--out', str(bench_output),
                           '--pairs', str(args.pairs), '--target-seconds', str(args.target_seconds)] + extra
                with (output / ('benchmark-' + group + '.log')).open('w') as log:
                    subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
                report['benchmark'][group] = {p.parent.parent.name: json.loads(p.read_text())
                    for p in sorted(bench_output.glob('*/acceptance/summary.json'))}
                write_json(output / 'results.json', report)
        else:
            report['benchmarks_skipped'] = True
        report['complete'] = True
        write_json(output / 'results.json', report)
        print(json.dumps({'complete': True, 'results': str(output / 'results.json'),
                          'benchmarks_skipped': args.skip_benchmarks}), flush=True)
    except Exception as error:
        report['error'] = str(error)
        write_json(output / 'results.json', report)
        raise


if __name__ == '__main__':
    main()
