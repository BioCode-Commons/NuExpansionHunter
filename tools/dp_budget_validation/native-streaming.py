#!/usr/bin/env python3
"""Native Linux AVX2 full-program streaming comparison on a fixed synthetic fixture.

Run run-native.py --require-avx2 --skip-benchmarks first. Its tested SIMD control
includes the candidate's identical mapped Graph::hasEdge fix. This runner builds
both complete applications with the same system dependencies and Clang flags,
then measures streaming only; compilation and cache preparation are not timed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import sysconfig

sys.dont_write_bytecode = True

GT = Path('ehunter/thirdparty/graph-tools-master-f421f4c')
GRAPH = GT / 'src/graphcore/Graph.cpp'
PINNED = GT / 'include/graphalign/PinnedDagAligner.hh'
BUDGET = GT / 'include/graphalign/dagAligner/PinnedBudgetAlignMatrix.hh'
EXPECTED_INPUT = {
    'catalog.json': 'aad4ac731a21cde309d33ea2bd69629da1fcb4363bb89f70ee777c23bdd6466e',
    'reads.sam': '5e50a36a857a7c9f94629cfcb70806004e29cc5056cb6a4985239ce7069ecbf3',
    'reference.fa': '37daf2afe74c6b31a8c5ac49016c659100d332d441329536474aa638a1fe49cc',
    'reference.fa.fai': 'b5a18c62a5293acaba23cf94fe997961e3bfe5b2b4bf5a3363b4cbb9f2558c89',
}
PACKAGES = ['clang', 'cmake', 'ninja-build', 'pkg-config', 'libboost-dev',
            'libboost-program-options-dev', 'libboost-filesystem-dev',
            'libboost-system-dev', 'libgtest-dev', 'libspdlog-dev', 'libabsl-dev',
            'libhts-dev', 'zlib1g-dev', 'libbz2-dev', 'liblzma-dev',
            'libcurl4-gnutls-dev']


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def capture(command):
    result = subprocess.run([str(v) for v in command], capture_output=True, text=True)
    return {'command': [str(v) for v in command], 'returncode': result.returncode,
            'stdout': result.stdout, 'stderr': result.stderr}


def source_manifest(tree):
    """Cover full-app C/C++ sources and build recipes, not only the DP harness."""
    return {str(p.relative_to(tree)): sha(p) for p in sorted((tree / 'ehunter').rglob('*'))
            if p.is_file() and (p.suffix in {'.c', '.cc', '.cpp', '.h', '.hh', '.hpp', '.cmake'}
                                or p.name == 'CMakeLists.txt')}


def check_manifest(tree, manifest):
    for name, digest in manifest.items():
        if sha(tree / name) != digest:
            raise RuntimeError('Source changed after standalone correctness: ' + str(tree / name))


def executable(name):
    found = shutil.which(name)
    if not found:
        raise RuntimeError('Required executable not found: ' + name)
    # Do not resolve the clang++ symlink into the C driver name.
    return str(Path(found).absolute())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', required=True, type=Path)
    parser.add_argument('--correctness-results', required=True, type=Path,
                        help='Completed run-native.py results.json on this same runner')
    parser.add_argument('--out', required=True, type=Path, help='New output directory')
    parser.add_argument('--cc', default='clang')
    parser.add_argument('--cxx', default='clang++')
    parser.add_argument('--jobs', type=int, default=2, help='Build concurrency only')
    parser.add_argument('--compute', type=int, default=min(7, max(1, (os.cpu_count() or 1) - 1)))
    parser.add_argument('--pairs', type=int, default=9)
    parser.add_argument('--system-libdir', type=Path,
                        default=Path('/usr/lib') / (sysconfig.get_config_var('MULTIARCH') or ''),
                        help='Directory with shared libhts.so and libcurl.so')
    parser.add_argument('--curl-include', type=Path,
                        default=Path('/usr/include') / (sysconfig.get_config_var('MULTIARCH') or ''),
                        help='Matching system directory containing curl/curl.h')
    args = parser.parse_args()
    if min(args.jobs, args.compute, args.pairs) < 1:
        parser.error('jobs, compute and pairs must be positive')
    if not sys.platform.startswith('linux') or platform.machine().lower() not in {'x86_64', 'amd64'}:
        parser.error('This full-app recipe requires native Linux x86_64; use run-native.py for ARM correctness')
    here = Path(__file__).resolve().parent
    candidate = args.candidate.resolve()
    correctness_path = args.correctness_results.resolve()
    correctness = json.loads(correctness_path.read_text())
    if not correctness.get('complete') or correctness.get('native_features') != {'arch': 'x86_64', 'avx2': True}:
        parser.error('Completed native AVX2 correctness is required')
    if Path(correctness['candidate']).resolve() != candidate:
        parser.error('Candidate does not match the standalone correctness run')
    baseline = Path(correctness['baseline']).resolve()
    control = Path(correctness['control']).resolve()
    output = args.out.resolve()
    for tree in (baseline, control, candidate):
        if output == tree or tree in output.parents:
            parser.error('--out must be outside source trees')
    output.mkdir(parents=True, exist_ok=False)
    report = {'complete': False, 'protocol': 'Native AVX2 full-program streaming only; fair SIMD+hasEdge control',
              'platform': platform.platform(), 'machine': platform.machine(),
              'logical_cpu_count': os.cpu_count(), 'compute_workers': args.compute,
              'preprocessing_workers': 1, 'hts_decode_threads': 0, 'hts_encode_threads': 0,
              'pairs': args.pairs, 'warmup_pairs': 1,
              'correctness_results': str(correctness_path), 'correctness_results_sha256': sha(correctness_path),
              'source_trees': {'published_simd': str(baseline), 'fair_control': str(control), 'candidate': str(candidate)},
              'scope': 'Synthetic covered10k; no extrapolation to clinical or whole-genome datasets, or another CPU',
              'time_scope': 'Whole streaming process including startup, input, alignment, genotyping and output; separate prepared catalog caches already exist. Compilation, input generation, cache preparation and comparisons excluded.',
              'memory_policy': 'Original no-growth RSS gates; same preprocessing count on both sides; no additional DP cache allowance',
              'commands': [], 'builds': {}}
    def save():
        write_json(output / 'results.json', report)

    def run(command, label):
        command = [str(v) for v in command]
        log = output / (label + '.log')
        print('Running ' + label, flush=True)
        with log.open('wb') as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        item = {'command': command, 'returncode': result.returncode, 'log': str(log)}
        report['commands'].append(item)
        save()
        if result.returncode:
            raise RuntimeError('Failed ' + label + '; see ' + str(log))
        return item

    save()
    try:
        # Recheck native execution on the current host, not just a copied report.
        features = capture([correctness_path.parent / 'native-features'])
        if features['returncode'] or json.loads(features['stdout']) != correctness['native_features']:
            raise RuntimeError('Native AVX2 feature witness does not match the current host')
        report['native_features_rechecked'] = features
        for kind, tree in [('baseline', baseline), ('candidate', candidate)]:
            check_manifest(tree, correctness['source_comparison'][kind])
        check_manifest(candidate, correctness['source_sha256'])
        if sha(control / GRAPH) != correctness['control_graph_sha256']:
            raise RuntimeError('Fair-control Graph.cpp changed after correctness')
        manifests = {kind: source_manifest(tree) for kind, tree in
                     [('published_simd', baseline), ('fair_control', control), ('candidate', candidate)]}
        expected_control = dict(manifests['published_simd'])
        expected_control[str(GRAPH)] = manifests['candidate'][str(GRAPH)]
        if manifests['fair_control'] != expected_control:
            raise RuntimeError('Full-app fair control differs from published SIMD beyond Graph::hasEdge')
        changes = sorted(path for path in set(manifests['fair_control']) | set(manifests['candidate'])
                         if manifests['fair_control'].get(path) != manifests['candidate'].get(path))
        if set(changes) != {str(PINNED), str(BUDGET)}:
            raise RuntimeError('Unexpected full-app source/build differences: ' + repr(changes))
        report['source_comparison'] = {'sha256': manifests, 'control_candidate_changed_paths': changes}
        cc, cxx, cmake = executable(args.cc), executable(args.cxx), executable('cmake')
        report['compiler'] = capture([cxx, '--version'])
        if report['compiler']['returncode'] or 'clang' not in report['compiler']['stdout'].lower():
            raise RuntimeError('Use Clang: GraphTools GNU CMake defaults add a global SSE4.2 flag')
        report['cmake'] = capture([cmake, '--version'])
        report['packages'] = capture(['dpkg-query', '-W', '-f=${Package}\t${Version}\n'] + PACKAGES)
        libdir = args.system_libdir.resolve()
        htslib, curl = libdir / 'libhts.so', libdir / 'libcurl.so'
        for library in (htslib, curl):
            if not library.is_file():
                raise RuntimeError('Shared dependency missing: ' + str(library))
        curl_include = args.curl_include.resolve()
        if not (curl_include / 'curl/curl.h').is_file():
            raise RuntimeError('Matching system CURL headers missing: ' + str(curl_include))
        report['shared_dependencies'] = {str(p): {'resolved': str(p.resolve()), 'sha256': sha(p)}
                                         for p in (htslib, curl)}
        flags = '-Wno-unknown-warning-option -Wno-error=deprecated-declarations -Wno-error=deprecated-literal-operator'
        for name, tree in [('baseline', control), ('candidate', candidate)]:
            build = output / ('build-' + name)
            configure = [cmake, '-S', tree / 'ehunter', '-B', build, '-G', 'Ninja',
                         '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_C_COMPILER=' + cc,
                         '-DCMAKE_CXX_COMPILER=' + cxx, '-DCMAKE_CXX_FLAGS=' + flags,
                         '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON', '-DBUILD_TESTS=OFF',
                         '-Dhtslib=' + str(htslib), '-DCURL_LIBRARY=' + str(curl),
                         '-DCURL_NO_CURL_CMAKE=ON', '-DCURL_INCLUDE_DIR=' + str(curl_include),
                         '-DCURL_USE_STATIC_LIBS=OFF']
            run(configure, 'configure-' + name)
            cache = {}
            for line in (build / 'CMakeCache.txt').read_text().splitlines():
                if line and not line.startswith(('#', '//')) and '=' in line:
                    key, value = line.split('=', 1)
                    cache[key.split(':', 1)[0]] = value
            shutil.copyfile(build / 'CMakeCache.txt', output / ('cmake-cache-' + name + '.log'))
            for key, expected in [('htslib', htslib), ('CURL_LIBRARY', curl), ('CURL_INCLUDE_DIR', curl_include)]:
                observed = Path(cache.get(key, '')).resolve()
                if observed != expected.resolve():
                    raise RuntimeError('CMake selected an unexpected dependency: %s; observed=%s; expected=%s'
                                       % (key, observed, expected.resolve()))
            compile_commands = json.loads((build / 'compile_commands.json').read_text())
            for entry in compile_commands:
                tokens = entry.get('arguments') or shlex.split(entry['command'])
                forbidden = [v for v in tokens if v.startswith(('-msse', '-mavx')) or v == '-march=native']
                if forbidden:
                    raise RuntimeError('Unexpected global SIMD compiler flags: ' + repr(forbidden))
            shutil.copyfile(build / 'compile_commands.json', output / ('compile-commands-' + name + '.json'))
            run([cmake, '--build', build, '--target', 'ExpansionHunter', '--parallel', args.jobs], 'build-' + name)
            binary = build / 'ExpansionHunter'
            report['builds'][name] = {'binary': str(binary), 'sha256': sha(binary),
                                      'dynamic_dependencies': capture(['ldd', binary]),
                                      'global_avx_sse_flags_absent': True}
            save()
        helper = output / 'sam_to_bam'
        run([cc, '-O2', here / 'sam_to_bam.c', htslib, '-o', helper], 'build-sam-helper')
        fixture = output / 'fixture'
        run([sys.executable, here / 'generate_catalog.py', '--out', fixture, '--loci', 10000,
             '--covered-loci', 10000, '--coverage-layout', 'even', '--pairs-per-locus', 8,
             '--background-pairs-per-locus', 0, '--read-length', 150, '--stride', 4096,
             '--loci-per-contig', 1000, '--motif', 'CAG', '--repeat-count', 10,
             '--quality', 'mixed', '--mismatch-every', 1, '--missing-quality-every', 11,
             '--seed', 'nueh-round5-v1', '--bam-helper', helper], 'generate-fixture')
        actual_hashes = {name: sha(fixture / name) for name in EXPECTED_INPUT}
        if actual_hashes != EXPECTED_INPUT:
            raise RuntimeError('Synthetic input differs from the ARM covered10k workload')
        report['fixture_uncompressed_sha256'] = actual_hashes
        run([sys.executable, here / 'verify_fixture.py', fixture], 'verify-fixture')
        report['fixture_verification'] = json.loads((fixture / 'independent-validation.json').read_text())
        shutil.copyfile(fixture / 'manifest.json', output / 'fixture-manifest.json')
        shutil.copyfile(fixture / 'independent-validation.json', output / 'fixture-validation.json')
        # No compilation/input generation remains. The harness prepares separate
        # schema-correct caches before warmup and performs all comparisons between,
        # never inside, the measured child processes.
        extra = json.dumps(['--read-decode-threads', '1'])
        run([sys.executable, here / 'paired_streaming.py',
             '--baseline', report['builds']['baseline']['binary'],
             '--candidate', report['builds']['candidate']['binary'],
             '--input', fixture, '--out', output / 'benchmark', '--pairs', args.pairs,
             '--compute', args.compute, '--hts-decode', 0, '--encode', 0,
             '--baseline-extra-json', extra, '--candidate-extra-json', extra,
             '--candidate-preprocess-workers', 1, '--cached'], 'streaming-benchmark')
        paired = json.loads((output / 'benchmark/results.json').read_text())
        if (len(paired['pairs']) != args.pairs or not paired['warmup']['artifacts_identical']
                or not paired['summary']['outputs_identical_under_user_policy']):
            raise RuntimeError('Incomplete streaming run or mismatched complete output artifacts')
        report['summary'] = paired['summary']
        report['output_policy'] = paired['output_policy']
        report['complete'] = True
        save()
        print(json.dumps(report['summary'], indent=2), flush=True)
    except Exception as error:
        report['error'] = str(error)
        save()
        raise


if __name__ == '__main__':
    main()
