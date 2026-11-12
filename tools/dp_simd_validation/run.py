#!/usr/bin/env python3
"""Portable native-build insertion differential and sanitizer verification."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parent
workspace = root.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out', required=True, help='New output directory; existing directories are not overwritten')
parser.add_argument('--baseline', default=str(workspace / 'baseline'), help='Baseline repository root')
parser.add_argument('--candidate', default=str(workspace / 'candidate'), help='Candidate repository root')
parser.add_argument('--boost-include', help='Optional Boost header directory; otherwise use compiler search paths')
parser.add_argument('--cxx', default='c++', help='C++ compiler executable')
parser.add_argument('--define', action='append', default=[], help='Candidate-only preprocessor definition; repeatable')
args = parser.parse_args()
out = Path(args.out).resolve()
out.mkdir(parents=True, exist_ok=False)
trees = {'baseline': Path(args.baseline).resolve(), 'candidate': Path(args.candidate).resolve()}
include = Path('ehunter/thirdparty/graph-tools-master-f421f4c/include')
boost = Path(args.boost_include).resolve() if args.boost_include else None
source = root / 'insertion-guard.cpp'
report = {'probe_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'commands': [], 'results': {},
          'baseline': str(trees['baseline']), 'candidate': str(trees['candidate']),
          'boost_include': str(boost) if boost else None, 'candidate_defines': args.define,
          'native_build': True}

def save():
    (out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

for name, sanitize in [('baseline', False), ('candidate', False), ('baseline', True), ('candidate', True)]:
    label = name + ('-sanitized' if sanitize else '')
    binary = out / label
    tree = trees[name]
    header = tree / include / 'graphalign/dagAligner/AffineAlignMatrixVectorized.hh'
    command = [args.cxx, '-std=c++11', '-O1' if sanitize else '-O3', '-DNDEBUG', '-I' + str(tree / include)]
    if boost:
        command += ['-isystem', str(boost)]
    if name == 'candidate':
        command += ['-D' + define for define in args.define]
    if sanitize:
        command += ['-g', '-fsanitize=address,undefined', '-fno-sanitize-recover=all', '-fno-omit-frame-pointer']
    command += [str(source), '-o', str(binary)]
    with (out / (label + '.build.log')).open('wb') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    report['commands'].append({'argv': command, 'returncode': result.returncode})
    save()
    if result.returncode: raise RuntimeError(label + ' build failed')
    expected = out / 'baseline.canonical.bin'
    env = dict(os.environ)
    if sanitize: env.update(ASAN_OPTIONS='detect_container_overflow=1:detect_leaks=0', UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    stderr = out / (label + '.stderr')
    if label == 'baseline':
        with expected.open('wb') as output, stderr.open('wb') as errors:
            result = subprocess.run([str(binary)], stdout=output, stderr=errors, env=env)
        if result.returncode: raise RuntimeError('baseline failed')
        digest = sha(expected)
        count = expected.stat().st_size
    else:
        digest = hashlib.sha256()
        count = 0
        with expected.open('rb') as wanted, stderr.open('wb') as errors:
            process = subprocess.Popen([str(binary)], stdout=subprocess.PIPE, stderr=errors, env=env)
            try:
                while True:
                    block = process.stdout.read(1024 * 1024)
                    reference = wanted.read(len(block)) if block else wanted.read(1)
                    if reference != block: raise RuntimeError(label + ' differs near byte ' + str(count))
                    if not block: break
                    digest.update(block)
                    count += len(block)
                if process.wait(): raise RuntimeError(label + ' failed')
            except BaseException:
                if process.poll() is None: process.kill()
                process.wait()
                raise
        digest = digest.hexdigest()
    report['results'][label] = {
        'header_sha256': sha(header),
        'binary_sha256': sha(binary),
        'canonical_bytes': count, 'canonical_sha256': digest,
        'summary': json.loads(stderr.read_text()), 'byte_comparison_pass': True,
        'sanitized': sanitize,
    }
    save()
    print(label + ': independent scalar oracle and canonical byte comparison passed', flush=True)
    print(json.dumps(report['results'][label]['summary']), flush=True)
report['complete'] = True
save()
