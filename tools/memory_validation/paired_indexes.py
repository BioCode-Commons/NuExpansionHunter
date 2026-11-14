#!/usr/bin/env python3
"""Exact index audit and alternating paired hot-loop measurements.

Build index_probe.cpp against each version's GraphTools library and pass the two
executables here. Both must use identical release compiler options and headers.
Timings are the probe's steady_clock interval; construction, archive writing,
query generation, and checks are excluded. Full-app tests remain necessary.
"""
import argparse
import hashlib
import json
import math
import pathlib
import platform
import statistics
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--baseline', required=True, type=pathlib.Path)
p.add_argument('--candidate', required=True, type=pathlib.Path)
p.add_argument('--out', required=True, type=pathlib.Path)
p.add_argument('--pairs', type=int, default=9)
p.add_argument('--target-seconds', type=float, default=0.2)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
binaries = {v: getattr(a, v).resolve() for v in ('baseline', 'candidate')}
report = {'platform': platform.platform(), 'machine': platform.machine(),
          'binaries': {v: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                       for v, path in binaries.items()}, 'audits': [], 'cases': []}

def save():
    (a.out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')

def run(version, k, kind, rounds, order):
    result = subprocess.run([str(binaries[version]), 'bench', str(k), kind, str(rounds), order],
                            capture_output=True, text=True, check=True)
    fields = result.stderr.split()
    return {'seconds': float(fields[1]), 'record_bytes': int(fields[3]), 'stdout': result.stdout.strip()}

for k in (10, 14):
    hashes = {}
    for version, binary in binaries.items():
        output = a.out / f'audit-{version}-{k}.txt'
        sizes = a.out / f'audit-{version}-{k}.sizes'
        with output.open('w') as out, sizes.open('w') as err:
            subprocess.run([str(binary), 'audit', str(k)], stdout=out, stderr=err, check=True)
        hashes[version] = hashlib.sha256(output.read_bytes()).hexdigest()
    if hashes['baseline'] != hashes['candidate']:
        raise RuntimeError(f'Full encode and ordered-path audit differs at k={k}')
    report['audits'].append({'k': k, 'sha256': hashes, 'equal': True})
    save()

cases = [(10, 'hit', 'shuffled'), (10, 'miss', 'shuffled'), (14, 'hit', 'shuffled'),
         (14, 'miss', 'shuffled'), (14, 'paths', 'shuffled'),
         (10, 'hit', 'rolling'), (14, 'hit', 'rolling'), (14, 'paths', 'rolling')]
for k, kind, order in cases:
    calibration = {v: run(v, k, kind, 100, order) for v in binaries}
    rounds = max(100, math.ceil(100 * a.target_seconds / min(r['seconds'] for r in calibration.values())))
    case = {'k': k, 'operation': kind, 'query_order': order, 'rounds': rounds,
            'calibration': calibration, 'pairs': []}
    report['cases'].append(case)
    for pair in range(a.pairs + 1):
        result = {'warmup': pair == 0, 'execution_order': ['baseline', 'candidate'] if pair % 2 == 0
                  else ['candidate', 'baseline']}
        for version in result['execution_order']:
            result[version] = run(version, k, kind, rounds, order)
        if result['baseline']['stdout'] != result['candidate']['stdout']:
            raise RuntimeError(f'Hot-loop checksums differ: {k} {kind} {order}')
        result['ratio'] = result['candidate']['seconds'] / result['baseline']['seconds']
        case['pairs'].append(result)
        save()
    measured = case['pairs'][1:]
    case['median_ratio'] = statistics.median(x['ratio'] for x in measured)
    case['baseline_seconds'] = statistics.median(x['baseline']['seconds'] for x in measured)
    case['candidate_seconds'] = statistics.median(x['candidate']['seconds'] for x in measured)
    case['faster_pairs'] = sum(x['ratio'] < 1 for x in measured)
    save()
    print(k, kind, order, case['baseline_seconds'], case['candidate_seconds'], case['median_ratio'], flush=True)
