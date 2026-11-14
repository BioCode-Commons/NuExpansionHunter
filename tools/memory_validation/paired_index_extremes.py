#!/usr/bin/env python3
"""Run independent exact audits and paired edge-regime index microbenchmarks."""
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
p.add_argument('--profiles', nargs='*')
p.add_argument('--pairs', type=int, default=5)
p.add_argument('--target-seconds', type=float, default=0.15)
p.add_argument('--audit-only', action='store_true')
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
binaries = {v: getattr(a, v).resolve() for v in ('baseline', 'candidate')}
profiles = a.profiles or subprocess.check_output([str(binaries['candidate']), 'list'], text=True).split()
report = {'platform': platform.platform(), 'machine': platform.machine(),
          'binaries': {v: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                       for v, path in binaries.items()}, 'profiles': []}

def save():
    (a.out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')

def run(version, profile, operation, rounds):
    result = subprocess.run([str(binaries[version]), 'bench', profile, operation, str(rounds)],
                            capture_output=True, text=True, check=True)
    value = json.loads(result.stderr)
    value['stdout'] = result.stdout.strip()
    return value

for profile in profiles:
    record = {'name': profile, 'audits': {}, 'cases': []}
    report['profiles'].append(record)
    hashes = {}
    for version, binary in binaries.items():
        output = a.out / f'audit-{version}-{profile}.txt'
        metadata = a.out / f'audit-{version}-{profile}.json'
        with output.open('w') as out, metadata.open('w') as err:
            subprocess.run([str(binary), 'audit', profile], stdout=out, stderr=err, check=True)
        record['audits'][version] = json.loads(metadata.read_text())
        hashes[version] = hashlib.sha256(output.read_bytes()).hexdigest()
    if hashes['baseline'] != hashes['candidate']:
        raise RuntimeError(f'Full encode and ordered-path audit differs: {profile}')
    candidate = record['audits']['candidate']
    if candidate['prefix'] != candidate['expected_prefix']:
        raise RuntimeError(f'Expected fast path / fallback was not exercised: {profile}')
    record['audit_sha256'] = hashes
    record['equal'] = True
    save()
    print('audit', profile, 'PASS', 'prefix' if candidate['prefix'] else 'fallback', flush=True)
    if a.audit_only:
        continue
    operations = ['miss'] if candidate['entries'] == 0 else ['hit']
    if candidate['entries'] and candidate['entries'] < (1 << (2 * candidate['k'])):
        operations.append('miss')
    if candidate['entries'] and not candidate['membership']:
        operations.append('paths')
    for operation in operations:
        calibration = {v: run(v, profile, operation, 100) for v in binaries}
        fastest = min(x['seconds'] for x in calibration.values())
        rounds = max(100, math.ceil(100 * a.target_seconds / max(fastest, 1e-9)))
        case = {'operation': operation, 'rounds': rounds, 'calibration': calibration, 'pairs': []}
        record['cases'].append(case)
        for pair in range(a.pairs + 1):
            result = {'warmup': pair == 0, 'execution_order': ['baseline', 'candidate'] if pair % 2 == 0
                      else ['candidate', 'baseline']}
            for version in result['execution_order']:
                result[version] = run(version, profile, operation, rounds)
            if result['baseline']['stdout'] != result['candidate']['stdout']:
                raise RuntimeError(f'Hot-loop checksums differ: {profile} {operation}')
            result['ratio'] = result['candidate']['seconds'] / result['baseline']['seconds']
            case['pairs'].append(result)
            save()
        measured = case['pairs'][1:]
        case['median_ratio'] = statistics.median(x['ratio'] for x in measured)
        case['faster_pairs'] = sum(x['ratio'] < 1 for x in measured)
        save()
        print(profile, operation, case['median_ratio'], flush=True)
