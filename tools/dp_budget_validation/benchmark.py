#!/usr/bin/env python3
"""Serialized full pinned-alignment comparison against the SIMD baseline."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from measure import files_identical, measure

CASES = [
    ('linear150', ['bench', 'owned', 'prefix', 'linear', 150, 162, 0, 100]),
    ('linear250', ['bench', 'owned', 'mixed', 'linear', 250, 262, 0, 100]),
    ('join150', ['bench', 'mapped', 'mixed', 'join', 150, 162, 0, 100]),
    ('repeat150', ['bench', 'owned', 'mixed', 'repeat', 150, 162, 0, 100]),
    ('repeat250', ['bench', 'mapped', 'mixed', 'repeat', 250, 262, 0, 100]),
    ('short3', ['bench', 'owned', 'prefix', 'linear', 3, 7, 0, 100]),
    ('short16', ['bench', 'owned', 'mixed', 'linear', 16, 24, 0, 100]),
    ('long-node150', ['bench', 'owned', 'prefix', 'linear', 150, 162, 2048, 100]),
    ('short-target150', ['bench', 'owned', 'prefix', 'linear', 150, 162, 22, 100]),
]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline', required=True, type=Path)
parser.add_argument('--candidate', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
parser.add_argument('--only', nargs='*')
parser.add_argument('--case-file', type=Path, help='Optional explicit list of named probe argument vectors')
parser.add_argument('--pairs', type=int, default=9)
parser.add_argument('--target-seconds', type=float, default=0.7)
args = parser.parse_args()
root = args.out.resolve()
root.mkdir(parents=True, exist_ok=False)
binaries = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
for name, original in (json.loads(args.case_file.read_text()) if args.case_file else CASES):
    if args.only is not None and name not in args.only:
        continue
    out = root / name
    out.mkdir()
    common = [str(v) for v in original]
    history = []
    for trial in range(6):
        records = {side: measure([str(binary)] + common, out / ('calibration-%d-%s' % (trial, side)))
                   for side, binary in binaries.items()}
        if not files_identical(records['baseline']['stdout_path'], records['candidate']['stdout_path']):
            raise RuntimeError('Calibration output mismatch: ' + name)
        history.append({'args': common[:], 'records': records})
        fastest = min(record['wall_seconds'] for record in records.values())
        if fastest >= args.target_seconds * 0.8:
            break
        common[-1] = str(max(int(common[-1]) + 1,
                            math.ceil(int(common[-1]) * args.target_seconds / max(fastest, 0.001))))
    else:
        raise RuntimeError('Calibration did not converge: ' + name)
    (out / 'case.json').write_text(json.dumps({'name': name, 'args': common, 'calibration': history,
         'scope': 'Full PinnedDagAligner call: extraction, unrolling, DP and ordered traceback',
         'binary_sha256': {side: hashlib.sha256(binary.read_bytes()).hexdigest()
                           for side, binary in binaries.items()}}, indent=2) + '\n')
    command = [sys.executable, str(Path(__file__).with_name('measure.py')),
               '--baseline', str(binaries['baseline']), '--candidate', str(binaries['candidate']),
               '--out', str(out / 'acceptance'), '--pairs', str(args.pairs), '--'] + common
    print(json.dumps({'started': name, 'args': common}), flush=True)
    with (out / 'run.log').open('w') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    summary = json.loads((out / 'acceptance/summary.json').read_text())
    print(json.dumps({'finished': name, 'wall_reduction_percent': summary['paired_speedup_percent'],
                      'cpu_reduction_percent': summary['paired_cpu_reduction_percent'],
                      'wall_ratio_interval': summary['median_ratio_bootstrap_95']}), flush=True)
