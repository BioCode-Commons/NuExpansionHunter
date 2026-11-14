#!/usr/bin/env python3
"""Run isolated streaming stress comparisons serially, retaining every failure.

Uses the existing complete VCF/JSON/BAM comparator and wait4 peak-RSS measurement.
Short workloads and confidence intervals crossing parity remain inconclusive;
this runner does not translate a microbenchmark into a general speed guarantee.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys


def classify(metric):
    lower, upper = metric['paired_ratio_bootstrap_95']
    return 'faster' if upper < 1 else 'slower' if lower > 1 else 'inconclusive'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--baseline-source', type=Path, required=True)
    parser.add_argument('--candidate-source', type=Path, required=True)
    parser.add_argument('--input-root', type=Path, action='append', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--select', nargs='+')
    parser.add_argument('--modes', nargs='+', choices=['owned', 'mapped'], default=['owned', 'mapped'])
    parser.add_argument('--threads', nargs='+', type=int, default=[1, 8])
    parser.add_argument('--pairs', type=int, default=3)
    parser.add_argument('--timeout', type=float, default=900, help='Seconds per fixture/mode, including preparation')
    args = parser.parse_args()
    if args.pairs < 1 or min(args.threads) < 1 or args.timeout <= 0:
        parser.error('Counts and timeout must be positive')
    fixtures = sorted(p for root in args.input_root for p in root.resolve().iterdir()
                      if p.is_dir() and (p/'manifest.json').is_file())
    if args.select:
        fixtures = [p for p in fixtures if p.name in args.select]
        if {p.name for p in fixtures} != set(args.select):
            parser.error('One or more selected fixtures are missing')
    if not fixtures or len({p.name for p in fixtures}) != len(fixtures):
        parser.error('Need at least one fixture and unique fixture names')
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    runner = Path(__file__).resolve().parents[1]/'frontend_validation/run_matrix.py'
    report = {'scope': 'Synthetic stress tests, not WGS or biological truth validation.',
              'policy': 'Exact VCF/JSON bytes; exact BAM header and record-byte multiset including duplicates.',
              'memory_target_ratio': 0.5, 'runs': [], 'complete': False}

    def save():
        (out/'results.json').write_text(json.dumps(report, indent=2)+'\n')

    for fixture in fixtures:
        for mode in args.modes:
            destination = out/(fixture.name+'-'+mode)
            # Eight colliding active loci exceed the two-worker queue capacity.
            threads = [2, 8] if fixture.name == 'stripe-collisions' else args.threads
            command = [sys.executable, str(runner), '--baseline', str(args.baseline.resolve()),
                       '--candidate', str(args.candidate.resolve()),
                       '--baseline-source', str(args.baseline_source.resolve()),
                       '--candidate-source', str(args.candidate_source.resolve()),
                       '--input', str(fixture), '--out', str(destination), '--threads',
                       *map(str, threads), '--preprocess', '0', '--decode', '2', '--encode', '1',
                       '--pairs', str(args.pairs)]
            if mode == 'owned':
                command.append('--uncached')
            record = {'fixture': fixture.name, 'mode': mode, 'command': command,
                      'log': str(destination)+'.log', 'status': 'running'}
            report['runs'].append(record)
            save()
            with open(record['log'], 'w') as log:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    record['returncode'] = process.wait(timeout=args.timeout)
                    record['status'] = 'passed' if record['returncode'] == 0 else 'failed'
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    record['status'] = 'timeout'
            result_file = destination/'results.json'
            if result_file.exists():
                result = json.loads(result_file.read_text())
                record['configurations'] = []
                for config in result['configurations']:
                    summary = config.get('summary')
                    if summary is None:
                        continue
                    rss = summary['rss']
                    record['configurations'].append({
                        'threads': config['compute_workers'], 'summary': summary,
                        'memory_target_met': (rss['candidate_median'] <= 0.5*rss['baseline_median']
                                              and rss['candidate_max'] <= 0.5*rss['baseline_max']),
                        'wall_evidence': classify(summary['wall']),
                        'read_processing_evidence': classify(summary['read_processing'])})
            save()
            print(json.dumps({'fixture': fixture.name, 'mode': mode, 'status': record['status']}), flush=True)
    report['complete'] = True
    report['all_processes_successful'] = all(r['status'] == 'passed' for r in report['runs'])
    save()
    return 0 if report['all_processes_successful'] else 1


if __name__ == '__main__':
    sys.exit(main())
