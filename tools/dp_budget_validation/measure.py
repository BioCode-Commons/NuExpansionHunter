#!/usr/bin/env python3
"""Alternating isolated-process wall/CPU/peak-RSS measurements on Linux and macOS.

Uses posix_spawn+wait4, avoiding inherited Python RSS and /usr/bin/time's
sandbox-blocked sysctl. Linux ru_maxrss is normalized from KiB. Correctness compares complete stdout for each pair.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import statistics
import sys
import time


def files_identical(left, right):
    """Compare every byte, without retaining complete output files in memory."""
    left, right = Path(left), Path(right)
    if left.stat().st_size != right.stat().st_size:
        return False
    with left.open('rb') as a, right.open('rb') as b:
        while True:
            first = a.read(1024 * 1024)
            second = b.read(1024 * 1024)
            if first != second:
                return False
            if not first:
                return True


def measure(command, prefix):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    out, err = str(prefix) + '.stdout', str(prefix) + '.stderr'
    actions = [(os.POSIX_SPAWN_OPEN, 1, out, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644),
               (os.POSIX_SPAWN_OPEN, 2, err, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)]
    start = time.perf_counter()
    pid = os.posix_spawn(str(command[0]), [str(x) for x in command], os.environ, file_actions=actions)
    _, status, usage = os.wait4(pid, 0)
    elapsed = time.perf_counter() - start
    result = dict(command=[str(x) for x in command], wall_seconds=elapsed,
                  user_seconds=usage.ru_utime, system_seconds=usage.ru_stime,
                  peak_rss_bytes=usage.ru_maxrss * (1024 if sys.platform.startswith('linux') else 1), wait_status=status,
                  stdout_sha256=hashlib.sha256(Path(out).read_bytes()).hexdigest(),
                  stdout_path=out, stderr_path=err)
    Path(str(prefix) + '.json').write_text(json.dumps(result, indent=2) + '\n')
    if status:
        raise RuntimeError(f'command failed ({status}): {command}; see {err}')
    return result


def cpu_summary(pairs):
    before = [p['baseline']['user_seconds'] + p['baseline']['system_seconds'] for p in pairs]
    after = [p['candidate']['user_seconds'] + p['candidate']['system_seconds'] for p in pairs]
    ratios = [c/b for b,c in zip(before, after)]
    rng = random.Random(20261110)
    samples = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(10000))
    ratio = statistics.median(ratios)
    return dict(baseline_median_cpu_seconds=statistics.median(before),
                candidate_median_cpu_seconds=statistics.median(after),
                paired_cpu_ratio=ratio, paired_cpu_reduction_percent=(1-ratio)*100,
                cpu_ratio_bootstrap_95=[samples[249], samples[9749]],
                cpu_time_gate_pass=len(pairs)>=9 and min(statistics.median(before),statistics.median(after))>=0.5
                    and ratio<0.98 and samples[9749]<1)


def summarize(pairs):
    ratios = [p['candidate']['wall_seconds'] / p['baseline']['wall_seconds'] for p in pairs]
    rng = random.Random(20261110)
    samples = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(10000))
    b = statistics.median(p['baseline']['wall_seconds'] for p in pairs)
    c = statistics.median(p['candidate']['wall_seconds'] for p in pairs)
    br = statistics.median(p['baseline']['peak_rss_bytes'] for p in pairs)
    cr = statistics.median(p['candidate']['peak_rss_bytes'] for p in pairs)
    brmax = max(p['baseline']['peak_rss_bytes'] for p in pairs)
    crmax = max(p['candidate']['peak_rss_bytes'] for p in pairs)
    median_rss_allowance = max(256*1024, br*0.01)
    max_rss_allowance = max(256*1024, brmax*0.01)
    ratio = statistics.median(ratios)
    upper = samples[9749]
    lower = samples[249]
    identity = all(p['baseline']['stdout_sha256'] == p['candidate']['stdout_sha256']
                   and files_identical(p['baseline']['stdout_path'], p['candidate']['stdout_path']) for p in pairs)
    artifact_identity = all(p.get('artifacts_identical', True) for p in pairs)
    return dict(cpu_summary(pairs), pairs=len(pairs), baseline_median_seconds=b, candidate_median_seconds=c,
                paired_median_ratio=ratio, paired_speedup_percent=(1-ratio)*100,
                median_ratio_bootstrap_95=[lower, upper], faster_pairs=sum(x<1 for x in ratios),
                baseline_median_peak_rss_bytes=br, candidate_median_peak_rss_bytes=cr,
                baseline_max_peak_rss_bytes=brmax, candidate_max_peak_rss_bytes=crmax,
                peak_rss_delta_bytes=cr-br, peak_rss_delta_percent=(cr/br-1)*100 if br else None,
                max_peak_rss_delta_bytes=crmax-brmax,
                stdout_identical=identity,
                artifacts_identical=artifact_identity,
                time_gate_pass=len(pairs) >= 9 and min(b,c) >= 0.5 and ratio < 0.98 and upper < 1,
                evidence_role='acceptance' if len(pairs) >= 9 else 'smoke_correctness_only',
                rss_noise_allowance_bytes=median_rss_allowance,
                rss_max_noise_allowance_bytes=max_rss_allowance,
                rss_median_gate_pass=cr <= br + median_rss_allowance,
                rss_max_gate_pass=crmax <= brmax + max_rss_allowance,
                rss_gate_pass=cr <= br + median_rss_allowance and crmax <= brmax + max_rss_allowance)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--candidate', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--pairs', type=int, default=9)
    parser.add_argument('--artifact', action='append', default=[],
                        help='Output overwritten by each process; retain and compare every byte outside timing')
    parser.add_argument('args', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    common = args.args[1:] if args.args[:1] == ['--'] else args.args
    root = Path(args.out).resolve()
    root.mkdir(parents=True, exist_ok=False)
    commands = {'baseline': [str(Path(args.baseline).resolve())] + common,
                'candidate': [str(Path(args.candidate).resolve())] + common}
    allpairs = []
    for pair in range(-1, args.pairs):
        order = ['baseline', 'candidate'] if pair % 2 == 0 else ['candidate', 'baseline']
        records = {'pair': pair, 'warmup': pair < 0, 'order': order}
        for name in order:
            for path in args.artifact:
                Path(path).unlink(missing_ok=True)
            records[name] = measure(commands[name], root / f'pair-{pair:02d}-{name}')
            records[name]['artifacts'] = []
            for index, path in enumerate(args.artifact):
                target = root / f'pair-{pair:02d}-{name}-artifact-{index}'
                shutil.copyfile(path, target)
                records[name]['artifacts'].append(str(target))
        if (records['baseline']['stdout_sha256'] != records['candidate']['stdout_sha256']
                or not files_identical(records['baseline']['stdout_path'], records['candidate']['stdout_path'])):
            raise RuntimeError(f'Output differs at pair {pair}')
        records['artifacts_identical'] = all(files_identical(left, right) for left, right in zip(
            records['baseline']['artifacts'], records['candidate']['artifacts']))
        if not records['artifacts_identical']:
            raise RuntimeError(f'Artifact bytes differ at pair {pair}')
        if pair >= 0:
            allpairs.append(records)
        (root / 'pairs.json').write_text(json.dumps(allpairs, indent=2) + '\n')
        print(json.dumps({'pair': pair, 'baseline': records['baseline']['wall_seconds'],
                          'candidate': records['candidate']['wall_seconds']}), flush=True)
    summary = summarize(allpairs)
    (root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
