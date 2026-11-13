#!/usr/bin/env python3
"""Serial three-way streaming benchmark; no concurrent builds or probes.

All variants use two compute workers, one reader and one BAM writer. This is a
four-active-thread budget, not physical-core affinity on macOS. One warmup is
excluded in advance, then a balanced order measures every variant each cycle.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import random
import re
import statistics
import sys

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS / 'frontend_validation'))
from run_matrix import capture, sha256, source_record, summary
from compare_outputs import compare_prefixes
from measure import measure


def paired_metric(cycles, baseline, metric):
    before = [metric(c['runs'][baseline]) for c in cycles]
    after = [metric(c['runs']['candidate']) for c in cycles]
    ratios = [a / b for a, b in zip(after, before)]
    rng = random.Random(20261113)
    boot = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(10000))
    return {'baseline_median': statistics.median(before), 'candidate_median': statistics.median(after),
            'paired_median_ratio': statistics.median(ratios), 'paired_ratio_bootstrap_95': [boot[249], boot[9749]],
            'candidate_lower_pairs': sum(a < b for a, b in zip(after, before))}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ('pr', 'happy', 'candidate'):
        ap.add_argument('--' + name, type=Path, required=True)
        ap.add_argument('--' + name + '-source', type=Path)
    ap.add_argument('--input', type=Path, action='append', required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--cycles', type=int, default=9)
    ap.add_argument('--target-speedup', type=float, default=1.2)
    ap.add_argument('--max-rss-ratio', type=float, default=1.05)
    args = ap.parse_args()
    if args.cycles < 1 or args.target_speedup <= 1 or not 1 <= args.max_rss_ratio < float('inf'):
        ap.error('Invalid cycles/target/RAM gate')
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    binaries = {k: getattr(args, k).resolve() for k in ('pr', 'happy', 'candidate')}
    inputs = [p.resolve() for p in args.input]
    if len({p.name for p in inputs}) != len(inputs):
        ap.error('Input directory basenames must be distinct')
    for binary in binaries.values():
        if not binary.is_file():
            ap.error('Missing binary: ' + str(binary))
    report = {
        'complete': False, 'scope': 'Synthetic 150bp streaming measurements; no real-WGS or hardware-general guarantee.',
        'configuration': {'compute_workers': 2, 'reader': 1, 'writer': 1, 'read_decode': 0, 'hts_decode': 0,
                          'hts_encode': 0, 'physical_core_affinity': None},
        'environment': {'platform': platform.platform(), 'machine': platform.machine(), 'logical_cpus': os.cpu_count()},
        'parameters': {'cycles': args.cycles, 'target_speedup': args.target_speedup, 'maximum_rss_ratio': args.max_rss_ratio},
        'binaries': {k: {'path': str(v), 'bytes': v.stat().st_size, 'sha256': sha256(v)} for k, v in binaries.items()},
        'sources': {k: source_record(getattr(args, k + '_source')) for k in binaries},
        'harness_sha256': {str(p): sha256(p) for p in (Path(__file__), TOOLS/'frontend_validation/run_matrix.py',
                            TOOLS/'streaming_validation/measure.py', TOOLS/'streaming_validation/compare_outputs.py')},
        'timing_policy': 'Serial processes, one excluded warmup triple, balanced cyclic execution order; preparation, hashes and output comparisons are outside measured runs. No page-cache clearing.',
        'output_policy': 'Exact VCF/JSON bytes and complete BAM header/raw-record multiset, retaining duplicates; record order and compression bytes may differ.',
        'inputs': {}, 'fixtures': []}
    if sys.platform == 'darwin':
        report['environment']['cpu'] = capture(['/usr/sbin/sysctl', 'hw.physicalcpu', 'hw.logicalcpu', 'machdep.cpu.brand_string'])
    orders = [('pr', 'happy', 'candidate'), ('happy', 'candidate', 'pr'), ('candidate', 'pr', 'happy'),
              ('candidate', 'happy', 'pr'), ('pr', 'candidate', 'happy'), ('happy', 'pr', 'candidate')]

    def save():
        (out/'results.json').write_text(json.dumps(report, indent=2) + '\n')

    def run(command, prefix):
        record = measure(command, prefix)
        log = Path(record['stdout_path']).read_text() + Path(record['stderr_path']).read_text()
        record['phases_seconds'] = {k: float(v) for k, v in re.findall(r'Phase (\w+_seconds)=([0-9.eE+-]+)', log)}
        record['reader_assisted_locus_jobs'] = [int(v) for v in re.findall(r'assisted_locus_jobs=(\d+)', log)]
        return record

    for inp in inputs:
        report['inputs'][inp.name] = {'path': str(inp), 'files': {name: {'bytes': (inp/name).stat().st_size,
            'sha256': sha256(inp/name)} for name in ('reference.fa', 'reference.fa.fai', 'catalog.json', 'reads.bam')},
            'manifest': json.loads((inp/'manifest.json').read_text()) if (inp/'manifest.json').is_file() else None}
        common = ['--reads', str(inp/'reads.bam'), '--reference', str(inp/'reference.fa'), '--variant-catalog',
                  str(inp/'catalog.json'), '--analysis-mode', 'streaming', '--threads', '2', '--read-decode-threads', '0',
                  '--decode-threads', '0', '--encode-threads', '0', '--region-extension-length', '1000', '--log-level', 'info']
        fixture_dir = out/inp.name; fixture_dir.mkdir()
        fixture = {'fixture': inp.name, 'preparations': {}, 'warmup': None, 'cycles': []}
        report['fixtures'].append(fixture); save()
        caches = {}
        for side, binary in binaries.items():
            caches[side] = fixture_dir/(side+'.catalog')
            command = [str(binary)] + common + ['--catalog-cache', str(caches[side]), '--prepare-catalog-cache',
                                              '--output-prefix', str(fixture_dir/('prepare-'+side))]
            fixture['preparations'][side] = run(command, fixture_dir/('prepare-'+side+'-measure'))
            fixture['preparations'][side]['cache_sha256'] = sha256(caches[side]); save()
        for cycle in range(-1, args.cycles):
            order = orders[cycle % len(orders)] if cycle >= 0 else ('candidate', 'happy', 'pr')
            entry = {'cycle': cycle, 'excluded_warmup': cycle < 0, 'order': list(order), 'runs': {}, 'comparisons': {}}
            prefixes = {}
            for side in order:
                directory = fixture_dir/(f'cycle-{cycle:02d}-'+side); directory.mkdir()
                prefixes[side] = directory/'result'
                command = [str(binaries[side])] + common + ['--catalog-cache', str(caches[side]), '--output-prefix', str(prefixes[side])]
                record = run(command, directory/'measure')
                required = ('read_processing_seconds', 'genotyping_seconds', 'output_seconds')
                if not all(k in record['phases_seconds'] for k in required):
                    raise RuntimeError('Missing phase telemetry: '+str(directory))
                record['outputs'] = {suffix: {'bytes': Path(str(prefixes[side])+suffix).stat().st_size,
                    'sha256': sha256(str(prefixes[side])+suffix)} for suffix in ('.vcf', '.json', '_realigned.bam')}
                entry['runs'][side] = record
            for baseline in ('pr', 'happy'):
                entry['comparisons'][baseline] = compare_prefixes(prefixes[baseline], prefixes['candidate'], out/'comparison-scratch')
            entry['outputs_identical'] = all(c['outputs_identical_under_user_policy'] for c in entry['comparisons'].values())
            if cycle < 0:
                fixture['warmup'] = entry
            else:
                fixture['cycles'].append(entry)
            save()
            print(json.dumps({'fixture': inp.name, 'cycle': cycle, 'outputs_identical': entry['outputs_identical'],
                'read_seconds': {k: v['phases_seconds']['read_processing_seconds'] for k, v in entry['runs'].items()}}), flush=True)
            if not entry['outputs_identical']:
                raise RuntimeError('Output mismatch; evidence retained in '+str(out))
        fixture['comparisons'] = {}
        for baseline in ('pr', 'happy'):
            pairs = [{'baseline': c['runs'][baseline], 'candidate': c['runs']['candidate'], 'outputs_identical': c['outputs_identical']} for c in fixture['cycles']]
            stats = summary(pairs, args.max_rss_ratio)
            stats['post_preprocessing'] = paired_metric(fixture['cycles'], baseline, lambda r: sum(r['phases_seconds'][k]
                for k in ('read_processing_seconds', 'genotyping_seconds', 'output_seconds')))
            stats['target_gate'] = {'minimum_speedup': args.target_speedup,
                'read_processing_passed': stats['read_processing']['paired_median_ratio'] <= 1/args.target_speedup,
                'read_and_output_passed': stats['read_and_output']['paired_median_ratio'] <= 1/args.target_speedup,
                'formal_cycle_count': args.cycles >= 9}
            fixture['comparisons'][baseline] = stats
        save()
        print(json.dumps({'fixture': inp.name, 'speedups': {b: 1/s['read_processing']['paired_median_ratio']
            for b, s in fixture['comparisons'].items()}, 'ram_gates': {b: s['ram_gate']['passed'] for b, s in fixture['comparisons'].items()}}), flush=True)
    report['complete'] = True
    report['all_outputs_identical'] = all(c['outputs_identical'] for f in report['fixtures'] for c in [f['warmup']]+f['cycles'])
    comparisons = [s for f in report['fixtures'] for s in f['comparisons'].values()]
    report['all_read_processing_targets_met'] = all(s['target_gate']['read_processing_passed'] for s in comparisons)
    report['all_read_and_output_targets_met'] = all(s['target_gate']['read_and_output_passed'] for s in comparisons)
    report['all_ram_gates_passed'] = all(s['ram_gate']['passed'] for s in comparisons)
    report['all_read_speedup_confidence_lower_bounds_meet_target'] = all(
        s['read_processing']['paired_ratio_bootstrap_95'][1] <= 1 / args.target_speedup for s in comparisons)
    save()


if __name__ == '__main__':
    main()
