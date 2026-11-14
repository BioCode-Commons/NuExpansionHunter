#!/usr/bin/env python3
"""Serial paired streaming benchmarks with exact full-output comparisons.

Counts 16 and 32 are compute workers, not simulated physical cores. Do not
overlap this program with builds, fixture generation, or other benchmarks.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import re
import statistics
import subprocess
import sys

HELPERS = Path(__file__).resolve().parents[1] / 'streaming_validation'
sys.path.insert(0, str(HELPERS))
from compare_outputs import compare_prefixes
from measure import measure


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def capture(command, cwd=None):
    result = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    return {'command': command, 'returncode': result.returncode,
            'stdout': result.stdout.strip(), 'stderr': result.stderr.strip()}


def source_record(path):
    if path is None:
        return None
    path = path.resolve()
    diff = capture(['git', 'diff', '--binary', 'HEAD'], path)
    return {'path': str(path), 'head': capture(['git', 'rev-parse', 'HEAD'], path),
            'status': capture(['git', 'status', '--short'], path),
            'tracked_diff_sha256': hashlib.sha256(diff['stdout'].encode()).hexdigest(),
            'tracked_diff_returncode': diff['returncode']}


def summary(pairs, max_rss_ratio=1.05):
    out = {'pairs': len(pairs), 'all_outputs_identical': all(p['outputs_identical'] for p in pairs)}
    for metric, getter in (
        ('wall', lambda r: r['wall_seconds']),
        ('read_processing', lambda r: r['phases_seconds']['read_processing_seconds']),
        ('read_and_output', lambda r: r['phases_seconds']['read_processing_seconds']
                                    + r['phases_seconds']['output_seconds']),
        ('cpu', lambda r: r['user_seconds'] + r['system_seconds']),
        ('rss', lambda r: r['peak_rss_bytes'])
    ):
        before = [getter(p['baseline']) for p in pairs]
        after = [getter(p['candidate']) for p in pairs]
        ratios = [c / b for b, c in zip(before, after)]
        rng = random.Random(20261112)
        boot = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(10000))
        ratio = statistics.median(ratios)
        out[metric] = {'baseline_median': statistics.median(before),
                       'candidate_median': statistics.median(after),
                       'baseline_max': max(before), 'candidate_max': max(after),
                       'paired_median_ratio': ratio, 'paired_reduction_percent': (1-ratio)*100,
                       'paired_ratio_bootstrap_95': [boot[249], boot[9749]],
                       'candidate_lower_pairs': sum(c < b for b, c in zip(before, after))}
    out['phase_median_seconds'] = {}
    for side in ('baseline', 'candidate'):
        keys = set().union(*(p[side]['phases_seconds'] for p in pairs))
        out['phase_median_seconds'][side] = {key: statistics.median(
            p[side]['phases_seconds'][key] for p in pairs if key in p[side]['phases_seconds']) for key in sorted(keys)}
    out['interpretation'] = ('Descriptive paired speed observations; bootstrap intervals with few pairs are not '
                             'a hardware-general speed guarantee. The RAM gate applies only to these runs.')
    rss = out['rss']
    out['ram_gate'] = {
        'maximum_allowed_ratio': max_rss_ratio,
        'median_rss_ratio': rss['candidate_median'] / rss['baseline_median'],
        'maximum_rss_ratio': rss['candidate_max'] / rss['baseline_max'],
        'passed': (rss['candidate_median'] <= max_rss_ratio * rss['baseline_median']
                   and rss['candidate_max'] <= max_rss_ratio * rss['baseline_max']),
        'scope': 'Observed process peak RSS across paired runs; not a bound for other workloads.'}
    return out


def main():
    def io_workers(value):
        if value == 'auto':
            return value
        count = int(value)
        if count < 0:
            raise argparse.ArgumentTypeError('I/O workers must be nonnegative or auto')
        return count

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--baseline-source', type=Path)
    parser.add_argument('--candidate-source', type=Path)
    parser.add_argument('--input', type=Path, action='append', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--threads', type=int, nargs='+', default=[16, 32])
    parser.add_argument('--preprocess', type=int, nargs='+', default=[0])
    parser.add_argument('--decode', type=io_workers, default=2)
    parser.add_argument('--encode', type=io_workers, default=1)
    parser.add_argument('--pairs', type=int, default=5)
    parser.add_argument('--uncached', action='store_true')
    parser.add_argument('--baseline-cache', type=Path)
    parser.add_argument('--candidate-cache', type=Path)
    parser.add_argument('--provenance-file', type=Path, action='append', default=[])
    parser.add_argument('--cpu-list', help='Linux CPU IDs (comma-separated); pin this harness and all children to distinct physical cores')
    parser.add_argument('--max-rss-ratio', type=float, default=1.05,
                        help='RAM gate for both median and maximum process peak RSS (default 1.05)')
    args = parser.parse_args()
    if not 1 <= args.max_rss_ratio < float('inf'):
        parser.error('--max-rss-ratio must be finite and at least 1')
    affinity = None
    if args.cpu_list is not None:
        if not sys.platform.startswith('linux') or not hasattr(os, 'sched_setaffinity'):
            parser.error('--cpu-list requires Linux; worker counts alone do not enforce a CPU limit')
        try:
            cpus = [int(cpu) for cpu in args.cpu_list.split(',')]
            if not cpus or len(cpus) != len(set(cpus)):
                raise ValueError('CPU IDs must be nonempty and unique')
            available = os.sched_getaffinity(0)
            if not set(cpus) <= available:
                raise ValueError('Requested CPUs are outside the current allowed affinity')
            topology = []
            for cpu in cpus:
                directory = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
                topology.append((int((directory/'physical_package_id').read_text()),
                                 int((directory/'core_id').read_text())))
            if len(set(topology)) != len(cpus):
                raise ValueError('CPU list includes SMT siblings; select distinct physical cores')
            os.sched_setaffinity(0, cpus)
            if os.sched_getaffinity(0) != set(cpus):
                raise ValueError('CPU affinity verification failed')
            affinity = {'cpu_ids': cpus, 'package_core_ids': topology,
                        'physical_cores': len(cpus), 'inherited_by_all_children': True}
        except (OSError, ValueError) as error:
            parser.error(str(error))
    if min(args.threads + [args.pairs]) < 1 or min(args.preprocess) < 0:
        parser.error('Invalid worker/pair counts')
    supplied_caches = {'baseline': args.baseline_cache, 'candidate': args.candidate_cache}
    if any(supplied_caches.values()) and not all(supplied_caches.values()):
        parser.error('Provide both --baseline-cache and --candidate-cache')
    if any(supplied_caches.values()) and args.uncached:
        parser.error('Supplied caches cannot be combined with --uncached')
    for side, cache in supplied_caches.items():
        if cache is not None:
            supplied_caches[side] = cache.resolve()
            if not supplied_caches[side].is_file():
                parser.error('Missing supplied cache: ' + str(supplied_caches[side]))
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binaries = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
    inputs = [p.resolve() for p in args.input]
    if len({p.name for p in inputs}) != len(inputs):
        parser.error('Input directory basenames must be unique')
    for path in binaries.values():
        if not path.is_file():
            parser.error('Missing binary: ' + str(path))
    for directory in inputs:
        for name in ('reference.fa', 'reference.fa.fai', 'catalog.json', 'reads.bam'):
            if not (directory / name).is_file():
                parser.error('Missing input: ' + str(directory / name))
    environment = {'platform': platform.platform(), 'machine': platform.machine(),
                   'python': sys.version, 'logical_cpus': os.cpu_count(),
                   'thread_scope': 'Compute worker counts; no physical-core simulation or CPU affinity applied.'}
    environment['cpu_affinity'] = affinity
    if affinity is not None:
        environment['thread_scope'] = 'All subprocess threads restricted to the recorded distinct physical cores.'
    if sys.platform == 'darwin':
        environment['cpu'] = capture(['/usr/sbin/sysctl', 'hw.physicalcpu', 'hw.logicalcpu', 'machdep.cpu.brand_string'])
    result = {'scope': 'Synthetic streaming measurements only; not real WGS. See recorded thread budget and CPU affinity.',
              'environment': environment,
              'parameters': {k: [str(v) for v in val] if isinstance(val, list) else str(val) if isinstance(val, Path) else val
                             for k, val in vars(args).items()},
              'binaries': {k: {'path': str(v), 'bytes': v.stat().st_size, 'sha256': sha256(v)} for k, v in binaries.items()},
              'sources': {'baseline': source_record(args.baseline_source), 'candidate': source_record(args.candidate_source)},
              'harness_sha256': {str(p): sha256(p) for p in [Path(__file__), HELPERS/'measure.py', HELPERS/'compare_outputs.py']},
              'provenance_files': {str(p.resolve()): sha256(p) for p in args.provenance_file},
              'inputs': {}, 'preparations': {}, 'configurations': [],
              'output_policy': 'Exact VCF/JSON bytes and BAM header/full-record-byte multiset, preserving duplicate multiplicity; only BAM record order ignored.',
              'timing_policy': 'One excluded warmup then alternating AB/BA pairs. Serial subprocesses; hashes, preparation, comparisons outside timings. Cache loading stays timed. No OS page-cache clearing.'}

    def save():
        (root / 'results.json').write_text(json.dumps(result, indent=2) + '\n')

    def run(command, prefix):
        record = measure(command, prefix)
        log = Path(record['stdout_path']).read_text() + Path(record['stderr_path']).read_text()
        record['phases_seconds'] = {name: float(value) for name, value in
            re.findall(r'Phase (\w+_seconds)=([0-9.eE+-]+)', log)}
        pools = re.findall(r'Thread pools: compute=(\d+), decode=(\d+), encode=(\d+)', log)
        if pools:
            record['resolved_thread_pools'] = dict(zip(('compute', 'decode', 'encode'), map(int, pools[-1])))
        telemetry = re.findall(r'Read preprocessing workers=(\d+) peak_worker_cache_bytes=(\d+) synchronous_large_records=(\d+)', log)
        if telemetry:
            record['preprocessing'] = dict(zip(('workers', 'peak_worker_cache_bytes', 'synchronous_large_records'), map(int, telemetry[-1])))
        routing = re.findall(r'Deferred decoding selected_records=(\d+) routed_pairs=(\d+) skipped_pairs=(\d+)', log)
        if routing:
            record['deferred_decoding'] = dict(zip(('selected_records', 'routed_pairs', 'skipped_pairs'), map(int, routing[-1])))
        assistance = re.findall(r'Streaming reader assisted_locus_jobs=(\d+)', log)
        if assistance:
            record['reader_assisted_locus_jobs'] = int(assistance[-1])
        return record

    for inputs_dir in inputs:
        fixture = inputs_dir.name
        result['inputs'][fixture] = {'directory': str(inputs_dir),
            'files': {name: {'sha256': sha256(inputs_dir/name), 'bytes': (inputs_dir/name).stat().st_size}
                      for name in ('reference.fa', 'reference.fa.fai', 'catalog.json', 'reads.bam')},
            'manifest': json.loads((inputs_dir/'manifest.json').read_text()) if (inputs_dir/'manifest.json').is_file() else None}
        common = ['--reads', str(inputs_dir/'reads.bam'), '--reference', str(inputs_dir/'reference.fa'),
                  '--variant-catalog', str(inputs_dir/'catalog.json'), '--analysis-mode', 'streaming',
                  '--decode-threads', str(args.decode), '--encode-threads', str(args.encode),
                  '--region-extension-length', '1000', '--log-level', 'info']
        cache_args = {side: [] for side in binaries}
        if not args.uncached:
            directory = root / fixture / 'preparation'
            directory.mkdir(parents=True)
            result['preparations'][fixture] = {}
            for side, binary in binaries.items():
                if supplied_caches[side] is not None:
                    cache = supplied_caches[side]
                    result['preparations'][fixture][side] = {
                        'reused_cache': True,
                        'cache': {'path': str(cache), 'sha256': sha256(cache), 'bytes': cache.stat().st_size},
                        'compatibility_policy': 'The application validates this archive on every timed run; the harness does not infer compatibility.'}
                    cache_args[side] = ['--catalog-cache', str(cache)]
                    save()
                    continue
                cache = directory / (side + '.catalog')
                command = [str(binary)] + common + ['--threads', '1', '--read-decode-threads', '0',
                    '--output-prefix', str(directory/side), '--catalog-cache', str(cache), '--prepare-catalog-cache']
                result['preparations'][fixture][side] = run(command, directory/(side+'-measure'))
                result['preparations'][fixture][side]['cache'] = {'path': str(cache), 'sha256': sha256(cache), 'bytes': cache.stat().st_size}
                cache_args[side] = ['--catalog-cache', str(cache)]
                save()
        reference_prefix = None
        for threads in args.threads:
            for preprocess in args.preprocess:
                config_dir = root / fixture / f't{threads}-p{preprocess}'
                config_dir.mkdir(parents=True)
                config = {'fixture': fixture, 'compute_workers': threads, 'preprocessing_workers': preprocess,
                          'hts_decode_workers': args.decode, 'hts_encode_workers': args.encode, 'pairs': [], 'warmup': None}
                result['configurations'].append(config)
                for pair in range(-1, args.pairs):
                    order = ['baseline', 'candidate'] if pair % 2 == 0 else ['candidate', 'baseline']
                    record = {'pair': pair, 'order': order, 'excluded_warmup': pair < 0}
                    prefixes = {}
                    for side in order:
                        directory = config_dir / f'pair-{pair:02d}-{side}'
                        directory.mkdir()
                        prefix = directory / 'result'
                        prefixes[side] = prefix
                        command = [str(binaries[side])] + common + ['--threads', str(threads),
                            '--read-decode-threads', str(preprocess)] + cache_args[side] + ['--output-prefix', str(prefix)]
                        record[side] = run(command, directory/'measure')
                        if record[side].get('preprocessing', {}).get('workers') != preprocess:
                            raise RuntimeError('Missing or unexpected preprocessing worker telemetry')
                        record[side]['outputs'] = {suffix: {'bytes': Path(str(prefix)+suffix).stat().st_size,
                            'sha256': sha256(str(prefix)+suffix)} for suffix in ('.vcf', '.json', '_realigned.bam')}
                    comparison = compare_prefixes(prefixes['baseline'], prefixes['candidate'], root/'comparison-scratch')
                    record['comparison'] = comparison
                    record['outputs_identical'] = comparison['outputs_identical_under_user_policy']
                    if pair < 0:
                        if reference_prefix is None:
                            reference_prefix = prefixes['baseline']
                        else:
                            record['cross_configuration_comparison'] = compare_prefixes(reference_prefix, prefixes['baseline'], root/'comparison-scratch')
                            record['outputs_identical'] &= record['cross_configuration_comparison']['outputs_identical_under_user_policy']
                        config['warmup'] = record
                    else:
                        config['pairs'].append(record)
                    save()
                    print(json.dumps({'fixture': fixture, 'threads': threads, 'preprocess': preprocess, 'pair': pair,
                        'baseline_wall': record['baseline']['wall_seconds'], 'candidate_wall': record['candidate']['wall_seconds'],
                        'outputs_identical': record['outputs_identical']}), flush=True)
                    if not record['outputs_identical']:
                        raise RuntimeError('Output mismatch; all evidence retained in ' + str(root/'results.json'))
                config['summary'] = summary(config['pairs'], args.max_rss_ratio)
                save()
                print(json.dumps({'fixture': fixture, 'threads': threads, 'preprocess': preprocess,
                                  'summary': config['summary']}), flush=True)
    result['complete'] = True
    save()


if __name__ == '__main__':
    main()
