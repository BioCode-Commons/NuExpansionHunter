#!/usr/bin/env python3
"""Same-host, paired streaming/cache measurements with complete output checks.

Python 3.8+, Linux or macOS. Run alone on the reserved host: no builds, fixture
generation or other benchmarks. Worker counts are not process-wide CPU caps.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import signal
import statistics
import subprocess
import sys
import traceback

HELPERS = Path(__file__).resolve().parents[1] / 'streaming_validation'
sys.path.insert(0, str(HELPERS))
from compare_outputs import compare_prefixes
from measure import measure

SIDES = ('baseline', 'candidate')
SUFFIXES = ('.vcf', '.json', '_realigned.bam')


def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def file_record(path):
    path = Path(path).resolve()
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha256(path)}


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def capture(command, cwd=None):
    try:
        result = subprocess.run(command, cwd=cwd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, check=False)
        return {'command': command, 'returncode': result.returncode,
                'stdout': result.stdout.strip(), 'stderr': result.stderr.strip()}
    except OSError as error:
        return {'command': command, 'error': str(error)}


def source_record(path, destination):
    if path is None:
        return None
    path = path.resolve()
    head = capture(['git', 'rev-parse', 'HEAD'], path)
    status = capture(['git', 'status', '--short'], path)
    diff = capture(['git', 'diff', '--binary', 'HEAD'], path)
    destination.write_text(diff.get('stdout', '') + '\n')
    untracked = capture(['git', 'ls-files', '--others', '--exclude-standard'], path)
    return {'path': str(path), 'head': head, 'status': status,
            'tracked_diff': file_record(destination),
            'tracked_diff_returncode': diff.get('returncode'),
            'untracked_files': [file_record(path / name) for name in untracked.get('stdout', '').splitlines()
                                if (path / name).is_file()]}


def read_text(path):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def cpu_topology(cpus):
    result = []
    for cpu in sorted(cpus):
        root = Path('/sys/devices/system/cpu/cpu{}'.format(cpu))
        result.append({'cpu': cpu,
                       'package': read_text(root / 'topology/physical_package_id'),
                       'core': read_text(root / 'topology/core_id'),
                       'thread_siblings': read_text(root / 'topology/thread_siblings_list'),
                       'numa_nodes': sorted(p.name for p in root.glob('node[0-9]*')),
                       'governor': read_text(root / 'cpufreq/scaling_governor'),
                       'llc_shared_cpu_lists': [read_text(p / 'shared_cpu_list')
                                                for p in sorted((root / 'cache').glob('index*'))
                                                if read_text(p / 'level') == '3']})
    return result


def parse_cpu_list(encoding):
    result = []
    for part in encoding.split(','):
        if not re.fullmatch(r'\d+(?:-\d+)?', part):
            raise ValueError('CPU IDs must be comma-separated integers or ascending ranges')
        ends = [int(x) for x in part.split('-')]
        if len(ends) == 2 and ends[1] < ends[0]:
            raise ValueError('CPU ranges must be ascending')
        result.extend(range(ends[0], ends[-1] + 1))
    if len(result) != len(set(result)):
        raise ValueError('CPU IDs must be unique')
    return sorted(result)


def configure_affinity(encoding, maximum_workers):
    if not sys.platform.startswith('linux'):
        if encoding is not None:
            raise ValueError('--cpu-list is supported only on Linux')
        return {'applied': False, 'reason': 'No portable CPU pinning on this platform'}
    available = sorted(os.sched_getaffinity(0))
    selected = parse_cpu_list(encoding) if encoding is not None else available
    topology = cpu_topology(selected)
    if encoding is not None:
        if not set(selected) <= set(available):
            raise ValueError('Requested CPUs are outside the current allocation/affinity')
        if len(selected) < maximum_workers:
            raise ValueError('--cpu-list must contain at least max(--threads) distinct physical cores')
        if any(row['package'] is None or row['core'] is None for row in topology):
            raise ValueError('Cannot validate physical cores from Linux sysfs')
        identities = {(row['package'], row['core']) for row in topology}
        if len(identities) != len(selected):
            raise ValueError('--cpu-list contains SMT siblings; choose one CPU per physical core')
        os.sched_setaffinity(0, selected)
    if sorted(os.sched_getaffinity(0)) != selected:
        raise ValueError('CPU affinity did not match the requested mask')
    return {'applied': encoding is not None, 'original_allowed_cpus': available,
            'effective_cpu_ids': selected, 'topology': topology,
            'same_mask_for_all_thread_counts_and_children': True,
            'physical_core_validation': encoding is not None}


def parse_telemetry(record):
    # Raw logs are always retained, including unknown/new fields.
    log = '\n'.join(Path(record[key]).read_text(errors='replace')
                    for key in ('stdout_path', 'stderr_path') if key in record and Path(record[key]).exists())
    record['phases_seconds'] = {key: float(value) for key, value in
        re.findall(r'Phase (\w+_seconds)=([0-9.eE+-]+)', log)}
    pools = re.findall(r'Thread pools: compute=(\d+), decode=(\d+), encode=(\d+)', log)
    record['resolved_thread_pools'] = (dict(zip(('compute', 'decode', 'encode'), map(int, pools[-1])))
                                        if pools else None)
    for name, pattern in (('streaming_batches', r'Streaming batches: ([^\n]+)'),
                          ('pending_mates', r'Streaming pending mates: ([^\n]+)'),
                          ('bam_sink', r'BAM sink: ([^\n]+)'),
                          ('bamlet_queue', r'BAMlet queue[^\n]*')):
        matches = re.findall(pattern, log, flags=re.IGNORECASE)
        record[name] = ({key: float(value) for key, value in re.findall(r'(\w+)=([0-9.eE+-]+)', matches[-1])}
                        if matches else None)
    assisted = re.findall(r'Streaming reader assisted_locus_jobs=(\d+)', log)
    record['reader_assisted_locus_jobs'] = int(assisted[-1]) if assisted else None
    hints = re.findall(r'Streaming preprocessing hint=(\d+)', log)
    record['read_decode_compatibility_hint'] = int(hints[-1]) if hints else None


def measure_driver(request_path):
    request = json.loads(Path(request_path).read_text())
    try:
        measure(request['command'], request['prefix'])
        return 0
    except BaseException:
        Path(request['prefix'] + '.driver-error.txt').write_text(traceback.format_exc())
        return 1


def paired_metric(before, after):
    ratios = [c / b for b, c in zip(before, after) if b > 0]
    result = {'baseline_median': statistics.median(before),
              'candidate_median': statistics.median(after),
              'baseline_max': max(before), 'candidate_max': max(after),
              'paired_deltas': [c - b for b, c in zip(before, after)], 'paired_ratios': ratios}
    if ratios:
        rng = random.Random(20261114)
        samples = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(2000))
        result.update(paired_median_ratio=statistics.median(ratios),
                      paired_ratio_bootstrap_95=[samples[49], samples[1949]])
    return result


def summarize(data, extra_bytes_per_worker):
    summaries = {}
    complete_pairs = [pair for pair in data['pairs'] if not pair['warmup'] and pair['status'] == 'verified']
    for workers in data['parameters']['threads']:
        pairs = [pair for pair in complete_pairs if pair['threads'] == workers]
        if not pairs:
            continue
        result = {'verified_pairs': len(pairs), 'metrics': {}}
        for name, get in (('wall_seconds', lambda r: r['wall_seconds']),
                          ('cpu_seconds', lambda r: r['user_seconds'] + r['system_seconds']),
                          ('peak_rss_bytes', lambda r: r['peak_rss_bytes'])):
            result['metrics'][name] = paired_metric([get(p['runs']['baseline']) for p in pairs],
                                                    [get(p['runs']['candidate']) for p in pairs])
        keys = set.intersection(*(set(p['runs'][side]['phases_seconds']) for p in pairs for side in SIDES))
        result['phases_seconds'] = {key: paired_metric(
            [p['runs']['baseline']['phases_seconds'][key] for p in pairs],
            [p['runs']['candidate']['phases_seconds'][key] for p in pairs]) for key in sorted(keys)}
        allowance = extra_bytes_per_worker * workers
        rss = result['metrics']['peak_rss_bytes']
        result['memory_allowance'] = {
            'extra_bytes_per_worker': extra_bytes_per_worker, 'allowed_extra_bytes': allowance,
            'paired_actual_delta_bytes': rss['paired_deltas'],
            'difference_of_medians_bytes': rss['candidate_median'] - rss['baseline_median'],
            'difference_of_maxima_bytes': rss['candidate_max'] - rss['baseline_max'],
            'median_paired_delta_bytes': statistics.median(rss['paired_deltas']),
            'largest_paired_delta_bytes': max(rss['paired_deltas']),
            'paired_candidate_limits_bytes': [p['runs']['baseline']['peak_rss_bytes'] + allowance for p in pairs],
            'all_pairs_within_allowance': all(delta <= allowance for delta in rss['paired_deltas']),
            'median_within_allowance': rss['candidate_median'] <= rss['baseline_median'] + allowance,
            'maximum_within_allowance': rss['candidate_max'] <= rss['baseline_max'] + allowance}
        wall = result['metrics']['wall_seconds']
        result['observed_median_wall_no_slower'] = wall['paired_median_ratio'] <= 1
        result['wall_bootstrap_upper95_no_slower'] = wall['paired_ratio_bootstrap_95'][1] <= 1
        summaries[str(workers)] = result
    scaling = {}
    reference_workers = min(data['parameters']['threads'])
    for side in SIDES:
        reference = {p['round']: p['runs'][side] for p in complete_pairs if p['threads'] == reference_workers}
        scaling[side] = {}
        for workers in data['parameters']['threads']:
            pairs = [p for p in complete_pairs if p['threads'] == workers and p['round'] in reference]
            if not pairs:
                continue
            speedups = [reference[p['round']]['wall_seconds'] / p['runs'][side]['wall_seconds'] for p in pairs]
            scaling[side][str(workers)] = {'reference_compute_workers': reference_workers,
                'same_round_pairs': len(pairs), 'wall_speedup_by_round': speedups,
                'median_wall_speedup': statistics.median(speedups),
                'median_compute_scaling_efficiency': statistics.median(speedups) / (workers / reference_workers),
                'note': 'Same allocation/mask; compute workers are not physical-core or total-thread counts.'}
    expected = len(data['parameters']['threads']) * data['parameters']['pairs']
    all_verified = len(complete_pairs) == expected and all(p['status'] == 'verified' for p in data['pairs'])
    return {'by_compute_workers': summaries, 'within_version_scaling': scaling,
            'all_requested_pairs_verified': all_verified,
            'all_observed_median_wall_no_slower': all_verified and all(
                item['observed_median_wall_no_slower'] for item in summaries.values()),
            'all_observed_peak_rss_within_allowance': all_verified and all(
                item['memory_allowance']['all_pairs_within_allowance']
                and item['memory_allowance']['median_within_allowance']
                and item['memory_allowance']['maximum_within_allowance'] for item in summaries.values()),
            'interpretation': 'Observed timings/RSS only. Few-pair bootstrap intervals are descriptive, not a '
                              'hardware-general guarantee. Exit 0 confirms complete measurements and exact outputs, '
                              'not a speed or memory claim; review the explicit performance fields.'}


class Runner:
    def __init__(self, args):
        self.args = args
        self.root = args.out.resolve()
        self.root.mkdir(parents=True, exist_ok=False)
        self.data = {'status': 'initializing', 'started_utc': timestamp(), 'complete': False,
                     'baseline_comparison_label': args.baseline_label,
                     'expected_baseline_source_commit': 'b00650d' if args.baseline_label == 'RAMOPT b00650d' else None,
                     'parameters': {key: str(value) if isinstance(value, Path) else value
                                    for key, value in vars(args).items()},
                     'preparations': {}, 'pairs': [], 'errors': [],
                     'output_policy': 'Every VCF/JSON byte and complete BAM header/full-record-byte multiset, '
                                      'including duplicate multiplicity. Only BAM record order may differ.',
                     'timing_policy': 'Cache preparation, hashes, warmup, comparison and cleanup are excluded. '
                                      'Archive loading remains timed. No OS page-cache clearing. '
                                      'All application processes run serially on the same inherited CPU mask.'}
        self.affinity = None
        self.extra_bytes_per_worker = int(args.extra_cache_mb_per_worker * 1000000)
        self.save()

    def save(self):
        atomic_json(self.root / 'results.json', self.data)

    def check_affinity(self):
        if self.affinity and 'effective_cpu_ids' in self.affinity:
            if sorted(os.sched_getaffinity(0)) != self.affinity['effective_cpu_ids']:
                raise RuntimeError('The CPU allocation/affinity changed during the benchmark')

    def run_process(self, command, directory, record):
        directory.mkdir(parents=True)
        prefix = directory / 'measure'
        request = directory / 'request.json'
        atomic_json(request, {'command': command, 'prefix': str(prefix)})
        record.update(status='running', command=command, started_utc=timestamp())
        self.check_affinity()
        self.save()
        driver_command = [sys.executable, str(Path(__file__).resolve()), '--measure-one', str(request)]
        with (directory / 'driver.stdout').open('wb') as out, (directory / 'driver.stderr').open('wb') as err:
            process = subprocess.Popen(driver_command, stdout=out, stderr=err, start_new_session=True)
            try:
                returncode = process.wait()
            except BaseException:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                except ProcessLookupError:
                    pass
                record.update(status='interrupted', finished_utc=timestamp())
                self.save()
                raise
        measurement = Path(str(prefix) + '.json')
        if measurement.is_file():
            record.update(json.loads(measurement.read_text()))
            parse_telemetry(record)
        record.update(driver_returncode=returncode, finished_utc=timestamp(),
                      status='completed' if returncode == 0 and measurement.is_file() else 'failed')
        self.save()
        self.check_affinity()
        if record['status'] != 'completed':
            raise RuntimeError('Application/preparation failed; see ' + str(directory))

    def command(self, side, workers, prefix, cache=None, prepare=False):
        args = self.args
        result = [str(getattr(args, side).resolve()), '--reads', str(args.reads.resolve()),
                  '--reference', str(args.reference.resolve()), '--variant-catalog', str(args.variant_catalog.resolve()),
                  '--analysis-mode', 'streaming', '--threads', str(workers),
                  '--decode-threads', args.decode, '--encode-threads', args.encode,
                  '--read-decode-threads', '0', '--region-extension-length', '1000', '--log-level', 'info',
                  '--output-prefix', str(prefix)]
        if args.bamlet_queue_mb is not None:
            result += ['--bamlet-queue-mb', args.bamlet_queue_mb]
        if cache is not None:
            result += ['--catalog-cache', str(cache)]
        if prepare:
            result.append('--prepare-catalog-cache')
        return result

    def discard(self, pair, protected_prefix=None):
        if not self.args.discard_equal_output or pair['status'] != 'verified':
            return
        # Caller has already durably saved exact comparison and complete hashes.
        for record in pair['runs'].values():
            if record['output_prefix'] == protected_prefix:
                continue
            for item in record['outputs'].values():
                if not item.get('removed', False):
                    Path(item['path']).unlink()
                    item['removed'] = True
                    self.save()

    def pair(self, workers, round_number, caches, reference_prefix=None):
        warmup = round_number < 0
        order = list(SIDES if round_number % 2 == 0 else reversed(SIDES))
        pair = {'threads': workers, 'round': round_number, 'warmup': warmup,
                'order': order, 'status': 'running', 'runs': {}}
        self.data['pairs'].append(pair)
        self.save()
        directory = self.root / ('warmup' if warmup else 'round-{:02d}'.format(round_number)) / ('t{}'.format(workers))
        for side in order:
            prefix = directory / side / 'result'
            record = {'output_prefix': str(prefix)}
            pair['runs'][side] = record
            self.run_process(self.command(side, workers, prefix, caches[side]), directory / side, record)
            record['outputs'] = {suffix: file_record(str(prefix) + suffix) for suffix in SUFFIXES}
            self.save()
        pair['comparison'] = compare_prefixes(pair['runs']['baseline']['output_prefix'],
                                             pair['runs']['candidate']['output_prefix'],
                                             self.root / 'comparison-scratch', self.args.compare_chunk_mib * 1024 * 1024)
        identical = pair['comparison']['outputs_identical_under_user_policy']
        if reference_prefix is not None:
            pair['cross_configuration_comparison'] = compare_prefixes(reference_prefix,
                pair['runs']['baseline']['output_prefix'], self.root / 'comparison-scratch',
                self.args.compare_chunk_mib * 1024 * 1024)
            identical &= pair['cross_configuration_comparison']['outputs_identical_under_user_policy']
        pair['status'] = 'verified' if identical else 'mismatch'
        self.save()
        print(json.dumps({'round': round_number, 'warmup': warmup, 'threads': workers,
                          'status': pair['status'], 'wall_seconds': {s: pair['runs'][s]['wall_seconds'] for s in SIDES}}),
              flush=True)
        if not identical:
            raise RuntimeError('Complete output comparison failed; all files retained for this pair')
        return pair

    def execute(self):
        args = self.args
        try:
            self.affinity = configure_affinity(args.cpu_list, max(args.threads))
            self.data['environment'] = {'platform': platform.platform(), 'machine': platform.machine(),
                'hostname': platform.node(), 'python': sys.version, 'logical_cpus': os.cpu_count(),
                'cpu_affinity': self.affinity, 'lscpu': capture(['lscpu', '-J']) if sys.platform.startswith('linux') else None,
                'numa_policy': capture(['numactl', '--show']) if sys.platform.startswith('linux') else None,
                'thread_scope': '--threads is compute workers. Codec pools and coordinator roles are additional. '
                                'A CPU mask limits eligible CPUs, not the number of threads.'}
            self.save()
            self.data['binaries'] = {side: file_record(getattr(args, side)) for side in SIDES}
            self.data['sources'] = {side: source_record(getattr(args, side + '_source'), self.root / (side + '-source.diff'))
                                    for side in SIDES}
            inputs = [args.reads.resolve(), args.reference.resolve(), args.variant_catalog.resolve()]
            inputs += [Path(str(path) + suffix) for path, suffixes in
                       ((args.reference.resolve(), ('.fai', '.gzi')), (args.reads.resolve(), ('.bai', '.csi', '.crai')))
                       for suffix in suffixes if Path(str(path) + suffix).is_file()]
            self.data['inputs'] = [file_record(path) for path in inputs]
            self.data['harness_files'] = [file_record(path) for path in
                (Path(__file__), Path(__file__).with_name('README.md'), HELPERS / 'measure.py', HELPERS / 'compare_outputs.py')]
            self.data['status'] = 'preparing'
            self.save()
            caches = {}
            for side in SIDES:
                supplied = getattr(args, side + '_cache')
                record = {'reused_cache': supplied is not None}
                self.data['preparations'][side] = record
                if supplied is None:
                    cache = self.root / 'preparation' / (side + '.catalog')
                    prefix = self.root / 'preparation' / side / 'result'
                    self.run_process(self.command(side, 1, prefix, cache, True), prefix.parent, record)
                else:
                    cache = supplied.resolve()
                record['cache'] = file_record(cache)
                record['compatibility'] = 'The matching application validates its archive; the harness does not infer compatibility.'
                caches[side] = cache
                self.save()
            self.data['status'] = 'warming'
            reference_prefix = None
            reference_pair = None
            for workers in reversed(args.threads):
                pair = self.pair(workers, -1, caches, reference_prefix)
                if reference_prefix is None:
                    reference_prefix = pair['runs']['baseline']['output_prefix']
                    reference_pair = pair
                self.discard(pair, reference_prefix)
            if reference_pair is not None:
                self.discard(reference_pair)
            self.data['status'] = 'measuring'
            self.data['round_configuration_orders'] = []
            for round_number in range(args.pairs):
                offset = (round_number // 2) % len(args.threads)
                order = args.threads[offset:] + args.threads[:offset]
                if round_number % 2:
                    order.reverse()
                self.data['round_configuration_orders'].append(order)
                self.save()
                for workers in order:
                    pair = self.pair(workers, round_number, caches)
                    self.discard(pair)
                    self.data['summary'] = summarize(self.data, self.extra_bytes_per_worker)
                    self.save()
            self.data.update(status='complete', complete=True, finished_utc=timestamp())
            self.data['summary'] = summarize(self.data, self.extra_bytes_per_worker)
            self.save()
            atomic_json(self.root / 'summary.json', self.data['summary'])
            return 0
        except BaseException as error:
            self.data.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'failed', complete=False,
                             finished_utc=timestamp())
            self.data['errors'].append({'type': type(error).__name__, 'message': str(error), 'traceback': traceback.format_exc()})
            for pair in self.data['pairs']:
                if pair['status'] == 'running':
                    pair['status'] = self.data['status']
            self.data['summary'] = summarize(self.data, self.extra_bytes_per_worker)
            self.data['summary']['all_requested_pairs_verified'] = False
            self.data['summary']['all_observed_median_wall_no_slower'] = False
            self.data['summary']['all_observed_peak_rss_within_allowance'] = False
            self.save()
            atomic_json(self.root / 'summary.json', self.data['summary'])
            print('Benchmark incomplete: {}. Evidence: {}'.format(error, self.root / 'results.json'), file=sys.stderr)
            return 130 if isinstance(error, KeyboardInterrupt) else 1


def nonnegative_or_auto(value):
    if value != 'auto' and (not value.isascii() or not value.isdecimal()):
        raise argparse.ArgumentTypeError('Expected auto or a non-negative integer')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'candidate', 'reads', 'reference', 'variant-catalog', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    for name in ('baseline-cache', 'candidate-cache', 'baseline-source', 'candidate-source'):
        parser.add_argument('--' + name, type=Path)
    parser.add_argument('--baseline-label', default='RAMOPT b00650d',
                        help='Comparison metadata only; binary hashes and optional source provenance remain authoritative')
    parser.add_argument('--threads', type=int, nargs='+', default=[16, 24, 32])
    parser.add_argument('--pairs', type=int, default=5)
    parser.add_argument('--decode', type=nonnegative_or_auto, default='2')
    parser.add_argument('--encode', type=nonnegative_or_auto, default='1')
    parser.add_argument('--bamlet-queue-mb', type=nonnegative_or_auto,
                        help='Optional common override; omitted preserves each binary\'s queue default')
    parser.add_argument('--extra-cache-mb-per-worker', type=float, default=20,
                        help='Allowed candidate peak RSS increase per compute worker, in decimal MB')
    parser.add_argument('--cpu-list', help='Linux CPU IDs/ranges, e.g. 0-15,32-47; one CPU per physical core')
    parser.add_argument('--discard-equal-output', action='store_true',
                        help='Remove raw outputs only after full equality and hashes are saved; keep logs/caches')
    parser.add_argument('--compare-chunk-mib', type=int, default=16)
    args = parser.parse_args()
    if not sys.platform.startswith('linux') and sys.platform != 'darwin':
        parser.error('This runner requires Linux or macOS posix_spawn/wait4 RSS accounting')
    if min(args.threads + [args.pairs, args.compare_chunk_mib]) < 1 or len(args.threads) != len(set(args.threads)):
        parser.error('Threads must be positive and unique; pairs and compare chunk size must be positive')
    args.threads.sort()
    if not math.isfinite(args.extra_cache_mb_per_worker) or args.extra_cache_mb_per_worker < 0:
        parser.error('Extra cache allowance must be finite and nonnegative')
    if args.extra_cache_mb_per_worker > sys.maxsize / (max(args.threads) * 1000000):
        parser.error('Extra cache allowance is too large for a byte-count budget')
    for name in ('baseline', 'candidate', 'reads', 'reference', 'variant_catalog'):
        path = getattr(args, name)
        if not path.is_file():
            parser.error('Missing file: ' + str(path))
    for name in SIDES:
        if not os.access(getattr(args, name), os.X_OK):
            parser.error('Binary is not executable: ' + str(getattr(args, name)))
    if bool(args.baseline_cache) != bool(args.candidate_cache):
        parser.error('Provide both --baseline-cache and --candidate-cache, or neither')
    if args.baseline_cache:
        if not args.baseline_cache.is_file() or not args.candidate_cache.is_file():
            parser.error('Both supplied caches must exist')
        if args.baseline_cache.resolve() == args.candidate_cache.resolve():
            parser.error('Baseline and candidate must use separate cache paths')
    if args.out.exists():
        parser.error('--out must be a new directory')
    def interrupt(signum, frame):
        raise KeyboardInterrupt('Received signal {}'.format(signum))
    signal.signal(signal.SIGTERM, interrupt)
    return Runner(args).execute()


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--measure-one':
        sys.exit(measure_driver(sys.argv[2]))
    sys.exit(main())
