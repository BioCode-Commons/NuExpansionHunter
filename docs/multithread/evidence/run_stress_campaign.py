#!/usr/bin/env python3
"""Sequential, resumable frozen-binary stress campaign. Default: print plan only.

Run --execute only after fixture generation and every other benchmark finish.
No case is concurrent with another. Failed/interrupted cases require explicit
--retry-failed CASE and use a fresh attempt directory; previous evidence stays.
Before such a retry, confirm that no orphan application is still running.
"""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import traceback

SIDES = ('baseline', 'candidate')
SUFFIXES = ('.vcf', '.json', '_realigned.bam')


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def save(path, value):
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def load(path):
    return json.loads(Path(path).read_text())


def exact(comparison):
    return (all(comparison.get(name) is True for name in
                ('vcf_identical', 'json_identical', 'outputs_identical_under_user_policy'))
            and all(comparison.get('bam', {}).get(name) is True for name in
                    ('header_identical', 'record_multiset_identical')))


class Campaign:
    def __init__(self, plan, plan_path, options):
        self.plan, self.options = plan, options
        self.workspace = Path(plan['workspace_root']).resolve()
        self.root = self.resolve(options.out or plan['output_root'])
        self.runner = self.resolve(plan['runner']['path'])
        self.plan_path = plan_path.resolve()
        self.policy = {'plan': plan, 'cpu_list_override': options.cpu_list,
                       'allow_unavailable_process_scan': options.allow_unavailable_process_scan}
        self.fingerprint = hashlib.sha256(json.dumps(self.policy, sort_keys=True).encode()).hexdigest()
        self.state_path = self.root / 'campaign.json'
        self.state = None
        self.lock = None
        self.scan_history = []

    def resolve(self, path):
        path = Path(path)
        return path.resolve() if path.is_absolute() else (self.workspace / path).resolve()

    def command(self, case, out):
        policy = self.plan['policy']
        command = [sys.executable, str(self.runner), '--baseline', str(self.resolve(self.plan['binaries']['baseline']['path'])),
                   '--candidate', str(self.resolve(self.plan['binaries']['candidate']['path'])),
                   '--reads', str(self.resolve(case['inputs']['reads'])),
                   '--reference', str(self.resolve(case['inputs']['reference'])),
                   '--variant-catalog', str(self.resolve(case['inputs']['variant_catalog'])),
                   '--out', str(out), '--threads'] + list(map(str, case['threads']))
        command += ['--pairs', str(policy['pairs']), '--decode', str(case.get('decode', policy['decode'])),
                    '--encode', str(case.get('encode', policy['encode'])), '--compare-chunk-mib', str(policy['compare_chunk_mib']),
                    '--extra-cache-mb-per-worker', str(policy['extra_cache_mb_per_worker'])]
        if policy['discard_equal_output']:
            command.append('--discard-equal-output')
        if policy.get('bamlet_queue_mb') is not None:
            command += ['--bamlet-queue-mb', str(policy['bamlet_queue_mb'])]
        if self.options.cpu_list or policy.get('cpu_list'):
            command += ['--cpu-list', self.options.cpu_list or policy['cpu_list']]
        for side in SIDES:
            if case.get('caches'):
                command += ['--' + side + '-cache', str(self.resolve(case['caches'][side]['path']))]
        return command

    def persist(self):
        self.state['updated_utc'] = now()
        save(self.state_path, self.state)

    def host(self):
        return {'hostname': platform.node(), 'machine': platform.machine(), 'platform': platform.platform(),
                'effective_cpu_ids': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None}

    def require_no_other_run(self):
        # Includes a still-running benchmark harness during comparison/preparation
        # gaps, and possible orphan applications after an ungraceful termination.
        try:
            result = subprocess.run(['ps', '-axww', '-o', 'pid=,ppid=,command='], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, text=True, check=True)
        except (OSError, subprocess.CalledProcessError) as error:
            observation = {'time': now(), 'status': 'unavailable',
                           'reason': (getattr(error, 'stderr', None) or str(error)).strip(),
                           'explicit_override': self.options.allow_unavailable_process_scan,
                           'limitation': 'Cannot inspect other processes; exclusive campaign lock only. '
                                         'External applications may run. Machine exclusivity is not established.'}
            self.record_scan(observation)
            if not self.options.allow_unavailable_process_scan:
                raise RuntimeError('Process scan unavailable; refusing by default. If owned jobs are confirmed '
                                   'closed, explicitly use --allow-unavailable-process-scan. Reason: ' + observation['reason'])
            print('Process scan unavailable (explicit override); machine exclusivity is NOT established.',
                  file=sys.stderr, flush=True)
            return
        needles = [str(self.runner), self.plan['runner']['path']]
        for binary in self.plan['binaries'].values():
            needles += [str(self.resolve(binary['path'])), binary['path']]
        found = []
        for line in result.stdout.splitlines():
            parts = line.strip().split(None, 2)
            if len(parts) != 3 or int(parts[0]) in (os.getpid(), os.getppid()):
                continue
            if any(needle in parts[2] for needle in needles):
                found.append(line.strip())
        if found:
            self.record_scan({'time': now(), 'status': 'matching_processes_found', 'matches': found})
            raise RuntimeError('Another matching benchmark/application is active; refusing overlap:\n' + '\n'.join(found))
        self.record_scan({'time': now(), 'status': 'no_matching_processes',
                          'limitation': 'Only matching application/harness commands checked; not proof of machine exclusivity.'})

    def record_scan(self, observation):
        self.scan_history.append(observation)
        if self.state is not None:
            self.state.setdefault('process_scan_observations', []).append(observation)
            self.persist()
        else:
            save(self.root / 'preflight.json', {'policy_sha256': self.fingerprint,
                 'process_scan_observations': self.scan_history,
                 'exclusive_campaign_lock_held': self.lock is not None})

    def verify_frozen_executables(self):
        # Small executable/harness reads only, outside application timing. Never
        # accept a replaced binary simply because its filename stayed the same.
        for binary in self.plan['binaries'].values():
            if sha256(self.resolve(binary['path'])) != binary['sha256']:
                raise RuntimeError('A frozen binary changed; explicitly update hashes/paths in a NEW campaign plan')
        if sha256(self.runner) != self.plan['runner']['sha256']:
            raise RuntimeError('The runner changed after planning; refresh the plan before starting a new campaign')

    def validate_case_inputs(self, case):
        # Deliberately deferred until this case starts: mixed100k may be generated
        # after this plan is written, before the campaign begins.
        for path in case['inputs'].values():
            if not self.resolve(path).is_file():
                raise RuntimeError('Missing case input: ' + str(self.resolve(path)))
        manifest = None
        if case.get('manifest'):
            path = self.resolve(case['manifest'])
            manifest = load(path)
            if manifest.get('coordinate_sorted') is not True:
                raise RuntimeError('Fixture is not declared coordinate sorted: ' + str(path))
            if case.get('manifest_sha256') and sha256(path) != case['manifest_sha256']:
                raise RuntimeError('Fixture manifest changed after planning: ' + str(path))
            for key in ('loci', 'read_records'):
                if case.get(key) is not None and manifest.get(key) != case[key]:
                    raise RuntimeError('Unexpected fixture {} for {}'.format(key, case['id']))
            for filename, expected in case.get('required_manifest_file_hashes', {}).items():
                if manifest.get('files', {}).get(filename, {}).get('sha256') != expected:
                    raise RuntimeError('Reused-cache reference/catalog mismatch in manifest: ' + filename)
        for cache in case.get('caches', {}).values():
            if not self.resolve(cache['path']).is_file():
                raise RuntimeError('Missing supplied cache: ' + str(self.resolve(cache['path'])))
        return manifest

    def validate_result(self, case, directory, manifest):
        data = load(directory / 'results.json')
        if data.get('status') != 'complete' or data.get('complete') is not True or data.get('errors'):
            raise RuntimeError('Runner did not finish cleanly: ' + str(directory))
        parameters = data.get('parameters', {})
        if sorted(parameters.get('threads', [])) != sorted(case['threads']) or parameters.get('pairs') != self.plan['policy']['pairs']:
            raise RuntimeError('Recorded thread/pair configuration differs from campaign plan')
        for key in ('decode', 'encode'):
            if str(parameters.get(key)) != str(case.get(key, self.plan['policy'][key])):
                raise RuntimeError('Recorded I/O configuration differs from campaign plan')
        for side in SIDES:
            if data.get('binaries', {}).get(side, {}).get('sha256') != self.plan['binaries'][side]['sha256']:
                raise RuntimeError('Measured binary differs from frozen campaign binary: ' + side)
            if case.get('caches'):
                actual = data.get('preparations', {}).get(side, {}).get('cache', {})
                if actual.get('sha256') != case['caches'][side]['sha256']:
                    raise RuntimeError('Supplied cache differs from matching preparation: ' + side)
        expected_files = manifest.get('files', {}) if manifest else {}
        inputs = data.get('inputs', [])
        if not isinstance(inputs, list):
            raise RuntimeError('Expected the new runner inputs list; legacy matrix metadata is unsupported here')
        observed_files = {str(Path(item['path']).resolve()): item for item in inputs}
        for path in case['inputs'].values():
            absolute = self.resolve(path)
            expected_hash = expected_files.get(absolute.name, {}).get('sha256')
            observed = observed_files.get(str(absolute), {})
            if not observed.get('sha256') or (expected_hash and observed['sha256'] != expected_hash):
                raise RuntimeError('Measured input missing or differs from its manifest: ' + str(absolute))
        pairs = data.get('pairs', [])
        if len(pairs) != len(case['threads']) * (self.plan['policy']['pairs'] + 1):
            raise RuntimeError('Missing or extra measured/warmup pairs')
        for workers in case['threads']:
            group = [pair for pair in pairs if pair.get('threads') == workers]
            warmups = [pair for pair in group if pair.get('warmup') is True]
            measured = [pair for pair in group if pair.get('warmup') is False]
            if len(warmups) != 1 or sorted(pair.get('round', -99) for pair in measured) != list(range(self.plan['policy']['pairs'])):
                raise RuntimeError('Warmup or measured round coverage is incomplete')
            for pair in group:
                if pair.get('status') != 'verified' or not exact(pair.get('comparison', {})):
                    raise RuntimeError('Complete output equality was not verified')
                if 'cross_configuration_comparison' in pair and not exact(pair['cross_configuration_comparison']):
                    raise RuntimeError('Cross-configuration output equality failed')
                for side in SIDES:
                    run = pair.get('runs', {}).get(side, {})
                    if run.get('status') != 'completed' or run.get('wait_status') != 0 or run.get('driver_returncode') != 0:
                        raise RuntimeError('Application run failed or has missing exit status')
                    for suffix in SUFFIXES:
                        output = run.get('outputs', {}).get(suffix, {})
                        if len(output.get('sha256', '')) != 64 or not isinstance(output.get('bytes'), int):
                            raise RuntimeError('Missing full-output digest/size evidence')
        # Performance flags are observations, deliberately not stop conditions.
        return {'results_json': str(directory / 'results.json'),
                'summary': data.get('summary', {}), 'inputs': data.get('inputs', []),
                'exact_outputs_verified': True}

    def execute_case(self, case, record):
        manifest = self.validate_case_inputs(case)
        self.require_no_other_run()
        attempts = record.setdefault('attempts', [])
        directory = self.root / case['id'] / 'attempt-{:03d}'.format(len(attempts) + 1)
        directory.parent.mkdir(parents=True, exist_ok=True)
        if directory.exists():
            raise RuntimeError('Fresh attempt directory already exists: ' + str(directory))
        command = self.command(case, directory)
        attempt = {'status': 'starting', 'started_utc': now(), 'output': str(directory), 'command': command}
        attempts.append(attempt)
        record['status'] = 'running'
        self.persist()
        log = directory.parent / (directory.name + '.launcher.log')
        print(json.dumps({'case': case['id'], 'status': 'starting', 'threads': case['threads'], 'output': str(directory)}), flush=True)
        with log.open('wb') as stream:
            process = subprocess.Popen(command, cwd=str(self.workspace), stdout=stream, stderr=subprocess.STDOUT,
                                       start_new_session=True, pass_fds=(self.lock.fileno(),))
            attempt.update(status='running', runner_pid=process.pid, launcher_log=str(log))
            self.persist()
            try:
                returncode = process.wait()
            except BaseException:
                # The runner handles SIGTERM by terminating its measurement
                # driver/application process group and persisting partial data.
                process.send_signal(signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    attempt['forced_runner_kill_requires_orphan_check'] = True
                attempt.update(status='interrupted', finished_utc=now())
                record['status'] = 'interrupted'
                self.persist()
                raise
        attempt.update(returncode=returncode, finished_utc=now())
        if returncode:
            attempt['status'] = 'failed'
            record['status'] = 'failed'
            self.persist()
            raise RuntimeError('Runner/correctness failure for {}; see {}'.format(case['id'], log))
        try:
            evidence = self.validate_result(case, directory, manifest)
        except BaseException as error:
            attempt.update(status='failed', validation_error=str(error))
            record['status'] = 'failed'
            self.persist()
            raise
        attempt.update(status='complete', evidence=evidence)
        record.update(status='complete', current_results_json=evidence['results_json'])
        self.persist()
        print(json.dumps({'case': case['id'], 'status': 'complete',
                          'performance_flags': {key: evidence['summary'].get(key) for key in
                              ('all_observed_median_wall_no_slower', 'all_observed_peak_rss_within_allowance')}}), flush=True)

    def execute(self):
        if self.plan.get('execution_enabled') is not True:
            raise RuntimeError('This is a provisional plan. Update the final binary/cache paths and expected hashes, '
                               'then explicitly set execution_enabled=true before execution')
        self.root.mkdir(parents=True, exist_ok=True)
        lock_path = self.resolve(self.plan.get('exclusive_lock', 'work/parallel-redesign/.stress-campaign.lock'))
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        self.lock = lock_path.open('a+')
        try:
            fcntl.flock(self.lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another campaign process still holds this output lock')
        self.require_no_other_run()
        for prerequisite in self.plan.get('prerequisites', []):
            data = load(self.resolve(prerequisite['results_json']))
            if data.get('complete') is not True or data.get('status') != 'complete':
                raise RuntimeError('Prerequisite benchmark is not complete: ' + prerequisite['results_json'])
        self.verify_frozen_executables()
        if self.state_path.exists():
            self.state = load(self.state_path)
            if self.state.get('policy_sha256') != self.fingerprint or self.state.get('host') != self.host():
                raise RuntimeError('Resume policy or host/allocation changed; use a separate campaign output')
            self.state.setdefault('process_scan_observations', []).extend(self.scan_history)
            self.persist()
        else:
            self.state = {'policy_sha256': self.fingerprint, 'plan': self.plan, 'plan_path': str(self.plan_path),
                          'started_utc': now(), 'host': self.host(), 'status': 'running',
                          'cases': {case['id']: {'status': 'pending', 'attempts': []} for case in self.plan['cases']},
                          'process_scan_observations': list(self.scan_history),
                          'measurement_limitation': 'Serial owned jobs and an exclusive campaign lock do not '
                                                    'establish exclusive machine use or stable OS memory pressure.'}
            self.persist()
        try:
            for case in self.plan['cases']:
                if self.options.only and case['id'] not in self.options.only:
                    continue
                self.verify_frozen_executables()
                record = self.state['cases'][case['id']]
                if record['status'] == 'complete':
                    self.validate_result(case, Path(record['current_results_json']).parent, self.validate_case_inputs(case))
                    print(json.dumps({'case': case['id'], 'status': 'already complete'}), flush=True)
                    continue
                if record['status'] != 'pending' and case['id'] not in self.options.retry_failed:
                    raise RuntimeError('Case {} needs review and explicit --retry-failed {}'.format(case['id'], case['id']))
                self.execute_case(case, record)
            all_complete = all(record['status'] == 'complete' for record in self.state['cases'].values())
            self.state.update(status='complete' if all_complete else 'partial', complete=all_complete)
            self.persist()
            save(self.root / 'report-inputs.json', {key: value['current_results_json'] for key, value in self.state['cases'].items()
                                                 if value['status'] == 'complete'})
            return 0
        except BaseException as error:
            self.state.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'stopped', complete=False)
            self.state.setdefault('errors', []).append({'time': now(), 'error': str(error), 'traceback': traceback.format_exc()})
            self.persist()
            raise
        finally:
            self.lock.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=Path(__file__).with_name('stress-campaign-v6-plan.json'))
    parser.add_argument('--out', help='Optional separate campaign output directory')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--execute', action='store_true', help='Actually run; default only prints commands')
    modes.add_argument('--dry-run', action='store_true', help='Print commands only; no process scan, executable or input hashing')
    parser.add_argument('--only', action='append', default=[], metavar='CASE')
    parser.add_argument('--retry-failed', action='append', default=[], metavar='CASE')
    parser.add_argument('--cpu-list', help='Optional one Linux physical-core CPU mask for every runner invocation')
    parser.add_argument('--allow-unavailable-process-scan', action='store_true',
                        help='Explicitly allow a blocked ps scan, record its reason and retain the campaign lock; '
                             'does not assert machine exclusivity')
    args = parser.parse_args()
    plan = load(args.plan)
    ids = {case['id'] for case in plan['cases']}
    if set(args.only + args.retry_failed) - ids:
        parser.error('Unknown --only/--retry-failed case ID')
    campaign = Campaign(plan, args.plan, args)
    if not args.execute:
        print(json.dumps({'mode': 'plan only; no application launched', 'counts': plan['counts'],
                          'publication_branch': plan['candidate_publication_branch'],
                          'commands': [{'case': case['id'], 'command': campaign.command(case, campaign.root / case['id'] / 'attempt-001')}
                                       for case in plan['cases'] if not args.only or case['id'] in args.only]}, indent=2))
        return 0
    def interrupt(signum, frame):
        raise KeyboardInterrupt('Received signal {}'.format(signum))
    signal.signal(signal.SIGTERM, interrupt)
    try:
        return campaign.execute()
    except BaseException as error:
        print('Campaign stopped: ' + str(error), file=sys.stderr)
        return 130 if isinstance(error, KeyboardInterrupt) else 1


if __name__ == '__main__':
    sys.exit(main())
