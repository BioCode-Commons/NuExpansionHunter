#!/usr/bin/env python3
"""Tiny fake-application checks; no genomic benchmark or performance assertions."""
import json
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


RUNNER = Path(__file__).with_name('run.py').resolve()
SPEC = importlib.util.spec_from_file_location('parallel_runner', RUNNER)
RUNNER_MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER_MODULE)
MOCK = r'''#!/usr/bin/env python3
import argparse, gzip, json, struct, sys
from pathlib import Path
p = argparse.ArgumentParser()
for name in ('reads', 'reference', 'variant-catalog', 'analysis-mode', 'threads',
             'decode-threads', 'encode-threads', 'read-decode-threads',
             'region-extension-length', 'log-level', 'output-prefix', 'catalog-cache'):
    p.add_argument('--' + name)
p.add_argument('--prepare-catalog-cache', action='store_true')
a = p.parse_args()
side = 'baseline' if 'baseline' in Path(__file__).name else 'candidate'
cache = Path(a.catalog_cache)
if a.prepare_catalog_cache:
    cache.write_text(side)
    sys.exit(0)
if cache.read_text() != side:
    sys.exit(6)
if 'failure' in Path(__file__).name:
    sys.exit(7)
prefix = a.output_prefix
Path(prefix + '.vcf').write_text('##fileformat=VCFv4.2\n')
Path(prefix + '.json').write_text('{"same": true}\n')
def record(name):
    core = struct.pack('<iiIIiiii', -1, -1, 2, 4 << 16, 1, -1, -1, 0)
    body = core + name.encode() + b'\0' + bytes([0x10, 40])
    return struct.pack('<i', len(body)) + body
x, y = record('x'), record('y')
records = [x, x, y] if side == 'baseline' else [y, x, x]
if 'unequal' in Path(__file__).name:
    records.pop()
with gzip.open(prefix + '_realigned.bam', 'wb') as out:
    out.write(b'BAM\1' + struct.pack('<ii', 0, 0) + b''.join(records))
print('Thread pools: compute={}, decode=2, encode=1'.format(a.threads))
print('Phase read_processing_seconds=0.01')
print('Phase output_seconds=0.001')
print('Streaming batches: pairs=3 chunks=1 alignment_worker_seconds=0.005 reduction_worker_seconds=0.001')
'''


class RunnerTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='parallel-runner-test-')
        self.root = Path(self.temp.name)
        for name in ('reads.bam', 'reference.fa', 'catalog.json'):
            (self.root / name).write_text('mock input\n')
        for name in ('baseline', 'candidate', 'candidate-unequal', 'candidate-failure'):
            path = self.root / name
            path.write_text(MOCK)
            path.chmod(0o755)

    def tearDown(self):
        self.temp.cleanup()

    def run_matrix(self, candidate='candidate', discard=False, baseline_label=None):
        command = [sys.executable, str(RUNNER), '--baseline', str(self.root / 'baseline'),
                   '--candidate', str(self.root / candidate), '--reads', str(self.root / 'reads.bam'),
                   '--reference', str(self.root / 'reference.fa'),
                   '--variant-catalog', str(self.root / 'catalog.json'),
                   '--threads', '1', '2', '3', '--pairs', '2', '--out', str(self.root / 'result')]
        if discard:
            command.append('--discard-equal-output')
        if baseline_label is not None:
            command += ['--baseline-label', baseline_label]
        process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=60)
        self.assertTrue((self.root / 'result/results.json').is_file(), process.stderr)
        return process, json.loads((self.root / 'result/results.json').read_text())

    def test_complete_counterbalanced_matrix_and_safe_cleanup(self):
        process, result = self.run_matrix(discard=True)
        self.assertEqual(0, process.returncode, process.stderr)
        self.assertTrue(result['complete'])
        self.assertEqual('RAMOPT b00650d', result['baseline_comparison_label'])
        self.assertEqual('b00650d', result['expected_baseline_source_commit'])
        self.assertTrue(result['summary']['all_requested_pairs_verified'])
        self.assertEqual([[1, 2, 3], [3, 2, 1]], result['round_configuration_orders'])
        self.assertEqual(9, len(result['pairs']))
        for pair in result['pairs']:
            self.assertEqual('verified', pair['status'])
            self.assertEqual([3, 3], pair['comparison']['bam']['record_counts'])
            self.assertTrue(pair['comparison']['bam']['record_multiset_identical'])
            self.assertFalse(pair['comparison']['bam']['compressed_bytes_identical'])
            if not pair['warmup']:
                self.assertEqual(['baseline', 'candidate'] if pair['round'] == 0 else ['candidate', 'baseline'], pair['order'])
            for record in pair['runs'].values():
                self.assertEqual(pair['threads'], record['resolved_thread_pools']['compute'])
                self.assertEqual(0.005, record['streaming_batches']['alignment_worker_seconds'])
                self.assertTrue(Path(record['stderr_path']).is_file())
                for output in record['outputs'].values():
                    self.assertTrue(output['removed'])
                    self.assertFalse(Path(output['path']).exists())
                    self.assertEqual(64, len(output['sha256']))
        for workers, summary in result['summary']['by_compute_workers'].items():
            self.assertEqual(2, summary['verified_pairs'])
            self.assertEqual(int(workers) * 20000000, summary['memory_allowance']['allowed_extra_bytes'])
        self.assertNotEqual(result['preparations']['baseline']['cache']['sha256'],
                            result['preparations']['candidate']['cache']['sha256'])

    def test_duplicate_mismatch_fails_and_keeps_outputs(self):
        process, result = self.run_matrix(candidate='candidate-unequal', discard=True)
        self.assertEqual(1, process.returncode)
        self.assertFalse(result['complete'])
        self.assertFalse(result['summary']['all_requested_pairs_verified'])
        pair = result['pairs'][0]
        self.assertEqual('mismatch', pair['status'])
        self.assertEqual([3, 2], pair['comparison']['bam']['record_counts'])
        for record in pair['runs'].values():
            for output in record['outputs'].values():
                self.assertTrue(Path(output['path']).is_file())
                self.assertNotIn('removed', output)

    def test_process_failure_preserves_measurement_and_partial_status(self):
        process, result = self.run_matrix(candidate='candidate-failure', baseline_label='prototype5 memory control')
        self.assertEqual(1, process.returncode)
        self.assertFalse(result['complete'])
        self.assertEqual('prototype5 memory control', result['baseline_comparison_label'])
        self.assertIsNone(result['expected_baseline_source_commit'])
        self.assertFalse(result['summary']['all_requested_pairs_verified'])
        self.assertEqual('failed', result['pairs'][0]['status'])
        failure = result['pairs'][0]['runs']['candidate']
        self.assertEqual('failed', failure['status'])
        self.assertNotEqual(0, failure['wait_status'])
        self.assertTrue(Path(failure['stderr_path']).is_file())
        self.assertTrue(result['errors'])


class AffinityTest(unittest.TestCase):
    def test_ranges_and_duplicate_rejection(self):
        self.assertEqual([0, 1, 4, 5], RUNNER_MODULE.parse_cpu_list('0-1,4-5'))
        for invalid in ('', '1-0', '0,0', '0-2,2', '-1'):
            with self.assertRaises(ValueError):
                RUNNER_MODULE.parse_cpu_list(invalid)

    def test_linux_rejects_smt_siblings_without_changing_affinity(self):
        topology = [{'package': '0', 'core': '0'}, {'package': '0', 'core': '0'}]
        with mock.patch.object(RUNNER_MODULE.sys, 'platform', 'linux'), \
             mock.patch.object(RUNNER_MODULE.os, 'sched_getaffinity', return_value={0, 1}, create=True), \
             mock.patch.object(RUNNER_MODULE.os, 'sched_setaffinity', create=True) as setter, \
             mock.patch.object(RUNNER_MODULE, 'cpu_topology', return_value=topology):
            with self.assertRaisesRegex(ValueError, 'SMT siblings'):
                RUNNER_MODULE.configure_affinity('0,1', 2)
            setter.assert_not_called()

    def test_linux_verifies_one_shared_mask(self):
        topology = [{'package': '0', 'core': str(core)} for core in range(4)]
        with mock.patch.object(RUNNER_MODULE.sys, 'platform', 'linux'), \
             mock.patch.object(RUNNER_MODULE.os, 'sched_getaffinity', side_effect=[set(range(8)), set(range(4))], create=True), \
             mock.patch.object(RUNNER_MODULE.os, 'sched_setaffinity', create=True) as setter, \
             mock.patch.object(RUNNER_MODULE, 'cpu_topology', return_value=topology):
            result = RUNNER_MODULE.configure_affinity('0-3', 4)
            setter.assert_called_once_with(0, [0, 1, 2, 3])
            self.assertTrue(result['same_mask_for_all_thread_counts_and_children'])
            self.assertEqual([0, 1, 2, 3], result['effective_cpu_ids'])


if __name__ == '__main__':
    unittest.main()
