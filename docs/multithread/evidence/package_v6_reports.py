#!/usr/bin/env python3
"""Package completed v6 metadata only; never run analyses or read genomic inputs.

Run only on the root agent's explicit instruction after the complete campaign
and T8 recheck. The destination must not exist. All source evidence is retained.
"""
import argparse
import copy
import datetime
import hashlib
import html
import importlib.util
import json
from pathlib import Path
import re
import statistics
import subprocess
import sys
import tempfile

MIB = 1048576
MAX_SMALL_FILE = 16 * 1024 * 1024
READ_NOTE = (
    'Read-phase boundaries differ: b00650d logs the phase before unmatched-read '
    'and input-stream destruction; v6 logs it after the input scope closes. '
    'These diagnostic elapsed times are not strictly equivalent intervals. '
    'Pending-mate cleanup may contribute to an observed difference; no entire '
    'delta is attributed to cleanup. End-to-end process wall time is the primary comparison.')
SELECTION_NOTE = (
    'The original genotyped10k-fixed T8 configuration has an objective timing '
    'inconsistency in baseline formal round 2: process wall and application/UTC '
    'elapsed times disagree. Its entire original result file and all runs remain '
    'in evidence. Main tables use the original T1/T4 configurations and the '
    'complete independent T8 recheck (three formal pairs plus one warmup), '
    'rather than selecting individual favorable repetitions. No measurements '
    'are changed. Original T8 is excluded only from selected main statistics.')


def load(path):
    path = Path(path)
    if path.stat().st_size > MAX_SMALL_FILE:
        raise ValueError('Refusing unexpectedly large metadata file: ' + str(path))
    return json.loads(path.read_text())


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def import_report(path):
    spec = importlib.util.spec_from_file_location('v6_publication_report', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_result(data, report, label):
    if data.get('complete') is not True or data.get('status') != 'complete' or data.get('errors'):
        raise ValueError(label + ': incomplete or failed result')
    rows = report.summarize_case(data)
    if not rows or any(row['state'] != 'Complete' for row in rows):
        raise ValueError(label + ': missing, duplicate, failed, or unequal pairs')
    if data['parameters']['pairs'] != 3:
        raise ValueError(label + ': expected three formal pairs')
    return rows


class Package:
    def __init__(self, workspace, out):
        self.workspace, self.out = workspace, out
        self.artifacts = []
        self.hostnames = set()

    def clean_text(self, text):
        text = text.replace(str(self.workspace), '<WORKSPACE>')
        text = re.sub(r'<USER_HOME>/\s"\x27<>]+', '<USER_HOME>', text)
        text = re.sub(r'<USER_HOME>/\s"\x27<>]+', '<USER_HOME>', text)
        text = text.replace('<APPLE_COMMAND_LINE_TOOLS>', '<APPLE_COMMAND_LINE_TOOLS>')
        for name in sorted(self.hostnames, key=len, reverse=True):
            if name:
                text = text.replace(name, '<LOCAL_HOST>')
        return text

    def clean(self, value):
        if isinstance(value, str):
            return self.clean_text(value)
        if isinstance(value, list):
            return [self.clean(v) for v in value]
        if isinstance(value, dict):
            return {self.clean_text(k): ('<LOCAL_HOST>' if k in ('hostname', 'node')
                    and isinstance(v, str) else self.clean(v)) for k, v in value.items()}
        return value

    def write_json(self, relative, value):
        path = self.out / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.clean(value), indent=2, allow_nan=False) + '\n')
        return path

    def copy_small(self, source, relative, replacements=()):
        source = Path(source)
        if source.stat().st_size > MAX_SMALL_FILE:
            raise ValueError('Refusing large source artifact: ' + str(source))
        raw = source.read_bytes()
        text = self.clean_text(raw.decode())
        for before, after in replacements:
            text = text.replace(before, after)
        dest = self.out / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text)
        self.artifacts.append(dict(source=self.clean_text(str(source)), published=relative,
            original_sha256=sha(raw), published_sha256=sha(text.encode()),
            original_bytes=len(raw), published_bytes=len(text.encode()),
            redaction='Workspace/user/host identifiers only; numeric measurements preserved.'))
        return dest


def md_table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(str(x) for x in row) + ' |' for row in rows])


def html_table(headers, rows):
    return '<div class="scroll"><table><thead><tr>' + ''.join(
        '<th>' + html.escape(str(x)) + '</th>' for x in headers) + '</tr></thead><tbody>' + ''.join(
        '<tr>' + ''.join('<td>' + html.escape(str(x)) + '</td>' for x in row) + '</tr>'
        for row in rows) + '</tbody></table></div>'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    base = workspace / 'work/parallel-redesign'
    out = args.out.resolve() if args.out else base / 'publication-v6'
    if out.exists():
        raise ValueError('Destination must not exist: ' + str(out))
    report_path = workspace / 'work/NuExpansionHunter/tools/parallel_validation/report.py'
    report = import_report(report_path)
    candidate_manifest = load(base / 'prototype6-manifest.json')
    baseline_manifest = load(workspace / 'work/frozen-v1/manifest.json')
    campaign = load(base / 'final-campaign-v6/campaign.json')
    plan = load(base / 'stress-campaign-v6-plan.json')
    if campaign.get('status') != 'complete' or campaign.get('complete') is not True:
        raise ValueError('Final stress campaign is not complete; no package created')
    if len(plan['cases']) != 13 or len(campaign['cases']) != 13:
        raise ValueError('Unexpected stress-case count')
    primary_spec = dict(id='primary100k', fixture_id='genotyped100k', threads=[8],
        manifest='work/fixtures/genotyped100k/manifest.json', loci=100000, read_records=20000000)
    originals = [(primary_spec, base / 'primary100k-v6/results.json')]
    for case in plan['cases']:
        record = campaign['cases'][case['id']]
        if record.get('status') != 'complete' or len(record.get('attempts', [])) != 1:
            raise ValueError('Unexpected incomplete/retried case; review all attempts: ' + case['id'])
        path = Path(record['current_results_json'])
        originals.append((case, path if path.is_absolute() else workspace / path))
    recheck_spec = dict(id='genotyped10k-t8-recheck', fixture_id='genotyped10k', threads=[8],
        manifest='work/fixtures/genotyped10k/manifest.json', loci=10000, read_records=2000000)
    originals.append((recheck_spec, base / 'genotyped10k-t8-recheck-v6/results.json'))
    loaded = []
    for spec, path in originals:
        data = load(path)
        validate_result(data, report, spec['id'])
        for role, expected in (('candidate', candidate_manifest['binary_sha256']),
                               ('baseline', baseline_manifest['binary_sha256'])):
            if data.get('binaries', {}).get(role, {}).get('sha256') != expected:
                raise ValueError(spec['id'] + ': unexpected ' + role + ' executable identity')
        if data['parameters']['threads'] != spec['threads']:
            raise ValueError(spec['id'] + ': unexpected thread configurations')
        loaded.append((spec, path, data))
    formal = sum(not p['warmup'] for _, _, d in loaded for p in d['pairs'])
    warmups = sum(p['warmup'] for _, _, d in loaded for p in d['pairs'])
    if (formal, warmups) != (84, 28):
        raise ValueError('Unexpected executed pair counts: ' + repr((formal, warmups)))

    # All completion/equality prerequisites above precede destination creation.
    out.mkdir(parents=True)
    package = Package(workspace, out)
    for _, _, data in loaded:
        environment = data.get('environment', {})
        for key in ('hostname', 'node'):
            if isinstance(environment.get(key), str):
                package.hostnames.add(environment[key])
    staged = base / 'report-stage/evidence'
    for path in sorted(staged.rglob('*')):
        if path.is_file():
            relative = path.relative_to(staged).as_posix()
            if relative == 'staging-manifest.json':
                relative = 'initial-staging-manifest.json'
            replacements = [('"artifact_manifest": "staging-manifest.json"',
                             '"artifact_manifest": "initial-staging-manifest.json"')]
            package.copy_small(path, 'evidence/' + relative, replacements)
    copies = [
        ('stress-campaign-v6-plan.json', 'evidence/final-v6-stress-plan.json'),
        ('final-campaign-v6/campaign.json', 'evidence/final-v6-stress-results.json'),
        ('prototype6-manifest.json', 'evidence/frozen-v6-source-manifest.json'),
        ('owner-ablation/results/results.json', 'evidence/owner-ablation-results.json'),
        ('owner-ablation/manifest.json', 'evidence/owner-ablation-manifest.json'),
        ('genotyped10k-t8-recheck-command.json', 'evidence/genotyped10k-t8-recheck-command.json'),
        ('audit_timing_integrity.py', 'evidence/audit_timing_integrity.py'),
        ('run_stress_campaign.py', 'evidence/run_stress_campaign.py'),
        ('diagnosis-and-design.md', 'diagnosis-and-design.md'),
        ('final-report-draft.md', 'report-draft.md'),
        ('release-review-v6.md', 'evidence/release-review-v6.md'),
        ('publication-audit.md', 'evidence/publication-audit.md'),
    ]
    for source, target in copies:
        package.copy_small(base / source, target)
    for path in sorted((workspace / 'work/frozen-v1').glob('*manifest*.json')):
        package.copy_small(path, 'evidence/baseline-' + path.name)
    package.copy_small(Path(__file__), 'evidence/package_v6_reports.py')
    package.copy_small(report_path, 'evidence/report-generator.py')

    selected = []
    selection = []
    for spec, path, data in loaded:
        label = spec['id']
        raw_relative = 'evidence/results/' + label + '.json'
        package.copy_small(path, raw_relative)
        selected_data = copy.deepcopy(data)
        if label == 'genotyped10k-fixed':
            selected_data['parameters']['threads'] = [1, 4]
            selected_data['pairs'] = [p for p in selected_data['pairs'] if p['threads'] in (1, 4)]
            # Keep matched confidence intervals; do not carry aggregate T8 flags.
            by_t = selected_data.get('summary', {}).get('by_compute_workers', {})
            selected_data['summary'] = {'by_compute_workers': {k: v for k, v in by_t.items() if k in ('1', '4')}}
            selected_data['publication_selection'] = dict(reason=SELECTION_NOTE,
                original_result=raw_relative, original_file_sha256=sha(path.read_bytes()),
                retained_threads=[1, 4], excluded_threads=[8], measurement_values_modified=False)
        validate_result(selected_data, report, label + ' selected')
        input_path = package.write_json('evidence/report-inputs/' + label + '.json', selected_data)
        manifest_path = None
        if spec.get('manifest'):
            manifest_path = package.copy_small(workspace / spec['manifest'],
                'evidence/fixtures/' + label + '.manifest.json')
        selection.append(dict(label=label, logical_case='genotyped10k-fixed' if label == recheck_spec['id'] else label,
            original_result=raw_relative, selected_input=str(input_path.relative_to(out)),
            selected_threads=selected_data['parameters']['threads']))
        selected.append((spec, selected_data, input_path, manifest_path))
    configs = sum(len(d['parameters']['threads']) for _, d, _, _ in selected)
    selected_formal = sum(not p['warmup'] for _, d, _, _ in selected for p in d['pairs'])
    selected_warm = sum(p['warmup'] for _, d, _, _ in selected for p in d['pairs'])
    if (configs, selected_formal, selected_warm) != (27, 81, 27):
        raise ValueError('Unexpected selected counts')

    # Refresh the timing audit from small metadata, preserving original anomaly.
    with tempfile.TemporaryDirectory(prefix='v6-publication-audit-', dir=str(base)) as temporary:
        audit_path = Path(temporary) / 'audit.json'
        subprocess.run([sys.executable, str(base / 'audit_timing_integrity.py'),
            '--campaign', str(base / 'final-campaign-v6'), '--result',
            'primary100k=' + str(base / 'primary100k-v6/results.json'), '--result',
            'genotyped10k-t8-recheck=' + str(base / 'genotyped10k-t8-recheck-v6/results.json'),
            '--out', str(audit_path)], check=True)
        audit = load(audit_path)
        if audit['sample_analysis_runs_checked'] != 224:
            raise ValueError('Final timing audit does not cover all 224 analyses')
        package.copy_small(audit_path, 'evidence/timing-integrity-audit.json')
        unexpected_flags = [r for r in audit['flagged_runs'] if not (
            r.get('case') == 'genotyped10k-fixed/attempt-001'
            and r.get('kind') == 'sample_analysis' and r.get('role') == 'baseline'
            and r.get('threads') == 8 and r.get('round') == 2 and r.get('warmup') is False)]
        if unexpected_flags:
            raise ValueError('New timing-integrity flags require explicit review; partial package retained')
        for index, flagged in enumerate(audit['flagged_runs'], 1):
            for field in ('stdout_path', 'stderr_path'):
                if flagged.get(field):
                    package.copy_small(Path(flagged[field]),
                        'evidence/timing-anomalies/{:02d}-{}'.format(index, Path(flagged[field]).name))

    command = [sys.executable, str(report_path), '--details', '--include-raw-metadata',
               '--out', str(out / 'measurements')]
    for spec, _, source, manifest in selected:
        label = spec['id']
        command += ['--result', label + '=' + str(source)]
        if manifest:
            command += ['--manifest', label + '=' + str(manifest)]
        for key, flag in (('loci', '--loci'), ('read_records', '--records')):
            if spec.get(key) is not None:
                command += [flag, label + '=' + str(spec[key])]
    subprocess.run(command, check=True)
    package.write_json('evidence/report-command.json', dict(command=command,
        note='Only sanitized publication inputs are passed to report.py; no analyses are run.'))
    # Place the boundary disclosure beside every generated read-phase table.
    p = out / 'measurements/report.md'
    p.write_text(p.read_text().replace('### CPU and read-phase details\n',
        '### CPU and read-phase details\n\n' + READ_NOTE + '\n'))
    p = out / 'measurements/index.html'
    p.write_text(p.read_text().replace('<h3>CPU and read-phase details</h3>',
        '<h3>CPU and read-phase details</h3><p>' + html.escape(READ_NOTE) + '</p>'))

    headers = ['Case/settings', 'T', 'E2E median B/C (s)', 'Paired E2E change',
               'Read median B/C (s)', 'Paired read change', 'Peak RSS max B/C (MiB)',
               'RSS maxima change', 'Largest paired RSS delta', 'Memory reference']
    rows, machine_rows = [], []
    for spec, data, _, _ in selected:
        for row in report.summarize_case(data):
            wall, read, rss = [row[k] for k in ('wall_seconds', 'read_processing_seconds', 'peak_rss_bytes')]
            metric_pairs = [p for p in data['pairs'] if not p['warmup'] and p['threads'] == row['threads']]
            rss_delta_percent = max((p['runs']['candidate']['peak_rss_bytes'] /
                p['runs']['baseline']['peak_rss_bytes'] - 1) * 100 for p in metric_pairs)
            seconds = lambda value: ('{:.6f}' if abs(value) < 0.01 else '{:.3f}').format(value)
            fmt = lambda m: seconds(m['baseline_median']) + ' / ' + seconds(m['candidate_median']) if m else 'Unavailable'
            pct = lambda m: '{:+.2f}%'.format((m['paired_median_ratio'] - 1) * 100) if m else 'Unavailable'
            cells = [spec['id'], row['threads'], fmt(wall), pct(wall), fmt(read), pct(read),
                '{:.3f} / {:.3f}'.format(rss['baseline_max'] / MIB, rss['candidate_max'] / MIB),
                '{:+.2f}%'.format((rss['candidate_max'] / rss['baseline_max'] - 1) * 100),
                '{:+.3f} MiB / {:+.2f}%'.format(rss['largest_paired_delta'] / MIB, rss_delta_percent),
                'Within' if row['observed_memory_within_allowance'] else 'Exceeded']
            rows.append(cells)
            machine_rows.append(dict(case=spec['id'], **row,
                largest_paired_rss_increase_percent=rss_delta_percent))
    count_note = ('14 fixture/settings cases (fixed and automatic I/O share the genotyped10k fixture); '
        '27 main configurations; 81 selected formal pairs plus 27 warmup pairs. '
        'Including the preserved original T8 and its supplement: 84 formal pairs plus 28 warmup pairs, '
        '224 analysis processes. Cache preparation processes are excluded from these counts.')
    explanation = ('B = frozen RAMOPT baseline; C = frozen v6. Medians/maxima use only formal runs. '
        'Paired time change is the median C/B ratio minus one. RSS maxima change compares the two '
        'formal maxima. Largest paired RSS delta reports separate maxima of absolute and percentage '
        'differences over matched pairs; these two maxima can come from different pairs. Negative '
        'values are retained. Memory reference means the original baseline peak RSS plus '
        '20,000,000 bytes per compute worker. The user later allowed a modest disclosed excess; '
        'this reference is not an acceptance gate. '
        'Every observed slowdown and memory-reference exceedance remains visible. Complete outputs '
        'do not establish performance acceptance. Three pairs on one ARM host are not a server-wide claim.')
    overview = '# Frozen v6: complete selected measurement matrix\n\n' + count_note + '\n\n' + SELECTION_NOTE + '\n\n' + READ_NOTE + '\n\n' + explanation + '\n\n' + md_table(headers, rows) + '\n'
    (out / 'measurement-overview.md').write_text(overview)

    run_headers = ['Case/settings', 'T', 'Round (zero-based)', 'Execution order', 'Role',
                   'E2E (s)', 'Read phase (s)', 'Peak RSS (MiB)', 'Main statistics']
    run_rows = []
    for spec, _, data in loaded:
        for pair in data['pairs']:
            if pair['warmup']:
                continue
            chosen = not (spec['id'] == 'genotyped10k-fixed' and pair['threads'] == 8)
            for role in pair['order']:
                run = pair['runs'][role]
                read = run.get('phases_seconds', {}).get('read_processing_seconds')
                run_rows.append([spec['id'], pair['threads'], pair['round'], '/'.join(pair['order']),
                    role, '{:.9f}'.format(run['wall_seconds']), 'Unavailable' if read is None else '{:.9f}'.format(read),
                    '{:.6f}'.format(run['peak_rss_bytes'] / MIB), 'Selected' if chosen else 'Preserved original T8; excluded'])
    (out / 'individual-formal-runs.md').write_text('# Every executed formal analysis\n\n' +
        '168 individual analyses from all 84 formal pairs, including original anomalous T8. '
        'All 56 warmup analyses remain in raw result JSON.\n\n' + READ_NOTE + '\n\n' + md_table(run_headers, run_rows) + '\n')
    html_body = '<h1>Frozen v6 measurements</h1>' + ''.join('<p>' + html.escape(t) + '</p>'
        for t in (count_note, SELECTION_NOTE, READ_NOTE, explanation)) + html_table(headers, rows)
    html_body += '<h2>Every executed formal analysis</h2><p>' + html.escape(READ_NOTE) + '</p>' + html_table(run_headers, run_rows)
    (out / 'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Frozen v6 measurements</title>'
        '<style>body{font:14px system-ui;margin:2em;line-height:1.5}.scroll{overflow:auto}table{border-collapse:collapse;white-space:nowrap}'
        'th,td{border:1px solid #ccd;padding:6px;text-align:right}th:first-child,td:first-child{text-align:left}th{background:#eef}</style>'
        '<body>' + html_body + '</body></html>\n')
    package.write_json('measurement-matrix.json', dict(counts=dict(cases=14, configurations=27,
        selected_formal_pairs=81, selected_warmup_pairs=27, executed_formal_pairs=84,
        executed_warmup_pairs=28, analysis_processes=224), read_boundary_note=READ_NOTE,
        selection_note=SELECTION_NOTE, rows=machine_rows))
    package.write_json('evidence/selection-provenance.json', dict(note=SELECTION_NOTE, inputs=selection))
    (out / 'methodology-and-selection.md').write_text('# Methodology and selection\n\n' +
        count_note + '\n\n' + SELECTION_NOTE + '\n\n' + READ_NOTE + '\n\n' + explanation + '\n\n'
        'Original complete results, source and fixture manifests, comparison evidence, and timing-audit '
        'flags are preserved under evidence/. Workspace and user/host identifiers are redacted; original '
        'and publication-file hashes are recorded in manifest.json. Historical v4 and owner-refcount '
        'ablation results are separate evidence, never pooled into the 27-configuration v6 matrix. '
        'The original baseline freeze manifest predates b00650d and retains its historical base field; '
        'it is not rewritten as contemporaneous commit provenance.\n')
    final_files = {}
    for path in sorted(out.rglob('*')):
        if path.is_file():
            raw = path.read_bytes()
            final_files[path.relative_to(out).as_posix()] = dict(sha256=sha(raw), bytes=len(raw))
    package.write_json('manifest.json', dict(schema_version=1, created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status='Completed measurement evidence; performance acceptance is not asserted; prose draft requires final review.',
        branch='Jiapeng-multithread', baseline_commit='b00650d297754f68d78778a441e2fcbf5096250b',
        frozen_v6_source_manifest='evidence/frozen-v6-source-manifest.json',
        source_manifest_status_note='Original experimental status is retained as historical provenance.',
        historical_staging_manifest='evidence/initial-staging-manifest.json',
        historical_staging_note='Describes earlier partial staging and its hashes only. Canonical final publication-file hashes are in this manifest files mapping.',
        counts=dict(cases=14, configurations=27, selected_formal_pairs=81, selected_warmup_pairs=27,
            executed_formal_pairs=84, executed_warmup_pairs=28, analysis_processes=224),
        selection='evidence/selection-provenance.json', artifacts=package.artifacts, files=final_files))
    print('Prepared publication metadata at ' + str(out))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
