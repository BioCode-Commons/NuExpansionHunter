#!/usr/bin/env python3
"""Read small result JSON files and flag timing inconsistencies without mutation.

No application, comparison, build, data hash, or large log read is performed.
UTC wrapper timestamps include driver/metadata overhead; their tolerance is
deliberately wider than the tolerance for nested application phase durations.
"""
import argparse
import datetime
import json
from pathlib import Path


def inspect_run(label, source, context, run, args):
    wall = run.get('wall_seconds')
    if run.get('status') != 'completed' or wall is None:
        return None
    phases = run.get('phases_seconds') or {}
    tolerance = max(args.phase_absolute_tolerance, abs(wall) * args.phase_relative_tolerance)
    phase_sum = sum(phases.values())
    flags = []
    if wall <= 0:
        flags.append('nonpositive_process_wall')
    for name, value in phases.items():
        if value < 0:
            flags.append('negative_phase:' + name)
        if value > wall + tolerance:
            flags.append('phase_exceeds_process_wall:' + name)
    if phase_sum > wall + tolerance:
        flags.append('phase_sum_exceeds_process_wall')
    span = None
    if run.get('started_utc') and run.get('finished_utc'):
        span = (datetime.datetime.fromisoformat(run['finished_utc'])
                - datetime.datetime.fromisoformat(run['started_utc'])).total_seconds()
        utc_tolerance = max(args.utc_absolute_tolerance, abs(wall) * args.utc_relative_tolerance)
        if abs(span - wall) > utc_tolerance:
            flags.append('utc_wrapper_span_differs_from_process_wall')
    else:
        utc_tolerance = None
    return dict(case=label, source_results=str(source), **context,
                wall_seconds=wall, phases_seconds=phases, phase_sum_seconds=phase_sum,
                phase_tolerance_seconds=tolerance, started_utc=run.get('started_utc'),
                finished_utc=run.get('finished_utc'), utc_wrapper_span_seconds=span,
                utc_minus_process_wall_seconds=None if span is None else span - wall,
                utc_tolerance_seconds=utc_tolerance, flags=flags,
                stdout_path=run.get('stdout_path'), stderr_path=run.get('stderr_path'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path)
    parser.add_argument('--result', action='append', default=[], help='LABEL=path/to/results.json')
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--phase-absolute-tolerance', type=float, default=0.05)
    parser.add_argument('--phase-relative-tolerance', type=float, default=0.02)
    parser.add_argument('--utc-absolute-tolerance', type=float, default=2.0)
    parser.add_argument('--utc-relative-tolerance', type=float, default=0.05)
    args = parser.parse_args()
    inputs = []
    if args.campaign:
        inputs.extend((p.parent.parent.name + '/' + p.parent.name, p)
                      for p in sorted(args.campaign.glob('*/attempt-*/results.json')))
    for specification in args.result:
        label, separator, path = specification.partition('=')
        if not separator:
            parser.error('--result requires LABEL=path')
        inputs.append((label, Path(path)))
    if not inputs:
        parser.error('Provide --campaign or --result')
    files, runs = [], []
    for label, path in inputs:
        # Read each changing result once. A result is an observation at this
        # instant; incomplete cases are never promoted to completed evidence.
        data = json.loads(path.read_text())
        files.append(dict(label=label, path=str(path), status=data.get('status'),
                          complete=data.get('complete'), bytes=path.stat().st_size))
        for role, run in data.get('preparations', {}).items():
            context = dict(kind='preparation', role=role, threads=None, round=None, warmup=None)
            record = inspect_run(label, path, context, run, args)
            if record is not None:
                runs.append(record)
        for pair in data.get('pairs', []):
            for role, run in pair.get('runs', {}).items():
                context = dict(kind='sample_analysis', role=role, threads=pair.get('threads'),
                               round=pair.get('round'), warmup=pair.get('warmup'),
                               pair_status=pair.get('status'))
                record = inspect_run(label, path, context, run, args)
                if record is not None:
                    runs.append(record)
    flagged = [run for run in runs if run['flags']]
    result = dict(
        captured_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        purpose='Timing integrity audit only; no results are changed and no performance acceptance is inferred.',
        policy={key: getattr(args, key) for key in (
            'phase_absolute_tolerance', 'phase_relative_tolerance',
            'utc_absolute_tolerance', 'utc_relative_tolerance')},
        definitions={
            'process_wall': 'Original runner time.perf_counter duration around posix_spawn plus wait4.',
            'phase': 'Application logged elapsed seconds, normally a subset of process runtime.',
            'utc_wrapper_span': 'Runner UTC timestamps surrounding the separate measurement driver; includes additional overhead.',
            'tolerance': 'max(absolute_seconds, relative_fraction * process_wall_seconds).',
            'limitations': 'Flags identify inconsistent clocks/boundaries or suspicious metadata; absence of flags does not establish comparable host conditions. Completed runs in incomplete cases are included only for integrity checks.'},
        files=files, completed_runs_checked=len(runs),
        sample_analysis_runs_checked=sum(r['kind'] == 'sample_analysis' for r in runs),
        preparation_runs_checked=sum(r['kind'] == 'preparation' for r in runs),
        flagged_run_count=len(flagged), flagged_runs=flagged, all_checked_runs=runs)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('completed_runs_checked',
          'sample_analysis_runs_checked', 'preparation_runs_checked', 'flagged_run_count')}, indent=2))


if __name__ == '__main__':
    main()
