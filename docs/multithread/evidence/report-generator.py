#!/usr/bin/env python3
"""Produce English Markdown/standalone HTML from parallel_validation JSON only.

No executables, BAMs, references, catalogs or logs are opened. Separate input
labels remain separate comparisons, including ablations with another baseline.
"""
import argparse
import collections
import datetime
import hashlib
import html
import json
import math
from pathlib import Path
import re
import statistics

MIB = 1024 * 1024
SIDES = ('baseline', 'candidate')


def number(value, positive=False):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and (value > 0 if positive else value >= 0))


def integer(value, positive=False):
    return isinstance(value, int) and not isinstance(value, bool) and (value > 0 if positive else value >= 0)


def json_snapshot(path):
    raw = Path(path).read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError('Expected a JSON object: ' + str(path))
    return value, raw, hashlib.sha256(raw).hexdigest()


def comparison_verified(comparison):
    if not isinstance(comparison, dict):
        return False
    bam = comparison.get('bam', {})
    return all(comparison.get(key) is True for key in ('vcf_identical', 'json_identical',
                                                     'outputs_identical_under_user_policy')) and all(
        bam.get(key) is True for key in ('header_identical', 'record_multiset_identical'))


def pair_verified(pair):
    return (pair.get('status') == 'verified' and comparison_verified(pair.get('comparison'))
            and ('cross_configuration_comparison' not in pair
                 or comparison_verified(pair['cross_configuration_comparison'])))


def pair_mismatch(pair):
    if pair.get('status') == 'mismatch':
        return True
    for key in ('comparison', 'cross_configuration_comparison'):
        comparison = pair.get(key, {})
        if any(comparison.get(name) is False for name in ('vcf_identical', 'json_identical',
                                                         'outputs_identical_under_user_policy')):
            return True
        if any(comparison.get('bam', {}).get(name) is False for name in ('header_identical', 'record_multiset_identical')):
            return True
    return False


def measured_pair(pair):
    for side in SIDES:
        record = pair.get('runs', {}).get(side, {})
        if record.get('status') != 'completed' or record.get('wait_status') != 0 or record.get('driver_returncode') != 0:
            return False
        if not number(record.get('wall_seconds'), positive=True):
            return False
        if not all(number(record.get(key)) for key in ('user_seconds', 'system_seconds', 'peak_rss_bytes')):
            return False
    return True


def metric(pairs, getter):
    values = []
    for pair in pairs:
        try:
            before, after = [getter(pair['runs'][side]) for side in SIDES]
        except (KeyError, TypeError):
            continue
        if number(before) and number(after):
            values.append((before, after))
    if not values:
        return None
    ratios = [after / before for before, after in values if before > 0]
    return {'pairs': len(values),
            'baseline_median': statistics.median(before for before, _ in values),
            'candidate_median': statistics.median(after for _, after in values),
            'baseline_max': max(before for before, _ in values), 'candidate_max': max(after for _, after in values),
            'paired_ratios': ratios, 'paired_median_ratio': statistics.median(ratios) if ratios else None,
            'paired_deltas': [after - before for before, after in values],
            'paired_median_delta': statistics.median(after - before for before, after in values),
            'largest_paired_delta': max(after - before for before, after in values)}


def saved_interval(data, workers, ratios):
    saved = data.get('summary', {}).get('by_compute_workers', {}).get(str(workers), {}).get('metrics', {}).get('wall_seconds', {})
    old = saved.get('paired_ratios', [])
    interval = saved.get('paired_ratio_bootstrap_95')
    if (ratios and len(old) == len(ratios) and all(number(a) and math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
                                                 for a, b in zip(old, ratios))
            and isinstance(interval, list) and len(interval) == 2 and all(number(x) for x in interval)):
        return interval
    return None


def summarize_case(data):
    parameters = data.get('parameters', {})
    requested = parameters.get('threads', [])
    requested = [value for value in requested if integer(value, positive=True)] if isinstance(requested, list) else []
    expected = parameters.get('pairs')
    expected = expected if integer(expected, positive=True) else None
    pairs = data.get('pairs', [])
    pairs = pairs if isinstance(pairs, list) else []
    pairs = [pair for pair in pairs if isinstance(pair, dict)]
    observed = [pair.get('threads') for pair in pairs if integer(pair.get('threads'), positive=True)]
    workers_values = sorted(set(requested + observed)) or [None]
    allowance_mb = parameters.get('extra_cache_mb_per_worker', 20)
    allowance_mb = allowance_mb if number(allowance_mb) else None
    declared_complete = data.get('complete') is True and data.get('status') == 'complete'
    global_failure = data.get('status') in ('failed', 'interrupted') or bool(data.get('errors'))
    rows = []
    for workers in workers_values:
        records = [pair for pair in pairs if pair.get('threads') == workers]
        warmups = [pair for pair in records if pair.get('warmup') is True]
        measurements = [pair for pair in records if pair.get('warmup') is False]
        grouped = collections.defaultdict(list)
        issues = []
        for pair in measurements:
            index = pair.get('round')
            if integer(index) and (expected is None or index < expected):
                grouped[index].append(pair)
            else:
                issues.append('Unexpected or missing measured round number')
        duplicates = sorted(index for index, values in grouped.items() if len(values) != 1)
        if duplicates:
            issues.append('Duplicate rounds: ' + ', '.join(map(str, duplicates)))
        unique = [values[0] for index, values in sorted(grouped.items()) if len(values) == 1]
        timed = [pair for pair in unique if measured_pair(pair)]
        exact = [pair for pair in unique if measured_pair(pair) and pair_verified(pair)]
        missing = sorted(set(range(expected)) - {pair['round'] for pair in timed}) if expected else None
        if missing:
            issues.append('Missing completed measured rounds: ' + ', '.join(map(str, missing)))
        if workers not in requested or len(requested) != len(set(requested)):
            issues.append('Thread configuration metadata is missing, duplicated or unexpected')
        if expected is None:
            issues.append('Requested pair count is not recorded')
        warmup_verified = len(warmups) == 1 and measured_pair(warmups[0]) and pair_verified(warmups[0])
        if not warmup_verified:
            issues.append('Exactly one fully verified warmup pair is required')
        if not declared_complete:
            issues.append('Input run is not marked complete')
        mismatch = any(pair_mismatch(pair) for pair in records)
        failed = global_failure or any(pair.get('status') in ('failed', 'interrupted') for pair in records)
        complete = (declared_complete and not failed and not mismatch and not issues
                    and expected is not None and len(exact) == expected)
        state = 'Output mismatch' if mismatch else 'Failed/interrupted' if failed else 'Complete' if complete else 'Incomplete'
        wall = metric(timed, lambda record: record['wall_seconds'])
        rss = metric(timed, lambda record: record['peak_rss_bytes'])
        allowance = int(allowance_mb * 1000000) * workers if allowance_mb is not None and workers is not None else None
        within = (all(delta <= allowance for delta in rss['paired_deltas'])
                  and rss['candidate_median'] <= rss['baseline_median'] + allowance
                  and rss['candidate_max'] <= rss['baseline_max'] + allowance) if rss and allowance is not None else None
        rows.append({'threads': workers, 'state': state, 'requested_pairs': expected, 'measured_pairs': len(timed),
                     'exact_verified_pairs': len(exact), 'warmup_pairs': len(warmups), 'warmup_verified': warmup_verified,
                     'missing_measured_rounds': missing, 'duplicate_rounds': duplicates, 'issues': issues,
                     'outputs_mismatch': mismatch, 'wall_seconds': wall,
                     'wall_ratio_bootstrap_95': saved_interval(data, workers, wall['paired_ratios']) if wall else None,
                     'peak_rss_bytes': rss, 'allowed_extra_bytes': allowance,
                     'extra_decimal_mb_per_worker': allowance_mb, 'observed_memory_within_allowance': within,
                     'cpu_seconds': metric(timed, lambda record: record['user_seconds'] + record['system_seconds']),
                     'read_processing_seconds': metric(timed, lambda record: record['phases_seconds']['read_processing_seconds'])})
    return rows


def manifest_counts(manifest):
    result = {}
    for name, keys in (('records', ('read_records', 'records', 'bam_records', 'total_records')),
                       ('loci', ('loci', 'num_loci', 'locus_count'))):
        for key in keys:
            if integer(manifest.get(key)):
                result[name] = manifest[key]
                break
    return result


def clean_binary(record):
    return {'filename': Path(record.get('path', 'not-recorded')).name,
            'sha256': record.get('sha256'), 'bytes': record.get('bytes')}


def clean_sources(data):
    result = {}
    for side in SIDES:
        source = data.get('sources', {}).get(side) or {}
        result[side] = {'head': source.get('head', {}).get('stdout'),
                        'dirty': bool(source.get('status', {}).get('stdout')) if source else None,
                        'tracked_diff_sha256': source.get('tracked_diff', {}).get('sha256')}
    return result


def mapping(values, name):
    result = {}
    for value in values:
        if '=' not in value:
            raise ValueError(name + ' requires LABEL=VALUE')
        key, item = value.split('=', 1)
        if not key or not item or key in result:
            raise ValueError(name + ' requires nonempty, unique labels and values')
        result[key] = item
    return result


def fmt(value, digits=2, signed=False):
    if value is None:
        return 'Not available'
    return ('{:+,.%df}' % digits if signed else '{:,.%df}' % digits).format(value)


def percent(metric_value, expose_slowdown=False):
    if not metric_value or metric_value['paired_median_ratio'] is None:
        return 'Not available'
    delta = (metric_value['paired_median_ratio'] - 1) * 100
    return fmt(delta, signed=True) + '%' + (' slower' if expose_slowdown and delta > 0 else '')


def median_pair(value, divisor=1):
    return (fmt(value['baseline_median'] / divisor) + ' / ' + fmt(value['candidate_median'] / divisor)) if value else 'Not available'


def row_cells(case, row):
    wall, rss = row['wall_seconds'], row['peak_rss_bytes']
    interval = row['wall_ratio_bootstrap_95']
    ci = ('[' + fmt((interval[0] - 1) * 100, signed=True) + '%, '
          + fmt((interval[1] - 1) * 100, signed=True) + '%]') if interval else 'Not available'
    rss_cells = [(fmt(rss[side + '_median'] / MIB) + ' / ' + fmt(rss[side + '_max'] / MIB))
                 if rss else 'Not available' for side in SIDES]
    delta = (fmt(rss['paired_median_delta'] / MIB, signed=True) + ' / '
             + fmt(rss['largest_paired_delta'] / MIB, signed=True)) if rss else 'Not available'
    allowance = row['allowed_extra_bytes']
    memory = ('Not available' if row['observed_memory_within_allowance'] is None else
              'Within observed allowance' if row['observed_memory_within_allowance'] else 'Allowance exceeded')
    if allowance is not None:
        memory += ' (+' + fmt(allowance / MIB) + ' MiB)'
    if row['state'] != 'Complete' and row['observed_memory_within_allowance'] is not None:
        memory += '; partial'
    expected = str(row['requested_pairs']) if row['requested_pairs'] is not None else '?'
    exact = ('Mismatch' if row['outputs_mismatch'] else
             '{}/{} verified; warmup {}'.format(row['exact_verified_pairs'], expected,
                                                'verified' if row['warmup_verified'] else 'pending/invalid'))
    return [case['display_label'], '{:,}'.format(case['fixture']['records']) if case['fixture'].get('records') is not None else 'Not recorded',
            '{:,}'.format(case['fixture']['loci']) if case['fixture'].get('loci') is not None else 'Not recorded',
            str(row['threads']) if row['threads'] is not None else '?',
            '{}/{}'.format(row['measured_pairs'], expected), row['state'], median_pair(wall),
            percent(wall, True), ci, rss_cells[0], rss_cells[1], delta, memory, exact]


HEADERS = ['Fixture', 'Records', 'Loci', 'T', 'Pairs measured/requested', 'State', 'Wall B / C (s)',
           'Paired wall change', '95% CI of change', 'RSS B median / max (MiB)', 'RSS C median / max (MiB)',
           'Paired RSS delta median / largest (MiB)', 'RAM allowance', 'Exact outputs']


def md_escape(value):
    return html.escape(str(value), quote=False).replace('|', '\\|').replace('\n', ' ').replace('[', '\\[').replace(']', '\\]')


def table_markdown(headers, rows):
    return '\n'.join(['| ' + ' | '.join(md_escape(cell) for cell in headers) + ' |',
                      '| ' + ' | '.join('---' for _ in headers) + ' |']
                     + ['| ' + ' | '.join(md_escape(cell) for cell in row) + ' |' for row in rows])


def table_html(headers, rows):
    head = ''.join('<th scope="col">' + html.escape(str(cell)) + '</th>' for cell in headers)
    body = []
    for row in rows:
        cells = []
        for cell in row:
            text = str(cell)
            flagged = any(marker in text for marker in ('slower', 'exceeded', 'mismatch', 'Mismatch', 'Failed', 'Incomplete'))
            cells.append('<td' + (' class="flag"' if flagged else '') + '>' + html.escape(text) + '</td>')
        body.append('<tr>' + ''.join(cells) + '</tr>')
    return '<div class="table-wrap"><table><thead><tr>' + head + '</tr></thead><tbody>' + ''.join(body) + '</tbody></table></div>'


def render(cases, details, generated):
    intro = ('This is an observational report, not an acceptance declaration. B = baseline; C = candidate; '
             'T = compute workers, not a process-wide CPU cap. Positive paired wall changes are actual observed '
             'median slowdowns. One warmup pair per configuration is excluded from timing/RSS statistics but '
             'must be verified. Failed application runs are excluded from timing statistics; unverified completed '
             'pairs can appear and are identified separately. Comparisons with different baselines remain separate.')
    method = ('Seconds are end-to-end process wall time. Changes are the median of paired C/B ratios, not the '
              'ratio of medians. RSS is process peak memory; MiB = 1,048,576 bytes. Each memory allowance is '
              'baseline peak RSS + T × the recorded decimal-MB allowance (default 20,000,000 bytes per worker); '
              'both paired and median/maximum differences are checked. A partial result cannot establish a gate. '
              'Exact outputs require every VCF/JSON byte and complete BAM header/full-record multiset, preserving '
              'duplicates; only BAM record order is ignored. These checks are read from the runner, not repeated here.')
    ci_note = ('Confidence intervals are copied only when the runner\'s stored paired ratios match the displayed '
               'pair set. Missing or stale intervals are shown as unavailable. Few-pair bootstrap intervals are '
               'descriptive and do not establish a hardware-general absence of regression.')
    md = ['# Streaming parallel comparison report', '', 'Generated: ' + generated, '', intro, '', method, '', ci_note, '']
    body = ['<h1>Streaming parallel comparison report</h1>', '<p class="muted">Generated: ' + html.escape(generated) + '</p>']
    body += ['<p>' + html.escape(text) + '</p>' for text in (intro, method, ci_note)]
    for case in cases:
        label = case['display_label']
        status = 'Input run status: {}; declared complete: {}.'.format(case['run_status'], str(case['declared_complete']).lower())
        md += ['## ' + md_escape(label), '', status, '', table_markdown(HEADERS, [row_cells(case, row) for row in case['rows']]), '']
        body += ['<section><h2>' + html.escape(label) + '</h2><p class="status">' + html.escape(status) + '</p>',
                 table_html(HEADERS, [row_cells(case, row) for row in case['rows']])]
        links = [('Compact case JSON', case['compact_link'])]
        if case.get('raw_link'):
            links.append(('Full raw metadata', case['raw_link']))
        if case.get('manifest_link'):
            links.append(('Fixture manifest', case['manifest_link']))
        md += ['; '.join('[' + title + '](' + link + ')' for title, link in links) + '.', '']
        body.append('<p>' + '; '.join('<a href="' + html.escape(link, quote=True) + '">' + title + '</a>' for title, link in links) + '.</p>')
        for row in case['rows']:
            if row['issues']:
                text = 'T={}: {}'.format(row['threads'], '; '.join(row['issues']))
                md += ['- ' + md_escape(text)]
                body.append('<p class="flag">' + html.escape(text) + '</p>')
        md.append('')
        for side in SIDES:
            binary = case['binaries'][side]
            text = side.title() + ' SHA-256: ' + (binary.get('sha256') or 'Not recorded')
            source = case['sources'][side]
            if source.get('head'):
                text += '; source HEAD: ' + source['head'] + '; dirty tree: ' + str(source['dirty']).lower()
            md += [text + '.', '']
            body.append('<p class="provenance">' + html.escape(text) + '.</p>')
        settings = case['settings']
        io_policy = 'decode={}; encode={}; BAM queue override={}.'.format(
            settings.get('decode', 'Not recorded'), settings.get('encode', 'Not recorded'),
            settings.get('bamlet_queue_mb') if settings.get('bamlet_queue_mb') is not None else 'version defaults')
        for name, text in (('Timing policy', case['timing_policy']), ('Requested I/O settings', io_policy),
                           ('Fixture counts', case['fixture']['count_provenance']),
                           ('Fixture scope', case['fixture'].get('scope', 'Not recorded'))):
            md += [name + ': ' + md_escape(text), '']
            body.append('<p><strong>' + name + ':</strong> ' + html.escape(text) + '</p>')
        if details:
            headers = ['T', 'CPU pairs', 'CPU B / C (s)', 'Paired CPU change', 'Read-phase pairs',
                       'Read phase B / C (s)', 'Paired read-phase change']
            rows = [[row['threads'], row['cpu_seconds']['pairs'] if row['cpu_seconds'] else 0,
                     median_pair(row['cpu_seconds']), percent(row['cpu_seconds']),
                     row['read_processing_seconds']['pairs'] if row['read_processing_seconds'] else 0,
                     median_pair(row['read_processing_seconds']), percent(row['read_processing_seconds'])] for row in case['rows']]
            note = ('CPU = process user + system time. Read phase = logged elapsed wall time. Summed worker '
                    'alignment/reduction elapsed counters are not CPU time and are not substituted here.')
            md += ['### CPU and read-phase details', '', note, '', table_markdown(headers, rows), '']
            body += ['<h3>CPU and read-phase details</h3><p>' + html.escape(note) + '</p>', table_html(headers, rows)]
        body.append('</section>')
    css = ('body{font:15px/1.5 system-ui,sans-serif;color:#172332;background:#fff;margin:32px auto;padding:0 24px;max-width:1680px}'
           'h1,h2,h3{line-height:1.2}h2{margin-top:36px}.muted{color:#526173}.table-wrap{overflow-x:auto;margin:18px 0}'
           'table{border-collapse:collapse;font-size:13px;width:100%}th,td{padding:9px;border:1px solid #d9e0e8;text-align:right;vertical-align:top}'
           'th{background:#edf2f7}th:first-child,td:first-child{text-align:left}.flag{color:#a12c22;font-weight:600}'
           '.status{font-weight:600}.provenance{font:12px/1.6 ui-monospace,monospace;overflow-wrap:anywhere}a{color:#125f9a}'
           '@media print{body{margin:0;padding:0;font-size:10pt}.table-wrap{overflow:visible}table{font-size:7pt}section{break-inside:avoid}}')
    document = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
    document += '<title>Streaming parallel comparison report</title><style>' + css + '</style></head><body>' + ''.join(body) + '</body></html>\n'
    return '\n'.join(md), document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', action='append', required=True, metavar='LABEL=RESULTS_JSON')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--label', action='append', default=[], metavar='LABEL=DISPLAY_NAME')
    parser.add_argument('--manifest', action='append', default=[], metavar='LABEL=MANIFEST_JSON')
    parser.add_argument('--records', action='append', default=[], metavar='LABEL=COUNT')
    parser.add_argument('--loci', action='append', default=[], metavar='LABEL=COUNT')
    parser.add_argument('--details', action='store_true', help='Include CPU/read-phase appendix tables')
    parser.add_argument('--include-raw-metadata', action='store_true',
                        help='Copy complete input JSON/manifests, which can include private absolute paths')
    args = parser.parse_args()
    try:
        sources = mapping(args.result, '--result')
        overrides = {name: mapping(getattr(args, name), '--' + name) for name in ('label', 'manifest', 'records', 'loci')}
        if any(set(values) - set(sources) for values in overrides.values()):
            raise ValueError('Every override label must match a --result label')
        for name in ('records', 'loci'):
            for key, value in overrides[name].items():
                if not re.fullmatch(r'\d+', value):
                    raise ValueError('--' + name + ' counts must be non-negative integers')
                overrides[name][key] = int(value)
        if args.out.exists():
            raise ValueError('--out must be a new directory to avoid stale report files')
        cases = []
        raw_copies = {}
        for index, (label, path) in enumerate(sources.items(), 1):
            data, raw, digest = json_snapshot(path)
            slug = '{:02d}-'.format(index) + (re.sub(r'[^A-Za-z0-9._-]+', '-', label).strip('.-') or 'case')
            manifest = data.get('manifest') if isinstance(data.get('manifest'), dict) else {}
            manifest_raw = None
            manifest_path = overrides['manifest'].get(label)
            if manifest_path is None and not manifest:
                reads = data.get('parameters', {}).get('reads')
                adjacent = Path(reads).parent / 'manifest.json' if isinstance(reads, str) else None
                if adjacent is not None and adjacent.is_file():
                    manifest_path = str(adjacent)
            if manifest_path:
                manifest, manifest_raw, _ = json_snapshot(manifest_path)
            counts = manifest_counts(manifest)
            if isinstance(manifest.get('scope'), str):
                counts['scope'] = manifest['scope']
            provenance = ['Fixture manifest' if manifest else 'Counts unavailable; BAM/catalog files were not scanned']
            for name in ('records', 'loci'):
                if label in overrides[name]:
                    counts[name] = overrides[name][label]
                    provenance.append(name + ' supplied with --' + name)
            counts['count_provenance'] = '; '.join(provenance)
            case = {'label': label, 'display_label': overrides['label'].get(label, label),
                    'source_result_filename': Path(path).name, 'result_snapshot_sha256': digest,
                    'run_status': data.get('status', 'Not recorded'), 'declared_complete': data.get('complete') is True,
                    'fixture': counts, 'rows': summarize_case(data), 'sources': clean_sources(data),
                    'settings': {key: data.get('parameters', {}).get(key) for key in
                                 ('threads', 'pairs', 'decode', 'encode', 'bamlet_queue_mb', 'extra_cache_mb_per_worker')},
                    'binaries': {side: clean_binary(data.get('binaries', {}).get(side, {})) for side in SIDES},
                    'timing_policy': data.get('timing_policy', 'Not recorded'),
                    'compact_link': 'data/' + slug + '.json'}
            if args.include_raw_metadata:
                case['raw_link'] = 'data/' + slug + '.raw.json'
                raw_copies[case['raw_link']] = raw
                if manifest_raw is not None:
                    case['manifest_link'] = 'data/' + slug + '.manifest.json'
                    raw_copies[case['manifest_link']] = manifest_raw
            cases.append(case)
        generated = datetime.datetime.now(datetime.timezone.utc).isoformat()
        markdown, document = render(cases, args.details, generated)
        (args.out / 'data').mkdir(parents=True)
        for case in cases:
            (args.out / case['compact_link']).write_text(json.dumps(case, indent=2, allow_nan=False) + '\n')
        for name, raw in raw_copies.items():
            (args.out / name).write_bytes(raw)
        (args.out / 'report.md').write_text(markdown)
        (args.out / 'index.html').write_text(document)
        print('Wrote English report: ' + str(args.out / 'index.html'))
        return 0
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    raise SystemExit(main())
