#!/usr/bin/env python3
"""Collect complete application measurements into JSON and a static HTML table."""
import argparse
import html
import json
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', type=Path, action='append', default=[])
    p.add_argument('--extremes', type=Path, action='append', default=[])
    p.add_argument('--input-inventory', type=Path, action='append', default=[],
                   help='Optional independently checked fixture metadata when the original run had no manifest')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    inventories = [json.loads(path.read_text()) for path in a.input_inventory]
    inventory = {item['fixture']: item for item in inventories}
    paths = list(a.matrix)
    failures = []
    for root in a.extremes:
        outer = json.loads((root/'results.json').read_text())
        if not outer.get('complete'):
            raise ValueError('Incomplete extreme matrix: ' + str(root))
        for run in outer['runs']:
            if run['status'] != 'passed':
                failures.append({k: run[k] for k in ('fixture', 'mode', 'status', 'log')})
            path = root/(run['fixture']+'-'+run['mode'])/'results.json'
            if path.exists():
                paths.append(path)
    records, preparations, raw = [], [], []
    for path in paths:
        result = json.loads(path.read_text())
        if not result.get('complete'):
            raise ValueError('Incomplete application matrix: ' + str(path))
        raw.append({'source': str(path.resolve()), 'result': result})
        mode = 'owned' if result['parameters'].get('uncached') else 'mapped'
        for fixture, sides in result['preparations'].items():
            preparations.append({'fixture': fixture, 'sides': sides})
        for config in result['configurations']:
            if 'summary' not in config:
                continue
            summary = config['summary']
            manifest = result['inputs'][config['fixture']].get('manifest') or inventory.get(config['fixture'], {})
            records.append({'fixture': config['fixture'], 'mode': mode,
                            'threads': config['compute_workers'], 'pairs': len(config['pairs']),
                            'decode': config['hts_decode_workers'], 'encode': config['hts_encode_workers'],
                            'output_equal': summary['all_outputs_identical'], 'summary': summary,
                            'coordinate_sorted': manifest.get('coordinate_sorted'),
                            'loci': manifest.get('loci'), 'read_records': manifest.get('read_records'),
                            'candidate_sha256': result['binaries']['candidate']['sha256'],
                            'resolved_thread_pools': config['pairs'][0]['candidate'].get('resolved_thread_pools')
                                                     if config['pairs'] else None,
                            'wall_range_seconds': {side: [min(p[side]['wall_seconds'] for p in config['pairs']),
                                                         max(p[side]['wall_seconds'] for p in config['pairs'])]
                                                   for side in ('baseline', 'candidate')},
                            'source': str(path.resolve())})
    records.sort(key=lambda r: (r['mode'] != 'mapped', r['coordinate_sorted'] is False,
                               r['fixture'], r['threads'], str(r['decode'])))
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    data = {'rows': records, 'preparations': preparations, 'failures': failures, 'raw_matrices': raw,
            'supplemental_input_inventories': inventories,
            'notes': ['All application runs use streaming mode and preserve BAM output.',
                      'Times are process start through exit. Cache preparation is separate.',
                      'RAM is the median of individual process peak RSS observations; maximum also retained.',
                      'Percent changes are medians of paired candidate/baseline ratios.',
                      'Primary acceptance scope: coordinate-sorted BAM, streaming with catalog cache.',
                      'Unsorted inputs retained as supplementary evidence; no further unsorted tuning requested.',
                      'Synthetic diagnostics on macOS ARM64, not a WGS/general hardware guarantee.']}
    (out/'benchmark-data.json').write_text(json.dumps(data, indent=2)+'\n')

    def esc(value):
        return html.escape(str(value))

    def change(metric):
        delta = -metric['paired_reduction_percent']
        lo, hi = metric['paired_ratio_bootstrap_95']
        color = 'bad' if lo > 1 else 'good' if hi < 1 else 'neutral'
        return f'<td class="{color}">{delta:+.1f}%</td>'

    rows = []
    for r in records:
        s = r['summary']; w = s['wall']; m = s['rss']
        pools = r.get('resolved_thread_pools') or {}
        decode = str(r['decode']) if r['decode'] != 'auto' else str(pools.get('decode', '?'))+'(auto)'
        encode = str(r['encode']) if r['encode'] != 'auto' else str(pools.get('encode', '?'))+'(auto)'
        rows.append('<tr><td>'+esc(r['fixture'])+'</td><td>'+('Cached' if r['mode']=='mapped' else 'Uncached')
                    +'</td><td>'+('Sorted' if r['coordinate_sorted'] is True else 'Unsorted' if r['coordinate_sorted'] is False else 'Not recorded')
                    +f'</td><td>{r["threads"]}/{decode}/{encode}</td><td>{r["pairs"]}</td>'
                    +f'<td>{w["baseline_median"]:.3f}</td><td>{w["candidate_median"]:.3f}</td>'+change(w)
                    +f'<td>{m["baseline_median"]/1048576:.1f}</td><td>{m["candidate_median"]/1048576:.1f}</td>'+change(m)
                    +'<td>'+('Identical' if r['output_equal'] else 'Different')+'</td></tr>')
    details = []
    for r in records:
        s = r['summary']; w = s['wall']; rp = s['read_processing']; m = s['rss']
        interval = w['paired_ratio_bootstrap_95']
        conclusion = 'Faster' if interval[1] < 1 else 'Slower' if interval[0] > 1 else 'Inconclusive'
        ranges = r['wall_range_seconds']
        details.append('<tr><td>'+esc(r['fixture'])+'</td><td>'+esc(r['mode'])+f'</td><td>{r["threads"]}</td>'
                       +f'<td>{rp["baseline_median"]:.3f} → {rp["candidate_median"]:.3f}</td>'
                       +f'<td>{interval[0]:.3f}–{interval[1]:.3f}</td><td>{conclusion}</td>'
                       +f'<td>{ranges["baseline"][0]:.2f}–{ranges["baseline"][1]:.2f} / '
                        f'{ranges["candidate"][0]:.2f}–{ranges["candidate"][1]:.2f}</td>'
                       +f'<td>{m["baseline_max"]/1048576:.1f} → {m["candidate_max"]/1048576:.1f}</td></tr>')
    prep = []
    for r in preparations:
        for side, v in r['sides'].items():
            if 'wall_seconds' not in v:
                continue
            prep.append('<tr><td>'+esc(r['fixture'])+'</td><td>'+esc(side)
                        +f'</td><td>{v["wall_seconds"]:.3f}</td><td>{v["peak_rss_bytes"]/1048576:.1f}</td>'
                        +f'<td>{v["cache"]["bytes"]/1048576:.1f}</td></tr>')
    fixtures = {}
    for r in records:
        fixtures.setdefault(r['fixture'], r)
    inputs = []
    for name, r in sorted(fixtures.items()):
        inputs.append('<tr><td>'+esc(name)+'</td><td>'+esc(r['loci'] if r['loci'] is not None else 'Not recorded')
                      +'</td><td>'+esc(r['read_records'] if r['read_records'] is not None else 'Not recorded')
                      +'</td><td>'+('Sorted' if r['coordinate_sorted'] is True else 'Unsorted' if r['coordinate_sorted'] is False else 'Not recorded')+'</td></tr>')
    page = '''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NuExpansionHunter streaming memory comparison</title>
<style>body{font:15px/1.6 system-ui,sans-serif;margin:32px;color:#17232c;background:#fafbfd}
h1{font-size:26px}h2{font-size:20px;margin-top:30px}table{border-collapse:collapse;background:white;width:100%;font-variant-numeric:tabular-nums}
th,td{padding:8px 10px;border-bottom:1px solid #dce2e8;text-align:right;white-space:nowrap}th{background:#edf2f6;position:sticky;top:0}
td:first-child,th:first-child{text-align:left}.wrap{overflow:auto}.good{color:#116342}.bad{color:#9c3229}p{max-width:1000px}.note{background:#fff3d8;padding:12px;border-radius:6px}a{color:#126397}</style>
<h1>NuExpansionHunter: streaming memory comparison</h1>
<p>Baseline: main at f3b1256. Candidate: exact index compression, on-demand queues and earlier streaming cleanup.
The primary scope is coordinate-sorted BAM, streaming, and a prepared catalog cache.
All tests retain VCF, JSON and BAM output. Only BAM record order may differ; complete fields and duplicate multiplicity are compared.
Additional sequence-preprocessing workers are zero in all configurations. Counts below are compute / HTS decode / BAM encode workers, not CPU-core limits.</p>
<p>Total time includes process startup, index loading or construction, read processing, genotyping, output, and exit. Cache preparation is separate.
Time and peak RAM are medians across measured runs. RAM is in MiB; maximum observed process peak RSS is retained below.
Changes use the median of paired ratios; negative means lower and can differ from the ratio of the two marginal medians.</p>
<p>Changes are green/red only when the 95% bootstrap interval is entirely below/above parity. Intervals crossing parity use plain text and are inconclusive.</p>
<p class="note">These are local macOS ARM64 synthetic stress tests. Sub-second runs and intervals based on only three pairs require caution.
An interval crossing parity cannot establish speedup or regression and does not guarantee no regression. covered10k is approximately 8x LowDepth;
genotyped10k additionally checks actual Genotype fields. Existing unsorted results are supplementary and outside the user's primary acceptance scope.
Small catalogs, long reads, distant mates and other cases do not necessarily achieve 50% lower memory.</p>
<p><a href="benchmark-data.json">Complete measurement and comparison records (JSON)</a></p>
<h2>End-to-end time and peak RAM</h2><div class="wrap"><table><thead><tr>
<th>Dataset</th><th>Mode</th><th>BAM order</th><th>Compute/decode/encode</th><th>Pairs</th><th>Baseline s</th><th>Candidate s</th><th>Time change</th>
<th>Baseline peak MiB</th><th>Candidate peak MiB</th><th>RAM change</th><th>Outputs</th></tr></thead><tbody>'''+''.join(rows)+'''</tbody></table></div>
<h2>Phases, timing uncertainty and maximum peaks</h2><div class="wrap"><table><thead><tr><th>Dataset</th><th>Mode</th><th>Compute workers</th>
<th>Read-processing s: baseline → candidate</th><th>Wall-time paired ratio: 95% bootstrap interval</th><th>Direction</th>
<th>Wall-time range s: baseline / candidate</th><th>Maximum peak MiB: baseline → candidate</th>
</tr></thead><tbody>'''+''.join(details)+'''</tbody></table></div>
<h2>One-time cache preparation</h2><div class="wrap"><table><thead><tr><th>Dataset</th><th>Version</th><th>Preparation s</th><th>Peak MiB</th><th>Cache file MiB</th>
</tr></thead><tbody>'''+''.join(prep)+'''</tbody></table></div>'''
    page += '''<h2>Input sizes</h2><div class="wrap"><table><thead><tr><th>Dataset</th><th>Catalog loci</th>
<th>Input reads</th><th>Coordinate order</th></tr></thead><tbody>'''+''.join(inputs)+'''</tbody></table></div>'''
    if failures:
        page += '<h2>Incomplete or failed cases</h2><pre>'+esc(json.dumps(failures, indent=2, ensure_ascii=False))+'</pre>'
    page += '</html>\n'
    (out/'benchmark-comparison.html').write_text(page)
    print(json.dumps({'rows': len(records), 'failures': len(failures), 'out': str(out)}))


if __name__ == '__main__':
    main()
