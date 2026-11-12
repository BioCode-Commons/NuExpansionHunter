#!/usr/bin/env python3
"""Serial alternating-process streaming E2E measurements and exact output checks.

Do not run this while building, generating fixtures, hashing large files, or
running another benchmark. A 32-thread run is not a claim about 32 physical cores.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path
import sys

from compare_outputs import compare_prefixes
from measure import measure, summarize


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline', required=True, type=Path)
    p.add_argument('--candidate', required=True, type=Path)
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--pairs', type=int, default=9)
    p.add_argument('--compute', type=int, default=8)
    p.add_argument('--candidate-compute', type=int)
    p.add_argument('--hts-decode', type=int, default=0)
    p.add_argument('--encode', type=int, default=0)
    p.add_argument('--baseline-extra-json', default='[]')
    p.add_argument('--candidate-extra-json', default='[]')
    p.add_argument('--candidate-preprocess-workers', type=int, default=0)
    p.add_argument('--cache-bytes-per-preprocess-worker', type=int, default=20000000)
    p.add_argument('--cached', action='store_true')
    p.add_argument('--baseline-cache',type=Path)
    p.add_argument('--candidate-cache',type=Path)
    args = p.parse_args()
    if args.candidate_compute is None:
        args.candidate_compute = args.compute
    if min(args.pairs,args.compute,args.candidate_compute) < 1 or min(args.hts_decode,args.encode,args.candidate_preprocess_workers,args.cache_bytes_per_preprocess_worker) < 0:
        p.error('Invalid counts')
    extra = {'baseline': json.loads(args.baseline_extra_json), 'candidate': json.loads(args.candidate_extra_json)}
    if not all(isinstance(a,list) and all(isinstance(v,str) for v in a) for a in extra.values()):
        p.error('Extra arguments must be JSON string arrays')
    root = args.out.resolve()
    root.mkdir(parents=True,exist_ok=False)
    inputs = args.input.resolve()
    for name in ('reference.fa','reference.fa.fai','catalog.json','reads.bam'):
        if not (inputs/name).is_file():
            p.error('Missing input: '+str(inputs/name))
    binaries = {'baseline': str(args.baseline.resolve()),'candidate': str(args.candidate.resolve())}
    common = ['--reads',str(inputs/'reads.bam'),'--reference',str(inputs/'reference.fa'),
              '--variant-catalog',str(inputs/'catalog.json'),'--analysis-mode','streaming',
              '--decode-threads',str(args.hts_decode),'--encode-threads',str(args.encode),'--log-level','info']
    compute = {'baseline':args.compute,'candidate':args.candidate_compute}
    cache_args = {name: [] for name in binaries}
    supplied_cache = {'baseline': args.baseline_cache, 'candidate': args.candidate_cache}
    if any(supplied_cache.values()) and not all(supplied_cache.values()):
        p.error('Provide both cache paths')
    if any(supplied_cache.values()) and not args.cached:
        p.error('Reusable cache paths require --cached')
    result = {'scope':'Synthetic independent-coordinate loci; performance applies only to this workload and host',
              'binary_sha256':{name:file_sha(path) for name,path in binaries.items()},
              'input_manifest':json.loads((inputs/'manifest.json').read_text()) if (inputs/'manifest.json').is_file() else None,
              'parameters':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
              'preparations':{},'pairs':[],'warmup':None,
              'output_policy':'Exact VCF/JSON and BAM header/full-record-byte multiset including duplicate multiplicity; BAM record order may differ',
              'performance_policy':'9 AB/BA pairs plus one discarded warmup; paired ratio bootstrap. No automatic retention decision.',
              'ram_policy':'Report original no-growth gate and explicit preprocessing-cache envelope separately. Process RSS is not proof of per-worker cache accounting.'}
    def save():
        (root/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    if args.cached:
        for name,binary in binaries.items():
            if supplied_cache[name]:
                cache = supplied_cache[name].resolve()
                if not cache.is_file():
                    p.error('Missing reusable cache: '+str(cache))
                result['preparations'][name] = {'reused_cache':str(cache),'bytes':cache.stat().st_size,'sha256':file_sha(cache)}
                cache_args[name] = ['--catalog-cache',str(cache)]
                save()
                continue
            cache = root/(name+'.cache')
            cmd = [binary]+common+['--threads',str(compute[name])]+extra[name]+[
                '--output-prefix',str(root/('prepare-'+name)), '--catalog-cache',str(cache),'--prepare-catalog-cache']
            result['preparations'][name] = measure(cmd,root/('prepare-'+name+'-measure'))
            cache_args[name] = ['--catalog-cache',str(cache)]
            save()
    for pair in range(-1,args.pairs):
        order = ['baseline','candidate'] if pair % 2 == 0 else ['candidate','baseline']
        record = {'pair':pair,'warmup':pair<0,'order':order}
        prefixes = {}
        for name in order:
            directory = root/f'pair-{pair:02d}-{name}'
            directory.mkdir()
            prefix = directory/'result'
            cmd = [binaries[name]]+common+['--threads',str(compute[name])]+extra[name]+cache_args[name]+['--output-prefix',str(prefix)]
            record[name] = measure(cmd,directory/'measure')
            if name == 'candidate':
                log = Path(record[name]['stdout_path']).read_text() + Path(record[name]['stderr_path']).read_text()
                cache_matches = re.findall(r'Read preprocessing workers=(\d+) peak_worker_cache_bytes=(\d+)', log)
                if not cache_matches:
                    raise RuntimeError('Missing candidate preprocessing cache telemetry')
                workers, cache_peak = map(int, cache_matches[-1])
                record[name]['preprocessing_cache'] = {'workers':workers,'peak_worker_cache_bytes':cache_peak}
                if workers != args.candidate_preprocess_workers or cache_peak > args.cache_bytes_per_preprocess_worker:
                    raise RuntimeError('Candidate preprocessing worker count or cache budget mismatch')
            prefixes[name] = prefix
            record[name]['output_files'] = {suffix: {'bytes':Path(str(prefix)+suffix).stat().st_size,
                'sha256':file_sha(str(prefix)+suffix)} for suffix in ('.vcf','.json','_realigned.bam')}
        # All parsing/hashing/sorting is outside process timings and before next run.
        record['comparison'] = compare_prefixes(prefixes['baseline'],prefixes['candidate'],root/'comparison-scratch')
        record['artifacts_identical'] = record['comparison']['outputs_identical_under_user_policy']
        if pair < 0:
            result['warmup'] = record
        else:
            result['pairs'].append(record)
        save()
        print(json.dumps({'pair':pair,'outputs_identical':record['artifacts_identical'],
                          'baseline_wall':record['baseline']['wall_seconds'],'candidate_wall':record['candidate']['wall_seconds']}),flush=True)
        if not record['artifacts_identical']:
            raise RuntimeError(f'Output mismatch at pair {pair}; see results.json')
    summary = summarize(result['pairs'])
    extra_budget = args.candidate_preprocess_workers * args.cache_bytes_per_preprocess_worker
    summary['outputs_identical_under_user_policy'] = all(r['artifacts_identical'] for r in result['pairs'])
    summary['preprocessing_cache_envelope_bytes'] = extra_budget
    summary['rss_median_within_cache_envelope'] = summary['candidate_median_peak_rss_bytes'] <= summary['baseline_median_peak_rss_bytes'] + extra_budget + summary['rss_noise_allowance_bytes']
    summary['rss_max_within_cache_envelope'] = summary['candidate_max_peak_rss_bytes'] <= summary['baseline_max_peak_rss_bytes'] + extra_budget + summary['rss_max_noise_allowance_bytes']
    summary['automatic_retention_decision'] = None
    result['summary'] = summary
    save()
    print(json.dumps(summary),flush=True)


if __name__ == '__main__':
    main()
