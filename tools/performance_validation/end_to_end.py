#!/usr/bin/env python3
"""Full executable comparison; VCF, JSON and compressed BAM must match bytes."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from measure import measure, summarize, files_identical


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--candidate', required=True)
    parser.add_argument('--input', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--copies', type=int, default=32)
    parser.add_argument('--pairs', type=int, default=9)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--modes', nargs='+', default=['seeking', 'streaming'])
    parser.add_argument('--cached', action='store_true')
    args = parser.parse_args()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    source = Path(args.input).resolve()
    inputs = out / 'input'
    inputs.mkdir(exist_ok=True)
    for name in ['variants.bam', 'variants.bam.bai', 'reference.fa']:
        shutil.copyfile(source / name, inputs / name)
    catalog = json.loads((source / 'variants.json').read_text())
    records = []
    for i in range(args.copies):
        for original in catalog:
            item = dict(original)
            item['LocusId'] = f"{original['LocusId']}_{i:03d}"
            if 'VariantId' in item:
                value = item['VariantId']
                item['VariantId'] = [f'{v}_{i:03d}' for v in value] if isinstance(value,list) else f'{value}_{i:03d}'
            records.append(item)
    catalog_path = inputs / 'catalog.json'
    catalog_path.write_text(json.dumps(records, indent=2) + '\n')
    common = ['--reads', str(inputs / 'variants.bam'), '--reference', str(inputs / 'reference.fa'),
              '--variant-catalog', str(catalog_path), '--threads', str(args.threads),
              '--encode-threads', '0', '--log-level', 'error']
    binaries = {'baseline': str(Path(args.baseline).resolve()), 'candidate': str(Path(args.candidate).resolve())}
    cache_args = {name: [] for name in binaries}
    if args.cached:
        for name, binary in binaries.items():
            archive = out / f'{name}.cache'
            measure([binary] + common + ['--analysis-mode', 'streaming', '--decode-threads', '0',
                    '--output-prefix', str(out / f'prepare-{name}'),
                    '--catalog-cache', str(archive), '--prepare-catalog-cache'], out / f'prepare-{name}')
            cache_args[name] = ['--catalog-cache', str(archive)]
    output = {'scope': 'Repository synthetic example; repeated identical-coordinate loci, not WGS',
              'copies': args.copies, 'threads': args.threads, 'cached': args.cached,
              'evidence_role': 'acceptance' if args.pairs >= 9 else 'smoke_correctness_only',
              'identity_policy': 'VCF, JSON and compressed BAM byte-for-byte, including BAM order',
              'input_sha256': {p.name: sha(p) for p in inputs.iterdir() if p.is_file()},
              'binary_sha256': {name: sha(Path(binary)) for name,binary in binaries.items()}, 'modes': {}}
    for mode in args.modes:
        pairs = []
        for pair in range(-1, args.pairs):
            order = ['baseline', 'candidate'] if pair % 2 == 0 else ['candidate', 'baseline']
            record = {'pair': pair, 'warmup': pair < 0, 'order': order}
            paths_by_name = {}
            for name in order:
                prefix = out / mode / f'pair-{pair:02d}-{name}' / 'result'
                prefix.parent.mkdir(parents=True, exist_ok=True)
                # Seeking rejects an explicit --decode-threads setting.
                mode_args = ['--decode-threads', '0'] if mode == 'streaming' else []
                command = [binaries[name]] + common + ['--analysis-mode', mode, '--output-prefix', str(prefix)] + mode_args + cache_args[name]
                # Delete only the three files this run owns. Otherwise a successful
                # run that omits an output could inherit a stale file and pass.
                files = [Path(str(prefix) + suffix) for suffix in ('.vcf', '.json', '_realigned.bam')]
                for path in files:
                    if path.exists():
                        path.unlink()
                record[name] = measure(command, prefix.parent / 'measure')
                if not all(path.is_file() for path in files):
                    raise RuntimeError(f'Missing output files: {[str(p) for p in files if not p.is_file()]}')
                record[name]['outputs'] = {p.name: {'sha256': sha(p), 'bytes': p.stat().st_size} for p in files}
                paths_by_name[name] = {p.name: p for p in files}
            record['hashes_identical'] = record['baseline']['outputs'] == record['candidate']['outputs']
            record['outputs_identical'] = record['hashes_identical'] and all(
                files_identical(path, paths_by_name['candidate'][name])
                for name, path in paths_by_name['baseline'].items())
            if pair >= 0:
                pairs.append(record)
            output['modes'][mode] = {'pairs': pairs}
            (out / 'results.json').write_text(json.dumps(output,indent=2)+'\n')
            print(json.dumps({'mode': mode, 'pair': pair, 'identity': record['outputs_identical'],
                              'baseline': record['baseline']['wall_seconds'], 'candidate': record['candidate']['wall_seconds']}),flush=True)
            if not record['outputs_identical']:
                raise RuntimeError(f'Full output bytes differ: {mode} pair {pair}')
        output['modes'][mode]['summary'] = summarize(pairs)
        output['modes'][mode]['summary']['outputs_identical'] = all(p['outputs_identical'] for p in pairs)
        (out / 'results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({m:v['summary'] for m,v in output['modes'].items()}),flush=True)


if __name__ == '__main__':
    main()
