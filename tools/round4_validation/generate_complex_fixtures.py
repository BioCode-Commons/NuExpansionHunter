#!/usr/bin/env python3
"""150 bp, 10%-target synthetic regression diagnostics; no EH timing."""
from pathlib import Path
from collections import Counter
from types import SimpleNamespace
import argparse, hashlib, importlib.util, json, re, shutil, subprocess, sys

tools_directory = Path(__file__).resolve().parents[1]
generator_path = tools_directory / 'streaming_validation/generate_catalog.py'
spec = importlib.util.spec_from_file_location('base_generator', generator_path)
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
seed = 'nueh-frontend-20261112'
quality_table = bytes(53 + i % 21 for i in range(256))

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def spans(cigar):
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    assert ''.join(n + op for n, op in ops) == cigar
    return (sum(int(n) for n, op in ops if op in 'MIS=X'),
            sum(int(n) for n, op in ops if op in 'MDN=X'))

def apply_cigar(reference, cigar, key):
    position = 0
    query = []
    inserted = hashlib.shake_256(key.encode()).digest(32).translate(generator.TRANSLATE_DNA).decode()
    insertion = 0
    for n, op in re.findall(r'(\d+)([MID])', cigar):
        n = int(n)
        if op == 'M':
            query.append(reference[position:position+n]); position += n
        elif op == 'D':
            position += n
        else:
            query.append(inserted[insertion:insertion+n]); insertion += n
    assert position == 150
    return ''.join(query)

def mutate_bases(sequence, offsets):
    result = list(sequence)
    for offset in offsets:
        result[offset] = 'ACGT'[('ACGT'.index(result[offset]) + 1) % 4]
    return ''.join(result)

def build(kind):
    directory = out / kind
    directory.mkdir(exist_ok=False)
    catalog = json.loads((source / 'catalog.json').read_text())
    assert len(catalog) == 1000
    ref_lines = (source / 'reference.fa').read_text().splitlines()
    contig = ref_lines[0][1:]
    assert contig == 'chrS0001' and all(not x.startswith('>') for x in ref_lines[1:])
    reference = ''.join(ref_lines[1:])
    original_reference = reference
    complex_graph = kind == 'complex-indels'
    if complex_graph:
        modified = bytearray(reference.encode())
        for index, entry in enumerate(catalog):
            begin = index * 4096 + 2048
            modified[begin:begin+30] = b'CAG' * 4 + b'T' + b'CGG' * 3 + b'AA' + b'CT' * 3
            entry['LocusStructure'] = '(CAG)*T(CGG)*AA(CT)*'
            entry['ReferenceRegion'] = [f'{contig}:{begin}-{begin+12}',
                f'{contig}:{begin+13}-{begin+22}', f'{contig}:{begin+24}-{begin+30}']
            entry['VariantType'] = ['Repeat'] * 3
            entry['VariantId'] = [entry['LocusId'] + suffix for suffix in ('_CAG', '_CGG', '_CT')]
        reference = modified.decode()
        fasta = generator.FastaWriter(directory / 'reference.fa')
        fasta.begin(contig, len(reference)); fasta.append(reference.encode()); fasta.end(); fasta.close()
        (directory / 'catalog.json').write_text(json.dumps(catalog, indent=2) + '\n')
    else:
        for name in ('reference.fa', 'reference.fa.fai', 'catalog.json'):
            shutil.copy2(source / name, directory / name)
    counts = Counter()
    per_locus = []
    config = SimpleNamespace(stride=4096, motif='CAG', repeat_count=10, read_length=150)
    with (directory / 'reads.sam').open('w') as f:
        f.write(f'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{contig}\tLN:{len(reference)}\n@RG\tID:synthetic\tSM:synthetic\n')
        for index in range(1000):
            base = index * 4096
            # Confirm the reused reference was generated with the recorded helper/seed.
            assert original_reference[base:base+4096].encode() == generator.locus_block(seed, index, 4096, 'CAG', 10)
            heavy = kind == 'skewed-heavy-light' and index % 50 == 0
            light_after_heavy = kind == 'skewed-heavy-light' and index % 50 == 1
            targets = 500 if heavy else 1 if light_after_heavy else 10 if kind == 'skewed-heavy-light' else 20
            near, outside = targets, targets * 8
            per_locus.append({'index': index, 'target_pairs': targets, 'heavy': heavy,
                              'light_after_heavy': light_after_heavy})
            rows = []
            near_first, near_last = generator.background_window(config, base)
            for category, pairs in [('target', targets), ('near-background', near), ('outside-background', outside)]:
                for pair in range(pairs):
                    if category == 'target':
                        left, right = 1873 + pair % 31 - 15, 2023 + pair % 31 - 15
                        prefix = 'p'
                    elif category == 'near-background':
                        left = near_first + pair % (near_last - near_first + 1)
                        right, prefix = left + 170, 'b'
                    else:
                        left = 32 + pair % 256
                        right, prefix = left + 170, 'o'
                    name = f'{prefix}{index:07d}_{pair:07d}'
                    for mate, start, other, flag in [(0, left, right, 99), (1, right, left, 147)]:
                        key = f'round4-regression:{kind}:{name}:{mate}'
                        sequence = reference[base+start:base+start+150]
                        cigar = '150M'
                        if category == 'target':
                            quality = 'I' * 150
                            if kind == 'late-mismatches':
                                offsets = [12 + pair % 9, 136 + pair % 7]
                                if (pair + mate) % 2:
                                    offsets.append(146 + pair % 3)
                                sequence = mutate_bases(sequence, offsets)
                                counts['target_with_' + str(len(offsets)) + '_high_quality_mismatches'] += 1
                            elif complex_graph or heavy:
                                cigar = ('28M2I37M3D25M1I57M', '35M3I40M3D72M', '24M2D49M2I75M')[(pair + mate) % 3]
                                sequence = apply_cigar(sequence, cigar, key)
                                if (pair + mate) % 3 == 0:
                                    sequence = mutate_bases(sequence, [16, 138])
                                counts['target_with_indels'] += 1
                            else:
                                counts['target_exact_match'] += 1
                        else:
                            quality = '*' if (pair * 2 + mate + 1) % 17 == 0 else hashlib.shake_256(key.encode()).digest(150).translate(quality_table).decode()
                        assert len(sequence) == 150 and spans(cigar) == (150, 150)
                        length = right + 150 - left
                        tlen = length if mate == 0 else -length
                        line = f'{name}\t{flag}\t{contig}\t{base+start+1}\t60\t{cigar}\t=\t{base+other+1}\t{tlen}\t{sequence}\t{quality}\tRG:Z:synthetic\n'
                        rows.append((start, name, mate, line))
                        counts[category] += 1
                        counts['records'] += 1
            for _, _, _, line in sorted(rows):
                f.write(line)
    assert counts['target'] * 10 == counts['records']
    # Independent SAM pass checks every row, pair, coordinate, span and category.
    pending = {}
    checked = Counter()
    last = -1
    max_pending = 0
    with (directory / 'reads.sam').open() as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.rstrip('\n').split('\t')
            name, flag, pos, cigar, mpos, tlen = fields[0], int(fields[1]), int(fields[3])-1, fields[5], int(fields[7])-1, int(fields[8])
            assert fields[2] == contig and fields[6] == '=' and flag in (99, 147)
            assert last <= pos and 0 <= pos and pos + 150 <= len(reference)
            last = pos
            assert len(fields[9]) == 150 and set(fields[9]) <= set('ACGT') and spans(cigar) == (150, 150)
            assert fields[10] == '*' or len(fields[10]) == 150
            block = pos // 4096
            local = pos % 4096
            if name.startswith('p'):
                assert 1048 <= local and local + 150 <= 3078 and fields[10] == 'I'*150
                checked['target'] += 1
            elif name.startswith('b'):
                assert local > 3078 and local >> 10 == 3
                checked['near-background'] += 1
            else:
                assert name.startswith('o') and local + 150 < 1024
                checked['outside-background'] += 1
            assert block == mpos // 4096
            if name in pending:
                old_flag, old_pos, old_mpos, old_tlen = pending.pop(name)
                assert old_flag != flag and old_pos == mpos and old_mpos == pos and old_tlen == -tlen
                assert abs(tlen) == abs(pos-mpos) + 150
            else:
                pending[name] = (flag, pos, mpos, tlen)
                max_pending = max(max_pending, len(pending))
            checked['records'] += 1
    assert not pending and all(checked[k] == counts[k] for k in checked)
    command = [str(helper), str(directory / 'reads.sam'), str(directory / 'reads.bam')]
    result = subprocess.run(command, capture_output=True, text=True)
    (directory / 'bam-conversion.log').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stderr
    validation = {'valid_sam_field_and_pair_checks': True, 'all_query_lengths_150': True,
                  'all_reference_spans_150': True, 'coordinate_sorted': True,
                  'max_pending_pairs_in_independent_validator': max_pending,
                  'counts': dict(checked), 'bam_conversion_and_index_exit': result.returncode}
    (directory / 'independent-validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    manifest = {'scope': 'Synthetic regression diagnostics only; no guaranteed fast-path activation, fallback, helper stall or runtime regression.',
        'fixture_kind': kind, 'read_length': 150, 'read_records': counts['records'],
        'target_read_records': counts['target'], 'background_read_records': counts['records']-counts['target'],
        'target_read_fraction': 0.1, 'loci': 1000, 'contigs': 1, 'counts': dict(counts),
        'source_fixture': str(source), 'source_manifest_sha256': sha(source/'manifest.json'),
        'source_reference_sha256': sha(source/'reference.fa'), 'base_generator_sha256': sha(generator_path),
        'generator_sha256': sha(Path(__file__)), 'generation_command': [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        'conversion_command': command, 'required_region_extension_length': 1000,
        'target_qualities': 'All Phred40; background deterministic per-read Phred20..40 with periodic missing QUAL.',
        'target_mutations': {'late-mismatches': 'Every target read has substitutions at offsets12..20 and136..142; half also146..148.',
            'complex-indels': 'Three repeat nodes (CAG)*T(CGG)*AA(CT)*; all target reads have balanced multi-base insertion/deletion, one third additionally two substitutions.',
            'skewed-heavy-light': 'Every50th locus has500 indel-heavy pairs, following locus1 exact pair, remaining loci10 exact pairs; background counts scale9:1 locally.'}[kind],
        'per_locus_load': per_locus,
        'limitations': ['Read-level perturbations are not consistent diploid variant truth.',
            'No 250bp/long reads, remote mates, supplementary records or structural variant truth.',
            '10% target is coordinate routing fraction; acceptance/realignment rates remain algorithm-dependent.'],
        'files': {p.name: {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in sorted(directory.iterdir()) if p.is_file()}}
    (directory / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'fixture': kind, 'records': counts['records'], 'target': counts['target'], 'valid': True}), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    kinds = ['late-mismatches', 'complex-indels', 'skewed-heavy-light']
    parser.add_argument('--only', choices=kinds)
    parser.add_argument('--source-fixture', required=True, type=Path,
                        help='Existing1000-locus fixture with the original nueh-frontend-20261112 reference')
    parser.add_argument('--out', required=True, type=Path,
                        help='Output parent directory; each fixture is created as a new child directory')
    parser.add_argument('--sam-to-bam', required=True, type=Path,
                        help='SAM-to-BAM converter built from tools/streaming_validation/sam_to_bam.c')
    args = parser.parse_args()
    source = args.source_fixture.resolve()
    out = args.out.resolve()
    helper = args.sam_to_bam.resolve()
    for name in ('reference.fa', 'reference.fa.fai', 'catalog.json', 'manifest.json'):
        if not (source / name).is_file():
            parser.error('Missing source fixture file: ' + str(source / name))
    if not helper.is_file():
        parser.error('SAM-to-BAM converter does not exist: ' + str(helper))
    out.mkdir(parents=True, exist_ok=True)
    for kind in [args.only] if args.only else kinds:
        build(kind)
