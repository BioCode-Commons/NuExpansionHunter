#!/usr/bin/env python3
"""Prepare a 10%-target diagnostic with deterministic harder target read sequences."""
from pathlib import Path
from collections import Counter
from itertools import zip_longest
import hashlib
import json
import argparse
import shutil
import sys
import re
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-fixture', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
parser.add_argument('--sam-to-bam', type=Path, required=True)
args = parser.parse_args()
source = args.source_fixture.resolve()
target = args.out.resolve()
target.mkdir(parents=True, exist_ok=False)
for name in ('reference.fa', 'reference.fa.fai', 'catalog.json'):
    shutil.copy2(source / name, target / name)

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def spans(cigar):
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    assert ''.join(n + op for n, op in ops) == cigar
    return (sum(int(n) for n, op in ops if op in 'MIS=X'),
            sum(int(n) for n, op in ops if op in 'MDN=X'))

counts = Counter()
with (source / 'reads.sam').open() as original, (target / 'reads.sam').open('w') as changed:
    for line in original:
        if line.startswith('@'):
            changed.write(line)
            continue
        fields = line.rstrip('\n').split('\t')
        counts['records'] += 1
        if fields[0].startswith('p'):
            group = int(fields[0].rsplit('_', 1)[1]) % 3
            assert fields[5] == '150M' and len(fields[9]) == 150
            counts['target'] += 1
            counts['target_group_' + str(group)] += 1
            if group == 1:
                sequence = list(fields[9])
                for offset in (20, 120):
                    sequence[offset] = 'ACGT'[('ACGT'.index(sequence[offset]) + 1) % 4]
                fields[9] = ''.join(sequence)
            elif group == 2:
                old = fields[9]
                fields[9] = old[:50] + 'G' + old[50:100] + old[101:]
                if fields[10] != '*':
                    old_qualities = fields[10]
                    fields[10] = old_qualities[:50] + 'I' + old_qualities[50:100] + old_qualities[101:]
                fields[5] = '50M1I50M1D49M'
            line = '\t'.join(fields) + '\n'
        else:
            counts['background'] += 1
        changed.write(line)

# Independent byte/field pass over both ordered SAM streams. Verify every
# background record byte-for-byte and all non-sequence target metadata unchanged.
background_original = hashlib.sha256()
background_changed = hashlib.sha256()
verified = Counter()
last_coordinate = (-1, -1)
contigs = {}
with (source / 'reads.sam').open() as original, (target / 'reads.sam').open() as changed:
    for left, right in zip_longest(original, changed):
        assert left is not None and right is not None
        if left.startswith('@'):
            assert right == left
            if left.startswith('@SQ'):
                name = next(v[3:] for v in left.rstrip().split('\t') if v.startswith('SN:'))
                contigs[name] = len(contigs)
            continue
        old = left.rstrip('\n').split('\t')
        new = right.rstrip('\n').split('\t')
        assert len(old) == len(new)
        coordinate = (contigs[new[2]], int(new[3]))
        assert coordinate >= last_coordinate
        last_coordinate = coordinate
        assert int(new[1]) in (99, 147)
        assert spans(new[5]) == (len(new[9]), 150)
        assert new[10] == '*' or len(new[10]) == len(new[9])
        assert len(new[9]) == 150
        for i in range(len(old)):
            if i not in (5, 9, 10):
                assert old[i] == new[i]
        if not old[0].startswith('p'):
            assert right == left
            background_original.update(left.encode())
            background_changed.update(right.encode())
            verified['background'] += 1
            continue
        verified['target'] += 1
        group = int(old[0].rsplit('_', 1)[1]) % 3
        if group == 0:
            assert right == left
        elif group == 1:
            assert [i for i, (a,b) in enumerate(zip(old[9],new[9])) if a != b] == [20,120]
            assert old[5] == new[5] and old[10] == new[10]
        else:
            assert new[5] == '50M1I50M1D49M'
            assert new[9][:50] == old[9][:50]
            assert new[9][50] == 'G'
            assert new[9][51:101] == old[9][50:100]
            assert new[9][101:] == old[9][101:]
            if old[10] == '*':
                assert new[10] == '*'
            else:
                assert new[10][:50] == old[10][:50]
                assert new[10][50] == 'I'
                assert new[10][51:101] == old[10][50:100]
                assert new[10][101:] == old[10][101:]
assert verified['target'] == counts['target']
assert verified['background'] == counts['background']
assert counts['target'] * 10 == counts['records']
assert background_original.hexdigest() == background_changed.hexdigest()

command = [str(args.sam_to_bam.resolve()), str(target / 'reads.sam'), str(target / 'reads.bam')]
subprocess.run(command, check=True)
manifest = json.loads((source / 'manifest.json').read_text())
manifest.update(fixture_kind='mixed-indels', source_fixture=str(source),
                source_manifest_sha256=sha(source/'manifest.json'),
                conversion_command=command, hard_fixture_generator_sha256=sha(Path(__file__)),
                target_read_fraction=counts['target']/counts['records'])
manifest['parameters']['out'] = str(target)
manifest['generation_command'] = [sys.executable] + sys.argv
manifest['target_mutations'] = {
    'grouping': 'Target pair index modulo 3; same group for both mates',
    'counts': dict(counts),
    'group_0': 'Original target sequence with its original one substitution',
    'group_1': 'Original target sequence plus substitutions at query offsets 20 and 120',
    'group_2': "Original sequence[:50] + G + sequence[50:100] + sequence[101:]; CIGAR 50M1I50M1D49M",
    'indel_qualities': 'Same transformation; inserted base Phred 40; missing qualities remain *',
    'scope': 'Synthetic diagnostic exercising harder read alignment; no guarantee a particular DP budget or fallback is selected.'
}
manifest['limitations'] += [
    'The target variants are read-level simulations; they are not a consistent diploid genotype truth set.',
    'Higher DP budgets/fallback are diagnostic motivations, not measured or guaranteed properties of this fixture.'
]
verification = {
    'valid': True, 'read_records': counts['records'], 'target_read_records': counts['target'],
    'target_read_fraction': 0.1, 'background_read_records': counts['background'],
    'background_sam_sha256': background_original.hexdigest(),
    'all_background_records_byte_identical': True,
    'coordinates_flags_mate_fields_names_tags_unchanged': True,
    'all_cigar_query_and_reference_spans_150': True,
    'coordinate_order_unchanged': True, 'target_groups': {str(k): counts['target_group_'+str(k)] for k in range(3)}
}
(target/'independent-validation.json').write_text(json.dumps(verification, indent=2)+'\n')
manifest['files'] = {p.name: {'bytes':p.stat().st_size, 'sha256':sha(p)}
                     for p in sorted(target.iterdir()) if p.is_file() and p.name != 'manifest.json'}
(target/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(verification, indent=2))
