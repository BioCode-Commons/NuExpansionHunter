#!/usr/bin/env python3
"""Deterministic bounded streaming stress inputs, not representative WGS truth.

Each output directory contains a FASTA/FAI, catalog, SAM, unindexed BAM and
manifest. The converter must preserve input order without attempting to index
unsorted BAM; build sam_to_bam_streaming.c against the benchmark HTSlib.
No ExpansionHunter run or output filtering is performed here.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import random
import struct
import subprocess
import sys
from types import SimpleNamespace

HELPERS = Path(__file__).resolve().parents[1] / 'streaming_validation'
sys.path.insert(0, str(HELPERS))
import generate_catalog as generator
from compare_outputs import bam_header, bam_records

KINDS = ('shuffled', 'far-apart', 'missing-mates', 'dense-fanout',
         'hot-locus', 'stripe-collisions', 'empty-input')
STRIDE = 4096
READ_LENGTH = 150
CONTIG = 'chrS0001'


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def pair_rows(blocks, locus, index, mate_locus=None):
    """SAM SEQ is in reference orientation, including reverse-strand records."""
    mate_locus = locus if mate_locus is None else mate_locus
    jitter = index % 31 - 15
    left_local, right_local = 1873 + jitter, 2023 + jitter
    left = locus * STRIDE + left_local
    right = mate_locus * STRIDE + right_local
    span = abs(right - left) + READ_LENGTH
    tlen = span if left < right else -span
    name = f'e{locus:05d}_{index:08d}'
    proper = 2 if locus == mate_locus else 0
    for mate, block, local, pos, other, flag, length in (
        (0, locus, left_local, left, right, 97 | proper, tlen),
        (1, mate_locus, right_local, right, left, 145 | proper, -tlen)
    ):
        sequence = blocks[block][local:local + READ_LENGTH].decode()
        # A fixed fraction has one mismatch, keeping the workload from being
        # exclusively perfect-match shortcuts. M CIGAR allows substitutions.
        if (index + mate) % 3 == 0:
            chars = list(sequence)
            chars[50] = 'ACGT'[('ACGT'.index(chars[50]) + 1) % 4]
            sequence = ''.join(chars)
        quality = '*' if (index * 2 + mate + 1) % 17 == 0 else 'I' * READ_LENGTH
        yield [name, str(flag), CONTIG, str(pos + 1), '60', '150M', '=',
               str(other + 1), str(length), sequence, quality, 'RG:Z:synthetic']


def validate_sam(rows, reference_bases, expected_singletons, sorted_input):
    """Independent field/pair pass, including the exact pending mate high water."""
    pending, seen = {}, set()
    completed = peak = backwards = far_pairs = 0
    last_pos = -1
    for fields in rows:
        name, flag, pos, other, tlen = (fields[0], int(fields[1]), int(fields[3]),
                                      int(fields[7]), int(fields[8]))
        mate = 1 if flag & 0x80 else 0
        key = (name, mate)
        if key in seen:
            raise ValueError('Duplicate record identity: ' + str(key))
        seen.add(key)
        if flag not in (99, 147, 97, 145) or fields[2] != CONTIG or fields[6] != '=':
            raise ValueError('Invalid fixture pairing fields')
        if fields[5] != '150M' or len(fields[9]) != READ_LENGTH or set(fields[9]) - set('ACGT'):
            raise ValueError('Invalid sequence/CIGAR')
        if fields[10] != '*' and len(fields[10]) != READ_LENGTH:
            raise ValueError('Invalid QUAL length')
        if not (1 <= pos <= reference_bases - READ_LENGTH + 1
                and 1 <= other <= reference_bases - READ_LENGTH + 1):
            raise ValueError('Read or mate outside reference')
        backwards += pos < last_pos
        last_pos = pos
        if name in pending:
            old_flag, old_pos, old_other, old_tlen = pending.pop(name)
            if not (bool(old_flag & 0x80) != bool(flag & 0x80)
                    and old_pos == other and old_other == pos and old_tlen == -tlen
                    and bool(old_flag & 2) == bool(flag & 2)
                    and abs(tlen) == abs(pos - other) + READ_LENGTH):
                raise ValueError('Nonreciprocal mate fields: ' + name)
            completed += 1
            far_pairs += abs(pos - other) >= 1000
        else:
            pending[name] = (flag, pos, other, tlen)
            peak = max(peak, len(pending))
    if len(pending) != expected_singletons or (sorted_input and backwards):
        raise ValueError('Unexpected singleton count or sorting violation')
    return {'records': len(rows), 'completed_pairs': completed,
            'unmatched_records_at_eof': len(pending), 'peak_pending_mates': peak,
            'backwards_coordinate_steps': backwards, 'far_apart_pairs': far_pairs,
            'all_complete_pair_fields_reciprocal': True,
            'singletons_scope': 'Mate alignments intentionally absent from this input subset; mate fields remain valid.'}


def validate_bam(path, rows, sort_order, reference_bases):
    """Read every resulting record, checking core metadata, order and count."""
    with gzip.open(path, 'rb') as source:
        header = bam_header(source)
        if (f'SO:{sort_order}'.encode() not in header
                or b'@RG\tID:synthetic\tSM:synthetic' not in header
                or f'LN:{reference_bases}'.encode() not in header):
            raise ValueError('Unexpected BAM header')
        count = 0
        for count, record in enumerate(bam_records(source), 1):
            if count > len(rows):
                raise ValueError('Converter added records')
            fields = rows[count - 1]
            refid, pos, bin_mq_name, flag_cigar, length, mate_ref, mate_pos, tlen = struct.unpack_from(
                '<iiIIiiii', record, 4)
            name_size = bin_mq_name & 255
            name = record[36:36 + name_size].split(b'\0', 1)[0].decode()
            observed = (name, flag_cigar >> 16, pos + 1, mate_pos + 1, tlen)
            expected = (fields[0], int(fields[1]), int(fields[3]), int(fields[7]), int(fields[8]))
            if observed != expected or (refid, mate_ref, length) != (0, 0, READ_LENGTH):
                raise ValueError('BAM conversion changed order or record metadata')
        if count != len(rows):
            raise ValueError('Converter dropped records')
    return {'all_records_read': True, 'records': count, 'input_order_preserved': True,
            'core_pair_metadata_preserved': True, 'index_required': False}


def build(kind, args):
    directory = args.out / kind
    reference_loci = 1 if kind == 'dense-fanout' else 9 if kind == 'stripe-collisions' else 128
    config = SimpleNamespace(out=directory, loci=reference_loci, covered_loci=0,
        coverage_layout='prefix', pairs_per_locus=0, background_pairs_per_locus=0,
        read_length=READ_LENGTH, stride=STRIDE, loci_per_contig=reference_loci,
        motif='CAG', repeat_count=10, quality='high', mismatch_every=0,
        missing_quality_every=0, seed=args.seed, bam_helper=None)
    generator.write_fixture(config)
    blocks = [generator.locus_block(args.seed, i, STRIDE, 'CAG', 10) for i in range(reference_loci)]
    catalog = json.loads((directory / 'catalog.json').read_text())
    rows, expected_singletons, intended_deliveries = [], 0, None
    sorted_input = kind not in ('shuffled', 'stripe-collisions')
    description = ''
    if kind == 'dense-fanout':
        catalog = [dict(catalog[0], LocusId=f'FANOUT_{i:03d}') for i in range(24)]
        for pair in range(5000 * args.scale):
            rows.extend(pair_rows(blocks, 0, pair))
        intended_deliveries = 24 * 5000 * args.scale
        description = '24 distinct locus IDs share one coordinate; each complete pair routes to all24.'
    elif kind == 'hot-locus':
        for pair in range(20000 * args.scale):
            rows.extend(pair_rows(blocks, 0, pair))
        description = 'All20000*scale pairs at one locus; other127 catalog loci have no reads.'
    elif kind == 'stripe-collisions':
        source_catalog = catalog
        catalog = [dict(source_catalog[i // 4096 if i % 4096 == 0 else 8],
                        LocusId=f'STRIPE_{i:07d}') for i in range(7 * 4096 + 1)]
        for pair in range(1024 * args.scale):
            for locus in range(8):
                rows.extend(pair_rows(blocks, locus, pair))
        description = ('Interleaved pairs at catalog slots0,4096,...,28672 share a queue mutex stripe. '
                       'Other slots alias an uncovered ninth locus. Use2 workers to make8 active loci exceed '
                       'the7-queue limit; scheduling-dependent actual simultaneous activity is not guaranteed.')
    elif kind != 'empty-input':
        pairs_per_locus = (128 if kind == 'missing-mates' else 64) * args.scale
        for locus in range(reference_loci):
            for pair in range(pairs_per_locus):
                mate_locus = ((locus + reference_loci // 2) % reference_loci
                              if kind == 'far-apart' and pair % 8 else locus)
                pair_data = list(pair_rows(blocks, locus, pair, mate_locus))
                rows.extend(pair_data[:1] if kind == 'missing-mates' else pair_data)
        if kind == 'missing-mates':
            expected_singletons = len(rows)
            description = 'Only first-mate records retained in the INPUT fixture; no partner exists at EOF.'
        elif kind == 'far-apart':
            description = ('Seven eighths of pairs have reciprocal mate coordinates about262kb apart; '
                           'proper-pair bit cleared for remote pairs. One eighth local-pair controls retained.')
        else:
            random.Random(args.seed + ':shuffled-v1').shuffle(rows)
            description = 'All individual records shuffled deterministically; mates can arrive in either order.'
    else:
        description = 'Valid BAM header, reference dictionary and sample read group with zero records.'
    if sorted_input:
        rows.sort(key=lambda fields: (int(fields[3]), fields[0], int(fields[1])))
    sort_order = 'coordinate' if sorted_input else 'unsorted'
    reference_bases = reference_loci * STRIDE
    validation = validate_sam(rows, reference_bases, expected_singletons, sorted_input)
    if kind in ('shuffled', 'stripe-collisions') and not validation['backwards_coordinate_steps']:
        raise ValueError('Unsorted scenario failed to generate backwards coordinates')
    (directory / 'catalog.json').write_text(json.dumps(catalog, indent=2) + '\n')
    with (directory / 'reads.sam').open('w') as sam:
        sam.write(f'@HD\tVN:1.6\tSO:{sort_order}\n@SQ\tSN:{CONTIG}\tLN:{reference_bases}\n'
                  '@RG\tID:synthetic\tSM:synthetic\n')
        for fields in rows:
            sam.write('\t'.join(fields) + '\n')
    command = [str(args.sam_to_bam), str(directory / 'reads.sam'), str(directory / 'reads.bam')]
    converted = subprocess.run(command, capture_output=True, text=True)
    (directory / 'bam-conversion.log').write_text(converted.stdout + converted.stderr)
    if converted.returncode:
        raise RuntimeError(f'BAM converter failed: {converted.stderr}')
    validation['bam'] = validate_bam(directory / 'reads.bam', rows, sort_order, reference_bases)
    (directory / 'independent-validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    manifest = {
        'scope': 'Bounded synthetic streaming extreme diagnostics; not real WGS or genotype truth.',
        'fixture_kind': kind, 'description': description, 'loci': len(catalog),
        'reference_blocks': reference_loci, 'reference_bases': reference_bases,
        'read_records': len(rows), 'read_length': READ_LENGTH,
        'coordinate_sorted': sorted_input, 'sample_id': 'synthetic',
        'expected_unmatched_records': expected_singletons,
        'intended_locus_pair_deliveries': intended_deliveries,
        'recommended_compute_workers': 2 if kind == 'stripe-collisions' else 8,
        'required_region_extension_length': 1000,
        'parameters': {'seed': args.seed, 'scale': args.scale},
        'generation_command': [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        'conversion_command': command, 'generator_sha256': sha256(Path(__file__)),
        'base_generator_sha256': sha256(HELPERS / 'generate_catalog.py'),
        'converter_sha256': sha256(args.sam_to_bam), 'validation': validation,
        'limitations': ['No read drops or caps in the application benchmark; fixture scale is fixed before runs.',
                       'No guarantee of particular scheduling, DP fast-path or nonempty BAM outputs.',
                       'Missing-mates and empty-input deliberately have no complete read pair.',
                       'Single reference contig,150bp reads,10-CAG allele; no clinical/expanded-allele claim.'],
        'files': {p.name: {'bytes': p.stat().st_size, 'sha256': sha256(p)}
                  for p in sorted(directory.iterdir()) if p.is_file() and p.name != 'manifest.json'}
    }
    (directory / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'kind': kind, 'loci': len(catalog), 'records': len(rows),
                      'peak_pending_mates': validation['peak_pending_mates'], 'out': str(directory)}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--sam-to-bam', type=Path, required=True)
    parser.add_argument('--only', choices=KINDS, nargs='+', default=list(KINDS))
    parser.add_argument('--scale', type=int, default=1, help='Positive multiplier for read counts; catalogs stay fixed')
    parser.add_argument('--seed', default='nueh-streaming-extremes-20261114-v1')
    args = parser.parse_args()
    args.out = args.out.resolve()
    args.sam_to_bam = args.sam_to_bam.resolve()
    if args.scale < 1 or not args.sam_to_bam.is_file():
        parser.error('Require scale>=1 and an existing streaming SAM-to-BAM converter')
    if len(set(args.only)) != len(args.only):
        parser.error('--only kinds must be unique')
    args.out.mkdir(parents=True, exist_ok=True)
    for kind in args.only:
        if (args.out / kind).exists():
            parser.error('Refusing to overwrite existing fixture: ' + str(args.out / kind))
    for kind in args.only:
        build(kind, args)


if __name__ == '__main__':
    main()
