#!/usr/bin/env python3
"""Build deterministic streaming diagnostics, not representative WGS or allele truth."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

HELPERS = Path(__file__).resolve().parents[1] / 'streaming_validation'
sys.path.insert(0, str(HELPERS))
import generate_catalog as generator


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def add_outside_mask_reads(args, config):
    """Insert paired records in each block's unmarked first 1 KiB bin.

    The existing default 4 KiB block has its TR at 2048 and a 1000 bp
    extraction extension, so only bins 1..3 are marked. Both 150 bp
    mates fit completely in bin 0. Original and inserted rows remain sorted.
    """
    source = args.out / 'reads.sam'
    target = args.out / 'reads.with-background.sam'
    quality = generator.quality_sequence(150, config.quality)
    global_index = 0
    with source.open() as old, target.open('w') as new:
        pending = old.readline()
        while pending.startswith('@'):
            new.write(pending)
            pending = old.readline()
        remaining = config.loci
        contig_index = 0
        while remaining:
            count = min(config.loci_per_contig, remaining)
            contig_index += 1
            contig = f'chrS{contig_index:04d}'
            for local_index in range(count):
                base = local_index * config.stride
                block = generator.locus_block(config.seed, global_index, config.stride,
                                              config.motif, config.repeat_count)
                rows = []
                for pair in range(args.background_pairs):
                    left = 32 + pair % 256
                    right = left + 170
                    assert right + 150 < 1024
                    for mate, start, mate_start, flag, tlen in (
                        (0, left, right, 99, 320), (1, right, left, 147, -320)
                    ):
                        name = f'o{global_index:07d}_{pair:07d}'
                        read = block[start:start + 150].decode()
                        qual = '*' if (pair * 2 + mate + 1) % 17 == 0 else quality
                        line = (f'{name}\t{flag}\t{contig}\t{base+start+1}\t60\t150M\t=\t'
                                f'{base+mate_start+1}\t{tlen}\t{read}\t{qual}\tRG:Z:synthetic\n')
                        rows.append((base + start + 1, line))
                while pending:
                    fields = pending.split('\t', 4)
                    if fields[2] != contig or int(fields[3]) > base + config.stride:
                        break
                    rows.append((int(fields[3]), pending))
                    pending = old.readline()
                for _, line in sorted(rows, key=lambda row: row[0]):
                    new.write(line)
                global_index += 1
            remaining -= count
        if pending:
            raise RuntimeError('Original SAM was not fully consumed')
    target.replace(source)


def replace_candidate_sequences(args):
    """Replace both reads in pairs with pair_index % 10 != 0.

    Coordinates, flags, quality, CIGAR and all other SAM fields stay unchanged.
    The unrelated sequences test candidate rejection; chance k-mer matches or
    accepted alignments are possible and are not excluded from this fixture.
    """
    source = args.out / 'reads.sam'
    target = args.out / 'reads.with-unrelated-sequences.sam'
    changed_records = 0
    selected_pairs = 0
    total_records = 0
    with source.open() as old, target.open('w') as new:
        for line in old:
            if line.startswith('@'):
                new.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            total_records += 1
            if not fields[0].startswith('p'):
                raise ValueError('Expected only target-layout pairs in candidate-negative fixture')
            pair_index = int(fields[0].rsplit('_', 1)[1])
            if pair_index % 10 != 0:
                flag = int(fields[1])
                mate = 1 if flag & 0x80 else 0
                if flag & 0x40:
                    selected_pairs += 1
                key = f'{args.seed}:{fields[0]}:{mate}:candidate-negative'.encode()
                sequence = hashlib.shake_256(key).digest(len(fields[9])).translate(generator.TRANSLATE_DNA).decode()
                changed_records += sequence != fields[9]
                fields[9] = sequence
                line = '\t'.join(fields) + '\n'
            new.write(line)
    target.replace(source)
    return {
        'pair_selection': 'pair_index % 10 != 0; index restarts at zero at each covered locus',
        'nominal_replaced_pair_fraction': 0.9,
        'selected_pairs': selected_pairs,
        'sequence_changed_records': changed_records,
        'unchanged_positive_control_records': total_records - 2 * selected_pairs,
        'target_layout_read_records': total_records,
        'actual_selected_pair_fraction': (2 * selected_pairs / total_records) if total_records else None,
        'sequence_generation': 'Independent SHAKE-256 deterministic ACGT sequence for each selected read',
        'scope': 'Unrelated sequence at target coordinates; not guaranteed to fail every orientation/seed/alignment test.'
    }


def vary_quality_scores(args):
    """Vary QUAL per read without changing sequences, records, or missing QUAL."""
    source = args.out / 'reads.sam'
    target = args.out / 'reads.with-variable-qualities.sam'
    translation = bytes(33 + 20 + byte % 21 for byte in range(256))
    records = missing = changed = quality_bases = 0
    with source.open() as old, target.open('w') as new:
        for line in old:
            if line.startswith('@'):
                new.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            records += 1
            if fields[10] == '*':
                missing += 1
            else:
                mate = 1 if int(fields[1]) & 0x80 else 0
                key = f'{args.seed}:{fields[0]}:{mate}:variable-qualities-v1'.encode()
                quality = hashlib.shake_256(key).digest(len(fields[9])).translate(translation).decode('ascii')
                changed += quality != fields[10]
                quality_bases += len(quality)
                fields[10] = quality
                line = '\t'.join(fields) + '\n'
            new.write(line)
    target.replace(source)
    return {
        'model': 'SHAKE-256 per read; Q = 20 + (digest byte % 21); SAM Phred+33',
        'key_utf8': '{seed}:{qname}:{mate}:variable-qualities-v1',
        'mate_encoding': '1 when SAM flag 0x80 is set, otherwise 0',
        'digest_length_bytes': 'SAM sequence length',
        'seed': args.seed,
        'minimum_phred': 20,
        'maximum_phred': 40,
        'scope': 'Synthetic independent per-base qualities; modulo mapping is slightly nonuniform, not an empirical sequencing-quality model.',
        'read_records': records,
        'missing_quality_records_preserved': missing,
        'quality_changed_records': changed,
        'quality_bases': quality_bases
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind', choices=('frontend', 'target', 'outside', 'mixed', 'candidate-negative'), required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--bam-helper', type=Path, required=True)
    parser.add_argument('--loci', type=int, default=10000)
    parser.add_argument('--covered-loci', type=int, default=1000)
    parser.add_argument('--target-pairs', type=int, default=50)
    parser.add_argument('--background-pairs', type=int, default=100)
    parser.add_argument('--near-background-pairs', type=int, default=0,
                        help='For mixed: background pairs per locus inside the coarse mask; remaining background pairs are outside')
    parser.add_argument('--seed', default='nueh-frontend-20261112')
    parser.add_argument('--variable-qualities', action='store_true',
                        help='Use reproducible per-read SHAKE-256 Phred20..40 qualities; preserve missing QUAL')
    args = parser.parse_args()
    args.out = args.out.resolve()
    if args.loci < 1 or not 0 <= args.covered_loci <= args.loci:
        parser.error('Require loci > 0 and 0 <= covered-loci <= loci')
    if min(args.target_pairs, args.background_pairs) < 0:
        parser.error('Pair counts must be nonnegative')
    if not 0 <= args.near_background_pairs <= args.background_pairs:
        parser.error('Require 0 <= near-background-pairs <= background-pairs')
    if args.kind != 'mixed' and args.near_background_pairs:
        parser.error('--near-background-pairs requires --kind mixed')
    if not args.bam_helper.is_file():
        parser.error('BAM helper does not exist')
    config = SimpleNamespace(out=args.out, loci=args.loci, covered_loci=args.covered_loci,
        coverage_layout='even', pairs_per_locus=args.target_pairs,
        background_pairs_per_locus=(args.background_pairs if args.kind == 'frontend'
                                    else args.near_background_pairs if args.kind == 'mixed' else 0),
        read_length=150, stride=4096, loci_per_contig=1000, motif='CAG', repeat_count=10,
        quality='mixed', mismatch_every=1, missing_quality_every=17, seed=args.seed,
        bam_helper=None)
    generator.write_fixture(config)
    negative_metadata = None
    if args.kind in ('outside', 'mixed'):
        # Both background layouts use disjoint name prefixes and remain sorted.
        all_background_pairs = args.background_pairs
        if args.kind == 'mixed':
            args.background_pairs -= args.near_background_pairs
        add_outside_mask_reads(args, config)
        args.background_pairs = all_background_pairs
    elif args.kind == 'candidate-negative':
        negative_metadata = replace_candidate_sequences(args)
    quality_metadata = vary_quality_scores(args) if args.variable_qualities else None
    command = [str(args.bam_helper.resolve()), str(args.out / 'reads.sam'), str(args.out / 'reads.bam')]
    subprocess.run(command, check=True)
    manifest = json.loads((args.out / 'manifest.json').read_text())
    targets = 2 * args.covered_loci * args.target_pairs
    background = 2 * args.loci * args.background_pairs if args.kind in ('frontend', 'outside', 'mixed') else 0
    manifest.update(fixture_kind=args.kind, target_read_records=targets,
                    background_read_records=background, read_records=targets+background,
                    conversion_command=command,
                    frontend_generator_sha256=sha256(Path(__file__)),
                    base_generator_sha256=sha256(HELPERS / 'generate_catalog.py'))
    manifest['background_scope'] = {
        'frontend': 'All background pairs pass the 1 KiB mask but miss exact extraction regions.',
        'outside': 'All background pairs have both mates wholly outside the coarse mask.',
        'mixed': 'Background pairs split between coarse-mask/exact-region misses and wholly outside-mask reads.',
        'target': 'No background reads; target-heavy control.',
        'candidate-negative': 'No coordinate-background reads; selected target-layout pairs have unrelated sequences.'
    }[args.kind]
    manifest['target_read_fraction'] = targets / (targets + background) if targets + background else None
    if args.kind == 'mixed':
        manifest['near_background_read_records'] = 2 * args.loci * args.near_background_pairs
        manifest['outside_background_read_records'] = background - manifest['near_background_read_records']
    manifest['limitations'] = [
        'Synthetic diagnostics only: not real WGS, clinical validation, or expanded-allele truth.',
        'Requires region-extension-length 1000; 150 bp reads and fixed 4096 bp blocks.',
        'Reference 10-CAG allele only; no remote mates, unmapped reads, or structural variants.',
        'Large catalogs have many empty loci; no claim that a WGS population has these proportions.'
    ]
    if negative_metadata is not None:
        manifest['candidate_negative'] = negative_metadata
        manifest['limitations'].append(
            'Selected reads intentionally contradict their reference-aligned coordinates; diagnostic only. '
            'Unrelated sequence may still pass k-mer filters by chance; this is not a measured 90% rejection rate.')
        manifest['pair_layout'] += '; pair_index % 10 != 0 has both sequences replaced; other SAM fields unchanged'
    if quality_metadata is not None:
        manifest['variable_qualities'] = quality_metadata
        manifest['parameters']['quality'] = 'variable-phred20-40'
        manifest['frontend_generation_command'] = [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]
    manifest['files'] = {path.name: {'bytes': path.stat().st_size, 'sha256': sha256(path)}
                         for path in sorted(args.out.iterdir())
                         if path.is_file() and path.name != 'manifest.json'}
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'kind': args.kind, 'loci': args.loci, 'target_records': targets,
                      'background_records': background, 'out': str(args.out)}), flush=True)


if __name__ == '__main__':
    main()
