#!/usr/bin/env python3
"""Deterministic independent-coordinate STR catalog/reference/SAM fixtures.

Synthetic workload, not WGS truth or a substitute for real clinical data.
Writes FASTA with an explicit FAI, streaming catalog JSON and coordinate-sorted SAM.
No BAM conversion is performed unless a separately built helper is supplied.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

TRANSLATE_DNA = bytes(b'ACGT'[i & 3] for i in range(256))


class FastaWriter:
    def __init__(self, path, line_width=60):
        self.handle = path.open('wb')
        self.index = Path(str(path) + '.fai').open('w')
        self.line_width = line_width
        self.column = 0

    def begin(self, name, length):
        assert self.column == 0
        self.handle.write(('>' + name + '\n').encode())
        self.index.write(f'{name}\t{length}\t{self.handle.tell()}\t{self.line_width}\t{self.line_width+1}\n')

    def append(self, data):
        offset = 0
        while offset < len(data):
            take = min(self.line_width - self.column, len(data) - offset)
            self.handle.write(data[offset:offset+take])
            offset += take
            self.column += take
            if self.column == self.line_width:
                self.handle.write(b'\n')
                self.column = 0

    def end(self):
        if self.column:
            self.handle.write(b'\n')
            self.column = 0

    def close(self):
        self.handle.close()
        self.index.close()


def locus_block(seed, index, stride, motif, repeat_count):
    data = bytearray(hashlib.shake_256(f'{seed}:{index}:reference'.encode()).digest(stride).translate(TRANSLATE_DNA))
    start = stride // 2
    repeat = motif.encode() * repeat_count
    data[start:start+len(repeat)] = repeat
    return bytes(data)


def quality_sequence(length, mode):
    if mode == 'high':
        return 'I' * length
    chars = ['I'] * length
    # Phred 20 is lowercase, 21 uppercase; both branches occur deterministically.
    for i in range(0, length, 47):
        chars[i] = '5'
    for i in range(1, length, 53):
        chars[i] = '6'
    for i in range(2, length, 101):
        chars[i] = '!'
    return ''.join(chars)


def is_covered(index, loci, covered, layout):
    if layout == 'prefix':
        return index < covered
    # Exactly covered selected indices, with no large set allocation.
    return ((index + 1) * covered // loci) != (index * covered // loci)


def background_window(args, base):
    """Pair starts pass the1KiB coarse mask but neither read hits exact target.

    LocusSpecDecoding.cpp extends this single reference region by1000 on each
    side at the default CLI setting. GenomeMask::addRegion marks start>>10
    through stop>>10 inclusive, including the partial final bin. Keep both
    reads beyond the exact stop but wholly inside that marked bin and block.
    """
    repeat_end = base + args.stride // 2 + len(args.motif) * args.repeat_count
    exact_stop = repeat_end + 1000
    first = exact_stop + 1
    limit = min(((exact_stop >> 10) + 1) << 10, base + args.stride)
    span = 2 * args.read_length + 20
    if first + span > limit:
        raise ValueError(
            'Cannot place background mate pair after exact target stop within '
            f'the same1KiB mask bin: block_base={base}, exact_stop={exact_stop}, '
            f'available={max(0, limit-first)}bp, required={span}bp. '
            'Change read-length/stride/repeat-count or use background-pairs-per-locus=0.')
    return first - base, limit - base - span


def write_fixture(args):
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    contigs = []
    remaining = args.loci
    while remaining:
        n = min(args.loci_per_contig, remaining)
        contigs.append((f'chrS{len(contigs)+1:04d}', n))
        remaining -= n
    fasta = FastaWriter(out / 'reference.fa')
    records = 0
    quality = quality_sequence(args.read_length, args.quality)
    catalog_path = out / 'catalog.json'
    with catalog_path.open('w') as catalog, (out / 'reads.sam').open('w') as sam:
        sam.write('@HD\tVN:1.6\tSO:coordinate\n')
        for name, count in contigs:
            sam.write(f'@SQ\tSN:{name}\tLN:{count * args.stride}\n')
        sam.write('@RG\tID:synthetic\tSM:synthetic\n')
        catalog.write('[\n')
        global_index = 0
        for contig, count in contigs:
            fasta.begin(contig, count * args.stride)
            for local_index in range(count):
                block = locus_block(args.seed, global_index, args.stride, args.motif, args.repeat_count)
                fasta.append(block)
                base = local_index * args.stride
                repeat_start = base + args.stride // 2
                repeat_end = repeat_start + len(args.motif) * args.repeat_count
                entry = {'LocusId': f'SYN_{global_index:07d}', 'LocusStructure': f'({args.motif})*',
                         'ReferenceRegion': f'{contig}:{repeat_start}-{repeat_end}', 'VariantType': 'Repeat'}
                if global_index:
                    catalog.write(',\n')
                catalog.write(json.dumps(entry, separators=(',', ':')))
                rows = []
                if is_covered(global_index, args.loci, args.covered_loci, args.coverage_layout):
                    for pair in range(args.pairs_per_locus):
                        jitter = pair % 31 - 15
                        left = args.stride // 2 - args.read_length - 25 + jitter
                        # At1000bp keep the whole left read within the default
                        # 1000bp extraction flank; permit slight mate overlap.
                        if args.read_length > 900:
                            left = args.stride // 2 - args.read_length + 20 + jitter
                        right = args.stride // 2 - 25 + jitter
                        fragment_length = right + args.read_length - left
                        for mate, start, mate_start, flag, tlen in [
                            (0, left, right, 99, fragment_length),
                            (1, right, left, 147, -fragment_length),
                        ]:
                            read = bytearray(block[start:start+args.read_length])
                            assert len(read) == args.read_length
                            if args.mismatch_every and (pair * 2 + mate + 1) % args.mismatch_every == 0:
                                offset = min(len(read)-1, len(read)//3)
                                read[offset] = b'ACGT'[(b'ACGT'.index(read[offset]) + 1) % 4]
                            qual = '*' if args.missing_quality_every and (pair * 2 + mate + 1) % args.missing_quality_every == 0 else quality
                            name = f'p{global_index:07d}_{pair:07d}'
                            line = f'{name}\t{flag}\t{contig}\t{base+start+1}\t60\t{args.read_length}M\t=\t{base+mate_start+1}\t{tlen}\t{read.decode()}\t{qual}\tRG:Z:synthetic\n'
                            rows.append((start, pair, mate, line))
                if args.background_pairs_per_locus:
                    first, last = background_window(args, base)
                    for pair in range(args.background_pairs_per_locus):
                        left = first + pair % (last - first + 1)
                        right = left + args.read_length + 20
                        fragment_length = right + args.read_length - left
                        for mate, start, mate_start, flag, tlen in [
                            (0, left, right, 99, fragment_length),
                            (1, right, left, 147, -fragment_length),
                        ]:
                            read = block[start:start+args.read_length]
                            assert len(read) == args.read_length
                            assert base + start > repeat_end + 1000
                            assert ((base + start) >> 10) == ((repeat_end + 1000) >> 10)
                            qual = '*' if args.missing_quality_every and (pair * 2 + mate + 1) % args.missing_quality_every == 0 else quality
                            name = f'b{global_index:07d}_{pair:07d}'
                            line = f'{name}\t{flag}\t{contig}\t{base+start+1}\t60\t{args.read_length}M\t=\t{base+mate_start+1}\t{tlen}\t{read.decode()}\t{qual}\tRG:Z:synthetic\n'
                            rows.append((start, args.pairs_per_locus + pair, mate, line))
                rows.sort(key=lambda row: row[:3])
                for _, _, _, line in rows:
                    sam.write(line)
                    records += 1
                global_index += 1
            fasta.end()
        catalog.write('\n]\n')
    fasta.close()
    manifest = {
        'scope': 'Synthetic, independent loci and distinct deterministic flanks; not a WGS performance claim',
        'parameters': {key: str(value) if isinstance(value, Path) else value for key,value in vars(args).items()},
        'loci': args.loci, 'covered_loci': args.covered_loci, 'contigs': len(contigs),
        'reference_bases': args.loci * args.stride, 'read_records': records,
        'reference_repeat_allele_only': True, 'coordinate_sorted': True,
        'pair_layout': 'Forward-strand aligned SAM SEQ for both records; flags99/147; right mate covers repeat; left mate covers flank',
        'files': {},
    }
    if args.background_pairs_per_locus:
        manifest['target_read_records'] = 2 * args.covered_loci * args.pairs_per_locus
        manifest['background_read_records'] = 2 * args.loci * args.background_pairs_per_locus
        manifest['background_scope'] = (
            'Synthetic front-end separation experiment, not typical WGS. Background pairs '
            'exist at every locus independent of covered-loci; they pass the1KiB coarse '
            'mask but have neither endpoint in any exact target extraction region. '
            'Requires the application default --region-extension-length1000.')
        manifest['pair_layout'] += '; additional b-prefixed background pairs lie beyond exact target stop inside its final coarse-mask bin'
    if args.bam_helper:
        command = [str(args.bam_helper.resolve()), str(out / 'reads.sam'), str(out / 'reads.bam')]
        subprocess.run(command, check=True)
        manifest['conversion_command'] = command
    for path in sorted(out.iterdir()):
        if path.is_file():
            h = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(1024*1024), b''):
                    h.update(block)
            manifest['files'][path.name] = {'bytes': path.stat().st_size, 'sha256': h.hexdigest()}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({key: manifest[key] for key in ('loci','covered_loci','contigs','reference_bases','read_records')}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--loci', type=int, default=1000)
    p.add_argument('--covered-loci', type=int)
    p.add_argument('--coverage-layout', choices=['even','prefix'], default='even')
    p.add_argument('--pairs-per-locus', type=int, default=16)
    p.add_argument('--background-pairs-per-locus', type=int, default=0,
                   help='Additional pairs per ALL loci: pass1KiB coarse mask, miss exact targets; requires extension1000')
    p.add_argument('--read-length', type=int, default=150)
    p.add_argument('--stride', type=int, default=4096)
    p.add_argument('--loci-per-contig', type=int, default=1000)
    p.add_argument('--motif', default='CAG')
    p.add_argument('--repeat-count', type=int, default=10)
    p.add_argument('--quality', choices=['high','mixed'], default='mixed')
    p.add_argument('--mismatch-every', type=int, default=0)
    p.add_argument('--missing-quality-every', type=int, default=0)
    p.add_argument('--seed', default='nueh-round5-v1')
    p.add_argument('--bam-helper', type=Path)
    args = p.parse_args()
    if args.covered_loci is None:
        args.covered_loci = args.loci
    if args.loci < 1 or not 0 <= args.covered_loci <= args.loci or min(args.pairs_per_locus,args.background_pairs_per_locus) < 0:
        p.error('Invalid locus/coverage/read counts')
    if not 1 <= args.loci_per_contig or args.loci_per_contig * args.stride >= (1 << 29):
        p.error('Each contig must be nonempty and shorter than the BAI 512Mb limit')
    if not 80 <= args.read_length <= 1000 or args.stride < 2 * max(1100, args.read_length + 80):
        p.error('Reads need80..1000bp and each block needs >=1000bp graph flanks plus read bounds')
    if not args.motif or set(args.motif) - set('ACGT') or args.repeat_count < 1:
        p.error('Invalid repeat motif/count')
    if len(args.motif) * args.repeat_count > min(100, args.read_length - 60):
        p.error('Repeat must fit into right reads with informative flanks')
    if args.mismatch_every < 0 or args.missing_quality_every < 0:
        p.error('Every-N options must be nonnegative')
    if args.background_pairs_per_locus:
        # Preflight every possible block alignment within a1KiB bin before
        # creating any output. Contig coordinates restart at0; at most1024
        # residues need checking even for hundreds of thousands of loci.
        period = 1024 // math.gcd(args.stride, 1024)
        try:
            for local_index in range(min(args.loci,args.loci_per_contig,period)):
                background_window(args, local_index * args.stride)
        except ValueError as error:
            p.error(str(error))
    write_fixture(args)


if __name__ == '__main__':
    main()
