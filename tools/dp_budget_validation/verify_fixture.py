#!/usr/bin/env python3
"""Independently inspect generated BAM bytes, reference matches and per-locus counts."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import struct
import sys
sys.dont_write_bytecode = True
from compare_outputs import bam_header, bam_records


def main():
    p = argparse.ArgumentParser()
    p.add_argument('directory', type=Path)
    args = p.parse_args()
    root = args.directory
    manifest = json.loads((root / 'manifest.json').read_text())
    config = manifest['parameters']
    refs = {}
    with (root / 'reference.fa').open() as stream:
        name, parts = None, []
        for line in stream:
            if line.startswith('>'):
                if name is not None:
                    refs[name] = ''.join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
        refs[name] = ''.join(parts)
    names = list(refs)
    catalog = json.loads((root / 'catalog.json').read_text())
    assert len(catalog) == config['loci'] == config['covered_loci']
    for i, entry in enumerate(catalog):
        contig = names[i // config['loci_per_contig']]
        start = (i % config['loci_per_contig']) * config['stride'] + config['stride'] // 2
        end = start + len(config['motif']) * config['repeat_count']
        assert entry['ReferenceRegion'] == '%s:%d-%d' % (contig, start, end)
        assert refs[contig][start:end] == config['motif'] * config['repeat_count']
    counts, pairs = Counter(), Counter()
    records = mismatch_records = missing_quality = 0
    previous = (-1, -1)
    bases = '=ACMGRSVTWYHKDBN'
    with gzip.open(root / 'reads.bam', 'rb') as stream:
        header = bam_header(stream)
        for record in bam_records(stream):
            tid, pos, bin_mq_nl, flag_nc, length, mtid, mpos, tlen = struct.unpack('<iiIIiiii', record[4:36])
            name_length = bin_mq_nl & 255
            cigar_count, flag = flag_nc & 65535, flag_nc >> 16
            read_name = record[36:36 + name_length - 1].decode()
            locus, pair = map(int, read_name[1:].split('_'))
            assert tid == locus // config['loci_per_contig'] and tid == mtid
            assert flag in (99, 147) and length == config['read_length'] and cigar_count == 1
            assert previous <= (tid, pos)
            previous = (tid, pos)
            assert struct.unpack('<I', record[36 + name_length:40 + name_length])[0] == (length << 4)
            start = 36 + name_length + 4 * cigar_count
            packed = record[start:start + (length + 1) // 2]
            query = ''.join(bases[(packed[i // 2] >> (4 if i % 2 == 0 else 0)) & 15] for i in range(length))
            reference = refs[names[tid]][pos:pos + length]
            differences = [i for i, (a, b) in enumerate(zip(query, reference)) if a != b]
            mate = 0 if flag == 99 else 1
            expect_mismatch = bool(config['mismatch_every'] and (pair * 2 + mate + 1) % config['mismatch_every'] == 0)
            assert differences == ([length // 3] if expect_mismatch else [])
            mismatch_records += bool(differences)
            qual = record[start + (length + 1) // 2:start + (length + 1) // 2 + length]
            expect_missing = bool(config['missing_quality_every'] and (pair * 2 + mate + 1) % config['missing_quality_every'] == 0)
            assert (qual == bytes([255]) * length) == expect_missing
            missing_quality += expect_missing
            counts[locus] += 1
            pairs[(locus, pair, mate)] += 1
            records += 1
    assert records == manifest['read_records'] == 2 * config['covered_loci'] * config['pairs_per_locus']
    assert len(counts) == config['covered_loci'] and set(counts.values()) == {2 * config['pairs_per_locus']}
    assert len(pairs) == records and set(pairs.values()) == {1}
    result = {'valid': True, 'records': records, 'covered_loci': len(counts),
              'pairs': records // 2, 'all_records_have_one_substitution': mismatch_records == records,
              'mismatch_records': mismatch_records, 'missing_quality_records': missing_quality,
              'coordinate_sorted': True, 'catalog_reference_intervals_verified': len(catalog),
              'bam_header_bytes': len(header), 'scope': 'BAM bytes and reference/catalog consistency; no DP timing or genotype truth claim'}
    (root / 'independent-validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
