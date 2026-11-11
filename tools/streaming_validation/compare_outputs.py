#!/usr/bin/env python3
"""Exact VCF/JSON and BAM header+full-record multiset comparison.

BAM records are externally sorted by their complete byte representation.
Multiplicity is preserved, and acceptance compares actual bytes, not hashes.
Raw compressed BAM byte identity is diagnostic because record order may differ.
"""
import argparse
import gzip
import hashlib
import heapq
import json
from pathlib import Path
import struct
import sys
import tempfile


def files_identical(left, right):
    left, right = Path(left), Path(right)
    if left.stat().st_size != right.stat().st_size:
        return False
    with left.open('rb') as a, right.open('rb') as b:
        while True:
            first, second = a.read(1024*1024), b.read(1024*1024)
            if first != second:
                return False
            if not first:
                return True


def exact_read(stream, size):
    result = stream.read(size)
    if len(result) != size:
        raise ValueError(f'Truncated BAM/sort run: expected {size}, got {len(result)}')
    return result


def bam_header(stream):
    magic = exact_read(stream, 4)
    if magic != b'BAM\x01':
        raise ValueError('Not a BAM stream')
    size_bytes = exact_read(stream, 4)
    size = struct.unpack('<i', size_bytes)[0]
    if size < 0:
        raise ValueError('Negative BAM text length')
    parts = [magic, size_bytes, exact_read(stream, size)]
    count_bytes = exact_read(stream, 4)
    count = struct.unpack('<i', count_bytes)[0]
    if count < 0:
        raise ValueError('Negative BAM reference count')
    parts.append(count_bytes)
    for _ in range(count):
        name_size = exact_read(stream, 4)
        size = struct.unpack('<i', name_size)[0]
        if size < 1:
            raise ValueError('Invalid BAM reference name length')
        parts.extend([name_size, exact_read(stream, size), exact_read(stream, 4)])
    return b''.join(parts)


def bam_records(stream):
    while True:
        length = stream.read(4)
        if not length:
            return
        if len(length) != 4:
            raise ValueError('Truncated BAM block length')
        size = struct.unpack('<i', length)[0]
        if size < 32:
            raise ValueError('Invalid BAM record block length')
        yield length + exact_read(stream, size)


def write_run(path, records):
    with path.open('wb') as out:
        for record in records:
            out.write(struct.pack('<I', len(record)))
            out.write(record)


def read_run(path):
    with path.open('rb') as stream:
        while True:
            length = stream.read(4)
            if not length:
                return
            if len(length) != 4:
                raise ValueError('Truncated sort run')
            yield exact_read(stream, struct.unpack('<I', length)[0])


def canonical_records(path, directory, chunk_bytes):
    paths, chunk, used, count = [], [], 0, 0
    with gzip.open(path, 'rb') as stream:
        header = bam_header(stream)
        for record in bam_records(stream):
            # Approximate Python cache budget includes bytes object and list pointer.
            # One unusually large record can exceed a chunk; only that record is held.
            charge = sys.getsizeof(record) + 8
            if chunk and used + charge > chunk_bytes:
                run = directory / f'run-{len(paths):06d}.bin'
                chunk.sort()
                write_run(run, chunk)
                paths.append(run)
                chunk, used = [], 0
            chunk.append(record)
            used += charge
            count += 1
        if chunk or not paths:
            run = directory / f'run-{len(paths):06d}.bin'
            chunk.sort()
            write_run(run, chunk)
            paths.append(run)
    # Limit simultaneously open file descriptors and merge heap size.
    generation = 0
    while len(paths) > 32:
        merged = []
        for start in range(0, len(paths), 32):
            group = paths[start:start+32]
            target = directory / f'merge-{generation:03d}-{len(merged):06d}.bin'
            write_run(target, heapq.merge(*(read_run(p) for p in group)))
            for old in group:
                old.unlink()
            merged.append(target)
        paths = merged
        generation += 1
    return header, count, heapq.merge(*(read_run(p) for p in paths))


def compare_bam(left, right, scratch, chunk_bytes=16*1024*1024):
    scratch = Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='bam-multiset-', dir=scratch) as temp:
        temp = Path(temp)
        a_dir, b_dir = temp/'a', temp/'b'
        a_dir.mkdir()
        b_dir.mkdir()
        ha, ca, a = canonical_records(Path(left), a_dir, chunk_bytes)
        hb, cb, b = canonical_records(Path(right), b_dir, chunk_bytes)
        same = ca == cb
        first_mismatch = None
        ah, bh = hashlib.sha256(), hashlib.sha256()
        sentinel = object()
        index = 0
        try:
            while True:
                ar, br = next(a, sentinel), next(b, sentinel)
                if ar is sentinel and br is sentinel:
                    break
                if ar is not sentinel:
                    ah.update(struct.pack('<Q',len(ar)))
                    ah.update(ar)
                if br is not sentinel:
                    bh.update(struct.pack('<Q',len(br)))
                    bh.update(br)
                if ar != br:
                    same = False
                    if first_mismatch is None:
                        first_mismatch = index
                index += 1
        finally:
            a.close()
            b.close()
        return {'header_identical': ha == hb, 'record_counts': [ca, cb],
                'record_multiset_identical': same, 'first_sorted_record_mismatch': first_mismatch,
                'canonical_diagnostic_sha256': [ah.hexdigest(),bh.hexdigest()],
                'compressed_bytes_identical': files_identical(left,right),
                'comparison': 'Exact complete BAM header and record bytes; external lexicographic sort preserves duplicates; SHA is diagnostic only'}


def compare_prefixes(left, right, scratch, chunk_bytes=16*1024*1024):
    left, right = str(left), str(right)
    result = {'vcf_identical': files_identical(left+'.vcf',right+'.vcf'),
              'json_identical': files_identical(left+'.json',right+'.json'),
              'bam': compare_bam(left+'_realigned.bam',right+'_realigned.bam',scratch,chunk_bytes)}
    result['outputs_identical_under_user_policy'] = (result['vcf_identical'] and result['json_identical']
        and result['bam']['header_identical'] and result['bam']['record_multiset_identical'])
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline-prefix', required=True)
    p.add_argument('--candidate-prefix', required=True)
    p.add_argument('--scratch', required=True)
    p.add_argument('--chunk-mib', type=int, default=16)
    p.add_argument('--json', type=Path)
    args = p.parse_args()
    if args.chunk_mib < 1:
        p.error('chunk-mib must be positive')
    result = compare_prefixes(args.baseline_prefix,args.candidate_prefix,args.scratch,args.chunk_mib*1024*1024)
    if args.json:
        args.json.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
    return 0 if result['outputs_identical_under_user_policy'] else 1


if __name__ == '__main__':
    sys.exit(main())
