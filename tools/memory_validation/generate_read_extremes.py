#!/usr/bin/env python3
"""Generate deterministic, coordinate-sorted streaming read stress fixtures.

Each selected case writes reference.fa/.fai, catalog.json, reads.sam/.bam/.bai,
independent-validation.json and manifest.json in a new child of --out. This
does not run ExpansionHunter or time anything. Run generation outside measured
benchmark intervals. Synthetic cases exercise input extremes, not WGS truth.

Example:
  python3 tools/memory_validation/generate_read_extremes.py \
    --out work/read-extremes --sam-to-bam work/sam_to_bam
"""

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys


KINDS = ("longread1000", "higherror-indels150", "iupac-missingquality150", "expandedrepeat150")
MOTIF = "CAG"
REFERENCE_COPIES = 10
STRIDE = 4096
REPEAT_START = STRIDE // 2
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
    "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}
BASE_GENERATOR = Path(__file__).resolve().parents[1] / "streaming_validation/generate_catalog.py"
spec = importlib.util.spec_from_file_location("read_extremes_base_generator", BASE_GENERATOR)
base_generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base_generator)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cigar_operations(cigar):
    operations = [(int(n), op) for n, op in re.findall(r"(\d+)([MID])", cigar)]
    require(operations and all(n > 0 for n, _ in operations)
            and "".join(str(n) + op for n, op in operations) == cigar, "Invalid CIGAR: " + cigar)
    return operations


def spans(cigar):
    operations = cigar_operations(cigar)
    return (sum(n for n, op in operations if op in "MI"),
            sum(n for n, op in operations if op in "MD"))


def concrete(sequence, choice):
    return "".join(IUPAC[base][choice % len(IUPAC[base])] for base in sequence)


def substitute(sequence, offsets):
    result = list(sequence)
    for offset in offsets:
        result[offset] = "ACGT"[("ACGT".index(result[offset]) + 1) % 4]
    return "".join(result)


def apply_cigar(reference, cigar, key):
    insertion = hashlib.shake_256(key.encode()).digest(32).translate(base_generator.TRANSLATE_DNA).decode()
    ref_offset = inserted = 0
    result = []
    for length, op in cigar_operations(cigar):
        if op == "M":
            result.append(reference[ref_offset:ref_offset + length])
            ref_offset += length
        elif op == "D":
            ref_offset += length
        else:
            result.append(insertion[inserted:inserted + length])
            inserted += length
    require(ref_offset == len(reference), "CIGAR does not consume supplied reference span")
    return "".join(result)


def expanded_read(block, sample_copies, start, length):
    """Project a substring of one expanded haplotype onto the reference.

    The extra copies are inserted after the reference's ten copies. Both
    selected reads have a nonempty flank anchor, so no all-insertion mapped
    record or fabricated reference span is needed.
    """
    anchor = REPEAT_START + REFERENCE_COPIES * len(MOTIF)
    extra = MOTIF * (sample_copies - REFERENCE_COPIES)
    insertion_end = anchor + len(extra)
    haplotype = block[:anchor] + extra + block[anchor:]
    end = start + length
    query = haplotype[start:end]
    operations = []
    before = max(0, min(end, anchor) - start)
    inserted = max(0, min(end, insertion_end) - max(start, anchor))
    after = max(0, end - max(start, insertion_end))
    for count, op in ((before, "M"), (inserted, "I"), (after, "M")):
        if count:
            operations.append(str(count) + op)
    reference_start = start if start < anchor else max(anchor, start - len(extra))
    cigar = "".join(operations)
    require(len(query) == length and spans(cigar)[1] > 0, "Expanded read requires a reference anchor")
    return reference_start, cigar, query


def pair_records(kind, block, locus, pair, sample_copies, counts, calibration=False):
    read_length = 1000 if kind == "longread1000" else 150
    jitter = pair % 11 - 5
    if kind == "longread1000":
        starts = (REPEAT_START - 970 + jitter, REPEAT_START - 20 + jitter)
    else:
        starts = (REPEAT_START - 175 + jitter, REPEAT_START - 30 + jitter)
    if kind == "expandedrepeat150":
        starts = (REPEAT_START - 55 + jitter, REPEAT_START + sample_copies * len(MOTIF) - 95 + jitter)
    if calibration:
        # Both complete reads lie on the same unique flank. This supplies a
        # measured 300bp fragment length without relying on repeat projection.
        starts = (REPEAT_START - 700 + jitter, REPEAT_START - 550 + jitter)
    records = []
    for mate, start in enumerate(starts):
        cigar = str(read_length) + "M"
        quality = "I" * read_length
        if kind == "expandedrepeat150" and not calibration:
            start, cigar, sequence = expanded_read(block, sample_copies, start, read_length)
            counts["expanded_haplotype_reads"] += 1
        else:
            sequence = concrete(block[start:start + read_length], locus + pair + mate)
        if calibration:
            counts["calibration_reads"] += 1
        if kind == "longread1000":
            sequence = substitute(sequence, (330 + pair % 7, 760 + mate * 3))
            counts["substitutions"] += 2
        elif kind == "higherror-indels150":
            if (pair + mate) % 2:
                cigar = "40M3I30M4D30M1I46M"
                sequence = apply_cigar(sequence, cigar, f"indels:{locus}:{pair}:{mate}")
                counts["reads_with_indels"] += 1
                offsets = tuple(range(50, 150, 17))
            else:
                offsets = tuple(range(43, 150, 13))
            sequence = substitute(sequence, offsets)
            counts["substitutions"] += len(offsets)
        elif kind == "iupac-missingquality150":
            query = list(sequence)
            ambiguous = "RYSWKMBDHVN"
            for number, offset in enumerate((42, 67, 92, 127)):
                query[offset] = ambiguous[(pair + mate + number) % len(ambiguous)]
            sequence = "".join(query)
            if (pair + mate) % 2:
                quality = "*"
            else:
                qualities = ["I"] * read_length
                for offset, char in ((10, "!"), (42, "5"), (67, "6"), (92, "!"), (127, "5")):
                    qualities[offset] = char
                quality = "".join(qualities)
            counts["ambiguous_query_bases"] += 4
        require(spans(cigar)[0] == len(sequence) == read_length, "Query/CIGAR length disagreement")
        counts["missing_quality_reads"] += quality == "*"
        counts["records"] += 1
        records.append((start, cigar, sequence, quality))
    outer_start = min(record[0] for record in records)
    outer_end = max(record[0] + spans(record[1])[1] for record in records)
    return records, outer_end - outer_start


def validate(directory, reference, kind, loci, pairs_per_locus, calibration_pairs_per_locus=0):
    """Read the written SAM independently; check fields and paired spans."""
    pending = {}
    counts = Counter()
    previous = -1
    max_pending = 0
    expected_length = 1000 if kind == "longread1000" else 150
    with (directory / "reads.sam").open() as stream:
        for line in stream:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            require(len(fields) == 12 and fields[-1] == "RG:Z:synthetic", "Unexpected SAM fields/tags")
            name, flag = fields[0], int(fields[1])
            start, mate_start, template_length = int(fields[3]) - 1, int(fields[7]) - 1, int(fields[8])
            sequence, quality = fields[9:11]
            q_span, r_span = spans(fields[5])
            require(fields[2] == "chrExtreme" and fields[6] == "=" and flag in (99, 147), "Invalid pair flags/contig")
            require(previous <= start and start >= 0 and 0 < r_span and start + r_span <= len(reference),
                    "Unsorted or out-of-reference alignment")
            previous = start
            require(q_span == len(sequence) == expected_length and set(sequence) <= set(IUPAC), "Invalid query")
            require(quality == "*" or (len(quality) == q_span and all(33 <= ord(x) <= 126 for x in quality)),
                    "Invalid quality length/encoding")
            require(start // STRIDE == mate_start // STRIDE, "Mates cross a locus block")
            # Expanded reads must be exact reference edits under their M/I
            # CIGAR. Other cases deliberately carry substitution/base-call errors.
            if kind == "expandedrepeat150":
                q = 0
                r = start
                for count, op in cigar_operations(fields[5]):
                    if op == "M":
                        require(sequence[q:q + count] == reference[r:r + count], "Expanded anchor disagrees with reference")
                        q += count
                        r += count
                    elif op == "I":
                        q += count
                    else:
                        raise ValueError("Unexpected deletion in expanded haplotype fixture")
                if name.startswith("c"):
                    require(fields[5] == "150M"
                            and start % STRIDE + r_span < REPEAT_START
                            and mate_start % STRIDE + 150 < REPEAT_START,
                            "Calibration pair must lie wholly in the left flank")
                    counts["calibration_reads"] += 1
                else:
                    require("I" in fields[5], "Expanded repeat read did not cross inserted copies")
                    counts["expanded_haplotype_reads"] += 1
            if name in pending:
                old_flag, old_start, old_mate, old_span, old_tlen = pending.pop(name)
                require(old_flag != flag and old_start == mate_start and old_mate == start, "Nonreciprocal mate fields")
                require(old_tlen == -template_length
                        and abs(template_length) == max(start + r_span, old_start + old_span) - min(start, old_start),
                        "Incorrect reference-coordinate template length")
                counts["pairs"] += 1
            else:
                pending[name] = (flag, start, mate_start, r_span, template_length)
                max_pending = max(max_pending, len(pending))
            counts["records"] += 1
            counts["missing_quality_reads"] += quality == "*"
    require(not pending and counts["pairs"] == loci * (pairs_per_locus + calibration_pairs_per_locus),
            "Missing/extra mate records")
    require(counts["calibration_reads"] == loci * calibration_pairs_per_locus * 2,
            "Calibration record count mismatch")
    catalog = json.loads((directory / "catalog.json").read_text())
    require(len(catalog) == loci, "Catalog locus count mismatch")
    for index, entry in enumerate(catalog):
        begin = index * STRIDE + REPEAT_START
        require(entry["ReferenceRegion"] == f"chrExtreme:{begin}-{begin + 30}"
                and reference[begin:begin + 30] == MOTIF * REFERENCE_COPIES, "Catalog/reference disagreement")
    fai = (directory / "reference.fa.fai").read_text().rstrip().split("\t")
    require(fai[0] == "chrExtreme" and int(fai[1]) == len(reference) and list(map(int, fai[3:])) == [60, 61],
            "FAI metadata disagreement")
    with (directory / "reference.fa").open("rb") as fasta:
        for position in (0, 59, 60, REPEAT_START, len(reference) - 1):
            fasta.seek(int(fai[2]) + position // 60 * 61 + position % 60)
            require(fasta.read(1).decode() == reference[position], "FAI offset disagreement")
    return {"valid": True, "coordinate_sorted": True, "reciprocal_mates_and_exact_tlen": True,
            "cigar_query_and_reference_bounds_valid": True, "fai_checked": True,
            "expanded_reference_edits_checked": kind == "expandedrepeat150",
            "read_length": expected_length, "max_pending_pairs": max_pending, "counts": dict(counts)}


def build(args, kind):
    directory = args.out.resolve() / kind
    directory.mkdir(parents=True, exist_ok=False)
    blocks = []
    for locus in range(args.loci):
        block = base_generator.locus_block(args.seed, locus, STRIDE, MOTIF, REFERENCE_COPIES).decode()
        if kind == "iupac-missingquality150":
            chars = list(block)
            for position, base in ((REPEAT_START - 90, "R"), (REPEAT_START + 110, "Y"), (REPEAT_START + 350, "N")):
                chars[position] = base
            block = "".join(chars)
        blocks.append(block)
    reference = "".join(blocks)
    fasta = base_generator.FastaWriter(directory / "reference.fa")
    fasta.begin("chrExtreme", len(reference))
    fasta.append(reference.encode())
    fasta.end()
    fasta.close()
    catalog = [{"LocusId": f"EXT_{locus:07d}", "LocusStructure": "(CAG)*", "VariantType": "Repeat",
                "ReferenceRegion": f"chrExtreme:{locus * STRIDE + REPEAT_START}-{locus * STRIDE + REPEAT_START + 30}"}
               for locus in range(args.loci)]
    (directory / "catalog.json").write_text(json.dumps(catalog, indent=2) + "\n")
    counts = Counter()
    with (directory / "reads.sam").open("w") as sam:
        sam.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrExtreme\tLN:{len(reference)}\n@RG\tID:synthetic\tSM:synthetic\n")
        for locus, block in enumerate(blocks):
            rows = []
            for pair in range(args.pairs_per_locus + args.calibration_pairs_per_locus):
                calibration = pair >= args.pairs_per_locus
                records, template_length = pair_records(
                    kind, block, locus, pair, args.expanded_copies, counts, calibration)
                for mate, (start, cigar, sequence, quality) in enumerate(records):
                    position = locus * STRIDE + start
                    mate_position = locus * STRIDE + records[1 - mate][0]
                    flag = 99 if mate == 0 else 147
                    tlen = template_length if mate == 0 else -template_length
                    name_prefix = "c" if calibration else "p"
                    line = (f"{name_prefix}{locus:07d}_{pair:07d}\t{flag}\tchrExtreme\t{position + 1}\t60\t{cigar}"
                            f"\t=\t{mate_position + 1}\t{tlen}\t{sequence}\t{quality}\tRG:Z:synthetic\n")
                    rows.append((position, pair, mate, line))
            for _, _, _, line in sorted(rows):
                sam.write(line)
    validation = validate(directory, reference, kind, args.loci, args.pairs_per_locus,
                          args.calibration_pairs_per_locus)
    require(validation["counts"]["records"] == counts["records"], "Generator/validator record disagreement")
    conversion = None
    if not args.sam_only:
        conversion = [str(args.sam_to_bam.resolve()), str(directory / "reads.sam"), str(directory / "reads.bam")]
        result = subprocess.run(conversion, capture_output=True, text=True)
        (directory / "bam-conversion.log").write_text(result.stdout + result.stderr)
        require(result.returncode == 0 and (directory / "reads.bam.bai").is_file(), "BAM conversion/index failed")
        validation["bam_conversion_and_index_exit"] = result.returncode
    (directory / "independent-validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    details = {
        "longread1000": "1000bp mates entirely within default 1000bp extraction flanks; two substitutions per read.",
        "higherror-indels150": "150bp reads retain an initial exact seed; alternate 9 substitutions or balanced 4I/4D plus6 substitutions.",
        "iupac-missingquality150": "Separated R/Y/N reference bases; four query IUPAC calls per read; half QUAL missing, remaining include Phred0/20/21/40.",
        "expandedrepeat150": f"150bp reads from one {args.expanded_copies}-CAG haplotype versus reference10; both mates cross inserted repeat copies and retain a flank anchor.",
    }
    if args.calibration_pairs_per_locus:
        details[kind] += (f" Plus {args.calibration_pairs_per_locus} exact 150bp calibration pairs per locus,"
                          " both reads wholly on the left flank, supplying a 300bp fragment-length observation.")
    manifest = {
        "scope": "Synthetic streaming input extremes; not representative WGS or a guaranteed alignment fast-path/fallback witness",
        "fixture_kind": kind, "loci": args.loci, "covered_loci": args.loci, "contigs": 1,
        "read_length": validation["read_length"], "read_records": counts["records"],
        "target_read_records": counts["records"], "target_read_fraction": 1.0,
        "reference_bases": len(reference), "reference_repeat_count": REFERENCE_COPIES,
        "sample_repeat_count": args.expanded_copies if kind == "expandedrepeat150" else REFERENCE_COPIES,
        "reference_repeat_allele_only": kind != "expandedrepeat150", "coordinate_sorted": True,
        "required_region_extension_length": 1000, "details": details[kind], "counts": dict(counts),
        "pair_layout": "Flags99/147; SAM SEQ stored in reference-aligned orientation for both mates; TLEN uses reference-consuming CIGAR spans.",
        "expectations": {"all_input_fields_valid": True, "all_pairs_complete": True,
                         "catalog_reference_agree": True, "same_input_for_both_binaries": True,
                         "output_comparison": "Exact VCF/JSON and BAM complete-header/full-record multiset, preserving duplicate multiplicity; only BAM record order may differ",
                         "genotype_truth_asserted": False},
        "limitations": ["Read acceptance, genotyping and specific DP-path activation must be observed, not assumed.",
                        "Mismatch/IUPAC cases model read errors and uncertain base calls; no diploid genotype truth is asserted.",
                        "Expanded-repeat SAM alignments project one explicit haplotype; repeat mapping is intrinsically nonunique."],
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "generator_sha256": sha(Path(__file__)), "base_generator_sha256": sha(BASE_GENERATOR),
        "generation_command": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "conversion_command": conversion,
        "files": {path.name: {"bytes": path.stat().st_size, "sha256": sha(path)}
                  for path in sorted(directory.iterdir()) if path.is_file()},
    }
    if args.calibration_pairs_per_locus:
        manifest["expectations"]["observed_output_checks_required"] = {
            "minimum_coverage": 10, "positive_fragment_length": True,
            "called_repeat_genotype_at_every_locus": True,
            "check_with": "tools/memory_validation/validate_genotype_outputs.py",
        }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"fixture": kind, "loci": args.loci, "records": counts["records"], "valid": True}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--only", choices=KINDS)
    parser.add_argument("--loci", type=int, default=200)
    parser.add_argument("--pairs-per-locus", type=int, default=8)
    parser.add_argument("--calibration-pairs-per-locus", type=int, default=0,
                        help="Extra exact pairs wholly on the left flank; only with --only expandedrepeat150. "
                             "Genotyping profile: --loci 200 --pairs-per-locus 256 --calibration-pairs-per-locus 8")
    parser.add_argument("--expanded-copies", type=int, default=80)
    parser.add_argument("--seed", default="nueh-memory-read-extremes-v1")
    parser.add_argument("--sam-to-bam", type=Path, help="Converter built from tools/streaming_validation/sam_to_bam.c")
    parser.add_argument("--sam-only", action="store_true", help="Skip BAM conversion, for generator-only checks")
    args = parser.parse_args()
    if args.loci < 1 or args.pairs_per_locus < 1 or args.loci * STRIDE >= 1 << 29:
        parser.error("Positive locus/pair counts and a contig shorter than the BAI512Mb limit are required")
    if not 50 <= args.expanded_copies <= 1000:
        parser.error("--expanded-copies must be50..1000 so both150bp reads cross inserted repeat copies")
    if args.calibration_pairs_per_locus < 0 or (args.calibration_pairs_per_locus
                                              and args.only != "expandedrepeat150"):
        parser.error("Nonnegative --calibration-pairs-per-locus requires --only expandedrepeat150 when nonzero")
    if not args.sam_only and (args.sam_to_bam is None or not args.sam_to_bam.is_file()):
        parser.error("Supply an existing --sam-to-bam converter or use --sam-only")
    for kind in [args.only] if args.only else KINDS:
        build(args, kind)


if __name__ == "__main__":
    main()
