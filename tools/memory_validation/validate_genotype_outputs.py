#!/usr/bin/env python3
"""Assert that supplemental benchmark outputs actually exercise genotyping.

Accept result.json files or run directories (recursively finding result.json).
Every locus must have coverage >= --min-coverage, positive read/fragment length,
and a called repeat genotype. Corresponding VCF GT/REPCN fields are also checked.
This validates observed computation, not biological genotype truth or equality
between binaries; run_matrix.py performs the complete output comparison.

Example:
  python3 tools/memory_validation/validate_genotype_outputs.py \
    --expected-loci 200 --require-evidence --out work/genotype-checks.json \
    work/expanded-positive-benchmark
"""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys


COUNT_FIELDS = ("CountsOfFlankingReads", "CountsOfInrepeatReads", "CountsOfSpanningReads")
GENOTYPE = re.compile(r"\d+(?:[/|]\d+)*\Z")


def numeric(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def distribution(values):
    return {"min": min(values), "median": statistics.median(values), "max": max(values)} if values else None


def validate(path, args):
    failures = Counter()
    examples = {}

    def check(condition, category, detail):
        if not condition:
            failures[category] += 1
            examples.setdefault(category, [])
            if len(examples[category]) < 5:
                examples[category].append(detail)

    payload = path.read_bytes()
    data = json.loads(payload)
    loci = data.get("LocusResults")
    if not isinstance(loci, dict):
        raise ValueError("Missing LocusResults object")
    check(len(loci) == args.expected_loci, "locus_count", len(loci))
    coverage, fragments, read_lengths = [], [], []
    genotypes = Counter()
    repeat_variants = {}
    called_loci = evidence_loci = 0
    for locus_id, locus in loci.items():
        for field, values, minimum in (("Coverage", coverage, args.min_coverage),
                                        ("FragmentLength", fragments, 1), ("ReadLength", read_lengths, 1)):
            value = locus.get(field)
            check(numeric(value) and value >= minimum, field, {"locus": locus_id, "value": value})
            if numeric(value):
                values.append(value)
        variants = locus.get("Variants", {})
        repeat_count = called_count = 0
        has_evidence = False
        for variant_id, variant in variants.items():
            if variant.get("VariantType") != "Repeat":
                continue
            repeat_count += 1
            genotype = variant.get("Genotype")
            called = isinstance(genotype, str) and GENOTYPE.fullmatch(genotype) is not None
            check(called, "missing_or_invalid_genotype", {"locus": locus_id, "variant": variant_id, "value": genotype})
            check(isinstance(variant.get("GenotypeConfidenceInterval"), str)
                  and bool(variant["GenotypeConfidenceInterval"]), "missing_confidence_interval", variant_id)
            check(variant_id not in repeat_variants, "duplicate_variant_id", variant_id)
            repeat_variants[variant_id] = genotype
            if called:
                called_count += 1
                genotypes[genotype] += 1
            has_evidence |= any(isinstance(variant.get(field), str)
                                and variant[field] not in ("", "()") for field in COUNT_FIELDS)
        check(repeat_count > 0, "no_repeat_variant", locus_id)
        called_loci += repeat_count > 0 and called_count == repeat_count
        evidence_loci += has_evidence
        if args.require_evidence:
            check(has_evidence, "no_reported_repeat_evidence", locus_id)

    vcf = path.with_suffix(".vcf")
    vcf_seen = set()
    vcf_called = 0
    with vcf.open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise ValueError("Expected a single-sample VCF with ten columns")
            info = dict(value.split("=", 1) for value in fields[7].split(";") if "=" in value)
            variant_id = info.get("VARID")
            if variant_id not in repeat_variants:
                check(False, "unexpected_vcf_variant", variant_id)
                continue
            check(variant_id not in vcf_seen, "duplicate_vcf_variant", variant_id)
            vcf_seen.add(variant_id)
            sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            gt = sample.get("GT", "")
            called = GENOTYPE.fullmatch(gt) is not None
            check(called, "missing_vcf_gt", variant_id)
            vcf_called += called
            check(sample.get("REPCN") == repeat_variants[variant_id], "vcf_json_genotype_mismatch", variant_id)
    for variant_id in repeat_variants.keys() - vcf_seen:
        check(False, "missing_vcf_variant", variant_id)
    return {"json": str(path), "vcf": str(vcf), "json_sha256": hashlib.sha256(payload).hexdigest(),
            "valid": not failures, "loci": len(loci), "called_loci": called_loci,
            "repeat_variants": len(repeat_variants), "called_vcf_records": vcf_called,
            "loci_with_reported_repeat_evidence": evidence_loci,
            "coverage": distribution(coverage), "fragment_length": distribution(fragments),
            "read_length": distribution(read_lengths), "distinct_genotypes": len(genotypes),
            "most_common_genotypes": genotypes.most_common(10),
            "failure_counts": dict(failures), "failure_examples": examples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="Output JSON files or benchmark run directories")
    parser.add_argument("--expected-loci", type=int, required=True)
    parser.add_argument("--min-coverage", type=float, default=10)
    parser.add_argument("--require-evidence", action="store_true")
    parser.add_argument("--out", type=Path, help="Optional complete JSON validation report")
    args = parser.parse_args()
    if args.expected_loci < 1 or not math.isfinite(args.min_coverage) or args.min_coverage < 0:
        parser.error("Expected locus count must be positive and minimum coverage finite and nonnegative")
    paths = set()
    for source in args.inputs:
        if source.is_dir():
            discovered = list(source.rglob("result.json"))
            if not discovered:
                parser.error(f"No result.json outputs under {source}")
            paths.update(path.resolve() for path in discovered)
        elif source.is_file():
            paths.add(source.resolve())
        else:
            parser.error(f"Missing input: {source}")
    report = {"scope": "Observed genotype computation; no biological truth assertion or binary-output comparison",
              "expected_loci_per_file": args.expected_loci, "minimum_coverage": args.min_coverage,
              "require_reported_repeat_evidence": args.require_evidence, "files": []}
    for path in sorted(paths):
        try:
            record = validate(path, args)
        except (OSError, ValueError, TypeError, KeyError, AttributeError) as error:
            record = {"json": str(path), "valid": False, "error": str(error)}
        report["files"].append(record)
    report["all_valid"] = all(record["valid"] for record in report["files"])
    report["files_checked"] = len(report["files"])
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report if not args.out else {key: report[key] for key in ("files_checked", "all_valid")}, indent=2))
    return 0 if report["all_valid"] else 1


if __name__ == "__main__":
    sys.exit(main())
