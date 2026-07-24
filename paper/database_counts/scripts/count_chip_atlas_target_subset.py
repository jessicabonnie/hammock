#!/usr/bin/env python3
"""Count target-matched ChIP-Atlas experiments and verified BED files.

The ChIP-Atlas experiment list may report multiple comma-separated assemblies
for one experiment. Those labels are expanded into candidate experiment/assembly
pairs, but a physical BED file is counted only after its URL is verified.

Physical-file summaries retain every verified assembly-specific BED. Pairwise
comparison summaries use exactly one designated reference assembly per species,
because conventional coordinate-overlap comparisons require a shared reference.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

EXPERIMENT_URL = "https://chip-atlas.org/data/ExperimentList.json"
DEFAULT_TARGET = "CTCF"
DEFAULT_AG_CLASS = "TFs and others"
DEFAULT_QVAL = "5"
DEFAULT_WORKERS = 16
COMPACT_FIELDS = (
    "srx", "sra", "geo", "genome", "agClass", "agSubClass", "clClass", "clSubClass"
)

ASSEMBLY_SPECIES = {
    "hg19": "Homo sapiens", "hg38": "Homo sapiens",
    "mm9": "Mus musculus", "mm10": "Mus musculus",
    "rn6": "Rattus norvegicus",
    "dm3": "Drosophila melanogaster", "dm6": "Drosophila melanogaster",
    "ce10": "Caenorhabditis elegans", "ce11": "Caenorhabditis elegans",
    "sacCer3": "Saccharomyces cerevisiae", "danRer10": "Danio rerio",
    "xenTro9": "Xenopus tropicalis", "galGal5": "Gallus gallus",
    "TAIR10": "Arabidopsis thaliana",
}

# Exactly one reference per species is used for pairwise comparison counts.
# Older references remain in the physical BED inventory but are excluded from
# these workload calculations.
DEFAULT_COMPARISON_REFERENCE = {
    "Homo sapiens": "hg38",
    "Mus musculus": "mm10",
    "Drosophila melanogaster": "dm6",
    "Rattus norvegicus": "rn6",
    "Caenorhabditis elegans": "ce11",
    "Saccharomyces cerevisiae": "sacCer3",
    "Danio rerio": "danRer10",
    "Xenopus tropicalis": "xenTro9",
    "Gallus gallus": "galGal5",
    "Arabidopsis thaliana": "TAIR10",
}


def normalize(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_") or "all"


def get_field(record: dict[str, Any], *names: str) -> str:
    lowered = {str(k).lower(): v for k, v in record.items()}
    for name in names:
        if name in record:
            return normalize(record[name])
        if name.lower() in lowered:
            return normalize(lowered[name.lower()])
    return ""


def compact_row_to_record(row: list[Any] | tuple[Any, ...]) -> dict[str, Any]:
    if len(row) < len(COMPACT_FIELDS):
        raise ValueError(f"Compact row has {len(row)} fields; expected at least 8")
    return {field: row[index] for index, field in enumerate(COMPACT_FIELDS)}


def load_records(url: str) -> list[dict[str, Any]]:
    request = urllib.request.Request(
        url, headers={"User-Agent": "hammock-database-counts/1.4"}
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        payload = json.load(response)

    if isinstance(payload, list):
        raw_records: Any = payload
    elif isinstance(payload, dict):
        raw_records = next(
            (
                payload[key]
                for key in ("experiments", "data", "results")
                if isinstance(payload.get(key), list)
            ),
            None,
        )
        if raw_records is None:
            raise ValueError(f"Unrecognized JSON keys: {sorted(payload)}")
    else:
        raise ValueError(f"Unexpected JSON type: {type(payload).__name__}")

    records: list[dict[str, Any]] = []
    for item in raw_records:
        if isinstance(item, dict):
            records.append(item)
        elif isinstance(item, (list, tuple)):
            records.append(compact_row_to_record(item))
        else:
            raise ValueError(f"Unsupported record type: {type(item).__name__}")
    return records


def canonical_record(record: dict[str, Any]) -> dict[str, str]:
    return {
        "experiment_id": get_field(record, "srx", "experiment_id", "expid", "id"),
        "sra_accession": get_field(record, "sra", "run", "sra_accession"),
        "geo_accession": get_field(record, "geo", "gsm", "geo_accession"),
        "source_assembly_set": get_field(record, "genome", "assembly", "genome_assembly"),
        "antigen_class": get_field(record, "agClass", "antigen_class"),
        "target_reported": get_field(record, "agSubClass", "antigen", "target"),
        "cell_type_class": get_field(record, "clClass", "cell_type_class"),
        "cell_type": get_field(record, "clSubClass", "cell_type"),
        "title": get_field(record, "title"),
    }


def split_assemblies(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def bed_url(assembly: str, experiment_id: str, qval: str) -> str:
    threshold = {"5": "05", "10": "10", "20": "20"}[qval]
    return (
        f"https://chip-atlas.dbcls.jp/data/{assembly}/eachData/"
        f"bed{threshold}/{experiment_id}.{threshold}.bed"
    )


def verify_bed(url: str, timeout: int = 45) -> dict[str, Any]:
    headers = {"User-Agent": "hammock-database-counts/1.4"}
    requests = [
        ("HEAD", urllib.request.Request(url, headers=headers, method="HEAD")),
        (
            "GET_RANGE",
            urllib.request.Request(
                url, headers={**headers, "Range": "bytes=0-0"}, method="GET"
            ),
        ),
    ]
    errors: list[str] = []
    for method, request in requests:
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = int(response.status)
                if status in (200, 206):
                    return {
                        "bed_verified": True,
                        "verification_method": method,
                        "http_status": status,
                        "content_length": response.headers.get("Content-Length", ""),
                        "verification_error": "",
                    }
                errors.append(f"HTTP {status}")
        except urllib.error.HTTPError as exc:
            errors.append(f"HTTP {exc.code}")
            if exc.code == 404:
                break
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            errors.append(f"{type(exc).__name__}: {exc}")
    return {
        "bed_verified": False,
        "verification_method": "HEAD+GET_RANGE",
        "http_status": "",
        "content_length": "",
        "verification_error": " | ".join(errors),
    }


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def pair_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=DEFAULT_TARGET)
    parser.add_argument("--antigen-class", default=DEFAULT_AG_CLASS)
    parser.add_argument("--species", nargs="+", help="Optional exact species names")
    parser.add_argument("--assemblies", nargs="+", help="Optional exact assembly labels")
    parser.add_argument("--qval", default=DEFAULT_QVAL, choices=("5", "10", "20"))
    parser.add_argument("--url", default=EXPERIMENT_URL)
    parser.add_argument("--verification-workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--skip-bed-verification",
        action="store_true",
        help="Retain candidates but do not infer physical or pairwise counts.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results",
    )
    args = parser.parse_args()

    target = normalize(args.target).casefold()
    antigen_class = normalize(args.antigen_class).casefold()
    species_filter = {normalize(value).casefold() for value in args.species or []}
    assembly_filter = {normalize(value) for value in args.assemblies or []}
    retrieved_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    records = load_records(args.url)
    candidates: list[dict[str, Any]] = []
    missing_ids = 0
    missing_assemblies = 0

    for raw in records:
        base = canonical_record(raw)
        if base["antigen_class"].casefold() != antigen_class:
            continue
        if base["target_reported"].casefold() != target:
            continue
        if not base["experiment_id"]:
            missing_ids += 1
            continue
        assemblies = split_assemblies(base["source_assembly_set"])
        if not assemblies:
            missing_assemblies += 1
            continue
        for assembly in assemblies:
            species = ASSEMBLY_SPECIES.get(assembly, "unknown")
            if species_filter and species.casefold() not in species_filter:
                continue
            if assembly_filter and assembly not in assembly_filter:
                continue
            row = dict(base)
            row.update(
                {
                    "assembly": assembly,
                    "species": species,
                    "assembly_count_reported": len(assemblies),
                    "is_multi_assembly_record": len(assemblies) > 1,
                    "target_normalized": args.target.upper(),
                    "assay": "ChIP-seq",
                    "bed_threshold": f"q < 1e-{args.qval}",
                    "candidate_bed_url": bed_url(
                        assembly, base["experiment_id"], args.qval
                    ),
                    "is_designated_comparison_reference": (
                        DEFAULT_COMPARISON_REFERENCE.get(species) == assembly
                    ),
                    "source_url": args.url,
                    "retrieved_at_utc": retrieved_at,
                }
            )
            candidates.append(row)

    candidates = list(
        {
            (row["experiment_id"], row["assembly"]): row
            for row in candidates
        }.values()
    )
    if not candidates:
        raise RuntimeError("No candidate records matched")

    if args.skip_bed_verification:
        for row in candidates:
            row.update(
                {
                    "bed_verified": "not_checked",
                    "verification_method": "not_checked",
                    "http_status": "",
                    "content_length": "",
                    "verification_error": "",
                }
            )
    else:
        with ThreadPoolExecutor(
            max_workers=max(1, args.verification_workers)
        ) as pool:
            futures = {
                pool.submit(verify_bed, row["candidate_bed_url"]): i
                for i, row in enumerate(candidates)
            }
            for future in as_completed(futures):
                candidates[futures[future]].update(future.result())

    candidates.sort(
        key=lambda row: (
            row["species"], row["assembly"], row["cell_type_class"],
            row["cell_type"], row["experiment_id"],
        )
    )
    verified = [row for row in candidates if row["bed_verified"] is True]
    comparison_set = [
        row for row in verified if row["is_designated_comparison_reference"]
    ]

    scope_parts = [slugify(args.target)]
    if args.species:
        scope_parts.append(
            "species_" + "_".join(slugify(value) for value in args.species)
        )
    if args.assemblies:
        scope_parts.append(
            "assemblies_" + "_".join(slugify(value) for value in args.assemblies)
        )
    if not args.species and not args.assemblies:
        scope_parts.append("all_references")
    stem = "chip_atlas_" + "_".join(scope_parts)

    candidate_path = args.output_dir / f"{stem}_candidate_manifest.tsv"
    verified_path = args.output_dir / f"{stem}_verified_bed_manifest.tsv"
    physical_summary_path = args.output_dir / f"{stem}_physical_file_summary.tsv"
    comparison_path = args.output_dir / f"{stem}_single_reference_manifest.tsv"
    comparison_summary_path = (
        args.output_dir / f"{stem}_single_reference_pairwise_summary.tsv"
    )
    provenance_path = args.output_dir / f"{stem}_provenance.json"

    manifest_fields = [
        "experiment_id", "sra_accession", "geo_accession", "source_assembly_set",
        "assembly", "species", "assembly_count_reported", "is_multi_assembly_record",
        "is_designated_comparison_reference", "assay", "antigen_class",
        "target_reported", "target_normalized", "cell_type_class", "cell_type",
        "title", "bed_threshold", "candidate_bed_url", "bed_verified",
        "verification_method", "http_status", "content_length", "verification_error",
        "source_url", "retrieved_at_utc",
    ]
    write_tsv(candidate_path, candidates, manifest_fields)
    write_tsv(verified_path, verified, manifest_fields)
    if not args.skip_bed_verification:
        write_tsv(comparison_path, comparison_set, manifest_fields)

    candidate_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    verified_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in candidates:
        candidate_groups[(row["species"], row["assembly"])].append(row)
    for row in verified:
        verified_groups[(row["species"], row["assembly"])].append(row)

    physical_rows: list[dict[str, Any]] = []
    for species, assembly in sorted(candidate_groups):
        c = candidate_groups[(species, assembly)]
        v = verified_groups.get((species, assembly), [])
        physical_rows.append(
            {
                "row_type": "assembly",
                "target": args.target.upper(),
                "species": species,
                "assembly": assembly,
                "candidate_experiment_count": len({r["experiment_id"] for r in c}),
                "verified_experiment_count": (
                    len({r["experiment_id"] for r in v})
                    if not args.skip_bed_verification else "not_checked"
                ),
                "verified_bed_file_count": (
                    len(v) if not args.skip_bed_verification else "not_checked"
                ),
                "failed_or_missing_candidate_count": (
                    len(c) - len(v)
                    if not args.skip_bed_verification else "not_checked"
                ),
                "cell_type_class_count": len(
                    {r["cell_type_class"] for r in v if r["cell_type_class"]}
                ),
                "cell_type_count": len({r["cell_type"] for r in v if r["cell_type"]}),
                "retrieved_at_utc": retrieved_at,
            }
        )

    unique_candidates = {r["experiment_id"] for r in candidates}
    unique_verified = {r["experiment_id"] for r in verified}
    physical_rows.append(
        {
            "row_type": "total",
            "target": args.target.upper(),
            "species": "all selected species",
            "assembly": "all selected assemblies",
            "candidate_experiment_count": len(unique_candidates),
            "verified_experiment_count": (
                len(unique_verified)
                if not args.skip_bed_verification else "not_checked"
            ),
            "verified_bed_file_count": (
                len(verified) if not args.skip_bed_verification else "not_checked"
            ),
            "failed_or_missing_candidate_count": (
                len(candidates) - len(verified)
                if not args.skip_bed_verification else "not_checked"
            ),
            "cell_type_class_count": len(
                {r["cell_type_class"] for r in verified if r["cell_type_class"]}
            ),
            "cell_type_count": len(
                {r["cell_type"] for r in verified if r["cell_type"]}
            ),
            "retrieved_at_utc": retrieved_at,
        }
    )
    write_tsv(physical_summary_path, physical_rows, list(physical_rows[-1].keys()))

    comparison_rows: list[dict[str, Any]] = []
    if not args.skip_bed_verification:
        comparison_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in comparison_set:
            comparison_groups[row["species"]].append(row)

        for species, rows in sorted(comparison_groups.items()):
            assembly = DEFAULT_COMPARISON_REFERENCE[species]
            n = len(rows)
            comparison_rows.append(
                {
                    "row_type": "species",
                    "target": args.target.upper(),
                    "species": species,
                    "comparison_assembly": assembly,
                    "bed_file_count": n,
                    "unique_experiment_count": len(
                        {row["experiment_id"] for row in rows}
                    ),
                    "cell_type_class_count": len(
                        {
                            row["cell_type_class"]
                            for row in rows
                            if row["cell_type_class"]
                        }
                    ),
                    "cell_type_count": len(
                        {row["cell_type"] for row in rows if row["cell_type"]}
                    ),
                    "within_reference_pairwise_comparisons": pair_count(n),
                    "retrieved_at_utc": retrieved_at,
                }
            )

        comparison_rows.append(
            {
                "row_type": "total",
                "target": args.target.upper(),
                "species": "all species; not pooled for comparison",
                "comparison_assembly": "one designated assembly per species",
                "bed_file_count": len(comparison_set),
                "unique_experiment_count": len(
                    {row["experiment_id"] for row in comparison_set}
                ),
                "cell_type_class_count": len(
                    {
                        row["cell_type_class"]
                        for row in comparison_set
                        if row["cell_type_class"]
                    }
                ),
                "cell_type_count": len(
                    {row["cell_type"] for row in comparison_set if row["cell_type"]}
                ),
                "within_reference_pairwise_comparisons": sum(
                    pair_count(len(rows)) for rows in comparison_groups.values()
                ),
                "retrieved_at_utc": retrieved_at,
            }
        )
        write_tsv(
            comparison_summary_path,
            comparison_rows,
            list(comparison_rows[-1].keys()),
        )

    provenance = {
        "source": "ChIP-Atlas",
        "source_url": args.url,
        "retrieved_at_utc": retrieved_at,
        "filters": {
            "species": args.species or "all",
            "assemblies": args.assemblies or "all",
            "antigen_class_exact": args.antigen_class,
            "target_exact_case_insensitive": args.target,
            "q_value_threshold": f"q < 1e-{args.qval}",
        },
        "counting_policy": {
            "physical_inventory": "all verified experiment/assembly BED files",
            "pairwise_comparison": (
                "only verified BED files on one designated reference per species; "
                "different references and different species are not paired"
            ),
            "designated_comparison_references": DEFAULT_COMPARISON_REFERENCE,
            "no_fallback_policy": (
                "an experiment missing from the designated reference is excluded "
                "from pairwise counts rather than replaced with another assembly"
            ),
        },
        "verification_performed": not args.skip_bed_verification,
        "verification_workers": args.verification_workers,
        "records_in_source": len(records),
        "unique_target_experiments": len(unique_candidates),
        "candidate_experiment_assembly_pairs": len(candidates),
        "verified_bed_files": len(verified) if not args.skip_bed_verification else None,
        "single_reference_comparison_files": (
            len(comparison_set) if not args.skip_bed_verification else None
        ),
        "records_excluded_for_missing_experiment_id": missing_ids,
        "records_excluded_for_missing_assembly": missing_assemblies,
        "verification_failure_summary": dict(
            Counter(
                row["verification_error"]
                for row in candidates
                if row["bed_verified"] is False
            )
        ),
        "outputs": [
            str(candidate_path), str(verified_path), str(physical_summary_path),
            str(comparison_path), str(comparison_summary_path),
        ],
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

    print(f"Unique target-matched experiments: {len(unique_candidates):,}")
    print(f"Candidate experiment/assembly pairs: {len(candidates):,}")
    if args.skip_bed_verification:
        print("BED verification skipped; physical and pairwise counts are not inferred.")
    else:
        print(f"Verified BED files across all references: {len(verified):,}")
        print(f"Single-reference comparison files: {len(comparison_set):,}")
        for row in comparison_rows:
            if row["row_type"] == "species":
                print(
                    f"{row['species']} ({row['comparison_assembly']}): "
                    f"{row['bed_file_count']:,} files, "
                    f"{row['within_reference_pairwise_comparisons']:,} pairs"
                )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
