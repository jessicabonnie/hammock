#!/usr/bin/env python3
"""Count target-matched ChIP-Atlas experiments and BED files.

Downloads the public ChIP-Atlas experiment list, filters an exact normalized
molecular target, and writes an auditable manifest plus assembly-level summary.

All assemblies are retained by default. Optional species and assembly filters
produce narrower derived views without changing the default global analysis.
One experiment contributes one experiment-level BED file at the selected
q-value threshold.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

EXPERIMENT_URL = "https://chip-atlas.org/data/ExperimentList.json"
DEFAULT_TARGET = "CTCF"
DEFAULT_AG_CLASS = "TFs and others"
DEFAULT_QVAL = "5"

# The current ChIP-Atlas endpoint returns compact positional rows under
# {"data": string[][]}. This order is documented in the ChIP-Atlas client.
COMPACT_FIELDS = (
    "srx",
    "sra",
    "geo",
    "genome",
    "agClass",
    "agSubClass",
    "clClass",
    "clSubClass",
)

# ChIP-Atlas assembly labels currently encountered or historically documented.
# Unknown labels are preserved and assigned species="unknown" rather than dropped.
ASSEMBLY_SPECIES = {
    "hg19": "Homo sapiens",
    "hg38": "Homo sapiens",
    "mm9": "Mus musculus",
    "mm10": "Mus musculus",
    "rn6": "Rattus norvegicus",
    "dm3": "Drosophila melanogaster",
    "dm6": "Drosophila melanogaster",
    "ce10": "Caenorhabditis elegans",
    "ce11": "Caenorhabditis elegans",
    "sacCer3": "Saccharomyces cerevisiae",
    "danRer10": "Danio rerio",
    "xenTro9": "Xenopus tropicalis",
    "galGal5": "Gallus gallus",
    "TAIR10": "Arabidopsis thaliana",
}


def normalize(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")
    return slug or "all"


def get_field(record: dict[str, Any], *names: str) -> str:
    lowered = {str(k).lower(): v for k, v in record.items()}
    for name in names:
        if name in record:
            return normalize(record[name])
        if name.lower() in lowered:
            return normalize(lowered[name.lower()])
    return ""


def compact_row_to_record(row: list[Any] | tuple[Any, ...]) -> dict[str, Any]:
    """Convert the positional ChIP-Atlas API representation to a mapping."""
    if len(row) < len(COMPACT_FIELDS):
        raise ValueError(
            "Compact experiment row has fewer than 8 fields: "
            f"received {len(row)} fields; first values={list(row)[:4]!r}"
        )
    return {field: row[index] for index, field in enumerate(COMPACT_FIELDS)}


def load_records(url: str) -> list[dict[str, Any]]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "hammock-database-counts/1.2"},
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
            raise ValueError(f"Unrecognized JSON object keys: {sorted(payload)}")
    else:
        raise ValueError(f"Unexpected JSON type: {type(payload).__name__}")

    records: list[dict[str, Any]] = []
    for index, item in enumerate(raw_records):
        if isinstance(item, dict):
            records.append(item)
        elif isinstance(item, (list, tuple)):
            records.append(compact_row_to_record(item))
        else:
            raise ValueError(
                "Experiment list contains an unsupported record type at "
                f"index {index}: {type(item).__name__}; value={item!r}"
            )
    return records


def infer_species(assembly: str) -> str:
    return ASSEMBLY_SPECIES.get(assembly, "unknown")


def canonical_record(record: dict[str, Any]) -> dict[str, str]:
    assembly = get_field(record, "genome", "assembly", "genome_assembly")
    return {
        "experiment_id": get_field(record, "srx", "experiment_id", "expid", "id"),
        "sra_accession": get_field(record, "sra", "run", "sra_accession"),
        "geo_accession": get_field(record, "geo", "gsm", "geo_accession"),
        "species": infer_species(assembly),
        "assembly": assembly,
        "antigen_class": get_field(record, "agClass", "antigen_class"),
        "target_reported": get_field(record, "agSubClass", "antigen", "target"),
        "cell_type_class": get_field(record, "clClass", "cell_type_class"),
        "cell_type": get_field(record, "clSubClass", "cell_type"),
        "title": get_field(record, "title"),
    }


def bed_url(assembly: str, experiment_id: str, qval: str) -> str:
    threshold = {"5": "05", "10": "10", "20": "20"}.get(qval, qval)
    return f"https://chip-atlas.dbcls.jp/data/{assembly}/eachData/bed{threshold}/{experiment_id}.{threshold}.bed"


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def pair_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=DEFAULT_TARGET)
    parser.add_argument("--antigen-class", default=DEFAULT_AG_CLASS)
    parser.add_argument("--species", nargs="+", help="Optional exact species names to retain")
    parser.add_argument("--assemblies", nargs="+", help="Optional exact assembly labels to retain")
    parser.add_argument("--qval", default=DEFAULT_QVAL, choices=("5", "10", "20"))
    parser.add_argument("--url", default=EXPERIMENT_URL)
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
    selected: list[dict[str, Any]] = []
    missing_ids = 0
    missing_assemblies = 0

    for raw in records:
        row = canonical_record(raw)
        if row["antigen_class"].casefold() != antigen_class:
            continue
        if row["target_reported"].casefold() != target:
            continue
        if not row["assembly"]:
            missing_assemblies += 1
            continue
        if species_filter and row["species"].casefold() not in species_filter:
            continue
        if assembly_filter and row["assembly"] not in assembly_filter:
            continue
        if not row["experiment_id"]:
            missing_ids += 1
            continue

        row.update(
            {
                "target_normalized": args.target.upper(),
                "assay": "ChIP-seq",
                "bed_threshold": f"q < 1e-{args.qval}",
                "bed_file_count": 1,
                "bed_url": bed_url(row["assembly"], row["experiment_id"], args.qval),
                "source_url": args.url,
                "retrieved_at_utc": retrieved_at,
            }
        )
        selected.append(row)

    selected.sort(
        key=lambda row: (
            row["species"], row["assembly"], row["cell_type_class"],
            row["cell_type"], row["experiment_id"],
        )
    )
    if not selected:
        raise RuntimeError("No records matched. Inspect API schema, target labels, and optional filters.")

    duplicate_ids = [key for key, count in Counter(row["experiment_id"] for row in selected).items() if count > 1]
    if duplicate_ids:
        raise RuntimeError(f"Duplicate experiment IDs after filtering: {duplicate_ids[:10]}")

    scope_parts = [slugify(args.target)]
    if args.species:
        scope_parts.append("species_" + "_".join(slugify(value) for value in args.species))
    if args.assemblies:
        scope_parts.append("assemblies_" + "_".join(slugify(value) for value in args.assemblies))
    if not args.species and not args.assemblies:
        scope_parts.append("all_references")
    stem = "chip_atlas_" + "_".join(scope_parts)

    manifest_path = args.output_dir / f"{stem}_manifest.tsv"
    summary_path = args.output_dir / f"{stem}_summary.tsv"
    provenance_path = args.output_dir / f"{stem}_provenance.json"

    manifest_fields = [
        "experiment_id", "sra_accession", "geo_accession", "species", "assembly",
        "assay", "antigen_class", "target_reported", "target_normalized",
        "cell_type_class", "cell_type", "title", "bed_threshold", "bed_file_count",
        "bed_url", "source_url", "retrieved_at_utc",
    ]
    write_tsv(manifest_path, selected, manifest_fields)

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in selected:
        grouped[(row["species"], row["assembly"])].append(row)

    summary_rows: list[dict[str, Any]] = []
    for (species, assembly), subset in sorted(grouped.items()):
        n = len(subset)
        summary_rows.append(
            {
                "row_type": "assembly",
                "target": args.target.upper(),
                "species": species,
                "assembly": assembly,
                "experiment_count": n,
                "bed_file_count": n,
                "cell_type_class_count": len({row["cell_type_class"] for row in subset if row["cell_type_class"]}),
                "cell_type_count": len({row["cell_type"] for row in subset if row["cell_type"]}),
                "within_assembly_pairwise_comparisons": pair_count(n),
                "all_pairwise_comparisons_if_coordinates_were_comparable": "",
                "cross_assembly_pairs": "",
                "retrieved_at_utc": retrieved_at,
            }
        )

    total = len(selected)
    within_assembly_pairs = sum(pair_count(len(subset)) for subset in grouped.values())
    all_pairs = pair_count(total)
    summary_rows.append(
        {
            "row_type": "total",
            "target": args.target.upper(),
            "species": "all selected species",
            "assembly": "all selected assemblies",
            "experiment_count": total,
            "bed_file_count": total,
            "cell_type_class_count": len({row["cell_type_class"] for row in selected if row["cell_type_class"]}),
            "cell_type_count": len({row["cell_type"] for row in selected if row["cell_type"]}),
            "within_assembly_pairwise_comparisons": within_assembly_pairs,
            "all_pairwise_comparisons_if_coordinates_were_comparable": all_pairs,
            "cross_assembly_pairs": all_pairs - within_assembly_pairs,
            "retrieved_at_utc": retrieved_at,
        }
    )
    summary_fields = list(summary_rows[-1].keys())
    write_tsv(summary_path, summary_rows, summary_fields)

    provenance = {
        "source": "ChIP-Atlas",
        "source_url": args.url,
        "source_json_shape": "object with data as compact positional rows",
        "compact_field_order": list(COMPACT_FIELDS),
        "retrieved_at_utc": retrieved_at,
        "filters": {
            "species": args.species or "all species inferred from assembly labels",
            "assemblies": args.assemblies or "all known non-empty assembly labels",
            "antigen_class_exact": args.antigen_class,
            "target_exact_case_insensitive": args.target,
            "q_value_threshold": f"q < 1e-{args.qval}",
        },
        "counting_unit": "one experiment-level peak BED per experiment at the selected threshold",
        "species_inference": ASSEMBLY_SPECIES,
        "unknown_assembly_species_policy": "retain the assembly and label species as unknown",
        "records_in_source": len(records),
        "records_selected": total,
        "assemblies_selected": sorted({row["assembly"] for row in selected}),
        "species_selected": sorted({row["species"] for row in selected}),
        "records_excluded_for_missing_experiment_id": missing_ids,
        "target_records_excluded_for_missing_assembly": missing_assemblies,
        "outputs": [str(manifest_path), str(summary_path)],
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"Selected {total:,} {args.target.upper()} experiments across {len(grouped):,} assembly groups")
    for row in summary_rows:
        print(row)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
