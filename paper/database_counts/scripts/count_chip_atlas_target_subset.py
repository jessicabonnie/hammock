#!/usr/bin/env python3
"""Count target-matched ChIP-Atlas experiments and BED files.

Downloads the public ChIP-Atlas experiment list, filters exact normalized targets,
and writes an auditable manifest plus assembly-level summary.

By default, the primary subset is human CTCF ChIP-seq on hg19 and hg38.
One experiment contributes one BED file at the selected q-value threshold.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

EXPERIMENT_URL = "https://chip-atlas.org/data/ExperimentList.json"
DEFAULT_TARGET = "CTCF"
DEFAULT_ASSEMBLIES = ("hg19", "hg38")
DEFAULT_AG_CLASS = "TFs and others"
DEFAULT_QVAL = "5"


def normalize(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def get_field(record: dict[str, Any], *names: str) -> str:
    lowered = {str(k).lower(): v for k, v in record.items()}
    for name in names:
        if name in record:
            return normalize(record[name])
        if name.lower() in lowered:
            return normalize(lowered[name.lower()])
    return ""


def load_records(url: str) -> list[dict[str, Any]]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "hammock-database-counts/1.0"},
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        payload = json.load(response)

    if isinstance(payload, list):
        records = payload
    elif isinstance(payload, dict):
        records = next(
            (payload[key] for key in ("experiments", "data", "results") if isinstance(payload.get(key), list)),
            None,
        )
        if records is None:
            raise ValueError(f"Unrecognized JSON object keys: {sorted(payload)}")
    else:
        raise ValueError(f"Unexpected JSON type: {type(payload).__name__}")

    if not all(isinstance(item, dict) for item in records):
        raise ValueError("Experiment list contains non-object records")
    return records


def canonical_record(record: dict[str, Any]) -> dict[str, str]:
    return {
        "experiment_id": get_field(record, "srx", "experiment_id", "expid", "id"),
        "sra_accession": get_field(record, "sra", "run", "sra_accession"),
        "geo_accession": get_field(record, "geo", "gsm", "geo_accession"),
        "assembly": get_field(record, "genome", "assembly", "genome_assembly"),
        "antigen_class": get_field(record, "agClass", "antigen_class"),
        "target_reported": get_field(record, "agSubClass", "antigen", "target"),
        "cell_type_class": get_field(record, "clClass", "cell_type_class"),
        "cell_type": get_field(record, "clSubClass", "cell_type"),
        "title": get_field(record, "title"),
    }


def bed_url(assembly: str, experiment_id: str, qval: str) -> str:
    # ChIP-Atlas publishes per-experiment BED4 files at multiple thresholds.
    # The exact URL is retained as a reproducible template rather than fetched.
    threshold = {"5": "05", "10": "10", "20": "20"}.get(qval, qval)
    return f"https://chip-atlas.dbcls.jp/data/{assembly}/eachData/bed{threshold}/{experiment_id}.{threshold}.bed"


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=DEFAULT_TARGET)
    parser.add_argument("--assemblies", nargs="+", default=list(DEFAULT_ASSEMBLIES))
    parser.add_argument("--antigen-class", default=DEFAULT_AG_CLASS)
    parser.add_argument("--qval", default=DEFAULT_QVAL, choices=("5", "10", "20"))
    parser.add_argument("--url", default=EXPERIMENT_URL)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results",
    )
    args = parser.parse_args()

    target = normalize(args.target).casefold()
    assemblies = set(args.assemblies)
    antigen_class = normalize(args.antigen_class).casefold()
    retrieved_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    records = load_records(args.url)
    selected: list[dict[str, Any]] = []
    missing_ids = 0

    for raw in records:
        row = canonical_record(raw)
        if row["assembly"] not in assemblies:
            continue
        if row["antigen_class"].casefold() != antigen_class:
            continue
        if row["target_reported"].casefold() != target:
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

    selected.sort(key=lambda row: (row["assembly"], row["cell_type_class"], row["cell_type"], row["experiment_id"]))
    if not selected:
        raise RuntimeError("No records matched. Inspect API schema and target labels before changing filters.")

    duplicate_ids = [key for key, count in Counter(row["experiment_id"] for row in selected).items() if count > 1]
    if duplicate_ids:
        raise RuntimeError(f"Duplicate experiment IDs after filtering: {duplicate_ids[:10]}")

    stem = f"chip_atlas_{args.target.lower()}_human"
    manifest_path = args.output_dir / f"{stem}_manifest.tsv"
    summary_path = args.output_dir / f"{stem}_summary.tsv"
    provenance_path = args.output_dir / f"{stem}_provenance.json"

    manifest_fields = [
        "experiment_id", "sra_accession", "geo_accession", "assembly", "assay",
        "antigen_class", "target_reported", "target_normalized", "cell_type_class",
        "cell_type", "title", "bed_threshold", "bed_file_count", "bed_url",
        "source_url", "retrieved_at_utc",
    ]
    write_tsv(manifest_path, selected, manifest_fields)

    summary_rows: list[dict[str, Any]] = []
    for assembly in sorted(assemblies):
        subset = [row for row in selected if row["assembly"] == assembly]
        n = len(subset)
        summary_rows.append(
            {
                "target": args.target.upper(),
                "species": "Homo sapiens",
                "assembly": assembly,
                "experiment_count": n,
                "bed_file_count": n,
                "cell_type_class_count": len({row["cell_type_class"] for row in subset if row["cell_type_class"]}),
                "cell_type_count": len({row["cell_type"] for row in subset if row["cell_type"]}),
                "within_assembly_pairwise_comparisons": math.comb(n, 2),
                "retrieved_at_utc": retrieved_at,
            }
        )

    total = len(selected)
    summary_rows.append(
        {
            "target": args.target.upper(),
            "species": "Homo sapiens",
            "assembly": "hg19+hg38",
            "experiment_count": total,
            "bed_file_count": total,
            "cell_type_class_count": len({row["cell_type_class"] for row in selected if row["cell_type_class"]}),
            "cell_type_count": len({row["cell_type"] for row in selected if row["cell_type"]}),
            "within_assembly_pairwise_comparisons": sum(math.comb(row["experiment_count"], 2) for row in summary_rows),
            "all_pairwise_comparisons_if_coordinates_were_comparable": math.comb(total, 2),
            "cross_assembly_pairs": math.prod([row["experiment_count"] for row in summary_rows]),
            "retrieved_at_utc": retrieved_at,
        }
    )
    summary_fields = list(summary_rows[-1].keys())
    write_tsv(summary_path, summary_rows, summary_fields)

    provenance = {
        "source": "ChIP-Atlas",
        "source_url": args.url,
        "retrieved_at_utc": retrieved_at,
        "filters": {
            "species": "Homo sapiens (inferred from hg19/hg38)",
            "assemblies": sorted(assemblies),
            "antigen_class_exact": args.antigen_class,
            "target_exact_case_insensitive": args.target,
            "q_value_threshold": f"q < 1e-{args.qval}",
        },
        "counting_unit": "one experiment-level peak BED per experiment at the selected threshold",
        "records_in_source": len(records),
        "records_selected": total,
        "records_excluded_for_missing_experiment_id": missing_ids,
        "outputs": [str(manifest_path), str(summary_path)],
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"Selected {total:,} {args.target.upper()} experiments")
    for row in summary_rows:
        print(row)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
