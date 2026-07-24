#!/usr/bin/env python3
"""Count native Roadmap hg19 and BLUEPRINT hg38 H3K27ac BED peak files.

Uses FILER's native-build metadata templates. Lifted records are intentionally
excluded so that cross-resource pairs represent a genuine reference mismatch.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import statistics
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HG19_METADATA_URL = "https://tf.lisanwanglab.org/GADB/metadata/filer.latest.hg19.template"
HG38_METADATA_URL = "https://tf.lisanwanglab.org/GADB/metadata/filer.latest.hg38.template"


def normalize(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def normalized_key(value: str) -> str:
    return "".join(ch for ch in normalize(value).casefold() if ch.isalnum())


def fetch_tsv(url: str) -> list[dict[str, str]]:
    req = urllib.request.Request(url, headers={"User-Agent": "hammock-database-counts/1.0"})
    with urllib.request.urlopen(req, timeout=180) as response:
        text = response.read().decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not reader.fieldnames:
        raise ValueError(f"No header found in {url}")
    return [{normalize(k): normalize(v) for k, v in row.items()} for row in reader]


def lookup(row: dict[str, str], *names: str) -> str:
    by_key = {normalized_key(k): v for k, v in row.items()}
    for name in names:
        value = by_key.get(normalized_key(name))
        if value is not None:
            return value
    return ""


def to_int(value: str) -> int | None:
    text = normalize(value).replace(",", "")
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def select(rows: list[dict[str, str]], source: str, assembly: str) -> list[dict[str, Any]]:
    selected: dict[str, dict[str, Any]] = {}
    for row in rows:
        data_source = lookup(row, "data_source", "data source")
        genome_build = lookup(row, "genome_build", "genome build")
        assay = lookup(row, "assay")
        antibody = lookup(row, "antibody", "target")
        output_type = lookup(row, "output_type", "output type")
        file_format = lookup(row, "file_format", "file format")
        download_url = lookup(row, "processed_file_download_url", "processed file download url")

        if data_source.casefold() != source.casefold():
            continue
        if genome_build.casefold() not in {assembly.casefold(), ("grch37" if assembly == "hg19" else "grch38")}:
            continue
        if normalized_key(assay) != "chipseq":
            continue
        if normalized_key(antibody) != "h3k27ac":
            continue
        if "peak" not in output_type.casefold():
            continue
        if "bed" not in file_format.casefold() and not download_url.casefold().endswith((".bed", ".bed.gz")):
            continue
        if not download_url:
            continue

        canonical = {
            "resource": source,
            "assembly": assembly,
            "identifier": lookup(row, "identifier"),
            "cell_type": lookup(row, "cell_type", "cell type"),
            "tissue_category": lookup(row, "tissue_category", "tissue category"),
            "biosample_type": lookup(row, "biosample_type", "biosample type"),
            "assay": assay,
            "target": "H3K27ac",
            "output_type": output_type,
            "file_format": file_format,
            "number_of_intervals": to_int(lookup(row, "number_of_intervals", "number of intervals")),
            "processed_file_download_url": download_url,
            "release_date": lookup(row, "release_date", "release date"),
        }
        selected[download_url] = canonical
    return sorted(selected.values(), key=lambda r: (r["tissue_category"], r["cell_type"], r["identifier"]))


def pair_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


def summarize(resource: str, assembly: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    intervals = [r["number_of_intervals"] for r in rows if isinstance(r["number_of_intervals"], int)]
    return {
        "row_type": "resource",
        "resource": resource,
        "target": "H3K27ac",
        "assembly": assembly,
        "bed_file_count": len(rows),
        "files_with_interval_counts": len(intervals),
        "total_intervals": sum(intervals),
        "median_intervals_per_file": int(statistics.median(intervals)) if intervals else "",
        "cell_type_count": len({r["cell_type"] for r in rows if r["cell_type"]}),
        "tissue_category_count": len({r["tissue_category"] for r in rows if r["tissue_category"]}),
        "within_resource_pairwise_comparisons": pair_count(len(rows)),
        "cross_resource_reference_mismatched_pairs": "",
    }


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hg19-metadata-url", default=HG19_METADATA_URL)
    parser.add_argument("--hg38-metadata-url", default=HG38_METADATA_URL)
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path(__file__).resolve().parents[1] / "results",
    )
    args = parser.parse_args()

    retrieved = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    hg19_rows = fetch_tsv(args.hg19_metadata_url)
    hg38_rows = fetch_tsv(args.hg38_metadata_url)

    roadmap = select(hg19_rows, "ROADMAP", "hg19")
    blueprint = select(hg38_rows, "Blueprint", "hg38")
    if not roadmap:
        raise RuntimeError("No native Roadmap hg19 H3K27ac BED peak records matched")
    if not blueprint:
        raise RuntimeError("No native BLUEPRINT hg38 H3K27ac BED peak records matched")

    manifest = roadmap + blueprint
    manifest_fields = [
        "resource", "assembly", "identifier", "cell_type", "tissue_category",
        "biosample_type", "assay", "target", "output_type", "file_format",
        "number_of_intervals", "processed_file_download_url", "release_date",
    ]
    manifest_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_manifest.tsv"
    write_tsv(manifest_path, manifest, manifest_fields)

    summary = [
        summarize("Roadmap", "hg19", roadmap),
        summarize("BLUEPRINT", "hg38", blueprint),
    ]
    blocked_pairs = len(roadmap) * len(blueprint)
    summary.append({
        "row_type": "cross_resource",
        "resource": "Roadmap × BLUEPRINT",
        "target": "H3K27ac",
        "assembly": "hg19 × hg38",
        "bed_file_count": len(roadmap) + len(blueprint),
        "files_with_interval_counts": sum(1 for r in manifest if isinstance(r["number_of_intervals"], int)),
        "total_intervals": sum(r["number_of_intervals"] for r in manifest if isinstance(r["number_of_intervals"], int)),
        "median_intervals_per_file": "",
        "cell_type_count": len({r["cell_type"] for r in manifest if r["cell_type"]}),
        "tissue_category_count": len({r["tissue_category"] for r in manifest if r["tissue_category"]}),
        "within_resource_pairwise_comparisons": sum(pair_count(len(rows)) for rows in (roadmap, blueprint)),
        "cross_resource_reference_mismatched_pairs": blocked_pairs,
    })
    summary_fields = list(summary[-1].keys())
    summary_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_summary.tsv"
    write_tsv(summary_path, summary, summary_fields)

    provenance = {
        "retrieved_at_utc": retrieved,
        "sources": {
            "native_hg19_metadata": args.hg19_metadata_url,
            "native_hg38_metadata": args.hg38_metadata_url,
        },
        "selection": {
            "Roadmap": "native hg19, ChIP-seq, exact H3K27ac, peak output, BED-like processed file",
            "BLUEPRINT": "native hg38, ChIP-seq, exact H3K27ac, peak output, BED-like processed file",
            "lifted_records": "excluded",
            "deduplication": "processed download URL",
        },
        "counts": {
            "roadmap_hg19_beds": len(roadmap),
            "blueprint_hg38_beds": len(blueprint),
            "cross_resource_reference_mismatched_pairs": blocked_pairs,
        },
        "outputs": [str(manifest_path), str(summary_path)],
    }
    provenance_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"Roadmap native hg19 H3K27ac BEDs: {len(roadmap):,}")
    print(f"BLUEPRINT native hg38 H3K27ac BEDs: {len(blueprint):,}")
    print(f"Cross-resource pairs blocked by reference mismatch: {blocked_pairs:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
