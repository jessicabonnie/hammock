#!/usr/bin/env python3
"""Count native Roadmap hg19 and BLUEPRINT hg38 H3K27ac BED peak files.

Uses FILER's native-build metadata templates. Lifted records are intentionally
excluded so that cross-resource pairs represent a genuine reference mismatch.
The FILER templates have changed column labels over time, so matching is based
on normalized aliases and emits diagnostics rather than silently guessing.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import statistics
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HG19_METADATA_URL = "https://tf.lisanwanglab.org/GADB/metadata/filer.latest.hg19.template"
HG38_METADATA_URL = "https://tf.lisanwanglab.org/GADB/metadata/filer.latest.hg38.template"

RESOURCE_ALIASES = {
    "roadmap": ("roadmap", "roadmap epigenomics", "roadmap epigenome"),
    "blueprint": ("blueprint", "blueprint epigenome", "blueprint epigenomics"),
}


def normalize(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def normalized_key(value: str) -> str:
    return "".join(ch for ch in normalize(value).casefold() if ch.isalnum())


def fetch_tsv(url: str) -> tuple[list[dict[str, str]], list[str]]:
    req = urllib.request.Request(url, headers={"User-Agent": "hammock-database-counts/1.1"})
    with urllib.request.urlopen(req, timeout=180) as response:
        text = response.read().decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not reader.fieldnames:
        raise ValueError(f"No header found in {url}")
    fields = [normalize(field) for field in reader.fieldnames]
    rows = [{normalize(k): normalize(v) for k, v in row.items()} for row in reader]
    return rows, fields


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


def contains_alias(value: str, aliases: tuple[str, ...]) -> bool:
    key = normalized_key(value)
    return any(normalized_key(alias) in key for alias in aliases)


def resource_label(row: dict[str, str]) -> str:
    return lookup(
        row,
        "data_source", "data source", "source", "resource", "project",
        "project_name", "project name", "consortium", "dataset_source",
    )


def target_label(row: dict[str, str]) -> str:
    return lookup(
        row,
        "antibody", "target", "target_of_assay", "target of assay",
        "feature", "epigenetic_mark", "epigenetic mark", "mark",
    )


def download_url(row: dict[str, str]) -> str:
    return lookup(
        row,
        "processed_file_download_url", "processed file download url",
        "file_download_url", "file download url", "download_url", "download url",
        "url",
    )


def select(
    rows: list[dict[str, str]], resource: str, assembly: str
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    selected: dict[str, dict[str, Any]] = {}
    aliases = RESOURCE_ALIASES[resource.casefold()]
    rejection = Counter()
    resource_examples = Counter()
    target_examples = Counter()

    for row in rows:
        source = resource_label(row)
        genome_build = lookup(row, "genome_build", "genome build", "assembly", "reference_genome")
        assay = lookup(row, "assay", "assay_type", "assay type", "data_type", "data type")
        target = target_label(row)
        output_type = lookup(row, "output_type", "output type", "data_category", "data category")
        file_format = lookup(row, "file_format", "file format", "format")
        url = download_url(row)

        if source:
            resource_examples[source] += 1
        if target:
            target_examples[target] += 1

        if not contains_alias(source, aliases):
            rejection["resource"] += 1
            continue

        # The template URL itself fixes the native build. Accept blank build cells,
        # but reject rows explicitly labeled as a different build.
        accepted_builds = {
            normalized_key(assembly),
            normalized_key("GRCh37" if assembly == "hg19" else "GRCh38"),
        }
        if genome_build and normalized_key(genome_build) not in accepted_builds:
            rejection["assembly"] += 1
            continue

        if normalized_key(assay) not in {"chipseq", "chipsequencing"}:
            rejection["assay"] += 1
            continue
        if normalized_key(target) != "h3k27ac":
            rejection["target"] += 1
            continue

        output_key = normalized_key(output_type)
        format_key = normalized_key(file_format)
        url_key = url.casefold()
        peak_like = "peak" in output_key or any(
            token in url_key for token in ("narrowpeak", "broadpeak", "gappedpeak", ".bed")
        )
        if not peak_like:
            rejection["output_type"] += 1
            continue

        bed_like = (
            "bed" in format_key
            or "peak" in format_key
            or url_key.endswith((".bed", ".bed.gz", ".narrowpeak", ".narrowpeak.gz", ".broadpeak", ".broadpeak.gz"))
            or "narrowpeak" in url_key
            or "broadpeak" in url_key
            or "gappedpeak" in url_key
        )
        if not bed_like:
            rejection["file_format"] += 1
            continue
        if not url:
            rejection["missing_url"] += 1
            continue

        canonical = {
            "resource": resource.upper() if resource.casefold() == "blueprint" else "Roadmap",
            "assembly": assembly,
            "identifier": lookup(row, "identifier", "file_id", "file id", "accession"),
            "cell_type": lookup(row, "cell_type", "cell type", "biosample", "biosample_name"),
            "tissue_category": lookup(row, "tissue_category", "tissue category", "tissue"),
            "biosample_type": lookup(row, "biosample_type", "biosample type"),
            "assay": assay,
            "target": "H3K27ac",
            "output_type": output_type,
            "file_format": file_format,
            "number_of_intervals": to_int(lookup(row, "number_of_intervals", "number of intervals")),
            "processed_file_download_url": url,
            "release_date": lookup(row, "release_date", "release date"),
            "source_label_raw": source,
            "genome_build_raw": genome_build,
            "target_label_raw": target,
        }
        selected[url] = canonical

    diagnostics = {
        "resource_requested": resource,
        "assembly_requested": assembly,
        "row_count": len(rows),
        "selected_count": len(selected),
        "rejection_counts": dict(rejection),
        "top_resource_labels": resource_examples.most_common(30),
        "top_h3k27ac_target_labels": [
            [label, count] for label, count in target_examples.most_common(100)
            if "h3k27ac" in normalized_key(label)
        ][:30],
    }
    result = sorted(selected.values(), key=lambda r: (r["tissue_category"], r["cell_type"], r["identifier"]))
    return result, diagnostics


def pair_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


def summarize(resource: str, assembly: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    intervals = [r["number_of_intervals"] for r in rows if isinstance(r["number_of_intervals"], int)]
    return {
        "row_type": "resource", "resource": resource, "target": "H3K27ac",
        "assembly": assembly, "bed_file_count": len(rows),
        "files_with_interval_counts": len(intervals), "total_intervals": sum(intervals),
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
    hg19_rows, hg19_fields = fetch_tsv(args.hg19_metadata_url)
    hg38_rows, hg38_fields = fetch_tsv(args.hg38_metadata_url)

    roadmap, roadmap_diag = select(hg19_rows, "roadmap", "hg19")
    blueprint, blueprint_diag = select(hg38_rows, "blueprint", "hg38")

    diagnostic_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_diagnostics.json"
    diagnostic_path.parent.mkdir(parents=True, exist_ok=True)
    diagnostic_path.write_text(json.dumps({
        "retrieved_at_utc": retrieved,
        "hg19_columns": hg19_fields,
        "hg38_columns": hg38_fields,
        "roadmap": roadmap_diag,
        "blueprint": blueprint_diag,
    }, indent=2) + "\n", encoding="utf-8")

    if not roadmap or not blueprint:
        missing = []
        if not roadmap:
            missing.append("Roadmap hg19")
        if not blueprint:
            missing.append("BLUEPRINT hg38")
        raise RuntimeError(
            "No matching records for " + " and ".join(missing) + ". "
            f"Inspect diagnostics: {diagnostic_path}"
        )

    manifest = roadmap + blueprint
    manifest_fields = [
        "resource", "assembly", "identifier", "cell_type", "tissue_category",
        "biosample_type", "assay", "target", "output_type", "file_format",
        "number_of_intervals", "processed_file_download_url", "release_date",
        "source_label_raw", "genome_build_raw", "target_label_raw",
    ]
    manifest_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_manifest.tsv"
    write_tsv(manifest_path, manifest, manifest_fields)

    summary = [summarize("Roadmap", "hg19", roadmap), summarize("BLUEPRINT", "hg38", blueprint)]
    blocked_pairs = len(roadmap) * len(blueprint)
    summary.append({
        "row_type": "cross_resource", "resource": "Roadmap × BLUEPRINT",
        "target": "H3K27ac", "assembly": "hg19 × hg38",
        "bed_file_count": len(manifest),
        "files_with_interval_counts": sum(isinstance(r["number_of_intervals"], int) for r in manifest),
        "total_intervals": sum(r["number_of_intervals"] for r in manifest if isinstance(r["number_of_intervals"], int)),
        "median_intervals_per_file": "",
        "cell_type_count": len({r["cell_type"] for r in manifest if r["cell_type"]}),
        "tissue_category_count": len({r["tissue_category"] for r in manifest if r["tissue_category"]}),
        "within_resource_pairwise_comparisons": pair_count(len(roadmap)) + pair_count(len(blueprint)),
        "cross_resource_reference_mismatched_pairs": blocked_pairs,
    })
    summary_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_summary.tsv"
    write_tsv(summary_path, summary, list(summary[-1].keys()))

    provenance = {
        "retrieved_at_utc": retrieved,
        "sources": {"native_hg19_metadata": args.hg19_metadata_url, "native_hg38_metadata": args.hg38_metadata_url},
        "selection": {
            "Roadmap": "native hg19, Roadmap resource alias, ChIP-seq, exact H3K27ac, peak BED-like output",
            "BLUEPRINT": "native hg38, BLUEPRINT resource alias, ChIP-seq, exact H3K27ac, peak BED-like output",
            "lifted_records": "excluded", "deduplication": "processed download URL",
        },
        "counts": {
            "roadmap_hg19_beds": len(roadmap), "blueprint_hg38_beds": len(blueprint),
            "cross_resource_reference_mismatched_pairs": blocked_pairs,
        },
        "outputs": [str(manifest_path), str(summary_path), str(diagnostic_path)],
    }
    provenance_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"Roadmap native hg19 H3K27ac BEDs: {len(roadmap):,}")
    print(f"BLUEPRINT native hg38 H3K27ac BEDs: {len(blueprint):,}")
    print(f"Cross-resource pairs blocked by reference mismatch: {blocked_pairs:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
