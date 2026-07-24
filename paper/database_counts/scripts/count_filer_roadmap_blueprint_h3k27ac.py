#!/usr/bin/env python3
"""Count native Roadmap hg19 and BLUEPRINT hg38 H3K27ac BED peak files.

Roadmap records are read from FILER's hg19 metadata template. The downloadable
FILER template currently lags the live FILER2 catalog and does not contain
BLUEPRINT, so BLUEPRINT records are collected from the live, filtered FILER2
browser. Lifted records are excluded.
"""

from __future__ import annotations

import argparse
import csv
import html
import io
import json
import math
import re
import statistics
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HG19_METADATA_URL = "https://tf.lisanwanglab.org/GADB/metadata/filer.latest.hg19.template"
FILER2_BROWSE_URL = "https://filer2.niagads.org/browse_faceted_with_search.php"
USER_AGENT = "hammock-database-counts/1.2"


def normalize(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def normalized_key(value: str) -> str:
    return "".join(ch for ch in normalize(value).casefold() if ch.isalnum())


def fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read().decode("utf-8-sig", errors="replace")


def fetch_tsv(url: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(fetch_text(url)), delimiter="\t")
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


def select_roadmap(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    selected: dict[str, dict[str, Any]] = {}
    for row in rows:
        source = lookup(row, "data_source", "data source")
        build = lookup(row, "genome_build", "genome build")
        assay = lookup(row, "assay")
        target = lookup(row, "antibody", "target")
        output_type = lookup(row, "output_type", "output type")
        file_format = lookup(row, "file_format", "file format")
        url = lookup(row, "processed_file_download_url", "processed file download url")

        if source.casefold() != "roadmap":
            continue
        if build.casefold() not in {"hg19", "grch37"}:
            continue
        if normalized_key(assay) != "chipseq":
            continue
        if normalized_key(target) != "h3k27ac":
            continue
        if "peak" not in output_type.casefold():
            continue
        if not url:
            continue
        if "bed" not in file_format.casefold() and not re.search(
            r"\.(?:bed|narrowpeak|broadpeak|gappedpeak)(?:\.gz)?$", url, re.I
        ):
            continue

        selected[url] = {
            "resource": "Roadmap",
            "assembly": "hg19",
            "identifier": lookup(row, "identifier"),
            "cell_type": lookup(row, "cell_type", "cell type"),
            "tissue_category": lookup(row, "tissue_category", "tissue category"),
            "biosample_type": lookup(row, "biosample_type", "biosample type"),
            "assay": assay,
            "target": "H3K27ac",
            "output_type": output_type,
            "file_format": file_format,
            "number_of_intervals": to_int(lookup(row, "number_of_intervals", "number of intervals")),
            "processed_file_download_url": url,
            "release_date": lookup(row, "release_date", "release date"),
            "metadata_source": "FILER hg19 template",
        }
    return sorted(selected.values(), key=lambda row: (row["cell_type"], row["identifier"]))


def strip_html(raw: str) -> str:
    raw = re.sub(r"<script\b[^>]*>.*?</script>", " ", raw, flags=re.I | re.S)
    raw = re.sub(r"<style\b[^>]*>.*?</style>", " ", raw, flags=re.I | re.S)
    raw = re.sub(r"<[^>]+>", "\n", raw)
    return html.unescape(raw)


def parse_blueprint_page(raw_html: str) -> tuple[list[dict[str, Any]], int | None]:
    text = strip_html(raw_html)
    total_match = re.search(r"out of\s+([\d,]+)", text, flags=re.I)
    total = int(total_match.group(1).replace(",", "")) if total_match else None

    url_pattern = re.compile(
        r"https?://[^\s<>\"']+/IHEC(?:/IHEC_Blueprint|/Blueprint)/ChIP-seq/[^\s<>\"']*?"
        r"hg38/[^\s<>\"']*?\.H3K27ac\.[^\s<>\"']*?\.bed\.gz",
        flags=re.I,
    )
    rows: dict[str, dict[str, Any]] = {}
    for match in url_pattern.finditer(raw_html):
        url = html.unescape(match.group(0))
        context = strip_html(raw_html[max(0, match.start() - 6000): match.end() + 1000])

        identifier_matches = re.findall(r"\bNGBLP[A-Z0-9]+\b", context)
        identifier = identifier_matches[-1] if identifier_matches else ""
        cell_matches = re.findall(
            r"Blueprint\s+(.+?)\s+ChIP-seq\s+H3K27ac(?:-histone-mark)?\s+peaks",
            context,
            flags=re.I,
        )
        cell_type = normalize(cell_matches[-1]) if cell_matches else ""
        interval_matches = re.findall(r"number[_ ]of[_ ]intervals\s+([\d,]+)", context, flags=re.I)
        interval_count = to_int(interval_matches[-1]) if interval_matches else None
        tissue_matches = re.findall(r"tissue[_ ]category\s+([^\n]+)", context, flags=re.I)
        tissue_category = normalize(tissue_matches[-1]) if tissue_matches else ""
        biosample_matches = re.findall(r"biosample[_ ]type\s+([^\n]+)", context, flags=re.I)
        biosample_type = normalize(biosample_matches[-1]) if biosample_matches else ""

        rows[url] = {
            "resource": "BLUEPRINT",
            "assembly": "hg38",
            "identifier": identifier,
            "cell_type": cell_type,
            "tissue_category": tissue_category,
            "biosample_type": biosample_type,
            "assay": "ChIP-seq",
            "target": "H3K27ac",
            "output_type": "peaks",
            "file_format": "bed6",
            "number_of_intervals": interval_count,
            "processed_file_download_url": url,
            "release_date": "",
            "metadata_source": "live FILER2 browser",
        }
    return list(rows.values()), total


def fetch_blueprint_h3k27ac(page_size: int = 500) -> list[dict[str, Any]]:
    selected: dict[str, dict[str, Any]] = {}
    start = 0
    total: int | None = None
    while total is None or start < total:
        query = urllib.parse.urlencode(
            {
                "count": page_size,
                "dataSource[]": "Blueprint",
                "genomeBuild[]": "hg38",
                "start": start,
            },
            doseq=True,
        )
        page_rows, page_total = parse_blueprint_page(fetch_text(f"{FILER2_BROWSE_URL}?{query}"))
        for row in page_rows:
            selected[row["processed_file_download_url"]] = row
        if total is None:
            total = page_total
        if page_total is None:
            break
        start += page_size
    return sorted(selected.values(), key=lambda row: (row["cell_type"], row["identifier"]))


def pair_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


def summarize(resource: str, assembly: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    intervals = [row["number_of_intervals"] for row in rows if isinstance(row["number_of_intervals"], int)]
    return {
        "row_type": "resource",
        "resource": resource,
        "target": "H3K27ac",
        "assembly": assembly,
        "bed_file_count": len(rows),
        "files_with_interval_counts": len(intervals),
        "total_intervals": sum(intervals),
        "median_intervals_per_file": int(statistics.median(intervals)) if intervals else "",
        "cell_type_count": len({row["cell_type"] for row in rows if row["cell_type"]}),
        "tissue_category_count": len({row["tissue_category"] for row in rows if row["tissue_category"]}),
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
    parser.add_argument("--blueprint-page-size", type=int, default=500)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results",
    )
    args = parser.parse_args()

    retrieved = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    roadmap = select_roadmap(fetch_tsv(args.hg19_metadata_url))
    blueprint = fetch_blueprint_h3k27ac(args.blueprint_page_size)

    if not roadmap:
        raise RuntimeError("No native Roadmap hg19 H3K27ac BED peak records matched")
    if not blueprint:
        raise RuntimeError(
            "No native BLUEPRINT hg38 H3K27ac BED records were parsed from the live "
            "FILER2 browser. Save the returned HTML and inspect its URL/field structure."
        )

    manifest = roadmap + blueprint
    manifest_fields = [
        "resource", "assembly", "identifier", "cell_type", "tissue_category",
        "biosample_type", "assay", "target", "output_type", "file_format",
        "number_of_intervals", "processed_file_download_url", "release_date",
        "metadata_source",
    ]
    manifest_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_manifest.tsv"
    write_tsv(manifest_path, manifest, manifest_fields)

    summary = [summarize("Roadmap", "hg19", roadmap), summarize("BLUEPRINT", "hg38", blueprint)]
    blocked_pairs = len(roadmap) * len(blueprint)
    summary.append({
        "row_type": "cross_resource",
        "resource": "Roadmap × BLUEPRINT",
        "target": "H3K27ac",
        "assembly": "hg19 × hg38",
        "bed_file_count": len(manifest),
        "files_with_interval_counts": sum(isinstance(row["number_of_intervals"], int) for row in manifest),
        "total_intervals": sum(
            row["number_of_intervals"] for row in manifest
            if isinstance(row["number_of_intervals"], int)
        ),
        "median_intervals_per_file": "",
        "cell_type_count": len({row["cell_type"] for row in manifest if row["cell_type"]}),
        "tissue_category_count": len({row["tissue_category"] for row in manifest if row["tissue_category"]}),
        "within_resource_pairwise_comparisons": pair_count(len(roadmap)) + pair_count(len(blueprint)),
        "cross_resource_reference_mismatched_pairs": blocked_pairs,
    })
    summary_path = args.output_dir / "filer_roadmap_blueprint_h3k27ac_summary.tsv"
    write_tsv(summary_path, summary, list(summary[-1]))

    provenance = {
        "retrieved_at_utc": retrieved,
        "sources": {
            "roadmap_hg19": args.hg19_metadata_url,
            "blueprint_hg38": FILER2_BROWSE_URL,
        },
        "important_source_note": (
            "The downloadable FILER template lacked BLUEPRINT although the live FILER2 "
            "browser contained BLUEPRINT tracks; BLUEPRINT was therefore parsed from "
            "the live filtered catalog."
        ),
        "selection": {
            "Roadmap": "native hg19, ChIP-seq, exact H3K27ac, peak BED",
            "BLUEPRINT": "native hg38 URL path, Blueprint collection, H3K27ac peak BED",
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
