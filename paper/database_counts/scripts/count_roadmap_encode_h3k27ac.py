#!/usr/bin/env python3
"""Count Roadmap hg19 and ENCODE GRCh38 H3K27ac peak BED files.

The output separates directly comparable within-resource pairs from cross-resource
pairs that cannot be compared with ordinary coordinate overlap because the two
collections use different human reference assemblies.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Iterable

ROADMAP_NARROWPEAK_URL = (
    "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/"
)
ENCODE_SEARCH_URL = "https://www.encodeproject.org/search/"
USER_AGENT = "hammock-database-counts/1.0"


def pair_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        return json.load(response)


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        href = dict(attrs).get("href")
        if href:
            self.links.append(href)


def count_roadmap(url: str, target: str) -> list[dict[str, Any]]:
    parser = LinkParser()
    parser.feed(fetch_text(url))
    pattern = re.compile(rf"^(E\d+)-{re.escape(target)}\.narrowPeak\.gz$", re.IGNORECASE)
    rows = []
    for href in sorted(set(parser.links)):
        name = urllib.parse.unquote(Path(href).name)
        match = pattern.match(name)
        if not match:
            continue
        rows.append(
            {
                "repository": "Roadmap Epigenomics",
                "dataset_id": match.group(1),
                "file_accession": "",
                "target": target,
                "assembly": "hg19",
                "file_format": "narrowPeak",
                "output_type": "consolidated narrow peaks",
                "download_url": urllib.parse.urljoin(url, href),
            }
        )
    return rows


def encode_search_url(target: str, assembly: str) -> str:
    params = [
        ("type", "File"),
        ("status", "released"),
        ("file_format", "bed"),
        ("assembly", assembly),
        ("assay_title", "Histone ChIP-seq"),
        ("target.label", target),
        ("limit", "all"),
        ("format", "json"),
        ("frame", "embedded"),
    ]
    return ENCODE_SEARCH_URL + "?" + urllib.parse.urlencode(params)


def embedded_value(record: dict[str, Any], path: str, default: str = "") -> str:
    value: Any = record
    for part in path.split("."):
        if not isinstance(value, dict):
            return default
        value = value.get(part)
    if isinstance(value, str):
        return value
    return default


def count_encode(target: str, assembly: str) -> list[dict[str, Any]]:
    url = encode_search_url(target, assembly)
    payload = fetch_json(url)
    candidates = []
    for record in payload.get("@graph", []):
        output_type = str(record.get("output_type", ""))
        file_type = str(record.get("file_type", ""))
        if "peak" not in output_type.casefold():
            continue
        if "bed" not in file_type.casefold() and record.get("file_format") != "bed":
            continue
        dataset_id = embedded_value(record, "dataset.accession")
        accession = str(record.get("accession", ""))
        href = str(record.get("href", ""))
        candidates.append(
            {
                "repository": "ENCODE",
                "dataset_id": dataset_id,
                "file_accession": accession,
                "target": target,
                "assembly": assembly,
                "file_format": file_type or "bed",
                "output_type": output_type,
                "download_url": urllib.parse.urljoin("https://www.encodeproject.org", href),
                "preferred_default": bool(record.get("preferred_default", False)),
                "date_created": str(record.get("date_created", "")),
            }
        )

    # Select one released peak BED per experiment. Prefer ENCODE's preferred_default,
    # then IDR/replicated outputs, then newest accession as a deterministic fallback.
    rank_terms = (
        "optimal IDR thresholded peaks",
        "conservative IDR thresholded peaks",
        "IDR thresholded peaks",
        "replicated peaks",
        "pseudoreplicated peaks",
        "peaks",
    )
    rank = {name.casefold(): index for index, name in enumerate(rank_terms)}
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in candidates:
        grouped.setdefault(row["dataset_id"] or row["file_accession"], []).append(row)

    selected = []
    for rows in grouped.values():
        rows.sort(
            key=lambda row: (
                not row["preferred_default"],
                rank.get(row["output_type"].casefold(), len(rank)),
                row["file_accession"],
            )
        )
        selected.append(rows[0])
    return sorted(selected, key=lambda row: (row["dataset_id"], row["file_accession"]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="H3K27ac")
    parser.add_argument("--roadmap-url", default=ROADMAP_NARROWPEAK_URL)
    parser.add_argument("--encode-assembly", default="GRCh38")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results",
    )
    args = parser.parse_args()

    retrieved_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    roadmap = count_roadmap(args.roadmap_url, args.target)
    encode = count_encode(args.target, args.encode_assembly)
    if not roadmap:
        raise RuntimeError("No Roadmap narrowPeak files matched the requested target")
    if not encode:
        raise RuntimeError("No ENCODE peak BED files matched the requested target")

    stem = f"roadmap_encode_{args.target.lower()}"
    manifest_path = args.output_dir / f"{stem}_manifest.tsv"
    summary_path = args.output_dir / f"{stem}_summary.tsv"
    provenance_path = args.output_dir / f"{stem}_provenance.json"

    manifest = roadmap + encode
    for row in manifest:
        row["retrieved_at_utc"] = retrieved_at
    fields = [
        "repository", "dataset_id", "file_accession", "target", "assembly",
        "file_format", "output_type", "download_url", "retrieved_at_utc",
    ]
    write_tsv(manifest_path, manifest, fields)

    n_roadmap = len(roadmap)
    n_encode = len(encode)
    summary = [
        {
            "row_type": "repository",
            "repository": "Roadmap Epigenomics",
            "target": args.target,
            "assembly": "hg19",
            "bed_file_count": n_roadmap,
            "within_reference_pairwise_comparisons": pair_count(n_roadmap),
            "cross_repository_incompatible_pairs": "",
            "retrieved_at_utc": retrieved_at,
        },
        {
            "row_type": "repository",
            "repository": "ENCODE",
            "target": args.target,
            "assembly": args.encode_assembly,
            "bed_file_count": n_encode,
            "within_reference_pairwise_comparisons": pair_count(n_encode),
            "cross_repository_incompatible_pairs": "",
            "retrieved_at_utc": retrieved_at,
        },
        {
            "row_type": "cross_repository",
            "repository": "Roadmap Epigenomics × ENCODE",
            "target": args.target,
            "assembly": f"hg19 × {args.encode_assembly}",
            "bed_file_count": n_roadmap + n_encode,
            "within_reference_pairwise_comparisons": pair_count(n_roadmap) + pair_count(n_encode),
            "cross_repository_incompatible_pairs": n_roadmap * n_encode,
            "retrieved_at_utc": retrieved_at,
        },
    ]
    write_tsv(summary_path, summary, list(summary[-1].keys()))

    provenance = {
        "retrieved_at_utc": retrieved_at,
        "target": args.target,
        "roadmap": {
            "source_url": args.roadmap_url,
            "assembly": "hg19",
            "selection": "consolidated narrowPeak filenames matching E###-TARGET.narrowPeak.gz",
            "file_count": n_roadmap,
        },
        "encode": {
            "search_url": encode_search_url(args.target, args.encode_assembly),
            "assembly": args.encode_assembly,
            "selection": "one released peak BED per experiment, preferring preferred_default and IDR/replicated outputs",
            "file_count": n_encode,
        },
        "cross_repository_pair_policy": (
            "Roadmap_count multiplied by ENCODE_count; these pairs share target and species "
            "but cannot be compared directly with coordinate-overlap methods because their assemblies differ"
        ),
        "outputs": [str(manifest_path), str(summary_path)],
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"Roadmap hg19 {args.target} narrowPeak files: {n_roadmap:,}")
    print(f"ENCODE {args.encode_assembly} {args.target} peak BED files: {n_encode:,}")
    print(f"Cross-resource incompatible pairs: {n_roadmap * n_encode:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
