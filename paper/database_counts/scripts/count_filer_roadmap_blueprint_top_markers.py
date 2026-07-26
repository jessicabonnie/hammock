#!/usr/bin/env python3
"""Count shared histone-mark BED files in Roadmap hg19 and BLUEPRINT hg38.

Roadmap records are read from FILER's native hg19 metadata template. BLUEPRINT
records are parsed from the live FILER2 browser because the downloadable hg38
template currently omits that collection.

For every histone mark represented in both resources, the script reports:
  * native Roadmap hg19 BED-file count
  * native BLUEPRINT hg38 BED-file count
  * all possible pairs in the combined collection
  * directly comparable within-reference pairs
  * cross-resource pairs blocked by the hg19/hg38 mismatch

Markers are ranked by blocked cross-reference comparisons (R * B), which directly
measures the comparison opportunity lost to incompatible coordinate systems.
"""

from __future__ import annotations

import argparse
import csv
import html
import io
import json
import math
import re
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

HG19_METADATA_URL = "https://tf.lisanwanglab.org/GADB/metadata/filer.latest.hg19.template"
FILER2_BROWSE_URL = "https://filer2.niagads.org/browse_faceted_with_search.php"
USER_AGENT = "hammock-database-counts/1.3"


def normalize(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def normalized_key(value: str) -> str:
    return "".join(ch for ch in normalize(value).casefold() if ch.isalnum())


def canonical_mark(value: str) -> str:
    """Return a consistent histone-mark label, or an empty string if not one."""
    compact = re.sub(r"[^A-Za-z0-9]", "", normalize(value))
    match = re.fullmatch(r"(?i)(H[1-4](?:K\d+|A|B|C|E|F|G|H|J|K|M|P|S|T|V|W|Y)[A-Za-z0-9]*)", compact)
    if not match:
        return ""
    mark = match.group(1)
    # Histone marks in these catalogs conventionally begin with uppercase H and K.
    mark = re.sub(r"^h", "H", mark, flags=re.I)
    mark = re.sub(r"k", "K", mark, flags=re.I)
    mark = re.sub(r"me", "me", mark, flags=re.I)
    mark = re.sub(r"ac", "ac", mark, flags=re.I)
    return mark


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


def is_peak_bed(output_type: str, file_format: str, url: str) -> bool:
    output_key = output_type.casefold()
    format_key = file_format.casefold()
    url_key = url.casefold()
    peak_like = "peak" in output_key or any(
        token in url_key for token in ("narrowpeak", "broadpeak", "gappedpeak", ".bed")
    )
    bed_like = (
        "bed" in format_key
        or "peak" in format_key
        or re.search(r"\.(?:bed|narrowpeak|broadpeak|gappedpeak)(?:\.gz)?(?:\?|$)", url_key)
        is not None
    )
    return peak_like and bed_like and bool(url)


def select_roadmap(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected: dict[str, dict[str, str]] = {}
    for row in rows:
        source = lookup(row, "data_source", "data source")
        build = lookup(row, "genome_build", "genome build")
        assay = lookup(row, "assay")
        mark = canonical_mark(lookup(row, "antibody", "target"))
        output_type = lookup(row, "output_type", "output type")
        file_format = lookup(row, "file_format", "file format")
        url = lookup(row, "processed_file_download_url", "processed file download url")

        if source.casefold() != "roadmap":
            continue
        if build.casefold() not in {"hg19", "grch37"}:
            continue
        if normalized_key(assay) != "chipseq":
            continue
        if not mark:
            continue
        if not is_peak_bed(output_type, file_format, url):
            continue

        selected[url] = {
            "resource": "Roadmap",
            "assembly": "hg19",
            "marker": mark,
            "identifier": lookup(row, "identifier"),
            "cell_type": lookup(row, "cell_type", "cell type"),
            "processed_file_download_url": url,
        }
    return list(selected.values())


def strip_html(raw: str) -> str:
    raw = re.sub(r"<script\b[^>]*>.*?</script>", " ", raw, flags=re.I | re.S)
    raw = re.sub(r"<style\b[^>]*>.*?</style>", " ", raw, flags=re.I | re.S)
    raw = re.sub(r"<[^>]+>", "\n", raw)
    return html.unescape(raw)


def mark_from_blueprint_context(context: str, url: str) -> str:
    patterns = (
        r"ChIP-seq\s+([A-Za-z0-9_-]+?)(?:-histone-mark)?\s+peaks",
        r"antibody\s+([A-Za-z0-9_-]+)",
        r"target\s+([A-Za-z0-9_-]+)",
    )
    for pattern in patterns:
        matches = re.findall(pattern, context, flags=re.I)
        for value in reversed(matches):
            mark = canonical_mark(value)
            if mark:
                return mark

    # FILER2 BLUEPRINT filenames commonly contain the mark as a dot-delimited token.
    for token in reversed(re.split(r"[./_\-]", urllib.parse.urlparse(url).path)):
        mark = canonical_mark(token)
        if mark:
            return mark
    return ""


def parse_blueprint_page(raw_html: str) -> tuple[list[dict[str, str]], int | None]:
    text = strip_html(raw_html)
    total_match = re.search(r"out of\s+([\d,]+)", text, flags=re.I)
    total = int(total_match.group(1).replace(",", "")) if total_match else None

    url_pattern = re.compile(
        r"https?://[^\s<>\"']+/IHEC(?:/IHEC_Blueprint|/Blueprint)/ChIP-seq/[^\s<>\"']*?"
        r"hg38/[^\s<>\"']*?\.bed\.gz",
        flags=re.I,
    )
    rows: dict[str, dict[str, str]] = {}
    for match in url_pattern.finditer(raw_html):
        url = html.unescape(match.group(0))
        context = strip_html(raw_html[max(0, match.start() - 6000): match.end() + 1200])
        mark = mark_from_blueprint_context(context, url)
        if not mark:
            continue
        identifiers = re.findall(r"\bNGBLP[A-Z0-9]+\b", context)
        cell_matches = re.findall(
            rf"Blueprint\s+(.+?)\s+ChIP-seq\s+{re.escape(mark)}(?:-histone-mark)?\s+peaks",
            context,
            flags=re.I,
        )
        rows[url] = {
            "resource": "BLUEPRINT",
            "assembly": "hg38",
            "marker": mark,
            "identifier": identifiers[-1] if identifiers else "",
            "cell_type": normalize(cell_matches[-1]) if cell_matches else "",
            "processed_file_download_url": url,
        }
    return list(rows.values()), total


def fetch_blueprint(page_size: int) -> list[dict[str, str]]:
    selected: dict[str, dict[str, str]] = {}
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
        rows, page_total = parse_blueprint_page(fetch_text(f"{FILER2_BROWSE_URL}?{query}"))
        for row in rows:
            selected[row["processed_file_download_url"]] = row
        if total is None:
            total = page_total
        if page_total is None:
            break
        start += page_size
    return list(selected.values())


def pair_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


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
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results",
    )
    args = parser.parse_args()

    retrieved = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    roadmap = select_roadmap(fetch_tsv(args.hg19_metadata_url))
    blueprint = fetch_blueprint(args.blueprint_page_size)
    if not roadmap:
        raise RuntimeError("No native Roadmap hg19 histone-mark BED peaks matched")
    if not blueprint:
        raise RuntimeError("No native BLUEPRINT hg38 histone-mark BED peaks matched")

    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in roadmap + blueprint:
        counts[row["marker"]][row["resource"]] += 1

    summary: list[dict[str, Any]] = []
    for marker, resource_counts in counts.items():
        r = resource_counts["Roadmap"]
        b = resource_counts["BLUEPRINT"]
        if r == 0 or b == 0:
            continue
        total_files = r + b
        blocked = r * b
        directly_comparable = pair_count(r) + pair_count(b)
        all_pairs = pair_count(total_files)
        summary.append(
            {
                "marker": marker,
                "roadmap_hg19_bed_count": r,
                "blueprint_hg38_bed_count": b,
                "combined_bed_count": total_files,
                "roadmap_within_reference_pairs": pair_count(r),
                "blueprint_within_reference_pairs": pair_count(b),
                "directly_comparable_pairs": directly_comparable,
                "blocked_cross_reference_pairs": blocked,
                "all_possible_pairs": all_pairs,
                "blocked_percent": 100 * blocked / all_pairs if all_pairs else 0,
                "retrieved_at_utc": retrieved,
            }
        )

    summary.sort(
        key=lambda row: (
            row["blocked_cross_reference_pairs"],
            row["combined_bed_count"],
            row["marker"],
        ),
        reverse=True,
    )
    if len(summary) < args.top_n:
        raise RuntimeError(
            f"Only {len(summary)} histone marks were found in both resources; "
            f"cannot select top {args.top_n}."
        )

    fields = list(summary[0])
    all_path = args.output_dir / "roadmap_blueprint_all_shared_markers.tsv"
    top_path = args.output_dir / f"roadmap_blueprint_top{args.top_n}_markers.tsv"
    manifest_path = args.output_dir / "roadmap_blueprint_shared_marker_manifest.tsv"
    provenance_path = args.output_dir / "roadmap_blueprint_top_markers_provenance.json"

    write_tsv(all_path, summary, fields)
    write_tsv(top_path, summary[: args.top_n], fields)
    write_tsv(
        manifest_path,
        sorted(roadmap + blueprint, key=lambda row: (row["marker"], row["resource"], row["identifier"])),
        ["resource", "assembly", "marker", "identifier", "cell_type", "processed_file_download_url"],
    )

    provenance = {
        "retrieved_at_utc": retrieved,
        "sources": {
            "roadmap_hg19": args.hg19_metadata_url,
            "blueprint_hg38": FILER2_BROWSE_URL,
        },
        "ranking": "blocked_cross_reference_pairs = Roadmap count * BLUEPRINT count",
        "selection": {
            "Roadmap": "native hg19 Roadmap ChIP-seq histone-mark peak BEDs",
            "BLUEPRINT": "native hg38 BLUEPRINT ChIP-seq histone-mark BEDs from live FILER2",
            "markers": "canonical histone-mark labels present in both resources",
            "deduplication": "processed BED URL",
        },
        "counts": {
            "roadmap_histone_mark_beds": len(roadmap),
            "blueprint_histone_mark_beds": len(blueprint),
            "shared_markers": len(summary),
            "top_n": args.top_n,
        },
        "outputs": [str(all_path), str(top_path), str(manifest_path)],
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"Roadmap native hg19 histone-mark BEDs: {len(roadmap):,}")
    print(f"BLUEPRINT native hg38 histone-mark BEDs: {len(blueprint):,}")
    print(f"Markers represented in both resources: {len(summary):,}")
    print(f"Top {args.top_n} markers by blocked cross-reference pairs:")
    for rank, row in enumerate(summary[: args.top_n], start=1):
        print(
            f"  {rank}. {row['marker']}: "
            f"{row['roadmap_hg19_bed_count']:,} hg19 + "
            f"{row['blueprint_hg38_bed_count']:,} hg38; "
            f"{row['blocked_cross_reference_pairs']:,} blocked pairs"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
