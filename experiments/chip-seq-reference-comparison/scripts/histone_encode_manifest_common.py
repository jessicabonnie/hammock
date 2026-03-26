"""
Shared helpers for histone ENCODE cart → per-file manifest TSVs.

Batch / provenance columns are taken from the ENCODE experiment cart report row.
"""

import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

ORGANISM_MOUSE = "Mus musculus"
ORGANISM_HUMAN = "Homo sapiens"

# ENCODE ``file.assembly`` values used to pick a reference FASTA (keep in sync with bedpaths script).
KNOWN_ENCODE_ASSEMBLY = frozenset({"mm10", "mm10-minimal", "mm9", "GRCh38", "hg19"})

# Map ``reference_build`` TSV cell / aliases → ENCODE assembly token.
REFERENCE_BUILD_TO_ENCODE_ASSEMBLY = {
    "GRCh38": "GRCh38",
    "hg19": "hg19",
    "GRCh37": "hg19",
    "mm10": "mm10",
    "GRCm38": "mm10",
    "mm9": "mm9",
    "NCBIM37": "mm9",
    "mm10-minimal": "mm10-minimal",
    "mm10_minimal": "mm10-minimal",
}

# Filename slugs for per-species outputs (samplesheets, etc.).
ORGANISM_SLUG = {
    ORGANISM_MOUSE: "Mus_musculus",
    ORGANISM_HUMAN: "Homo_sapiens",
}


def encode_assembly_to_bed_list_slug(encode_assembly: str) -> str:
    """Filesystem-safe token for ``histone_narrowpeak_bed_paths_<slug>.txt``."""
    a = (encode_assembly or "").strip()
    if a == "mm10-minimal":
        return "mm10_minimal"
    return a.replace(" ", "_")


def bed_list_slug_to_encode_assembly(slug: str) -> str:
    """Inverse of :func:`encode_assembly_to_bed_list_slug` for discovery of list files."""
    s = (slug or "").strip()
    if s == "mm10_minimal":
        return "mm10-minimal"
    if s in KNOWN_ENCODE_ASSEMBLY:
        return s
    return ""


def parse_files_cell(cell: str) -> List[str]:
    if not cell:
        return []
    return re.findall(r"/files/(ENCFF\w+)/", cell)


def normalize_organism_label(text: str) -> str:
    """Map cart ``Organism`` cell or free text to ORGANISM_MOUSE / ORGANISM_HUMAN / ''."""
    if not text:
        return ""
    raw = text.strip()
    if raw == ORGANISM_MOUSE or raw == ORGANISM_HUMAN:
        return raw
    s = raw.lower()
    if "mus musculus" in s or s in ("mouse", "mm10", "mm39"):
        return ORGANISM_MOUSE
    if "homo sapiens" in s or s in ("human", "hg19", "hg38", "grch37", "grch38"):
        return ORGANISM_HUMAN
    return ""


def batch_fields_from_experiment_row(meta: Dict[str, str]) -> Dict[str, str]:
    """Map cart report columns to manifest fields (for blocking / PCA)."""
    return {
        "lab": (meta.get("Lab") or "").strip(),
        "project": (meta.get("Project") or "").strip(),
        "library_construction_platform": (meta.get("Library construction platform") or "").strip(),
        "library_construction_method": (meta.get("Library construction method") or "").strip(),
    }


def fetch_encode_file_assembly(accession: str, timeout: int = 60) -> str:
    """Return ENCODE ``assembly`` for a file accession (e.g. ``mm10``, ``GRCh38``, ``hg19``)."""
    url = "https://www.encodeproject.org/files/{}/?format=json".format(accession)
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        data = json.load(resp)
    return (data.get("assembly") or "").strip()


def reference_build_for_manifest(encode_assembly: str, organism: str) -> str:
    """
    Stable ``reference_build`` label for the manifest TSV.

    When ENCODE assembly is known (``--fetch-assembly``), repeats that token. Otherwise falls back
    to mm10 / GRCh38 from cart ``organism`` (organism-only fallback when ENCODE assembly is unknown).
    """
    a = (encode_assembly or "").strip()
    if a:
        return a
    org = normalize_organism_label(organism or "")
    if org == ORGANISM_MOUSE:
        return "mm10"
    if org == ORGANISM_HUMAN:
        return "GRCh38"
    return ""


def effective_encode_assembly_for_row(row: Dict[str, str]) -> str:
    """
    Resolve ENCODE assembly for choosing a reference FASTA: ``assembly`` column, then
    ``reference_build`` (with aliases), then organism-based mm10/GRCh38.
    """
    asm = (row.get("assembly") or "").strip()
    if asm:
        return asm
    rb = (row.get("reference_build") or "").strip()
    if rb in REFERENCE_BUILD_TO_ENCODE_ASSEMBLY:
        return REFERENCE_BUILD_TO_ENCODE_ASSEMBLY[rb]
    if rb in KNOWN_ENCODE_ASSEMBLY:
        return rb
    return reference_build_for_manifest("", (row.get("organism") or "").strip())


def write_narrowpeak_bed_path_lists_by_reference_build(
    rows: List[Dict[str, str]],
    *,
    output_dir: Path,
    basename_prefix: str = "histone_narrowpeak_bed_paths",
) -> Tuple[List[str], int]:
    """
    Write ``{basename_prefix}_<build>.txt`` (one ``bed_path`` per line) grouped by resolved ENCODE
    assembly (``mm10``, ``GRCh38``, ``hg19``, ``mm9``, ``mm10_minimal``, …).

    Returns (list of written file paths as strings, count of rows with no resolvable reference build).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    by_build = defaultdict(list)  # type: Dict[str, List[str]]
    unknown = 0
    for row in rows:
        path = (row.get("bed_path") or "").strip()
        if not path:
            continue
        eff = effective_encode_assembly_for_row(row)
        if not eff:
            unknown += 1
            continue
        by_build[eff].append(path)

    written = []  # type: List[str]
    for build in sorted(by_build.keys()):
        paths = sorted(set(by_build[build]))
        slug = encode_assembly_to_bed_list_slug(build)
        dest = out_dir / "{}_{}.txt".format(basename_prefix, slug)
        dest.write_text("\n".join(paths) + ("\n" if paths else ""), encoding="utf-8")
        written.append(str(dest))

    return written, unknown
