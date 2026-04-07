#!/usr/bin/env python3
"""
For each ``.bed.gz`` path listed in the species BED list files, run ``bedtools getfasta`` and
write ``ENCFFxxxxxx.fa`` under ``fasta_from_bed/`` (same basename stem as the BED).

Default reference FASTAs live under ``/data/blangme2/jessica/mus_homo/references/`` (mm10, GRCh38,
hg19, mm9, mm10-minimal). Override with env vars or ``--ref-*``. Carts often mix assemblies; use a
manifest with ``--fetch-assembly`` (see ``build_histone_narrowpeak_bed_manifest.py``) so the
``assembly`` / ``reference_build`` columns match each BED.

Parallelism: use ``-j N`` (default ``min(8, CPU count)``) to run multiple ``bedtools getfasta``
processes; use ``-j 1`` for strictly sequential execution.
"""

import argparse
import csv
import multiprocessing
import os
import shlex
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

RowDict = Dict[str, str]

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from histone_encode_manifest_common import (  # noqa: E402
    KNOWN_ENCODE_ASSEMBLY,
    ORGANISM_SLUG,
    effective_encode_assembly_for_row,
    encode_assembly_to_bed_list_slug,
    normalize_organism_label,
)

PATH_LISTS_DIRNAME = "path_lists"
NARROWPEAK_LIST_PREFIX = "histone_narrowpeak_bed_paths_"
BROADPEAK_LIST_PREFIX = "histone_broadpeak_bed_paths_"
NARROWPEAK_FASTA_SUBDIR = "narrowpeak"
BROADPEAK_FASTA_SUBDIR = "broadpeak"

# Lab reference tree (see experiments/chip-seq-reference-comparison/README.md §5).
_REF_BASE = "/data/blangme2/jessica/mus_homo/references"
DEFAULT_MM10_FASTA = _REF_BASE + "/mm10/mm10_no_alt_analysis_set_ENCODE.fasta"
DEFAULT_GRCH38_FASTA = (
    _REF_BASE + "/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
)
DEFAULT_HG19_FASTA = _REF_BASE + "/hg19/hg19.fa"
DEFAULT_MM9_FASTA = _REF_BASE + "/mm9/mm9.fa"
DEFAULT_MM10_MINIMAL_FASTA = _REF_BASE + "/mm10_minimal/mm10_minimal_ENCODE.fasta"


def _ref_from_arg_or_env(arg_path: Optional[Path], env_key: str, default: Optional[str] = None) -> Path:
    if arg_path is not None:
        return Path(arg_path)
    ev = os.environ.get(env_key, "").strip()
    if ev:
        return Path(ev)
    if default:
        return Path(default)
    return Path("")


def build_ref_map(args: argparse.Namespace) -> Dict[str, Path]:
    return {
        "mm10": _ref_from_arg_or_env(args.ref_mm10, "MM10_FASTA", DEFAULT_MM10_FASTA),
        "mm10-minimal": _ref_from_arg_or_env(
            args.ref_mm10_minimal, "MM10_MINIMAL_FASTA", DEFAULT_MM10_MINIMAL_FASTA
        ),
        "mm9": _ref_from_arg_or_env(args.ref_mm9, "MM9_FASTA", DEFAULT_MM9_FASTA),
        "GRCh38": _ref_from_arg_or_env(args.ref_grch38, "GRCH38_FASTA", DEFAULT_GRCH38_FASTA),
        "hg19": _ref_from_arg_or_env(args.ref_hg19, "HG19_FASTA", DEFAULT_HG19_FASTA),
    }


def resolve_ref_for_assembly(assembly: str, ref_map: Dict[str, Path]) -> Path:
    asm = (assembly or "").strip()
    if asm not in KNOWN_ENCODE_ASSEMBLY:
        raise ValueError(
            "unknown ENCODE assembly {!r} (known: {})".format(asm, ", ".join(sorted(KNOWN_ENCODE_ASSEMBLY)))
        )
    p = ref_map.get(asm, Path(""))
    if not str(p) or not p.is_file():
        raise ValueError(
            "no valid FASTA for assembly {!r}; set {} or the matching env var".format(
                asm,
                {
                    "mm10": "--ref-mm10 / $MM10_FASTA",
                    "mm10-minimal": "--ref-mm10-minimal / $MM10_MINIMAL_FASTA",
                    "mm9": "--ref-mm9 / $MM9_FASTA",
                    "GRCh38": "--ref-grch38 / $GRCH38_FASTA",
                    "hg19": "--ref-hg19 / $HG19_FASTA",
                }[asm],
            )
        )
    return p


def validate_refs_for_assemblies(assemblies: Set[str], ref_map: Dict[str, Path]) -> List[str]:
    bad = []
    for asm in sorted(assemblies):
        try:
            resolve_ref_for_assembly(asm, ref_map)
        except ValueError as e:
            bad.append(str(e))
    return bad


def bed_to_fasta_stem(bed_path: Path) -> str:
    name = bed_path.name
    if name.endswith(".bed.gz"):
        return name[: -len(".bed.gz")]
    if name.endswith(".bed"):
        return name[: -len(".bed")]
    return bed_path.stem


def fasta_path_for_peaktype(out_dir: Path, peak_subdir: str, bed_path: Path) -> Path:
    return out_dir / peak_subdir / "{}.fa".format(bed_to_fasta_stem(bed_path))


def read_path_list(list_path: Path) -> List[Path]:
    paths = []  # type: List[Path]
    with list_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            paths.append(Path(line))
    return paths


def run_getfasta(bed_path: Path, ref_fasta: Path, out_fa: Path, dry_run: bool) -> None:
    if not bed_path.is_file():
        raise FileNotFoundError("BED not found: {}".format(bed_path))
    if not ref_fasta.is_file():
        raise FileNotFoundError("Reference FASTA not found: {}".format(ref_fasta))
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "bedtools",
        "getfasta",
        "-fi",
        str(ref_fasta),
        "-bed",
        str(bed_path),
        "-fo",
        str(out_fa),
    ]
    shell_cmd = (
        "module load bedtools/2.30.0-singularity >/dev/null 2>&1 && "
        "export SINGULARITY_BINDPATH=/vast,/data,/home && "
        "{}".format(" ".join(shlex.quote(x) for x in cmd))
    )
    if dry_run:
        print(shell_cmd)
        return
    subprocess.check_call(["bash", "-lc", shell_cmd])


def process_tasks(
    tasks: List[Tuple[Path, Path, Path]],
    *,
    dry_run: bool,
    skip_existing: bool,
    jobs: int,
) -> Tuple[List[Path], List[str]]:
    """Each task is (bed_path, ref_fasta, out_fa). Return (written paths, error strings)."""
    written = []  # type: List[Path]
    errors = []  # type: List[str]
    if jobs <= 1:
        for bed_path, ref_fasta, out_fa in tasks:
            if skip_existing and out_fa.is_file() and out_fa.stat().st_size > 0:
                written.append(out_fa.resolve())
                continue
            try:
                run_getfasta(bed_path, ref_fasta, out_fa, dry_run)
                written.append(out_fa.resolve())
            except Exception as e:  # noqa: BLE001
                errors.append("{}: {}".format(bed_path, e))
        return written, errors

    str_tasks = [
        (str(b.resolve()), str(r.resolve()), str(o), dry_run, skip_existing)
        for b, r, o in tasks
    ]
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futures = [ex.submit(_getfasta_worker, t) for t in str_tasks]
        for fut in as_completed(futures):
            ok, out_p, err = fut.result()
            if ok and out_p:
                written.append(Path(out_p))
            if err:
                errors.append(err)

    return written, errors


def _getfasta_worker(task: Tuple[str, str, str, bool, bool]) -> Tuple[bool, str, str]:
    """
    Run one bedtools job in a worker process.

    task = (bed_path, ref_fasta, out_fa, dry_run, skip_existing)
    Returns (ok, resolved_out_fa_or_empty, error_message_or_empty).
    """
    bed_s, ref_s, out_s, dry_run, skip_existing = task
    bed_path = Path(bed_s)
    ref_fasta = Path(ref_s)
    out_fa = Path(out_s)
    try:
        if skip_existing and out_fa.is_file() and out_fa.stat().st_size > 0:
            return (True, str(out_fa.resolve()), "")
        run_getfasta(bed_path, ref_fasta, out_fa, dry_run)
        return (True, str(out_fa.resolve()), "")
    except Exception as e:  # noqa: BLE001
        return (False, "", "{}: {}".format(bed_path, e))


def manifest_bed_tasks_from_rows(
    rows: List[RowDict],
    out_dir: Path,
    ref_map: Dict[str, Path],
) -> Tuple[List[Tuple[Path, Path, Path]], List[Tuple[str, str, str]], List[str]]:
    """Return (tasks, meta_for_path_lists, row_errors). Meta entries: (organism, stem, encode_assembly)."""
    tasks = []  # type: List[Tuple[Path, Path, Path]]
    meta = []  # type: List[Tuple[str, str, str]]
    row_errors = []  # type: List[str]
    for row in rows:
        bed_s = (row.get("bed_path") or "").strip()
        org = (row.get("organism") or "").strip()
        if not bed_s:
            continue
        bed = Path(bed_s)
        asm = effective_encode_assembly_for_row(row)
        if not asm:
            row_errors.append(
                "{}: cannot resolve reference (set assembly or reference_build, or a recognized "
                "organism)".format(bed)
            )
            continue
        try:
            ref = resolve_ref_for_assembly(asm, ref_map)
        except ValueError as e:
            row_errors.append("{}: {}".format(bed, e))
            continue
        stem = bed_to_fasta_stem(bed)
        out_fa = out_dir / NARROWPEAK_FASTA_SUBDIR / "{}.fa".format(stem)
        tasks.append((bed, ref, out_fa))
        meta.append((org, stem, asm))
    return tasks, meta, row_errors


def write_path_list(paths: List[Path], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as f:
        for p in sorted(paths):
            f.write("{}\n".format(p))


def write_species_fasta_path_lists(
    path_lists_dir: Path,
    out_dir: Path,
    meta: List[Tuple[str, str, str]],
    *,
    no_path_lists: bool,
) -> None:
    """meta entries are (organism, fasta_stem, encode_assembly); only organism/stem used here."""
    if no_path_lists:
        return
    by_slug = defaultdict(list)  # type: Dict[str, List[Path]]
    for org, stem, _asm in meta:
        norm = normalize_organism_label(org)
        if norm not in ORGANISM_SLUG:
            continue
        slug = ORGANISM_SLUG[norm]
        by_slug[slug].append(
            (out_dir / NARROWPEAK_FASTA_SUBDIR / "{}.fa".format(stem)).resolve()
        )
    for slug, paths in sorted(by_slug.items()):
        existing = sorted({p for p in paths if p.is_file()})
        if existing:
            dest = path_lists_dir / "histone_narrowpeak_fasta_paths_{}.txt".format(slug)
            write_path_list(existing, dest)
            print("Wrote {} ({} paths)".format(dest, len(existing)))


def write_reference_build_fasta_path_lists(
    path_lists_dir: Path,
    out_dir: Path,
    meta: List[Tuple[str, str, str]],
    *,
    no_path_lists: bool,
    basename_prefix: str,
) -> None:
    """meta entries are (organism, fasta_stem, encode_assembly)."""
    if no_path_lists:
        return
    by_build_slug = defaultdict(list)  # type: Dict[str, List[Path]]
    for _org, stem, asm in meta:
        bslug = encode_assembly_to_bed_list_slug(asm)
        by_build_slug[bslug].append(
            (out_dir / NARROWPEAK_FASTA_SUBDIR / "{}.fa".format(stem)).resolve()
        )
    for bslug, paths in sorted(by_build_slug.items()):
        existing = sorted({p for p in paths if p.is_file()})
        if existing:
            dest = path_lists_dir / "{}_{}.txt".format(basename_prefix, bslug)
            write_path_list(existing, dest)
            print("Wrote {} ({} paths)".format(dest, len(existing)))


def generate_tissue_breakout_png(rows: List[RowDict], output_png: Path) -> None:
    """
    Render grouped+stacked bar chart PNG:
    - for each tissue, one bar per organism
    - each bar segmented by reference build.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:  # noqa: BLE001
        raise RuntimeError("matplotlib/numpy unavailable for PNG plotting: {}".format(e))

    counts = defaultdict(lambda: defaultdict(int))  # type: Dict[Tuple[str, str], Dict[str, int]]
    tissues_set = set()  # type: Set[str]
    builds_set = set()  # type: Set[str]
    species_set = set()  # type: Set[str]
    for row in rows:
        tissue = (row.get("tissue_type") or "").strip() or "unknown"
        species = (row.get("organism") or "").strip() or "unknown"
        build = effective_encode_assembly_for_row(row) or "unknown"
        counts[(species, build)][tissue] += 1
        tissues_set.add(tissue)
        builds_set.add(build)
        species_set.add(species)

    tissues = sorted(tissues_set)
    species = sorted(species_set)
    # Keep common species first for readability.
    preferred = ["Mus musculus", "Homo sapiens"]
    species = [s for s in preferred if s in species] + [s for s in species if s not in preferred]
    builds = sorted(builds_set)

    x = np.arange(len(tissues))
    bar_w = 0.38 if len(species) <= 2 else max(0.2, 0.8 / max(1, len(species)))
    group_offsets = np.linspace(
        -((len(species) - 1) * bar_w) / 2.0, ((len(species) - 1) * bar_w) / 2.0, num=len(species)
    )

    fig, ax = plt.subplots(figsize=(14, 4.8))
    cmap = plt.get_cmap("tab20")
    build_color = {}  # type: Dict[str, Tuple[float, float, float, float]]
    for i, b in enumerate(builds):
        build_color[b] = cmap(i % 20)

    for si, sp in enumerate(species):
        xpos = x + group_offsets[si]
        bottom = np.zeros(len(tissues))
        for bi, b in enumerate(builds):
            vals = np.array([counts[(sp, b)].get(t, 0) for t in tissues], dtype=float)
            if not vals.any():
                continue
            label = "{} | {}".format(sp, b) if si == 0 else None
            ax.bar(
                xpos,
                vals,
                width=bar_w * 0.95,
                bottom=bottom,
                color=build_color[b],
                edgecolor="#333333",
                linewidth=0.3,
                label=label,
            )
            bottom += vals

    ax.set_title("Tissue types by species (stacked by reference build)")
    ax.set_xlabel("Tissue type")
    ax.set_ylabel("BED files")
    ax.set_xticks(x)
    ax.set_xticklabels(tissues, rotation=40, ha="right")
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)

    # Compact legend: one entry per build color.
    handles = []
    labels = []
    for b in builds:
        handles.append(plt.Rectangle((0, 0), 1, 1, color=build_color[b]))
        labels.append(b)
    ax.legend(handles, labels, title="Reference build", loc="upper right", frameon=True)

    # Species annotation under x-axis.
    if len(species) == 2:
        ax.text(0.01, -0.22, "left bar = {}".format(species[0]), transform=ax.transAxes, fontsize=9)
        ax.text(0.35, -0.22, "right bar = {}".format(species[1]), transform=ax.transAxes, fontsize=9)

    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_png), dpi=180)
    plt.close(fig)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="bedtools getfasta for histone narrowPeak BEDs (per-build lists, or manifest + assembly)."
    )
    p.add_argument(
        "--histone-marks-dir",
        type=Path,
        default=None,
        help="Directory containing bed/, fastq/, and path_lists/ (or use --manifest)",
    )
    p.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="histone_narrowpeak_bed_manifest.tsv (assembly and/or reference_build columns; see build script)",
    )
    p.add_argument(
        "--fasta-out-dir",
        type=Path,
        default=None,
        help="Output directory for .fa files (default: <histone-marks-dir>/fasta_from_bed)",
    )
    p.add_argument(
        "--ref-mm10",
        type=Path,
        default=None,
        help="mm10 reference FASTA (default: $MM10_FASTA or mus_homo/references/mm10/…)",
    )
    p.add_argument(
        "--ref-grch38",
        type=Path,
        default=None,
        help="GRCh38 reference FASTA (default: $GRCH38_FASTA or mus_homo/references/grch38/…)",
    )
    p.add_argument(
        "--ref-hg19",
        type=Path,
        default=None,
        help="hg19 FASTA (default: $HG19_FASTA or mus_homo/references/hg19/hg19.fa)",
    )
    p.add_argument(
        "--ref-mm9",
        type=Path,
        default=None,
        help="mm9 FASTA (default: $MM9_FASTA or mus_homo/references/mm9/mm9.fa)",
    )
    p.add_argument(
        "--ref-mm10-minimal",
        type=Path,
        default=None,
        dest="ref_mm10_minimal",
        help="mm10-minimal FASTA (default: $MM10_MINIMAL_FASTA or mus_homo/references/mm10_minimal/…)",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip getfasta if output .fa exists and is non-empty",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print bedtools commands only",
    )
    p.add_argument(
        "--no-path-lists",
        action="store_true",
        help="Do not write any histone_*_fasta_paths_*.txt lists",
    )
    p.add_argument(
        "--no-tissue-plot",
        action="store_true",
        help="Do not generate tissue_type_breakout_by_build_and_organism.png in --manifest mode",
    )
    p.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=None,
        metavar="N",
        help="Parallel bedtools processes (default: min(8, CPU count); use 1 to disable)",
    )
    args = p.parse_args(argv)

    ncpu = multiprocessing.cpu_count() or 1
    jobs = args.jobs if args.jobs is not None else min(8, max(1, ncpu))
    if jobs < 1:
        jobs = 1

    ref_map = build_ref_map(args)
    if args.manifest is None:
        p.error("--manifest is required (legacy non-manifest mode has been removed)")

    mp = args.manifest.resolve()
    if not mp.is_file():
        print("error: manifest not found: {}".format(mp), file=sys.stderr)
        return 2
    hm = (args.histone_marks_dir or mp.parent).resolve()
    out_dir = (args.fasta_out_dir or (hm / "fasta_from_bed")).resolve()
    path_lists_dir = (hm / PATH_LISTS_DIRNAME).resolve()

    with mp.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)  # type: List[RowDict]

    has_asm_col = "assembly" in fields
    nrows = len(rows)
    asm_nonempty = sum(1 for r in rows if (r.get("assembly") or "").strip())
    if has_asm_col and asm_nonempty not in (0, nrows):
        print(
            "error: manifest has mixed empty/non-empty assembly values; fix or rebuild with "
            "--fetch-assembly",
            file=sys.stderr,
        )
        return 2
    if (not has_asm_col) or asm_nonempty == 0:
        print(
            "warning: manifest has no usable per-file assembly; assuming mm10 (mouse) and "
            "GRCh38 (human). Rebuild with build_histone_narrowpeak_bed_manifest.py "
            "--fetch-assembly for correct hg19 / mm9 / mm10-minimal handling.",
            file=sys.stderr,
        )

    assemblies = set()  # type: Set[str]
    for r in rows:
        eff = effective_encode_assembly_for_row(r)
        if eff:
            assemblies.add(eff)

    bad_refs = validate_refs_for_assemblies(assemblies, ref_map)
    if bad_refs:
        for msg in bad_refs:
            print("error: {}".format(msg), file=sys.stderr)
        return 2

    tasks, meta, row_errors = manifest_bed_tasks_from_rows(rows, out_dir, ref_map)
    if row_errors:
        for msg in row_errors[:40]:
            print("error: {}".format(msg), file=sys.stderr)
        if len(row_errors) > 40:
            print("error: ... and {} more manifest row errors".format(len(row_errors) - 40), file=sys.stderr)
        return 2

    print(
        "Processing {} BEDs from manifest {} (jobs={})".format(len(tasks), mp, jobs)
    )
    all_written, errs_all = process_tasks(
        tasks,
        dry_run=args.dry_run,
        skip_existing=args.skip_existing,
        jobs=jobs,
    )
    for err in errs_all[:20]:
        print("error: {}".format(err), file=sys.stderr)
    if len(errs_all) > 20:
        print("error: ... and {} more".format(len(errs_all) - 20), file=sys.stderr)
    write_species_fasta_path_lists(
        path_lists_dir, out_dir, meta, no_path_lists=args.no_path_lists
    )
    write_reference_build_fasta_path_lists(
        path_lists_dir,
        out_dir,
        meta,
        no_path_lists=args.no_path_lists,
        basename_prefix="histone_narrowpeak_fasta_paths",
    )
    # Also process broadPeak lists (if present) from path_lists/ using build-specific refs.
    broad_lists = sorted(path_lists_dir.glob(BROADPEAK_LIST_PREFIX + "*.txt"))
    for lp in broad_lists:
        slug = lp.name[len(BROADPEAK_LIST_PREFIX) : -len(".txt")]
        asm = slug if slug != "mm10_minimal" else "mm10-minimal"
        if asm not in KNOWN_ENCODE_ASSEMBLY:
            print(
                "warning: skipping broadPeak list with unknown build slug {!r}: {}".format(
                    slug, lp
                ),
                file=sys.stderr,
            )
            continue
        try:
            ref = resolve_ref_for_assembly(asm, ref_map)
        except ValueError as e:
            print("error: {}: {}".format(lp, e), file=sys.stderr)
            return 2
        list_tasks = []  # type: List[Tuple[Path, Path, Path]]
        for bed in read_path_list(lp):
            out_fa = fasta_path_for_peaktype(out_dir, BROADPEAK_FASTA_SUBDIR, bed)
            list_tasks.append((bed, ref, out_fa))
        print(
            "Processing broadPeak BEDs for build {} from {} with {} (jobs={})".format(
                asm, lp, ref, jobs
            )
        )
        w2, e2 = process_tasks(
            list_tasks,
            dry_run=args.dry_run,
            skip_existing=args.skip_existing,
            jobs=jobs,
        )
        all_written.extend(w2)
        errs_all.extend(e2)
        if not args.no_path_lists:
            existing = sorted({p for p in w2 if p.is_file()})
            if existing:
                dest = path_lists_dir / "histone_broadpeak_fasta_paths_{}.txt".format(slug)
                write_path_list(existing, dest)
                print("Wrote {} ({} paths)".format(dest, len(existing)))
    if not args.no_tissue_plot and not args.dry_run:
        plot_png = hm / "tissue_type_breakout_by_build_and_organism.png"
        try:
            generate_tissue_breakout_png(rows, plot_png)
            print("Wrote {}".format(plot_png))
        except Exception as e:  # noqa: BLE001
            print(
                "warning: could not generate tissue breakout PNG: {}".format(e),
                file=sys.stderr,
            )
    print("FASTA directory: {} ({} outputs touched this run)".format(out_dir, len(all_written)))
    return 1 if errs_all else 0


if __name__ == "__main__":
    raise SystemExit(main())
