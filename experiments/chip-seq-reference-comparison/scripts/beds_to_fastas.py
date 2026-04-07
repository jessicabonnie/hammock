#!/usr/bin/env python3
"""
Run ``bedtools getfasta`` for a list of BED files against one reference build.

Input list format:
- One BED path per line
- Empty lines and lines starting with '#' are ignored

Output:
- Writes ``<bed-stem>.fa`` files in ``--output-dir``
- Prints one absolute ``.fa`` path per line on stdout (errors only on stderr)

Reference build selection:
- ``--reference`` must be one of: mm10, mm10-minimal, mm9, GRCh38, hg19
- Aliases: mm10_minimal -> mm10-minimal, GRCm38 -> mm10, GRCh37 -> hg19, NCBIM37 -> mm9
"""

import argparse
import multiprocessing
import os
import shlex
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Lab reference tree.
_REF_BASE = "/data/blangme2/jessica/mus_homo/references"
DEFAULT_MM10_FASTA = _REF_BASE + "/mm10/mm10_no_alt_analysis_set_ENCODE.fasta"
DEFAULT_GRCH38_FASTA = _REF_BASE + "/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
DEFAULT_HG19_FASTA = _REF_BASE + "/hg19/hg19.fa"
DEFAULT_MM9_FASTA = _REF_BASE + "/mm9/mm9.fa"
DEFAULT_MM10_MINIMAL_FASTA = _REF_BASE + "/mm10_minimal/mm10_minimal_ENCODE.fasta"

KNOWN_REFERENCE = frozenset({"mm10", "mm10-minimal", "mm9", "GRCh38", "hg19"})


def _ref_from_arg_or_env(arg_path: Optional[Path], env_key: str, default: str) -> Path:
    if arg_path is not None:
        return arg_path
    ev = os.environ.get(env_key, "").strip()
    if ev:
        return Path(ev)
    return Path(default)


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


def normalize_reference(reference: str) -> str:
    ref = (reference or "").strip()
    aliases = {
        "mm10_minimal": "mm10-minimal",
        "GRCm38": "mm10",
        "GRCh37": "hg19",
        "NCBIM37": "mm9",
    }
    ref = aliases.get(ref, ref)
    if ref not in KNOWN_REFERENCE:
        raise ValueError(
            "--reference {!r} is not recognized (known: {})".format(
                reference, ", ".join(sorted(KNOWN_REFERENCE))
            )
        )
    return ref


def resolve_ref_for_reference(reference: str, ref_map: Dict[str, Path]) -> Path:
    p = ref_map.get(reference, Path(""))
    if not str(p) or not p.is_file():
        raise ValueError(
            "no valid FASTA for reference {!r}; set {} or matching env var".format(
                reference,
                {
                    "mm10": "--ref-mm10 / $MM10_FASTA",
                    "mm10-minimal": "--ref-mm10-minimal / $MM10_MINIMAL_FASTA",
                    "mm9": "--ref-mm9 / $MM9_FASTA",
                    "GRCh38": "--ref-grch38 / $GRCH38_FASTA",
                    "hg19": "--ref-hg19 / $HG19_FASTA",
                }[reference],
            )
        )
    return p


def bed_to_fasta_stem(bed_path: Path) -> str:
    name = bed_path.name
    if name.endswith(".bed.gz"):
        return name[: -len(".bed.gz")]
    if name.endswith(".bed"):
        return name[: -len(".bed")]
    return bed_path.stem


def read_bed_list(list_path: Path) -> List[Path]:
    beds = []  # type: List[Path]
    seen = set()  # type: Set[str]
    with list_path.open(encoding="utf-8") as f:
        for line in f:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            p = Path(raw).expanduser()
            if not p.is_absolute():
                p = (list_path.parent / p).resolve()
            key = str(p)
            if key in seen:
                continue
            seen.add(key)
            beds.append(p)
    return beds


def _first_bed_chrom(bed_path: Path) -> str:
    with bed_path.open(encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            return s.split("\t", 1)[0].split(" ", 1)[0]
    return ""


def _first_fasta_contig(ref_fasta: Path) -> str:
    fai = Path(str(ref_fasta) + ".fai")
    if not fai.is_file():
        return ""
    with fai.open(encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            return s.split("\t", 1)[0]
    return ""


def _infer_chr_mode(bed_path: Path, ref_fasta: Path) -> str:
    bed_chrom = _first_bed_chrom(bed_path)
    fasta_contig = _first_fasta_contig(ref_fasta)
    if not bed_chrom or not fasta_contig:
        return "none"
    bed_has_chr = bed_chrom.startswith("chr")
    fasta_has_chr = fasta_contig.startswith("chr")
    if bed_has_chr == fasta_has_chr:
        return "none"
    return "add_chr" if fasta_has_chr else "strip_chr"


def _rewrite_bed_chrom_style(bed_path: Path, mode: str) -> Path:
    fd, tmp_path = tempfile.mkstemp(prefix=bed_path.stem + ".", suffix=".bed")
    os.close(fd)
    out = Path(tmp_path)
    with bed_path.open(encoding="utf-8") as fin, out.open("w", encoding="utf-8") as fout:
        for line in fin:
            s = line.rstrip("\n")
            if not s or s.startswith("#"):
                fout.write(line)
                continue
            parts = s.split("\t")
            if not parts:
                fout.write(line)
                continue
            chrom = parts[0]
            if mode == "add_chr" and not chrom.startswith("chr"):
                parts[0] = "chr{}".format(chrom)
            elif mode == "strip_chr" and chrom.startswith("chr"):
                parts[0] = chrom[3:]
            fout.write("\t".join(parts) + "\n")
    return out


def run_getfasta(bed_path: Path, ref_fasta: Path, out_fa: Path, dry_run: bool) -> None:
    if not bed_path.is_file():
        raise FileNotFoundError("BED not found: {}".format(bed_path))
    if not ref_fasta.is_file():
        raise FileNotFoundError("Reference FASTA not found: {}".format(ref_fasta))
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    temp_bed = None  # type: Optional[Path]
    bed_for_getfasta = bed_path
    mode = _infer_chr_mode(bed_path, ref_fasta)
    if mode != "none":
        temp_bed = _rewrite_bed_chrom_style(bed_path, mode)
        bed_for_getfasta = temp_bed
    if dry_run:
        if temp_bed is not None and temp_bed.exists():
            temp_bed.unlink()
        return
    cmd = [
        "bedtools",
        "getfasta",
        "-fi",
        str(ref_fasta),
        "-bed",
        str(bed_for_getfasta),
        "-fo",
        str(out_fa),
    ]
    shell_cmd = (
        "module load bedtools/2.30.0-singularity >/dev/null 2>&1 && "
        "export SINGULARITY_BINDPATH=/vast,/data,/home && "
        "{}".format(" ".join(shlex.quote(x) for x in cmd))
    )
    try:
        subprocess.check_call(["bash", "-lc", shell_cmd])
    finally:
        if temp_bed is not None and temp_bed.exists():
            temp_bed.unlink()


def _worker(task: Tuple[str, str, str, bool, bool]) -> Tuple[bool, str, str]:
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


def process_tasks(
    tasks: List[Tuple[Path, Path, Path]],
    *,
    dry_run: bool,
    skip_existing: bool,
    jobs: int,
) -> Tuple[List[Path], List[str]]:
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
        futures = [ex.submit(_worker, t) for t in str_tasks]
        for fut in as_completed(futures):
            ok, out_p, err = fut.result()
            if ok and out_p:
                written.append(Path(out_p))
            if err:
                errors.append(err)
    return written, errors


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="bedtools getfasta for a BED list and one reference build"
    )
    p.add_argument("--bed-list", type=Path, required=True, help="Text file with one BED path per line")
    p.add_argument("--output-dir", type=Path, required=True, help="Output directory for .fa files")
    p.add_argument(
        "--reference",
        type=str,
        required=True,
        help="Reference build: mm10, mm10-minimal, mm9, GRCh38, or hg19",
    )
    p.add_argument("--ref-mm10", type=Path, default=None, help="Override mm10 FASTA path")
    p.add_argument("--ref-grch38", type=Path, default=None, help="Override GRCh38 FASTA path")
    p.add_argument("--ref-hg19", type=Path, default=None, help="Override hg19 FASTA path")
    p.add_argument("--ref-mm9", type=Path, default=None, help="Override mm9 FASTA path")
    p.add_argument(
        "--ref-mm10-minimal",
        type=Path,
        default=None,
        dest="ref_mm10_minimal",
        help="Override mm10-minimal FASTA path",
    )
    p.add_argument("--skip-existing", action="store_true", help="Skip non-empty existing outputs")
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Do not run bedtools; still print intended output .fa paths on stdout",
    )
    p.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=None,
        metavar="N",
        help="Parallel bedtools processes (default: min(8, CPU count); use 1 for sequential)",
    )
    args = p.parse_args(argv)

    bed_list = args.bed_list.resolve()
    if not bed_list.is_file():
        print("error: --bed-list not found: {}".format(bed_list), file=sys.stderr)
        return 2

    try:
        reference = normalize_reference(args.reference)
    except ValueError as e:
        print("error: {}".format(e), file=sys.stderr)
        return 2

    ref_map = build_ref_map(args)
    try:
        ref_fasta = resolve_ref_for_reference(reference, ref_map)
    except ValueError as e:
        print("error: {}".format(e), file=sys.stderr)
        return 2

    beds = read_bed_list(bed_list)
    if not beds:
        print("error: no BED paths found in {}".format(bed_list), file=sys.stderr)
        return 2

    ncpu = multiprocessing.cpu_count() or 1
    jobs = args.jobs if args.jobs is not None else min(8, max(1, ncpu))
    if jobs < 1:
        jobs = 1

    out_dir = args.output_dir.resolve()
    tasks = []  # type: List[Tuple[Path, Path, Path]]
    for bed in beds:
        out_fa = out_dir / "{}.fa".format(bed_to_fasta_stem(bed))
        tasks.append((bed, ref_fasta, out_fa))

    written, errors = process_tasks(
        tasks,
        dry_run=args.dry_run,
        skip_existing=args.skip_existing,
        jobs=jobs,
    )

    for err in errors[:20]:
        print("error: {}".format(err), file=sys.stderr)
    if len(errors) > 20:
        print("error: ... and {} more".format(len(errors) - 20), file=sys.stderr)

    for p in sorted(written, key=lambda x: str(x)):
        print(p.resolve())
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
