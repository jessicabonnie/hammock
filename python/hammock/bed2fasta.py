"""BED → FASTA extraction for Mode D, via a single backend: ``bedtools getfasta``.

One backend, one contract: default strand handling (no reverse-complement) and
``chrom:start-end`` headers, so both sides of a cross-reference comparison see
comparable sequences. (``twoBitToFa -bed`` was rejected because it
reverse-complements minus-strand intervals and headers by BED-name by default,
which would silently diverge from ``getfasta`` and corrupt
``jaccard_similarity_with_ends``.)

Every conversion is validated: an empty output FASTA or a "chromosome not
found" warning from bedtools (the ``chr1`` vs ``1`` naming mismatch) is turned
into a hard error rather than flowing downstream as a fake ``Jaccard=0`` — the
project has a documented history of exactly that silent-zero failure.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from typing import List


class ConversionError(RuntimeError):
    """A BED→FASTA conversion failed or produced suspect output."""


def _require_bedtools() -> str:
    path = shutil.which("bedtools")
    if path is None:
        raise ConversionError(
            "bedtools not found on PATH. Load it first, e.g. `ml bedtools` "
            "(or `ml UCSC_Genome_Browser/2021`), then re-run.")
    return path


def bedtools_version() -> str:
    try:
        out = subprocess.run([_require_bedtools(), "--version"],
                             check=True, capture_output=True, text=True)
        return out.stdout.strip() or out.stderr.strip()
    except (ConversionError, subprocess.CalledProcessError):
        return "unknown"


def _count_bed_intervals(bed_path: str) -> int:
    n = 0
    opener = open
    if bed_path.lower().endswith(".gz"):
        import gzip
        opener = lambda p: gzip.open(p, "rt")  # noqa: E731
    with opener(bed_path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue
            n += 1
    return n


def _fasta_stats(path: str, scan_cap: int = 8 * 1024 * 1024) -> tuple[int, float]:
    """Return (record_count, n_fraction). ``n_fraction`` is estimated over at
    most ``scan_cap`` sequence bytes to stay cheap on large outputs."""
    records = 0
    n_bases = 0
    total_bases = 0
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                records += 1
                continue
            if total_bases < scan_cap:
                seq = line.strip()
                total_bases += len(seq)
                n_bases += seq.upper().count("N")
    frac = (n_bases / total_bases) if total_bases else 0.0
    return records, frac


def bed_to_fasta(bed_path: str, ref_fasta: str, out_fasta: str,
                 verbose: bool = False) -> None:
    """Extract sequences for ``bed_path`` from ``ref_fasta`` into ``out_fasta``.

    Raises ConversionError on a failed run, an empty result, or a chromosome
    naming mismatch (bedtools warns + exits 0 but skips the intervals).
    """
    bedtools = _require_bedtools()
    proc = subprocess.run(
        [bedtools, "getfasta", "-fi", ref_fasta, "-bed", bed_path, "-fo", out_fasta],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise ConversionError(
            f"bedtools getfasta failed for {bed_path} against {ref_fasta}\n"
            f"{proc.stderr.strip()}")

    # Match bedtools' specific chromosome-mismatch phrasing rather than any
    # "warning": getfasta prints a benign "Warning: the index file is older
    # than the FASTA file." on *successful* exit-0 runs (common on clusters
    # where the .fai mtime predates the .fa), which must NOT be treated as a
    # failure. The mismatch line is: "WARNING. chromosome (chr1) was not found
    # in the FASTA file. Skipping."
    stderr = proc.stderr or ""
    if "not found in the fasta file" in stderr.lower():
        raise ConversionError(
            f"bedtools getfasta could not find some chromosomes from {bed_path} "
            f"in {ref_fasta} — a chromosome-name mismatch (e.g. 'chr1' vs '1'). "
            f"Intervals were silently skipped; refusing to continue.\n"
            f"bedtools said: {stderr.strip()}")

    if not os.path.exists(out_fasta) or os.path.getsize(out_fasta) == 0:
        raise ConversionError(
            f"bedtools getfasta produced an empty FASTA for {bed_path} against "
            f"{ref_fasta}. Check that the BED chromosome names match the "
            f"reference (e.g. 'chr1' vs '1').")

    # getfasta emits exactly one record per BED interval; a shortfall means
    # intervals were skipped (chr mismatch or out-of-range). Treat as an error
    # rather than let a partial result flow into Mode D as a fake low Jaccard.
    records, n_frac = _fasta_stats(out_fasta)
    n_intervals = _count_bed_intervals(bed_path)
    if records < n_intervals:
        raise ConversionError(
            f"{os.path.basename(bed_path)}: extracted only {records} sequences "
            f"from {n_intervals} BED intervals against {ref_fasta} — some were "
            f"skipped (chromosome-name mismatch or out-of-range coordinates). "
            f"Refusing to continue.\nbedtools said: {stderr.strip()}")
    if n_frac > 0.5:
        print(f"[hammock] warning: {os.path.basename(bed_path)}: extracted "
              f"sequence is {n_frac:.0%} N (assembly-gap regions). N-runs create "
              f"spurious shared minimizers and can inflate Jaccard.",
              file=sys.stderr)
    if verbose:
        print(f"[hammock] {os.path.basename(bed_path)} → {out_fasta} "
              f"({records} records)", file=sys.stderr)


def _out_name(i: int, bed_path: str) -> str:
    """Unique-per-index output name so two BEDs sharing a basename (different
    dirs) never collide on the same output path."""
    base = os.path.basename(bed_path)
    for ext in (".gz", ".bed", ".narrowpeak", ".broadpeak", ".gappedpeak", ".bb"):
        if base.lower().endswith(ext):
            base = base[: -len(ext)]
    return f"{i:04d}_{base}.fa"


def convert_list(bed_paths: List[str], ref_fasta: str, out_dir: str,
                 threads: int = 1, verbose: bool = False) -> List[str]:
    """Convert every BED in ``bed_paths`` to a FASTA in ``out_dir`` (parallel).

    Returns the FASTA paths in the same order as ``bed_paths``.
    """
    from hammock.runner import _parallel_map  # lazy: avoid import cycle

    os.makedirs(out_dir, exist_ok=True)
    out_paths = [os.path.join(out_dir, _out_name(i, b))
                 for i, b in enumerate(bed_paths)]

    def _one(i: int) -> str:
        bed_to_fasta(bed_paths[i], ref_fasta, out_paths[i], verbose=verbose)
        return out_paths[i]

    return _parallel_map(len(bed_paths), _one, threads)
