"""Mode D (FASTA) unit tests.

These do NOT compare against hammock-orig — bioconda's `digest` package only
ships wheels for Python ≥ 3.9 and the orig pipx venv runs on 3.8, so the
runtime parity gate is deferred (see project_mixed_stride_rerun.md /
plan: phase 4 work).

Instead we verify:
  1. the pipeline runs without errors and produces a valid CSV header,
  2. self-comparison gives Jaccard = 1.0,
  3. repeated runs are byte-identical (determinism),
  4. the empty-minimizer fallback path is taken when seq is shorter
     than the digest window.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).parent
DATA = HERE / "data"

OURS = shutil.which("hammock")

pytestmark = pytest.mark.skipif(OURS is None, reason="hammock not on PATH")


def _files_list(tmp_path: Path, *names: str) -> Path:
    f = tmp_path / "files.txt"
    f.write_text("\n".join(str(DATA / n) for n in names) + "\n")
    return f


def _run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def test_mode_d_runs_and_self_jaccard_is_one(tmp_path: Path) -> None:
    files = _files_list(tmp_path, "tiny.fa")
    _run([OURS, str(files), str(files), "--mode", "D",
          "-p", "14", "-k", "8", "-w", "40",
          "-o", str(tmp_path / "out")], tmp_path)
    csv = next(tmp_path.glob("out*.csv")).read_text().splitlines()
    cont_block = ["containment_AB", "containment_BA",
                  "cosketch_geom", "cosketch_arith", "cosketch_max"]
    assert csv[0].split(",") == [
        "file1", "file2", "sketch_type", "mode",
        "precision", "num_hashes", "kmer_size", "window_size",
        "jaccard_similarity", "jaccard_similarity_ie", *cont_block,
        "ref1", "ref2",
    ]
    # self-pair must be Jaccard = 1.0 (and containment = 1.0).
    header = csv[0].split(",")
    self_row = csv[1].split(",")
    assert self_row[0] == "tiny.fa" and self_row[1] == "tiny.fa"
    # Plain FASTA Mode D (no reference flags): ref1/ref2 columns are "NA".
    assert self_row[header.index("ref1")] == "NA"
    assert self_row[header.index("ref2")] == "NA"
    # Exact 1.0, not approx, for the inclusion-exclusion columns too: for
    # identical registers the union estimate equals both cardinalities
    # bitwise, so inter = 2c - c = c exactly and 1/(1 + 1 - 1) == 1.0.
    for name in ("jaccard_similarity", "jaccard_similarity_ie",
                 "containment_AB", "containment_BA", "cosketch_geom"):
        assert float(self_row[header.index(name)]) == 1.0, name


def test_mode_d_is_deterministic(tmp_path: Path) -> None:
    files = _files_list(tmp_path, "tiny.fa", "tiny2.fa")
    common = [str(files), str(files), "--mode", "D",
              "-p", "14", "-k", "8", "-w", "40"]
    _run([OURS, *common, "-o", str(tmp_path / "run1")], tmp_path)
    _run([OURS, *common, "-o", str(tmp_path / "run2")], tmp_path)
    a = next(tmp_path.glob("run1*.csv")).read_text()
    b = next(tmp_path.glob("run2*.csv")).read_text()
    assert a == b


def test_mode_d_threads_match_sequential(tmp_path: Path) -> None:
    """--threads N must produce byte-equal output to single-threaded run."""
    files = _files_list(tmp_path, "tiny.fa", "tiny2.fa")
    common = [str(files), str(files), "--mode", "D",
              "-p", "14", "-k", "8", "-w", "40"]
    _run([OURS, *common, "-o", str(tmp_path / "seq")], tmp_path)
    _run([OURS, *common, "--threads", "4", "-o", str(tmp_path / "par")], tmp_path)
    seq = next(tmp_path.glob("seq*.csv")).read_text()
    par = next(tmp_path.glob("par*.csv")).read_text()
    assert seq == par


def test_mode_d_distinct_files_have_nontrivial_jaccard(tmp_path: Path) -> None:
    """tiny.fa vs tiny2.fa: Jaccard should be < 1 (they share some content but differ)."""
    files1 = _files_list(tmp_path, "tiny.fa")
    files2 = tmp_path / "primary.txt"
    files2.write_text(str(DATA / "tiny2.fa") + "\n")
    _run([OURS, str(files1), str(files2), "--mode", "D",
          "-p", "14", "-k", "8", "-w", "40",
          "-o", str(tmp_path / "cross")], tmp_path)
    lines = next(tmp_path.glob("cross*.csv")).read_text().splitlines()
    header = lines[0].split(",")
    row = lines[1].split(",")
    jac = float(row[header.index("jaccard_similarity")])
    assert 0.0 <= jac < 1.0, f"expected non-self Jaccard in [0, 1), got {jac}"


# Direct unit tests on the Python helpers -----------------------------------


def test_minimizer_sketch_empty_string_is_noop() -> None:
    from hammock.modes.sequence import MinimizerSketch
    s = MinimizerSketch(kmer_size=8, window_size=40, precision=10)
    s.add_string("")
    assert s.minimizer_hll.estimate_cardinality() == 0.0


def test_minimizer_sketch_short_seq_uses_empty_fallback() -> None:
    """A sequence shorter than k+w-1 yields no minimizers and falls through to
    the 'no minimizers' branch, which adds the whole seq as one element.

    That element is an *exact-match* indicator -- it only ever collides with a
    byte-identical record -- so sub-threshold records contribute identity, not
    graded similarity. See divergence #8 in CLAUDE.md.
    """
    from hammock.modes.sequence import MinimizerSketch
    s = MinimizerSketch(kmer_size=8, window_size=40, precision=10)
    s.add_string("ACGTACGT")  # len 8 < window 40 → digest returns []
    assert 0.5 < s.minimizer_hll.estimate_cardinality() < 2.0
