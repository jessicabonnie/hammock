"""Guard: passing an actual BED/FASTA file where a list-of-paths file is
expected must fail early with a clear message (not a confusing downstream error)."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).parent
DATA = HERE / "data"
OURS = shutil.which("hammock")

pytestmark = pytest.mark.skipif(OURS is None, reason="hammock not on PATH")


def _run(args, cwd):
    return subprocess.run([OURS, *args], cwd=cwd, capture_output=True, text=True)


def test_bed_file_as_positional_is_rejected(tmp_path: Path) -> None:
    # tiny_a.bed is a real BED file; passing it directly (not a list) must error.
    bed = str(DATA / "tiny_a.bed")
    r = _run([bed, bed, "-o", str(tmp_path / "x")], tmp_path)
    assert r.returncode != 0
    assert "looks like a BED file" in r.stderr
    assert "list of paths" in r.stderr


def test_fasta_file_as_positional_is_rejected(tmp_path: Path) -> None:
    fa = str(DATA / "tiny.fa")
    r = _run([fa, fa, "-o", str(tmp_path / "x")], tmp_path)
    assert r.returncode != 0
    assert "looks like a FASTA file" in r.stderr


def test_proper_list_file_is_accepted(tmp_path: Path) -> None:
    # A real list-of-paths file must NOT trip the guard.
    lst = tmp_path / "list.txt"
    lst.write_text(str(DATA / "tiny_a.bed") + "\n")
    r = _run([str(lst), str(lst), "-p", "12", "-o", str(tmp_path / "ok")], tmp_path)
    assert r.returncode == 0, r.stderr
    assert next(tmp_path.glob("ok*.csv")).exists()
