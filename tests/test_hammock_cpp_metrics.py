"""The standalone `hammock-cpp --metrics` block must match the Python CLI exactly.

Both programs reach the same quantities from different starting points -- the
Python side derives them from the containment arrays `pairwise_metrics_hll`
already returns, the C++ side computes a union per pair -- so this is the test
that keeps them from drifting. The derivation expression is deliberately
written the same way in both (`runner._jaccard_ie_from_containments` and
`jaccard_ie_from_containments` in hammock_cli.cpp), which is what makes exact
equality achievable rather than merely close agreement.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parent.parent
_METRIC_COLS = [
    "jaccard_similarity", "jaccard_similarity_ie",
    "containment_AB", "containment_BA",
    "cosketch_geom", "cosketch_arith", "cosketch_max",
]


def _find_hammock_cpp() -> Path | None:
    """The binary is built into build/<wheel-tag>/ and is NOT installed into
    the wheel (CMakeLists: users run `cmake --install` separately)."""
    hits = sorted(_REPO.glob("build/*/hammock-cpp"))
    return hits[0] if hits else None


_BIN = _find_hammock_cpp()
pytestmark = pytest.mark.skipif(
    _BIN is None, reason="hammock-cpp not built (see build/<tag>/hammock-cpp)")


@pytest.fixture
def corpus(tmp_path: Path):
    """Two BED files with substantial but unequal overlap against a reference."""
    import random

    def write(name: str, n: int, seed: int, shift: int) -> Path:
        r = random.Random(seed)
        p = tmp_path / name
        pos = 0
        lines = []
        for _ in range(n):
            pos += r.randint(500, 5000)
            lines.append(f"chr1\t{pos + shift}\t{pos + shift + r.randint(200, 2000)}\n")
        p.write_text("".join(lines))
        return p

    a = write("sample1.bed", 2000, 1, 0)
    b = write("sample2.bed", 800, 2, 0)
    ref = write("ref.bed", 2000, 1, 300)
    q = tmp_path / "q.txt"
    r = tmp_path / "r.txt"
    q.write_text(f"{a}\n{b}\n")
    r.write_text(f"{ref}\n")
    return q, r


def _run_python(q: Path, r: Path, out: Path) -> dict:
    subprocess.run(
        [sys.executable, "-m", "hammock.cli", str(q), str(r),
         "--mode", "B", "-p", "20", "-o", str(out)],
        check=True, capture_output=True, text=True, timeout=600)
    csv_path = next(out.parent.glob(f"{out.name}*.csv"))
    with csv_path.open() as f:
        return {(row["file1"], row["file2"]): row for row in csv.DictReader(f)}


def _run_cpp(q: Path, r: Path, out: Path, metrics: bool) -> dict:
    cmd = [str(_BIN), str(q), str(r), "--mode", "B", "-p", "20", "-o", str(out)]
    if metrics:
        cmd.append("--metrics")
    subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
    csv_path = next(out.parent.glob(f"{out.name}*.csv"))
    with csv_path.open() as f:
        return {(row["query"], row["reference"]): row
                for row in csv.DictReader(f, delimiter="\t")}


def test_metrics_block_matches_python_bit_for_bit(corpus, tmp_path: Path):
    q, r = corpus
    py = _run_python(q, r, tmp_path / "py")
    cpp = _run_cpp(q, r, tmp_path / "cpp", metrics=True)

    assert set(py) == set(cpp), (sorted(py), sorted(cpp))
    for key in py:
        for col in _METRIC_COLS:
            # Exact equality, not approx: same IEEE operations in the same
            # order on both sides. If this starts failing, the two derivations
            # have diverged -- fix that rather than loosening the assertion.
            assert float(py[key][col]) == float(cpp[key][col]), (key, col)


def test_default_output_shape_is_unchanged(corpus, tmp_path: Path):
    """Without --metrics the binary must still emit exactly three columns, so
    published benchmark timings and their parsers stay valid."""
    q, r = corpus
    out = tmp_path / "plain"
    _run_cpp(q, r, out, metrics=False)
    header = next(out.parent.glob(f"{out.name}*.csv")).read_text().splitlines()[0]
    assert header.split("\t") == ["query", "reference", "jaccard_similarity"]


def test_metrics_changes_the_output_filename(corpus, tmp_path: Path):
    """A 3-column and a 9-column file must not collide on one path."""
    q, r = corpus
    _run_cpp(q, r, tmp_path / "same", metrics=False)
    _run_cpp(q, r, tmp_path / "same", metrics=True)
    names = sorted(p.name for p in tmp_path.glob("same*.csv"))
    assert len(names) == 2, names
    assert any(n.endswith("_metrics.csv") for n in names), names


def test_metrics_rejected_with_peak_height(corpus, tmp_path: Path):
    """BagMinHash cardinality and intersection are on different scales, so
    inclusion-exclusion over them would silently return a constant 0."""
    q, r = corpus
    proc = subprocess.run(
        [str(_BIN), str(q), str(r), "--mode", "A", "-o", str(tmp_path / "x"),
         "--metrics", "--peak-height", "5"],
        capture_output=True, text=True, timeout=600)
    assert proc.returncode != 0
    assert "peak-height" in proc.stderr
