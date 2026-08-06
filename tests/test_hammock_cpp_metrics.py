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
import os
import subprocess
import sys
from pathlib import Path

import pytest

import hammock

_REPO = Path(__file__).resolve().parent.parent
_METRIC_COLS = [
    "jaccard_similarity", "jaccard_similarity_ie",
    "containment_AB", "containment_BA",
    "cosketch_geom", "cosketch_arith", "cosketch_max",
]


def _find_hammock_cpp() -> Path | None:
    """The binary is built into build/<wheel-tag>/. The wheel does install it
    (to <site-packages>/bin/), but that directory is not on PATH, so tests use
    the build tree copy.

    Pick the NEWEST, not the alphabetically first: with both cp310 and cp311
    build dirs present, sorted()[0] would pin a stale binary and "bit-for-bit"
    would then be asserted against old C++.
    """
    hits = list(_REPO.glob("build/*/hammock-cpp"))
    return max(hits, key=lambda p: p.stat().st_mtime) if hits else None


_BIN = _find_hammock_cpp()
# `build/` is gitignored and there is no CI, so a wheel install makes every
# test here vanish -- silently, as a bare `s`. This is the only cross-
# implementation check in the suite, so allow it to be made mandatory.
pytestmark = pytest.mark.skipif(
    _BIN is None and not os.environ.get("HAMMOCK_REQUIRE_CPP"),
    reason="hammock-cpp not built (see build/<tag>/hammock-cpp); "
           "set HAMMOCK_REQUIRE_CPP=1 to make this a failure instead")


def test_binary_is_not_stale():
    """A binary older than its source silently validates the wrong code."""
    assert _BIN is not None, "hammock-cpp not built (HAMMOCK_REQUIRE_CPP is set)"
    # Check the whole core, not just the CLI: a stale libhammock_core.a (i.e.
    # an old hll_sketch.cpp) with a freshly-compiled CLI would otherwise pass,
    # and "bit-for-bit" would be asserted against old estimator code.
    srcs = [_REPO / "cpp" / "app" / "hammock_cli.cpp"]
    srcs += list((_REPO / "cpp" / "src").glob("*.cpp"))
    srcs += list((_REPO / "cpp" / "include").rglob("*.hpp"))
    newest = max(srcs, key=lambda p: p.stat().st_mtime)
    # mtime-based, so a `git checkout`/`git merge` that only touches timestamps
    # trips it too. That direction is the safe one -- it fails loud and a
    # rebuild clears it -- but say so, or the message reads as a real drift.
    assert _BIN.stat().st_mtime >= newest.stat().st_mtime, (
        f"{_BIN} is older than {newest} -- rebuild before trusting these "
        f"results (a git checkout also refreshes source mtimes)")


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


def _run_cpp_path(q: Path, r: Path, out: Path, metrics: bool,
                  mode: str = "B", extra: list[str] | None = None) -> Path:
    """Run the binary and return the path it actually wrote.

    Reads that path off the `--verbose` "Wrote <path>" line rather than globbing
    the prefix. Globbing returns an arbitrary match once two shapes share a
    prefix, and reconstructing the name in Python would duplicate the suffix
    rules in outprefix_with_suffix -- which this release changes three ways.
    """
    cmd = [str(_BIN), str(q), str(r), "--mode", mode, "-p", "20", "-o", str(out),
           "--verbose", "--metrics" if metrics else "--no-metrics", *(extra or [])]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
    for line in proc.stderr.splitlines():
        if line.startswith("Wrote "):
            return Path(line[len("Wrote "):].strip())
    raise AssertionError(f"no 'Wrote <path>' line in stderr:\n{proc.stderr}")


def _run_cpp(q: Path, r: Path, out: Path, metrics: bool) -> dict:
    csv_path = _run_cpp_path(q, r, out, metrics)
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


# _METRIC_COLS already leads with jaccard_similarity, so the full header is
# just the two key columns plus the block.
_FULL_HEADER = ["query", "reference"] + _METRIC_COLS
_SLIM_HEADER = ["query", "reference", "jaccard_similarity"]


def test_default_output_carries_the_metrics_block(corpus, tmp_path: Path):
    """The default is the full block, on the untagged path.

    Since 0.7.0 the binary emits what the Python CLI emits. The untagged
    filename is deliberate: it is the path every pre-0.7.0 default run already
    used, so the *reduced* shape is the one that gets a tag.
    """
    q, r = corpus
    path = _run_cpp_path(q, r, tmp_path / "plain", metrics=True)
    assert path.read_text().splitlines()[0].split("\t") == _FULL_HEADER
    assert path.name == "plain_hll_p20_jaccB.csv"


def test_no_metrics_gives_three_columns(corpus, tmp_path: Path):
    """--no-metrics is the opt-out for timing runs, and it tags its output."""
    q, r = corpus
    path = _run_cpp_path(q, r, tmp_path / "slim", metrics=False)
    assert path.read_text().splitlines()[0].split("\t") == _SLIM_HEADER
    assert path.name == "slim_hll_p20_jaccB_j3.csv"


def test_no_metrics_changes_the_output_filename(corpus, tmp_path: Path):
    """A 3-column and a 9-column file must not collide on one path."""
    q, r = corpus
    _run_cpp(q, r, tmp_path / "same", metrics=True)
    _run_cpp(q, r, tmp_path / "same", metrics=False)
    names = sorted(p.name for p in tmp_path.glob("same*.csv"))
    assert len(names) == 2, names
    assert any(n.endswith("_j3.csv") for n in names), names


def test_jaccard_similarity_is_identical_across_shapes(corpus, tmp_path: Path):
    """The frozen column must not depend on which shape you asked for.

    `jaccard_similarity` is register-equality and every archived analysis is
    calibrated against it, so flipping the default must not perturb it. The
    metrics block adds columns to its right; it must not touch column 3.
    """
    q, r = corpus
    full = _run_cpp(q, r, tmp_path / "full", metrics=True)
    slim = _run_cpp(q, r, tmp_path / "slim", metrics=False)
    assert set(full) == set(slim)
    for key in full:
        assert full[key]["jaccard_similarity"] == slim[key]["jaccard_similarity"], key


def test_mode_defaults_to_b(corpus, tmp_path: Path):
    """No --mode means B, matching cli.py's autodetect for BED input.

    Until 0.7.0 the binary defaulted to A while the Python CLI chose B, so the
    same two lists through the two front-ends produced a column of the same name
    holding numbers from different algorithms.
    """
    q, r = corpus
    proc = subprocess.run(
        [str(_BIN), str(q), str(r), "-p", "20", "-o", str(tmp_path / "d"), "--verbose"],
        check=True, capture_output=True, text=True, timeout=600)
    assert "(mode B)" in proc.stderr, proc.stderr
    assert (tmp_path / "d_hll_p20_jaccB.csv").exists()


def test_suba_and_subb_reach_the_filename(corpus, tmp_path: Path):
    """Two runs differing only in --subA must not overwrite each other.

    The C++ filename builder had no subA branch before 0.7.0, and formatted subB
    with an unconditional %.2f, so `--subB 0.001` and `--subB 0.005` also landed
    on one path. Both are silent data loss: the second run wins.
    """
    q, r = corpus
    a = _run_cpp_path(q, r, tmp_path / "s", False, mode="C", extra=["--subA", "0.5"])
    b = _run_cpp_path(q, r, tmp_path / "s", False, mode="C", extra=["--subA", "0.25"])
    assert a != b, (a, b)

    lo = _run_cpp_path(q, r, tmp_path / "t", False, extra=["--subB", "0.001"])
    hi = _run_cpp_path(q, r, tmp_path / "t", False, extra=["--subB", "0.005"])
    assert lo != hi, (lo, hi)
    # The .4f rule is strict below 0.01, so the historical name is preserved.
    keep = _run_cpp_path(q, r, tmp_path / "u", False, extra=["--subB", "0.01"])
    assert keep.name == "u_hll_p20_jaccB_B0.01_j3.csv", keep.name


def test_version_is_on_stdout_and_matches_the_package():
    """A harness probes this to refuse a stale binary, so it must be stdout."""
    proc = subprocess.run([str(_BIN), "--version"],
                          check=True, capture_output=True, text=True, timeout=60)
    assert proc.stdout.strip() == f"hammock-cpp {hammock.__version__}"


def test_peak_height_is_gone(corpus, tmp_path: Path):
    """--peak-height and its BagMinHash backend were removed in 0.7.0.

    It named a count-weighted sketch that was never wired into either CLI, and
    it is now rejected by the generic unknown-argument path rather than by a
    special case -- so this asserts the flag name appears in the error, which
    only the unknown-argument branch produces.
    """
    q, r = corpus
    proc = subprocess.run(
        [str(_BIN), str(q), str(r), "--mode", "A", "-o", str(tmp_path / "x"),
         "--peak-height", "5"],
        capture_output=True, text=True, timeout=600)
    assert proc.returncode != 0
    assert "unknown argument '--peak-height'" in proc.stderr
