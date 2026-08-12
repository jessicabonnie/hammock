"""The standalone `hammock-cpp` must match the Python CLI exactly, in all
three output shapes (ie/re/full) -- docs/seed-metrics-column-restructure.md.

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

# Column lists per shape, matching runner.py's _metrics_shape exactly (order
# matters -- these double as the expected trailing-header slice).
_IE_COLS = ["jaccard_similarity_ie"]
_RE_COLS = ["jaccard_similarity", "register_equality_similarity"]
_FULL_COLS = [
    "jaccard_similarity", "jaccard_similarity_ie",
    "containment_AB", "containment_BA",
    "cosketch_geom", "cosketch_arith", "cosketch_max",
    "register_equality_similarity",
]
_IE_HEADER = ["query", "reference"] + _IE_COLS
_RE_HEADER = ["query", "reference"] + _RE_COLS
_FULL_HEADER = ["query", "reference"] + _FULL_COLS

# (shape name, CLI flags, filename tag, header)
_SHAPES = [
    ("ie", [], "ie", _IE_HEADER),
    ("re", ["--register-equality"], "re", _RE_HEADER),
    ("full", ["--metrics"], "full", _FULL_HEADER),
]
_SHAPE_FLAGS = {name: flags for name, flags, _tag, _hdr in _SHAPES}


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


def _run_python(q: Path, r: Path, out: Path, shape: str) -> dict:
    cmd = [sys.executable, "-m", "hammock.cli", str(q), str(r),
           "--mode", "B", "-p", "20", "-o", str(out), *_SHAPE_FLAGS[shape]]
    subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
    csv_path = next(out.parent.glob(f"{out.name}*.csv"))
    with csv_path.open() as f:
        return {(row["file1"], row["file2"]): row for row in csv.DictReader(f)}


def _run_cpp_path(q: Path, r: Path, out: Path, shape: str,
                  mode: str = "B", extra: list[str] | None = None) -> Path:
    """Run the binary and return the path it actually wrote.

    Reads that path off the `--verbose` "Wrote <path>" line rather than globbing
    the prefix. Globbing returns an arbitrary match once two shapes share a
    prefix, and reconstructing the name in Python would duplicate the suffix
    rules in outprefix_with_suffix -- which this release changes three ways.
    """
    cmd = [str(_BIN), str(q), str(r), "--mode", mode, "-p", "20", "-o", str(out),
           "--verbose", *_SHAPE_FLAGS[shape], *(extra or [])]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
    for line in proc.stderr.splitlines():
        if line.startswith("Wrote "):
            return Path(line[len("Wrote "):].strip())
    raise AssertionError(f"no 'Wrote <path>' line in stderr:\n{proc.stderr}")


def _run_cpp(q: Path, r: Path, out: Path, shape: str) -> dict:
    csv_path = _run_cpp_path(q, r, out, shape)
    with csv_path.open() as f:
        return {(row["query"], row["reference"]): row
                for row in csv.DictReader(f, delimiter="\t")}


@pytest.mark.parametrize("shape,cols", [("ie", _IE_COLS), ("re", _RE_COLS), ("full", _FULL_COLS)])
def test_shape_matches_python_bit_for_bit(corpus, tmp_path: Path, shape, cols):
    """Cross-tool bit-for-bit gate, for every shape, not just full.

    Exact equality, not approx: same IEEE operations in the same order on
    both sides. If this starts failing, the two derivations have diverged --
    fix that rather than loosening the assertion.
    """
    q, r = corpus
    py = _run_python(q, r, tmp_path / f"py_{shape}", shape)
    cpp = _run_cpp(q, r, tmp_path / f"cpp_{shape}", shape)

    assert set(py) == set(cpp), (sorted(py), sorted(cpp))
    for key in py:
        for col in cols:
            assert float(py[key][col]) == float(cpp[key][col]), (shape, key, col)


@pytest.mark.parametrize("shape,flags,tag,header", _SHAPES)
def test_shape_header_and_filename(corpus, tmp_path: Path, shape, flags, tag, header):
    """Every shape is tagged now -- none stays bare (pre-restructure, only the
    reduced shape was tagged '_j3' and the full block was the untagged
    default)."""
    q, r = corpus
    path = _run_cpp_path(q, r, tmp_path / shape, shape)
    assert path.read_text().splitlines()[0].split("\t") == header
    assert path.name == f"{shape}_hll_p20_jaccB_{tag}.csv"


def test_re_is_an_alias_for_register_equality(corpus, tmp_path: Path):
    """`_run_cpp_path`/`_SHAPE_FLAGS` always spell it `--register-equality`;
    this test is the one place `--re` itself gets exercised."""
    q, r = corpus
    cmd = [str(_BIN), str(q), str(r), "--mode", "B", "-p", "20",
           "-o", str(tmp_path / "aliased"), "--verbose", "--re"]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
    written = next(line for line in proc.stderr.splitlines() if line.startswith("Wrote "))
    written_path = Path(written[len("Wrote "):].strip())
    assert written_path.name.endswith("_re.csv")
    header = written_path.read_text().splitlines()[0].split("\t")
    assert header == _RE_HEADER


def test_shapes_produce_distinct_filenames(corpus, tmp_path: Path):
    """Three shapes at the same prefix must not collide on one path."""
    q, r = corpus
    for shape in ("ie", "re", "full"):
        _run_cpp(q, r, tmp_path / "same", shape)
    names = sorted(p.name for p in tmp_path.glob("same*.csv"))
    assert len(names) == 3, names
    assert {n.rsplit("_", 1)[-1] for n in names} == {"ie.csv", "re.csv", "full.csv"}, names


def test_jaccard_similarity_is_identical_across_re_and_full_shapes(corpus, tmp_path: Path):
    """The frozen column must not depend on which shape you asked for.

    `jaccard_similarity` is register-equality and every archived analysis is
    calibrated against it, so flipping shapes must not perturb it. "ie" has
    no jaccard_similarity column at all, so the comparable pair is re vs full.
    """
    q, r = corpus
    full = _run_cpp(q, r, tmp_path / "full", "full")
    re_ = _run_cpp(q, r, tmp_path / "slim", "re")
    assert set(full) == set(re_)
    for key in full:
        assert full[key]["jaccard_similarity"] == re_[key]["jaccard_similarity"], key


@pytest.mark.parametrize("shape", ["re", "full"])
def test_register_equality_similarity_duplicates_jaccard_similarity(corpus, tmp_path: Path, shape):
    q, r = corpus
    rows = _run_cpp(q, r, tmp_path / f"dup_{shape}", shape)
    for key, row in rows.items():
        assert row["jaccard_similarity"] == row["register_equality_similarity"], (shape, key, row)


def test_default_matches_metrics_jaccard_similarity_ie(corpus, tmp_path: Path):
    """jaccard_similarity_ie must be byte-identical between the bare default
    and --metrics."""
    q, r = corpus
    ie = _run_cpp(q, r, tmp_path / "ie", "ie")
    full = _run_cpp(q, r, tmp_path / "full", "full")
    assert set(ie) == set(full)
    for key in ie:
        assert ie[key]["jaccard_similarity_ie"] == full[key]["jaccard_similarity_ie"], key


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
    # No explicit shape flag -> the default MetricsMode::Ie, tagged '_ie'.
    # (Pre-restructure this asserted the untagged 'd_hll_p20_jaccB.csv'; that
    # path no longer exists -- every shape is tagged now, including this one.)
    assert (tmp_path / "d_hll_p20_jaccB_ie.csv").exists()


def test_suba_and_subb_reach_the_filename(corpus, tmp_path: Path):
    """Two runs differing only in --subA must not overwrite each other.

    The C++ filename builder had no subA branch before 0.7.0, and formatted subB
    with an unconditional %.2f, so `--subB 0.001` and `--subB 0.005` also landed
    on one path. Both are silent data loss: the second run wins.
    """
    q, r = corpus
    a = _run_cpp_path(q, r, tmp_path / "s", "re", mode="C", extra=["--subA", "0.5"])
    b = _run_cpp_path(q, r, tmp_path / "s", "re", mode="C", extra=["--subA", "0.25"])
    assert a != b, (a, b)

    lo = _run_cpp_path(q, r, tmp_path / "t", "re", extra=["--subB", "0.001"])
    hi = _run_cpp_path(q, r, tmp_path / "t", "re", extra=["--subB", "0.005"])
    assert lo != hi, (lo, hi)
    # The .4f rule is strict below 0.01, so the historical name is preserved
    # (only the trailing shape tag moved from '_j3' to '_re').
    keep = _run_cpp_path(q, r, tmp_path / "u", "re", extra=["--subB", "0.01"])
    assert keep.name == "u_hll_p20_jaccB_B0.01_re.csv", keep.name


def test_metrics_and_register_equality_are_mutually_exclusive(corpus, tmp_path: Path):
    q, r = corpus
    proc = subprocess.run(
        [str(_BIN), str(q), str(r), "--mode", "B", "-p", "20",
         "-o", str(tmp_path / "x"), "--metrics", "--register-equality"],
        capture_output=True, text=True, timeout=600)
    assert proc.returncode != 0
    assert "mutually exclusive" in proc.stderr, proc.stderr


def test_no_metrics_is_gone(corpus, tmp_path: Path):
    """--no-metrics was removed outright (not aliased) -- same pattern as the
    --peak-height removal (test_peak_height_is_gone, below): it now falls
    through to the generic unknown-argument error."""
    q, r = corpus
    proc = subprocess.run(
        [str(_BIN), str(q), str(r), "--mode", "B", "-p", "20",
         "-o", str(tmp_path / "y"), "--no-metrics"],
        capture_output=True, text=True, timeout=600)
    assert proc.returncode != 0
    assert "unknown argument '--no-metrics'" in proc.stderr, proc.stderr


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
