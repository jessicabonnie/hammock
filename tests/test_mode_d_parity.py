"""Mode D *structural* parity against the orig conda-env `hammock`.

We deliberately do NOT assert byte-equal numeric similarity here: divergence #6
in CLAUDE.md means our Mode D adds each minimizer's raw 64-bit hash via
`add_hash64`, whereas orig's conda env falls through to the Python slow path
that hashes the *decimal digits* of each minimizer hash as k-mers. The
minimizer sets are identical (same `digest`), but the HLL ingestion differs by
design, so `jaccard_similarity` diverges (e.g. 0.75 vs 0.7903 on tiny.fa). An
earlier byte-equal version of this test only ever "passed" because
`shutil.which("hammock")` resolved to the orig binary (conda `hammock` ahead of
our env on PATH), silently comparing orig-to-orig.

What we CAN and do assert:
  * structural columns match orig exactly (file labels, mode, params, order);
  * our similarity columns are well-formed: in [0,1], self-vs-self == 1.0,
    and symmetric.

The orig is installed in the user's `hammock` conda env (Python 3.12 + bioconda
`digest`); the pipx-installed `hammock-orig` runs Python 3.8 where bioconda
doesn't ship `digest`, so this test deliberately bypasses that one.

Skipped if the conda env hammock isn't on disk, or if our CLI on PATH IS the
orig binary (nothing to compare).
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).parent
DATA = HERE / "data"

CONDA_ORIG = Path("/home/jbonnie1/.conda/envs/hammock/bin/hammock")
OURS = shutil.which("hammock")

# If our CLI on PATH *is* the orig binary, there's nothing to compare — skip
# rather than trivially self-compare (the bug that masked this test before).
_OURS_IS_ORIG = OURS is not None and Path(OURS).resolve() == CONDA_ORIG.resolve()

pytestmark = pytest.mark.skipif(
    not CONDA_ORIG.exists() or OURS is None or _OURS_IS_ORIG,
    reason="Mode D parity needs the conda-env hammock AND a distinct refactor "
           "hammock on PATH (put the refactor env's bin first).",
)


def _run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def _files_list(tmp_path: Path, *names: str) -> Path:
    f = tmp_path / "files.txt"
    f.write_text("\n".join(str(DATA / n) for n in names) + "\n")
    return f


# Structural columns both orig and the refactor emit identically (labels,
# mode, params). The similarity columns diverge by design (divergence #6), so
# they are NOT part of the structural comparison.
_STRUCTURAL = ("file1", "file2", "sketch_type", "mode",
               "precision", "num_hashes", "kmer_size", "window_size")


def _rows(csv_text: str) -> tuple[list[str], list[list[str]]]:
    lines = csv_text.strip().split("\n")
    return lines[0].split(","), [ln.split(",") for ln in lines[1:]]


def _structural(csv_text: str) -> list[tuple]:
    header, rows = _rows(csv_text)
    idx = [header.index(c) for c in _STRUCTURAL]
    return [tuple(r[i] for i in idx) for r in rows]


@pytest.mark.parametrize("k,w,p", [
    (8, 40, 14),
    (8, 40, 12),
    (10, 30, 14),
    (5, 20, 14),
])
def test_mode_d_structural_parity(tmp_path: Path, k: int, w: int, p: int) -> None:
    files = _files_list(tmp_path, "tiny.fa", "tiny2.fa")
    common = [str(files), str(files), "--mode", "D",
              "-p", str(p), "-k", str(k), "-w", str(w)]
    _run([str(CONDA_ORIG), *common, "-o", str(tmp_path / "orig")], tmp_path)
    _run([OURS, *common, "-o", str(tmp_path / "ours")], tmp_path)

    orig_csv = next(tmp_path.glob("orig*.csv")).read_text()
    ours_csv = next(tmp_path.glob("ours*.csv")).read_text()

    # (a) Structural columns must match orig exactly — same labels, params, order.
    assert _structural(orig_csv) == _structural(ours_csv), (
        f"Mode D structural mismatch for k={k} w={w} p={p}\n"
        f"--- orig ---\n{orig_csv}\n--- ours ---\n{ours_csv}"
    )

    # (b) Our similarity columns must be well-formed (values diverge from orig
    #     by design — divergence #6 — so we sanity-check ours, not equality).
    header, rows = _rows(ours_csv)
    sim_cols = [i for i, c in enumerate(header)
                if c.startswith(("jaccard_similarity", "containment", "cosketch"))]
    assert sim_cols, "no similarity columns emitted"
    # containment_AB/BA are *antisymmetric* -- C_AB(a,b) == C_BA(b,a), not
    # C_AB(b,a) -- so they belong in the range and self-pair checks but must be
    # excluded from the symmetry check below.
    sym_cols = [i for i in sim_cols
                if header[i].startswith(("jaccard_similarity", "cosketch"))]
    by_pair = {}
    for r in rows:
        f1, f2 = r[header.index("file1")], r[header.index("file2")]
        for i in sim_cols:
            v = float(r[i])
            # containment_* is emitted unclamped and can exceed 1 by ~1 ulp
            # (inter is a difference of three Ertl estimates), so allow slack.
            assert -1e-9 <= v <= 1.0 + 1e-9, \
                f"{header[i]}={v} out of [0,1] ({f1},{f2})"
            if f1 == f2:
                assert v == pytest.approx(1.0), \
                    f"self-pair {header[i]}={v} != 1.0 for {f1}"
        by_pair[(f1, f2)] = {i: float(r[i]) for i in sym_cols}
    # (c) Symmetry: J(A,B) == J(B,A). Exact for the Jaccard columns (double
    #     addition is commutative and the union is register-wise max), so
    #     approx here is generous.
    for (a, b), vals in by_pair.items():
        if (b, a) in by_pair:
            for i, v in vals.items():
                assert v == pytest.approx(by_pair[(b, a)][i]), \
                    f"asymmetric {header[i]} for {a},{b}"
