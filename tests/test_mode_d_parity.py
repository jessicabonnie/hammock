"""Mode D parity: our `hammock` CLI must match the orig conda-env `hammock`
on FASTA fixtures, modulo the containment/cosketch columns hammock_claude
adds (orig has none).

The orig is installed in the user's `hammock` conda env (Python 3.12 + bioconda
`digest`); the pipx-installed `hammock-orig` runs Python 3.8 where bioconda
doesn't ship `digest`, so this test deliberately bypasses that one.

Skipped if the conda env hammock isn't on disk.
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

pytestmark = pytest.mark.skipif(
    not CONDA_ORIG.exists() or OURS is None,
    reason="Mode D parity requires the conda-env hammock and our hammock both available.",
)


def _run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)


def _files_list(tmp_path: Path, *names: str) -> Path:
    f = tmp_path / "files.txt"
    f.write_text("\n".join(str(DATA / n) for n in names) + "\n")
    return f


# Containment + cosketch columns are hammock_claude additions that orig 0.4.0
# does not emit. Drop them (in both forms — plain and `_with_ends`) before
# comparing CSVs.
_BASES = ("containment_AB", "containment_BA",
          "cosketch_geom", "cosketch_arith", "cosketch_max")
_PROJECTED_OUT = set(_BASES) | {b + "_with_ends" for b in _BASES}


def _projected(csv_text: str) -> list[tuple]:
    lines = csv_text.strip().split("\n")
    header = lines[0].split(",")
    keep = [i for i, name in enumerate(header) if name not in _PROJECTED_OUT]
    return [tuple(line.split(",")[i] for i in keep) for line in lines]


@pytest.mark.parametrize("k,w,p", [
    (8, 40, 14),
    (8, 40, 12),
    (10, 30, 14),
    (5, 20, 14),
])
def test_mode_d_byte_equal(tmp_path: Path, k: int, w: int, p: int) -> None:
    files = _files_list(tmp_path, "tiny.fa", "tiny2.fa")
    common = [str(files), str(files), "--mode", "D",
              "-p", str(p), "-k", str(k), "-w", str(w)]
    _run([str(CONDA_ORIG), *common, "-o", str(tmp_path / "orig")], tmp_path)
    _run([OURS, *common, "-o", str(tmp_path / "ours")], tmp_path)

    orig_csv = next(tmp_path.glob("orig*.csv")).read_text()
    ours_csv = next(tmp_path.glob("ours*.csv")).read_text()
    assert _projected(orig_csv) == _projected(ours_csv), (
        f"Mode D mismatch for k={k} w={w} p={p}\n"
        f"--- orig ---\n{orig_csv}\n"
        f"--- ours ---\n{ours_csv}"
    )
