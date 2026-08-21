"""Mode D process dispatch must work under spawn and preserve exact HLL state."""
from __future__ import annotations

import gzip
import multiprocessing as mp
from pathlib import Path
from types import SimpleNamespace

import pytest

from hammock import _core, cli, runner


DATA = Path(__file__).parent / "data"


def _args(threads: int) -> SimpleNamespace:
    return SimpleNamespace(
        mode="D",
        threads=threads,
        kmer_size=8,
        window_size=40,
        seed=42,
        precision=12,
        sketch_type="minimizer",
        verbose=False,
    )


def test_mode_d_spawn_processes_match_serial_and_keep_order() -> None:
    # Repetition makes completion order potentially differ from input order.
    paths = [str(DATA / "tiny2.fa"), str(DATA / "tiny.fa"), str(DATA / "tiny2.fa")]
    serial = runner._sketch_many(paths, _args(1), "query")
    spawned = runner._sketch_many(paths, _args(2), "query")

    assert [
        _core._hll_transport_state(sketch.minimizer_hll)
        for sketch in spawned
    ] == [
        _core._hll_transport_state(sketch.minimizer_hll)
        for sketch in serial
    ]


def test_mode_d_worker_failure_identifies_path_and_leaves_no_child(tmp_path: Path) -> None:
    missing = str(tmp_path / "missing.fa")
    with pytest.raises(runner.ModeDSketchWorkerError, match="missing.fa"):
        runner._sketch_many([missing, str(DATA / "tiny.fa")], _args(2), "query")

    assert not mp.active_children()


def test_cli_worker_failure_is_concise_and_publishes_no_csv(tmp_path: Path, capsys) -> None:
    missing = tmp_path / "missing.fa"
    paths = tmp_path / "paths.txt"
    paths.write_text(f"{missing}\n{DATA / 'tiny.fa'}\n")

    rc = cli.main([
        str(paths), str(paths), "--mode", "D", "--threads", "2",
        "-p", "12", "-k", "8", "-w", "40", "-o", str(tmp_path / "out"),
    ])

    captured = capsys.readouterr()
    assert rc == 1
    assert "missing.fa" in captured.err
    assert "Traceback" not in captured.err
    assert not list(tmp_path.glob("out*.csv"))
    assert not list(tmp_path.glob(".*.tmp"))


def test_mode_d_spawn_processes_support_gzipped_fasta(tmp_path: Path) -> None:
    gz_path = tmp_path / "tiny.fa.gz"
    with open(DATA / "tiny.fa", "rb") as source, gzip.open(gz_path, "wb") as target:
        target.write(source.read())

    serial = runner._sketch_many([str(gz_path)], _args(1), "query")
    # Two entries exercise both spawn dispatch and duplicate-path ordering.
    spawned = runner._sketch_many([str(gz_path), str(gz_path)], _args(2), "query")
    expected = _core._hll_transport_state(serial[0].minimizer_hll)
    assert [_core._hll_transport_state(s.minimizer_hll) for s in spawned] == [expected, expected]
