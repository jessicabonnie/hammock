"""Mode D process dispatch must work under spawn and preserve exact HLL state."""
from __future__ import annotations

import gzip
import multiprocessing as mp
from pathlib import Path
from types import SimpleNamespace

import pytest

from hammock import _core, cli, runner


DATA = Path(__file__).parent / "data"


def _block_forever() -> None:
    """Picklable spawn target for shutdown escalation coverage."""
    import time

    time.sleep(60)


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


def test_cli_process_csv_is_byte_identical_and_atomic(tmp_path: Path) -> None:
    paths = tmp_path / "paths.txt"
    paths.write_text(f"{DATA / 'tiny.fa'}\n{DATA / 'tiny2.fa'}\n")
    common = [str(paths), str(paths), "--mode", "D", "-p", "12", "-k", "8", "-w", "40"]

    assert cli.main([*common, "--threads", "1", "-o", str(tmp_path / "serial")]) == 0
    assert cli.main([*common, "--threads", "2", "-o", str(tmp_path / "process")]) == 0

    serial = tmp_path / "serial_mnmzr_p12_jaccD_k8_w40_ie.csv"
    process = tmp_path / "process_mnmzr_p12_jaccD_k8_w40_ie.csv"
    assert serial.read_bytes() == process.read_bytes()
    assert not list(tmp_path.glob(".*.tmp"))


def test_bounded_shutdown_terminates_a_blocked_spawn_worker() -> None:
    context = mp.get_context("spawn")
    tasks = context.Queue(maxsize=1)
    results = context.Queue(maxsize=1)
    worker = context.Process(target=_block_forever)
    worker.start()

    runner._stop_mode_d_workers([worker], tasks, results, graceful=False)

    assert not worker.is_alive()
    assert worker.exitcode is not None
