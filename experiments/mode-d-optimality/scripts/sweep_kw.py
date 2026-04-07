#!/usr/bin/env python3
"""
scripts/sweep_kw.py
===================
Evaluate minimizer sketch accuracy across (k, w) parameter combinations using
synthetic sequence pairs with known true k-mer Jaccard similarity.

Unlike the existing test_minimizer_similarity.py (which uses character-wise
similarity as ground truth), this script uses **exact k-mer set Jaccard** as
the reference — which is precisely what the minimizer sketch is trying to
estimate.

Ground truth for a pair (seq1, seq2) at kmer size k:
    J_true(k) = |kmers(seq1) ∩ kmers(seq2)| / |kmers(seq1) ∪ kmers(seq2)|

Outputs a CSV to results/ with one row per (k, w, mutation_rate, trial).
Use scripts/plot_sweep.py to visualize.

Usage:
    python scripts/sweep_kw.py
    python scripts/sweep_kw.py --config config/config.yaml
    python scripts/sweep_kw.py --out results/sweep.csv --trials 20
"""

import argparse
import csv
import os
import random
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from itertools import product

import yaml

sys.path.insert(0, "/home/jbonnie1/interval_sketch/hammock")
from hammock.lib.minimizer import MinimizerSketch


# ---------------------------------------------------------------------------
# Sequence generation
# ---------------------------------------------------------------------------

def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(seq: str, mutation_rate: float, rng: random.Random) -> str:
    """Return a copy of seq with each position independently mutated."""
    bases = "ACGT"
    out = []
    for b in seq:
        if rng.random() < mutation_rate:
            out.append(rng.choice([x for x in bases if x != b]))
        else:
            out.append(b)
    return "".join(out)


# ---------------------------------------------------------------------------
# Ground truth
# ---------------------------------------------------------------------------

def kmer_set(seq: str, k: int) -> set:
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}


def true_kmer_jaccard(seq1: str, seq2: str, k: int) -> float:
    s1 = kmer_set(seq1, k)
    s2 = kmer_set(seq2, k)
    union = s1 | s2
    if not union:
        return 0.0
    return len(s1 & s2) / len(union)


# ---------------------------------------------------------------------------
# Sketch similarity
# ---------------------------------------------------------------------------

def sketched_similarity(seq1: str, seq2: str, k: int, w: int, precision: int) -> dict:
    sk1 = MinimizerSketch(kmer_size=k, window_size=w, precision=precision)
    sk2 = MinimizerSketch(kmer_size=k, window_size=w, precision=precision)
    sk1.add_string(seq1)
    sk2.add_string(seq2)
    return sk1.similarity_values(sk2)


# ---------------------------------------------------------------------------
# Worker (runs in subprocess — must be importable at module level)
# ---------------------------------------------------------------------------

def _run_task(args):
    k, w, mut, trial, seq_length, precision, seed = args
    rng = random.Random(seed)
    seq1 = random_dna(seq_length, rng)
    seq2 = mutate(seq1, mut, rng)
    true_j = true_kmer_jaccard(seq1, seq2, k)
    try:
        sim = sketched_similarity(seq1, seq2, k, w, precision)
        sj  = sim["hash_similarity"]
        sje = sim["hash_with_ends_similarity"]
    except Exception as e:
        sj = sje = float("nan")
        print(f"  ERROR k={k} w={w} mut={mut} trial={trial}: {e}", file=sys.stderr)
    return (
        k, w, f"{mut:.3f}", trial,
        f"{true_j:.6f}",
        f"{sj:.6f}",
        f"{sje:.6f}",
        f"{sj  - true_j:.6f}" if sj == sj else "nan",
        f"{sje - true_j:.6f}" if sje == sje else "nan",
    )


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_sweep(k_values, w_values, mutation_rates, seq_length, trials,
              precision, seed, out_path, workers=1):
    # Filter degenerate combos
    kw_pairs = [(k, w) for k, w in product(k_values, w_values) if w >= k]
    skipped = len(k_values) * len(w_values) - len(kw_pairs)
    if skipped:
        print(f"Skipping {skipped} (k,w) pair(s) where w < k (degenerate sketches)")

    # Build task list; each task gets a unique seed derived from the master seed
    tasks = []
    master_rng = random.Random(seed)
    for k, w in kw_pairs:
        for mut in mutation_rates:
            for trial in range(trials):
                task_seed = master_rng.randint(0, 2**32 - 1)
                tasks.append((k, w, mut, trial + 1, seq_length, precision, task_seed))

    total = len(tasks)
    print(f"Total tasks: {total}  workers: {workers}")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    rows = []
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(_run_task, t): t for t in tasks}
            done = 0
            for fut in as_completed(futs):
                rows.append(fut.result())
                done += 1
                if done % 100 == 0:
                    print(f"  {done}/{total} done", flush=True)
    else:
        for i, t in enumerate(tasks):
            rows.append(_run_task(t))
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{total} done", flush=True)

    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "k", "w", "mutation_rate", "trial",
            "true_kmer_jaccard",
            "sketched_hash",
            "sketched_hash_with_ends",
            "error_hash",
            "error_hash_with_ends",
        ])
        for row in rows:
            writer.writerow(row)

    print(f"Results written to {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def main():
    p = argparse.ArgumentParser(description="k/w parameter sweep for minimizer Mode D accuracy")
    p.add_argument("--config", default="config/config.yaml",
                   help="YAML config file (default: config/config.yaml)")
    p.add_argument("--out", default=None,
                   help="Output CSV path (overrides config)")
    p.add_argument("--trials", type=int, default=None,
                   help="Trials per (k, w, mutation_rate) combination (overrides config)")
    p.add_argument("--workers", type=int, default=None,
                   help="Parallel workers (default: $SLURM_CPUS_PER_TASK if set, else 1)")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    cfg = load_config(args.config)

    k_values       = cfg["kmer_sizes"]
    w_values       = cfg["window_sizes"]
    mutation_rates = cfg["mutation_rates"]
    seq_length     = cfg.get("seq_length", 1000)
    trials         = args.trials or cfg.get("trials", 10)
    precision      = cfg.get("precision", 24)
    results_dir    = cfg.get("results_dir", "results")

    # Workers: CLI > SLURM env > 1
    if args.workers is not None:
        workers = args.workers
    elif os.environ.get("SLURM_CPUS_PER_TASK"):
        workers = int(os.environ["SLURM_CPUS_PER_TASK"])
    else:
        workers = 1

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = args.out or os.path.join(results_dir, f"sweep_{timestamp}.csv")

    print(f"k values:       {k_values}")
    print(f"w values:       {w_values}")
    print(f"mutation rates: {mutation_rates}")
    print(f"seq length:     {seq_length}")
    print(f"trials:         {trials}")
    print(f"precision:      {precision}")
    print(f"workers:        {workers}")
    print(f"output:         {out_path}")

    run_sweep(k_values, w_values, mutation_rates, seq_length, trials,
              precision, args.seed, out_path, workers=workers)


if __name__ == "__main__":
    main()
