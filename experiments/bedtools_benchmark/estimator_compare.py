#!/usr/bin/env python3
"""Compare hammock's two available Jaccard estimators against bedtools ground truth.

  J_re  = reported `jaccard_similarity` column  -> register-equality
  J_ie  = reconstructed from the containment columns via inclusion-exclusion:
            C_AB = I/|A|, C_BA = I/|B|
            J    = I/(|A|+|B|-I) = 1 / (1/C_AB + 1/C_BA - 1)

Ground truth is `bedtools jaccard` (exact bp set Jaccard).

Data is synthesized to span the whole J range (random BED pairs alone only
probe J ~ 0.01), by building each B file as a mixture of intervals copied
from its paired A file plus fresh random intervals.
"""
import argparse
import concurrent.futures
import csv
import os
import random
import subprocess
import sys

BEDTOOLS = "/data/apps/extern/BEDTools/2.30.0/bin/bedtools"
HAMMOCK = "/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock"

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
GENOME_LEN = 10_000_000

# Fraction of A's intervals reused in B. Spans J ~ 0 .. 1.
FRACS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 0.95, 1.0]


def rand_intervals(n, rng):
    out = []
    for _ in range(n):
        chrom = rng.choice(CHROMS)
        start = rng.randint(0, GENOME_LEN)
        end = rng.randint(start + 100, start + 10_000)
        out.append((chrom, start, end))
    return out


def write_bed(intervals, path):
    """Write sorted BED (bedtools jaccard requires sorted input)."""
    rows = sorted(intervals, key=lambda r: (r[0], r[1]))
    with open(path, "w") as f:
        for chrom, start, end in rows:
            f.write(f"{chrom}\t{start}\t{end}\n")


def make_data(tmp_dir, num_reps, num_intervals, seed):
    a_paths, b_paths = [], []
    for rep in range(num_reps):
        rng = random.Random(seed + rep)
        a_iv = rand_intervals(num_intervals, rng)
        a_path = os.path.join(tmp_dir, f"A_rep{rep}.bed")
        write_bed(a_iv, a_path)
        a_paths.append(a_path)
        for frac in FRACS:
            n_shared = int(round(frac * num_intervals))
            shared = rng.sample(a_iv, n_shared)
            novel = rand_intervals(num_intervals - n_shared, rng)
            b_path = os.path.join(tmp_dir, f"B_rep{rep}_f{frac:g}.bed")
            write_bed(shared + novel, b_path)
            b_paths.append(b_path)
    return a_paths, b_paths


def bedtools_jaccard(pair):
    a, b = pair
    res = subprocess.run([BEDTOOLS, "jaccard", "-a", a, "-b", b],
                         capture_output=True, text=True, check=True)
    lines = res.stdout.strip().split("\n")
    # header: intersection  union-intersection  jaccard  n_intersections
    fields = lines[1].split("\t")
    return (os.path.basename(a), os.path.basename(b)), float(fields[2])


def run_bedtools_all(a_paths, b_paths, workers):
    pairs = [(a, b) for a in a_paths for b in b_paths]
    truth = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        for key, j in ex.map(bedtools_jaccard, pairs):
            truth[key] = j
    return truth


def run_hammock(a_paths, b_paths, precision, tmp_dir, threads):
    l1 = os.path.join(tmp_dir, "list1.txt")
    l2 = os.path.join(tmp_dir, "list2.txt")
    with open(l1, "w") as f:
        f.write("\n".join(a_paths) + "\n")
    with open(l2, "w") as f:
        f.write("\n".join(b_paths) + "\n")
    prefix = os.path.join(tmp_dir, f"hm_p{precision}")
    cmd = [HAMMOCK, l1, l2, "--mode", "B", "--precision", str(precision),
           "-o", prefix, "--threads", str(threads)]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    # hammock decorates the prefix; find the CSV it actually wrote
    base = os.path.basename(prefix)
    cands = [f for f in os.listdir(tmp_dir)
             if f.startswith(base) and f.endswith(".csv")]
    if not cands:
        raise RuntimeError(f"no hammock CSV for p={precision} (prefix {base})")
    out = {}
    with open(os.path.join(tmp_dir, sorted(cands)[0])) as f:
        for row in csv.DictReader(f):
            out[(row["file1"], row["file2"])] = row
    return out


def ie_jaccard(c_ab, c_ba):
    """Inclusion-exclusion Jaccard from the two directional containments.

    J = I/(|A|+|B|-I); dividing through by I gives 1/(1/C_AB + 1/C_BA - 1).

    Returns (value, clamped). A zero containment is NOT "undefined" — it means
    `intersection_size` hit its `max(0.0, ...)` clamp (hll_sketch.cpp:183), so
    the estimator's actual output is J = 0 and that is what we score. Dropping
    these rows instead (an earlier version did) silently compares the two
    estimators on different row sets: the clamp fires only at low J, so it
    excises exactly the rows where inclusion-exclusion does worst, which
    inflates every rank statistic computed against it.
    """
    if c_ab <= 0.0 or c_ba <= 0.0:
        return 0.0, True
    denom = (1.0 / c_ab) + (1.0 / c_ba) - 1.0
    if denom <= 0.0:
        return 0.0, True
    return 1.0 / denom, False


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / (sxx * syy) ** 0.5


def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(ranks(xs), ranks(ys))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmp-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--intervals", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--precisions", default="12,16,20,24")
    args = ap.parse_args()

    precisions = [int(p) for p in args.precisions.split(",")]
    os.makedirs(args.tmp_dir, exist_ok=True)

    print("generating data...", file=sys.stderr)
    a_paths, b_paths = make_data(args.tmp_dir, args.reps, args.intervals, args.seed)
    print(f"  {len(a_paths)} A files x {len(b_paths)} B files "
          f"= {len(a_paths) * len(b_paths)} pairs", file=sys.stderr)

    print("bedtools ground truth...", file=sys.stderr)
    truth = run_bedtools_all(a_paths, b_paths, args.threads)

    rows = []
    for p in precisions:
        print(f"hammock p={p}...", file=sys.stderr)
        hm = run_hammock(a_paths, b_paths, p, args.tmp_dir, args.threads)
        for key, j_bt in truth.items():
            r = hm.get(key)
            if r is None:
                continue
            c_ab = float(r["containment_AB"])
            c_ba = float(r["containment_BA"])
            j_re = float(r["jaccard_similarity"])
            j_ie, clamped = ie_jaccard(c_ab, c_ba)
            rows.append({
                "precision": p,
                "file1": key[0],
                "file2": key[1],
                "bedtools_jaccard": j_bt,
                "j_register_equality": j_re,
                "containment_AB": c_ab,
                "containment_BA": c_ba,
                # cardinalities: the bias of the register-equality estimator is
                # governed by the load factor n/m and by |A|/|B|, so results are
                # uninterpretable without them (see RESULTS / the gap doc).
                "card_A": float(r.get("sketch_size_A", "nan") or "nan"),
                "card_B": float(r.get("sketch_size_B", "nan") or "nan"),
                "j_incl_excl": j_ie,
                "ie_clamped": int(clamped),
                "err_re": j_re - j_bt,
                "err_ie": j_ie - j_bt,
            })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out} ({len(rows)} rows)", file=sys.stderr)

    # ---- summary ----
    print()
    print("NOTE: MAE scores calibration only. The register-equality estimator "
          "is biased but\nlower-variance; per-bin error sd is the resolution "
          "comparison, and it does not\nagree with MAE. See "
          "docs/jaccard-definitional-gap.md before quoting these.\n")
    hdr = (f"{'p':>3} {'n':>4} {'MAE re':>9} {'MAE ie':>9} {'maxE re':>9} "
           f"{'maxE ie':>9} {'r re':>7} {'r ie':>7} {'rho re':>7} {'rho ie':>7} {'clamp':>6}")
    print(hdr)
    print("-" * len(hdr))
    # Every statistic below is computed on the SAME row set for both
    # estimators — clamped inclusion-exclusion rows score as J = 0, they are
    # not dropped.
    for p in precisions:
        sub = [r for r in rows if r["precision"] == p]
        bt = [r["bedtools_jaccard"] for r in sub]
        re_ = [r["j_register_equality"] for r in sub]
        ie_ = [r["j_incl_excl"] for r in sub]
        n_clamped = sum(r["ie_clamped"] for r in sub)
        mae_re = sum(abs(r["err_re"]) for r in sub) / len(sub)
        max_re = max(abs(r["err_re"]) for r in sub)
        mae_ie = sum(abs(r["err_ie"]) for r in sub) / len(sub)
        max_ie = max(abs(r["err_ie"]) for r in sub)
        print(f"{p:>3} {len(sub):>4} {mae_re:>9.5f} {mae_ie:>9.5f} "
              f"{max_re:>9.5f} {max_ie:>9.5f} "
              f"{pearson(bt, re_):>7.4f} {pearson(bt, ie_):>7.4f} "
              f"{spearman(bt, re_):>7.4f} {spearman(bt, ie_):>7.4f} {n_clamped:>6}")

    # ---- error vs true-J regime, at max precision ----
    pmax = max(precisions)
    print(f"\nerror by true-J bin (p={pmax}):")
    bins = [(0.0, 0.01), (0.01, 0.05), (0.05, 0.15), (0.15, 0.35),
            (0.35, 0.6), (0.6, 0.85), (0.85, 1.01)]
    print(f"{'J bin':>14} {'n':>4} {'mean J_bt':>10} {'mean re':>9} {'mean ie':>9} "
          f"{'MAE re':>9} {'MAE ie':>9}")
    for lo, hi in bins:
        sub = [r for r in rows
               if r["precision"] == pmax and lo <= r["bedtools_jaccard"] < hi]
        if not sub:
            continue
        print(f"{f'[{lo:g},{hi:g})':>14} {len(sub):>4} "
              f"{sum(r['bedtools_jaccard'] for r in sub) / len(sub):>10.5f} "
              f"{sum(r['j_register_equality'] for r in sub) / len(sub):>9.5f} "
              f"{sum(r['j_incl_excl'] for r in sub) / len(sub):>9.5f} "
              f"{sum(abs(r['err_re']) for r in sub) / len(sub):>9.5f} "
              f"{sum(abs(r['err_ie']) for r in sub) / len(sub):>9.5f}")

    # Resolution, not calibration: register-equality is a near-affine transform
    # of true J, so to compare its noise against inclusion-exclusion's the error
    # must first be rescaled into J units by the fitted slope (1 - c).
    def sd(vals):
        mu = sum(vals) / len(vals)
        return (sum((v - mu) ** 2 for v in vals) / len(vals)) ** 0.5

    print("\nresolution in J units (register-equality error sd rescaled by the")
    print("fitted slope 1-c, so both columns are in true-J units):")
    print(f"{'p':>3} {'c':>7} {'J bin':>14} {'n':>4} {'sd re/(1-c)':>12} "
          f"{'sd ie':>10} {'winner':>8}")
    for p in precisions:
        fit = [r for r in rows if r["precision"] == p and r["bedtools_jaccard"] < 0.9]
        c_hat = sum((r["j_register_equality"] - r["bedtools_jaccard"])
                    / (1 - r["bedtools_jaccard"]) for r in fit) / len(fit)
        for lo, hi in bins:
            sub = [r for r in rows
                   if r["precision"] == p and lo <= r["bedtools_jaccard"] < hi]
            if len(sub) < 3:
                continue
            sd_re = sd([r["err_re"] for r in sub]) / (1 - c_hat)
            sd_ie = sd([r["err_ie"] for r in sub])
            print(f"{p:>3} {c_hat:>7.4f} {f'[{lo:g},{hi:g})':>14} {len(sub):>4} "
                  f"{sd_re:>12.5f} {sd_ie:>10.5f} "
                  f"{('re' if sd_re < sd_ie else 'ie'):>8}")


if __name__ == "__main__":
    main()
