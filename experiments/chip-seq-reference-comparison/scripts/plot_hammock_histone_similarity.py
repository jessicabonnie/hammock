#!/usr/bin/env python3
"""
Join hammock pairwise similarity CSV to ``histone_fastq_manifest.tsv`` and plot summaries.

Typical inputs
--------------
* Hammock output CSV (from ``hammock ... --full-paths``): columns ``file1``, ``file2``, sketch
  metadata, and **``jaccard_similarity_with_ends``** (required unless ``--similarity-col``).
* Manifest: ``histone_fastq_manifest.tsv`` from ``build_histone_fastq_manifest.py``.

Optional: nf-core/chipseq samplesheet CSV (``nfcore_chipseq_samplesheet_Homo_sapiens.csv``) with
columns ``sample`` and ``encode_ip_file_accession`` to resolve peak/FASTA paths that contain the
nf-core sample id (``ENCSR..._H3K4me3``) but not the IP ``ENCFF`` accession.

Example (human **broadPeak** hammock CSV; names include ``p``, ``k``, ``w`` from your run)::

    python3 plot_hammock_histone_similarity.py \\
        homo_BROAD_mnmzr_p24_jaccD_k5_w40.csv \\
        --manifest ${HISTONE_MARKS}/histone_fastq_manifest.tsv \\
        --samplesheet ${HISTONE_MARKS}/nfcore_chipseq_samplesheet_Homo_sapiens.csv \\
        --organism \"Homo sapiens\" \\
        --path-filter fastas_from_broad \\
        --biological-replicate all \\
        --outdir figs-broad
"""

import argparse
from collections import Counter
import re
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    import seaborn as sns
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform
except ImportError as e:  # pragma: no cover
    print(
        "This script requires matplotlib, seaborn, and scipy: pip install matplotlib seaborn scipy",
        file=sys.stderr,
    )
    raise SystemExit(2) from e

ENCFF_RE = re.compile(r"\b(ENCFF[A-Za-z0-9]+)\b")
# nf-core peak / fasta basename: ENCSRxxx_H3K27ac_REP1_peaks.fa or ENCSRxxx_H3K27ac_peaks.fa
NF_PEAK_BASENAME_RE = re.compile(
    r"^(?P<sample>ENCSR[A-Za-z0-9]+_[A-Za-z0-9][A-Za-z0-9+-]*?)(?:_REP\d+)?_peaks\.fa(?:\.gz)?$",
    re.IGNORECASE,
)


def normalize_nfcore_sample_stem(stem: str) -> str:
    """Map path stem to nf-core ``sample`` column (strip ``_REP1`` before ``_peaks``)."""
    s = (stem or "").strip()
    return re.sub(r"_REP\d+$", "", s, flags=re.IGNORECASE)


def nf_sample_from_basename(path: str) -> str:
    """Return nf-core samplesheet ``sample`` id from peaks/fasta filename, or ``\"\"``."""
    name = Path(path.replace("\\", "/")).name
    m = NF_PEAK_BASENAME_RE.match(name)
    if not m:
        return ""
    return normalize_nfcore_sample_stem(m.group("sample"))


def _norm_histone(s: str) -> str:
    return (s or "").strip().replace(" ", "").replace("-", "").lower()


def heatmap_vmin_data_vmax_one(arr: np.ndarray, ignore_diagonal: bool = False) -> Tuple[float, float]:
    """
    Colormap: **vmin** = data minimum (optionally off-diagonal only for square matrices);
    **vmax** = 1 (fixed), so the top of the scale stays aligned with a natural upper bound.

    For clustermap distance matrices, ``ignore_diagonal=True`` avoids pinning vmin to 0 from the
    diagonal zeros.
    """
    a = np.asarray(arr, dtype=float)
    if ignore_diagonal and a.ndim == 2 and a.shape[0] == a.shape[1] and a.shape[0] > 1:
        a = a.copy()
        np.fill_diagonal(a, np.nan)
    v = a.ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0, 1.0
    vmin = float(np.min(v))
    vmax = 1.0
    if vmin >= vmax:
        vmin = max(0.0, vmax - 1e-9)
    return vmin, vmax


def reference_label(path: str) -> str:
    """Infer reference build from nf-core outdir or path tokens."""
    u = path.replace("\\", "/")
    # Order matters: more specific first
    if re.search(r"nf-core[-_]GRCh37|/GRCh37/|GRC[hH]37", u) and "GRCh38" not in u:
        return "GRCh37"
    if re.search(r"nf-core[-_]GRCh38|/GRCh38/", u):
        return "GRCh38"
    if re.search(r"nf-core[-_]hg19|/hg19/|\bhg19\b", u, re.I):
        return "hg19"
    if re.search(r"GRCh38|grch38", u):
        return "GRCh38"
    if re.search(r"GRCh37|grch37", u):
        return "GRCh37"
    if re.search(r"mm10|GRCm38", u, re.I):
        return "mm10"
    if re.search(r"\bmm9\b", u, re.I):
        return "mm9"
    return "unknown"


def load_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "file_accession" not in df.columns:
        raise ValueError(f"manifest missing file_accession column: {path}")
    # One row per ENCFF IP file
    df = df.drop_duplicates(subset=["file_accession"], keep="first")
    return df


def load_samplesheet_map(path: Optional[Path]) -> Dict[str, str]:
    if path is None or not path.is_file():
        return {}
    df = pd.read_csv(path, dtype=str).fillna("")
    if "sample" not in df.columns or "encode_ip_file_accession" not in df.columns:
        raise ValueError(
            f"samplesheet must have columns sample, encode_ip_file_accession: {path}"
        )
    return dict(zip(df["sample"].astype(str), df["encode_ip_file_accession"].astype(str)))


def match_path_to_manifest_row(
    path: str,
    manifest: pd.DataFrame,
    by_encff: Dict[str, pd.Series],
    samplesheet: Dict[str, str],
) -> Tuple[Optional[pd.Series], str]:
    """Return (row, method) or (None, reason)."""
    for encff in ENCFF_RE.findall(path):
        if encff in by_encff:
            return by_encff[encff], "encff"

    sample = nf_sample_from_basename(path)
    if sample:
        for cand in (sample, normalize_nfcore_sample_stem(sample)):
            if cand in samplesheet:
                encff = samplesheet[cand]
                if encff in by_encff:
                    return by_encff[encff], "sample_sheet"

    m2 = re.search(
        r"(?P<exp>ENCSR[A-Za-z0-9]+)_(?P<ab>[A-Za-z0-9][A-Za-z0-9+.-]*?)(?:_REP\d+)?(?:_peaks|\.fa|\.fasta|/|$)",
        Path(path.replace("\\", "/")).name,
    )
    if m2:
        exp = m2.group("exp")
        ab_raw = re.sub(r"_REP\d+$", "", m2.group("ab"), flags=re.IGNORECASE)
        ab = _norm_histone(ab_raw)
        sub = manifest[manifest["experiment_accession"].astype(str) == exp]
        if not sub.empty and "histone_target" in sub.columns:
            hit = sub[
                sub["histone_target"].map(_norm_histone).eq(ab)
                | sub["histone_target"].astype(str).str.contains(re.escape(ab_raw), case=False, na=False)
            ]
            if len(hit) == 1:
                return hit.iloc[0], "experiment_histone"
            if len(hit) > 1 and "biological_replicate" in hit.columns:
                # Disambiguate with rep in path e.g. ..._R1_
                rep_m = re.search(r"[._-]R(\d+)[._-]", path)
                if rep_m:
                    br = rep_m.group(1)
                    h2 = hit[hit["biological_replicate"].astype(str).str.strip() == br]
                    if len(h2) == 1:
                        return h2.iloc[0], "experiment_histone_rep"

    return None, "unmatched"


SIMILARITY_COL_DEFAULT = "jaccard_similarity_with_ends"


def pick_similarity_column(df: pd.DataFrame, override: Optional[str]) -> str:
    """All plots use ``jaccard_similarity_with_ends`` unless ``--similarity-col`` overrides."""
    if override:
        if override not in df.columns:
            raise ValueError(f"--similarity-col {override!r} not in CSV columns: {list(df.columns)}")
        return override
    if SIMILARITY_COL_DEFAULT not in df.columns:
        raise ValueError(
            "hammock CSV must contain column {!r} (minimizer end-kmer Jaccard). "
            "Re-run hammock with minimizer output, or pass --similarity-col. "
            "Columns found: {}".format(SIMILARITY_COL_DEFAULT, list(df.columns))
        )
    return SIMILARITY_COL_DEFAULT


def annotate_pairs(
    df: pd.DataFrame,
    manifest: pd.DataFrame,
    samplesheet: Dict[str, str],
    sim_col: str,
) -> pd.DataFrame:
    by_encff = {r["file_accession"]: r for _, r in manifest.iterrows()}

    rows = []
    unmatched = 0
    for _, row in df.iterrows():
        f1, f2 = str(row["file1"]), str(row["file2"])
        r1, how1 = match_path_to_manifest_row(f1, manifest, by_encff, samplesheet)
        r2, how2 = match_path_to_manifest_row(f2, manifest, by_encff, samplesheet)
        if r1 is None or r2 is None:
            unmatched += 1
        ref1, ref2 = reference_label(f1), reference_label(f2)

        def pack(r: Optional[pd.Series], path: str) -> Dict[str, Any]:
            if r is None:
                return {
                    "experiment_accession": "",
                    "histone_target": "",
                    "tissue_type": "",
                    "life_stage": "",
                    "biosample_accession": "",
                    "biological_replicate": "",
                    "file_accession": "",
                    "organism": "",
                }
            return {
                "experiment_accession": str(r.get("experiment_accession", "")),
                "histone_target": str(r.get("histone_target", "")),
                "tissue_type": str(r.get("tissue_type", "")),
                "life_stage": str(r.get("life_stage", "")),
                "biosample_accession": str(r.get("biosample_accession", "")),
                "biological_replicate": str(r.get("biological_replicate", "")),
                "file_accession": str(r.get("file_accession", "")),
                "organism": str(r.get("organism", "")),
            }

        p1, p2 = pack(r1, f1), pack(r2, f2)
        sim = float(row[sim_col])
        key = "|".join(
            [
                p1["experiment_accession"],
                p1["histone_target"],
                str(p1["biological_replicate"]),
            ]
        )
        rows.append(
            {
                sim_col: sim,
                "file1": f1,
                "file2": f2,
                "reference1": ref1,
                "reference2": ref2,
                "match1": how1 if r1 is not None else "unmatched",
                "match2": how2 if r2 is not None else "unmatched",
                "sample_key1": key,
                **{f"1_{k}": v for k, v in p1.items()},
                **{f"2_{k}": v for k, v in p2.items()},
            }
        )

    out = pd.DataFrame(rows)
    if unmatched:
        print(f"warning: {unmatched} / {len(out)} rows had missing manifest match on one or both paths", file=sys.stderr)
    return out


def plot_cross_reference_violin_faceted_histone(df: pd.DataFrame, sim_col: str, out: Path) -> None:
    same_bio = (
        (df["1_experiment_accession"] == df["2_experiment_accession"])
        & (df["1_experiment_accession"] != "")
        & (df["1_histone_target"] == df["2_histone_target"])
        & (df["1_tissue_type"] == df["2_tissue_type"])
        & (df["1_life_stage"] == df["2_life_stage"])
        & (df["reference1"] != df["reference2"])
        & (df["reference1"] != "unknown")
        & (df["reference2"] != "unknown")
    )
    sub = df.loc[same_bio].copy()
    if sub.empty:
        return
    sub["pair_label"] = sub["reference1"] + " vs " + sub["reference2"]
    agg = sub.groupby(
        [
            "1_experiment_accession",
            "1_histone_target",
            "1_tissue_type",
            "1_life_stage",
            "pair_label",
        ],
        as_index=False,
    )[sim_col].mean()
    marks = sorted(agg["1_histone_target"].astype(str).unique())
    if not marks:
        return
    n = min(len(marks), 6)
    ncol = 2
    nrow = int(np.ceil(n / ncol)) if n else 1
    fig, axes = plt.subplots(nrow, ncol, figsize=(12, 3.5 * nrow), squeeze=False)
    i = -1
    for i, mark in enumerate(marks[: n]):
        ax = axes.flat[i]
        part = agg[agg["1_histone_target"].astype(str) == mark]
        sns.violinplot(data=part, x="pair_label", y=sim_col, ax=ax, cut=0, inner="box")
        ax.set_title(mark)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=15)
    for j in range(i + 1, nrow * ncol):
        axes.flat[j].set_visible(False)
    fig.suptitle("Cross-reference similarity (faceted by histone mark)")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_cross_reference_violin_by_lifestage_pdf(df: pd.DataFrame, sim_col: str, out_pdf: Path) -> None:
    """
    One PDF page per life stage; each page has a single violin plot with one violin per histone.
    Reference pairs are canonicalized, so A vs B is identical to B vs A.
    """
    same_bio = (
        (df["1_experiment_accession"] == df["2_experiment_accession"])
        & (df["1_experiment_accession"] != "")
        & (df["1_histone_target"] == df["2_histone_target"])
        & (df["1_tissue_type"] == df["2_tissue_type"])
        & (df["1_life_stage"] == df["2_life_stage"])
        & (df["reference1"] != df["reference2"])
        & (df["reference1"] != "unknown")
        & (df["reference2"] != "unknown")
    )
    sub = df.loc[same_bio].copy()
    if sub.empty:
        print("warning: no same-experiment cross-reference pairs to plot (violin PDF)", file=sys.stderr)
        return

    def canon_pair(a: str, b: str) -> str:
        x, y = sorted([str(a), str(b)])
        return f"{x} vs {y}"

    sub["pair_label"] = [canon_pair(a, b) for a, b in zip(sub["reference1"], sub["reference2"])]
    agg = sub.groupby(
        [
            "1_experiment_accession",
            "1_histone_target",
            "1_tissue_type",
            "1_life_stage",
            "pair_label",
        ],
        as_index=False,
    )[sim_col].mean()

    ls_vals = sorted(x for x in agg["1_life_stage"].astype(str).str.strip().unique() if x and x.lower() != "nan")
    min_points_per_page = 5
    boxplot_threshold = 12
    pages = 0
    with PdfPages(out_pdf) as pdf:
        for ls in ls_vals:
            part = agg[agg["1_life_stage"].astype(str).str.strip() == ls].copy()
            if part.empty:
                continue
            n_obs = len(part)
            if n_obs < min_points_per_page:
                print(
                    f"warning: skip cross-reference page for life stage {ls!r}: only {n_obs} points (< {min_points_per_page})",
                    file=sys.stderr,
                )
                continue
            # In typical runs this is one canonical pair (e.g. GRCh37 vs GRCh38), but keep robust.
            pair_vals = sorted(part["pair_label"].astype(str).unique())
            pair_note = pair_vals[0] if len(pair_vals) == 1 else ", ".join(pair_vals[:2]) + (" ..." if len(pair_vals) > 2 else "")
            fig, ax = plt.subplots(figsize=(11, 5))
            x_order = sorted(part["1_histone_target"].astype(str).unique())
            if n_obs < boxplot_threshold:
                sns.boxplot(
                    data=part,
                    x="1_histone_target",
                    y=sim_col,
                    ax=ax,
                    order=x_order,
                    whis=1.5,
                    fliersize=0,
                    width=0.6,
                )
                sns.stripplot(
                    data=part,
                    x="1_histone_target",
                    y=sim_col,
                    ax=ax,
                    order=x_order,
                    color="#333333",
                    jitter=0.14,
                    size=4.2,
                    alpha=0.75,
                )
                plot_kind = f"box+strip (n={n_obs})"
            else:
                sns.violinplot(
                    data=part,
                    x="1_histone_target",
                    y=sim_col,
                    ax=ax,
                    cut=0,
                    inner="box",
                    order=x_order,
                )
                plot_kind = f"violin (n={n_obs})"
            ax.set_xlabel("Histone marker")
            ax.set_ylabel(sim_col)
            ax.set_title(
                f"Cross-reference similarity by histone — life stage: {ls}\n"
                f"Canonical reference pair(s): {pair_note}; plot: {plot_kind}"
            )
            ax.tick_params(axis="x", rotation=20)
            fig.tight_layout()
            pdf.savefig(fig, dpi=150, bbox_inches="tight")
            plt.close(fig)
            pages += 1
    if pages:
        print(f"wrote {out_pdf} ({pages} page(s))")
    else:
        print(f"warning: no pages written for cross-reference violin PDF: {out_pdf}", file=sys.stderr)


def build_similarity_matrix(
    df: pd.DataFrame,
    sim_col: str,
    reference: str,
    histone: Optional[str] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    """Mean similarity per unordered **biosample unit** (no biological replicate in key), both on ``reference``."""
    both = df[(df["reference1"] == reference) & (df["reference2"] == reference)].copy()
    if histone is not None:
        h = str(histone).strip()
        both = both[
            (both["1_histone_target"].astype(str).str.strip() == h)
            & (both["2_histone_target"].astype(str).str.strip() == h)
        ]
    if both.empty:
        return pd.DataFrame(), []

    def unit_key(row: pd.Series, side: str) -> str:
        return "|".join(
            [
                str(row[f"{side}_experiment_accession"]),
                str(row[f"{side}_histone_target"]),
                str(row[f"{side}_tissue_type"]),
                str(row[f"{side}_life_stage"]),
            ]
        )

    rows_agg = []
    for _, row in both.iterrows():
        u1 = unit_key(row, "1")
        u2 = unit_key(row, "2")
        if not u1.replace("|", "").strip() or not u2.replace("|", "").strip():
            continue
        a, b = sorted([u1, u2])
        rows_agg.append((a, b, float(row[sim_col])))

    if not rows_agg:
        return pd.DataFrame(), []

    g = pd.DataFrame(rows_agg, columns=["a", "b", sim_col]).groupby(["a", "b"], as_index=False)[sim_col].mean()
    labels = sorted(set(g["a"]) | set(g["b"]))
    mat = pd.DataFrame(np.eye(len(labels)), index=labels, columns=labels)
    for _, r in g.iterrows():
        i, j = r["a"], r["b"]
        v = float(r[sim_col])
        mat.loc[i, j] = v
        mat.loc[j, i] = v
    return mat, labels


def plot_clustermap(mat: pd.DataFrame, title: str, out: Path, sim_label: str = "similarity") -> None:
    if mat.empty or len(mat) < 2:
        print(f"warning: skip clustermap ({title}): matrix size {mat.shape}", file=sys.stderr)
        return
    dist = 1.0 - mat.clip(0, 1)
    np.fill_diagonal(dist.values, 0.0)
    vmin, vmax = heatmap_vmin_data_vmax_one(dist.values, ignore_diagonal=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cg = sns.clustermap(
            dist,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            xticklabels=True,
            yticklabels=True,
            figsize=(max(10, len(mat) * 0.35), max(10, len(mat) * 0.35)),
            dendrogram_ratio=0.12,
            cbar_kws={
                "label": "1 − %s (mean pairwise); vmin = data min (off-diag), vmax = 1"
                % sim_label
            },
        )
    cg.fig.suptitle(title)
    cg.savefig(out, dpi=150)
    plt.close(cg.fig)

    meta = {lab: parse_sample_uid(lab) for lab in mat.index}

    def tick_fn(lab: str) -> str:
        m = meta.get(lab, {}) or {}
        h = (m.get("histone_target") or "").strip()
        t = (m.get("tissue_type") or "").strip()
        ls = (m.get("life_stage") or "").strip()
        e = (m.get("experiment_accession") or "").strip()
        e_short = e[-10:] if len(e) > 10 else e
        parts = [h, t]
        if ls:
            parts.append(ls)
        parts.append(e_short)
        return "\n".join(parts)

    dend_out = out.parent / out.name.replace("clustermap_", "dendrogram_", 1)
    plot_dendrogram_standalone(
        dist,
        meta,
        "Dendrogram — %s" % title,
        dend_out,
        tick_fn,
        sim_label,
    )


def sample_uid_with_reference(row: pd.Series, side: str) -> str:
    """Biosample unit × reference: experiment | histone | tissue | life_stage | reference (replicates collapsed by mean)."""
    return "|".join(
        [
            str(row[f"{side}_experiment_accession"]),
            str(row[f"{side}_histone_target"]),
            str(row[f"{side}_tissue_type"]),
            str(row[f"{side}_life_stage"]),
            str(row[f"reference{side}"]),
        ]
    )


def parse_sample_uid(uid: str) -> Dict[str, str]:
    """Parse matrix row ids: ``exp|hist|tissue|life_stage`` or ``...|reference`` (no replicate in key)."""
    parts = (uid or "").split("|")
    out = {
        "experiment_accession": "",
        "histone_target": "",
        "biological_replicate": "",
        "tissue_type": "",
        "life_stage": "",
        "reference": "",
    }
    if len(parts) >= 5:
        out["experiment_accession"] = parts[0]
        out["histone_target"] = parts[1]
        out["tissue_type"] = parts[2]
        out["life_stage"] = parts[3]
        out["reference"] = parts[4]
    elif len(parts) >= 4:
        out["experiment_accession"] = parts[0]
        out["histone_target"] = parts[1]
        out["tissue_type"] = parts[2]
        out["life_stage"] = parts[3]
    return out


def build_similarity_matrix_all_tracks(
    df: pd.DataFrame,
    sim_col: str,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, str]]]:
    """
    Pairwise similarity between track ids ``exp|hist|tissue|life_stage|reference`` (mean-aggregated).
    Includes within-reference and cross-reference pairs. Unobserved pairs get similarity 0 (distance 1).
    """
    sub = df[
        (df["reference1"].astype(str) != "unknown")
        & (df["reference2"].astype(str) != "unknown")
    ].copy()
    if sub.empty:
        return pd.DataFrame(), {}

    rows_agg = []
    for _, row in sub.iterrows():
        u1 = sample_uid_with_reference(row, "1")
        u2 = sample_uid_with_reference(row, "2")
        if not u1.replace("|", "").strip() or not u2.replace("|", "").strip():
            continue
        a, b = sorted([u1, u2])
        rows_agg.append((a, b, float(row[sim_col])))

    if not rows_agg:
        return pd.DataFrame(), {}

    g = pd.DataFrame(rows_agg, columns=["a", "b", sim_col]).groupby(["a", "b"], as_index=False)[sim_col].mean()
    labels = sorted(set(g["a"]) | set(g["b"]))
    mat = pd.DataFrame(0.0, index=labels, columns=labels, dtype=float)
    np.fill_diagonal(mat.values, 1.0)
    for _, r in g.iterrows():
        i, j = r["a"], r["b"]
        v = float(r[sim_col])
        mat.loc[i, j] = v
        mat.loc[j, i] = v

    meta = {lab: parse_sample_uid(lab) for lab in labels}
    return mat, meta


# Shared 2D scatter aesthetics: face color = tissue, marker = histone, edge color = reference.
EMBED_SCATTER_ALPHA = 0.48
EMBED_POINT_SIZE = 135
EMBED_EDGE_WIDTH = 1.85


def pca_two_components_similarity_rows(S: np.ndarray) -> np.ndarray:
    """PCA (SVD) on rows of the similarity matrix, columns centered — first two PC scores per sample."""
    X = np.asarray(S, dtype=float)
    n = X.shape[0]
    if n < 2:
        raise ValueError("need at least 2 samples for PCA")
    X = X - np.mean(X, axis=0, keepdims=True)
    U, s, _Vt = np.linalg.svd(X, full_matrices=False)
    out = np.zeros((n, 2))
    out[:, 0] = U[:, 0] * s[0]
    if len(s) > 1:
        out[:, 1] = U[:, 1] * s[1]
    return out


def embedding_style_maps(
    labels: List[str],
    meta: Dict[str, Dict[str, str]],
) -> Tuple[Dict[str, Any], Dict[str, str], Dict[str, Any], List[str], List[str], List[str]]:
    """Maps for tissue (face), histone (marker), reference (edge); plus sorted category lists."""
    tissues = sorted({(meta.get(lab, {}) or {}).get("tissue_type", "") or "?" for lab in labels})
    histones = sorted({(meta.get(lab, {}) or {}).get("histone_target", "") or "?" for lab in labels})
    refs = sorted({(meta.get(lab, {}) or {}).get("reference", "") or "?" for lab in labels})
    n_t, n_h, n_r = max(len(tissues), 1), max(len(histones), 1), max(len(refs), 1)
    tissue_fc = dict(zip(tissues, sns.color_palette("husl", n_colors=n_t)))
    hist_mk = {
        h: ["o", "s", "^", "D", "v", "P", "X", "*", "h", "8"][i % 10] for i, h in enumerate(histones)
    }
    if n_r <= 8:
        ref_ec = dict(zip(refs, sns.color_palette("Dark2", n_colors=n_r)))
    else:
        ref_ec = dict(zip(refs, sns.color_palette("husl", n_colors=n_r)))
    return tissue_fc, hist_mk, ref_ec, tissues, histones, refs


def plot_embedding_scatter_figure(
    XY: np.ndarray,
    labels: List[str],
    meta: Dict[str, Dict[str, str]],
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    alpha: float = EMBED_SCATTER_ALPHA,
):
    """Scatter: face = tissue, marker = histone, edge = reference; semi-transparent points."""
    tissue_fc, hist_mk, ref_ec, tissues, histones, refs = embedding_style_maps(labels, meta)
    lab_i = {lab: i for i, lab in enumerate(labels)}

    fig, ax = plt.subplots(figsize=(12, 9))
    for lab in labels:
        m = meta.get(lab, {}) or {}
        t = m.get("tissue_type", "") or "?"
        h = m.get("histone_target", "") or "?"
        r = m.get("reference", "") or "?"
        i = lab_i[lab]
        fc = tissue_fc.get(t, "#888888")
        mk = hist_mk.get(h, "o")
        ec = ref_ec.get(r, "#222222")
        ax.scatter(
            XY[i, 0],
            XY[i, 1],
            c=[fc],
            marker=mk,
            s=EMBED_POINT_SIZE,
            edgecolors=[ec],
            linewidths=EMBED_EDGE_WIDTH,
            zorder=3,
            alpha=alpha,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.axhline(0, color="#dddddd", linewidth=0.7, zorder=0)
    ax.axvline(0, color="#dddddd", linewidth=0.7, zorder=0)

    h_tissue = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color="w",
            markerfacecolor=tissue_fc[x],
            markeredgecolor="#333333",
            markeredgewidth=0.4,
            markersize=8,
            label=x,
        )
        for x in tissues
    ]
    h_hist = [
        Line2D([0], [0], marker=hist_mk[x], linestyle="", color="gray", markersize=9, label=x)
        for x in histones
    ]
    h_ref = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color="w",
            markerfacecolor="white",
            markeredgecolor=ref_ec[x],
            markeredgewidth=2.0,
            markersize=9,
            label=x,
        )
        for x in refs
    ]

    leg1 = ax.legend(handles=h_tissue, title="Tissue", loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=7)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=h_hist, title="Histone", loc="upper left", bbox_to_anchor=(1.02, 0.72), fontsize=7)
    ax.add_artist(leg2)
    ax.legend(handles=h_ref, title="Reference (edge)", loc="upper left", bbox_to_anchor=(1.02, 0.38), fontsize=7)

    fig.tight_layout()
    return fig


def classical_mds_from_distance(D: np.ndarray, n_components: int = 2) -> np.ndarray:
    """Metric classical MDS from a Euclidean distance matrix (no sklearn)."""
    D = np.asarray(D, dtype=float)
    n = D.shape[0]
    if n < 3:
        raise ValueError("need at least 3 points")
    D = np.maximum(D, 0.0)
    J = np.eye(n) - np.ones((n, n)) / n
    D2 = D * D
    B = -0.5 * J @ D2 @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(-w)
    w = w[idx]
    V = V[:, idx]
    pos = w > 1e-10
    w = w[pos]
    V = V[:, pos]
    k = min(n_components, len(w))
    if k == 0:
        return np.zeros((n, n_components))
    s = np.sqrt(np.maximum(w[:k], 0.0))
    X = V[:, :k] * s
    if X.shape[1] < n_components:
        X = np.hstack([X, np.zeros((n, n_components - X.shape[1]))])
    return X[:, :n_components]


def plot_mds_lifestage_multipage_pdf(
    ann: pd.DataFrame,
    sim_col: str,
    out_pdf: Path,
    ls_vals: List[str],
) -> None:
    """
    One PDF page per life stage: classical MDS on distance = 1 − similarity for that stage's
    track × track matrix (same construction as clustermaps / PCA for that stage).
    """
    coord_rows: List[Dict[str, Any]] = []
    pages = 0
    with PdfPages(out_pdf) as pdf:
        for ls in ls_vals:
            mat, meta = build_similarity_matrix_lifestage(ann, sim_col, ls)
            if mat.empty or len(mat) < 3:
                print(
                    f"warning: skip MDS page for life stage {ls!r}: need >= 3 tracks",
                    file=sys.stderr,
                )
                continue
            dist = 1.0 - mat.clip(0, 1)
            np.fill_diagonal(dist.values, 0.0)
            D = dist.values.astype(float)
            try:
                XY = classical_mds_from_distance(D, n_components=2)
            except Exception as ex:  # pragma: no cover
                print(f"warning: MDS failed for life stage {ls!r}: {ex}", file=sys.stderr)
                continue
            labels = list(mat.index)
            for i, lab in enumerate(labels):
                m = meta.get(lab, {}) or {}
                coord_rows.append(
                    {
                        "life_stage": ls,
                        "track_id": lab,
                        "mds1": XY[i, 0],
                        "mds2": XY[i, 1],
                        "reference": m.get("reference", ""),
                        "histone_target": m.get("histone_target", ""),
                        "tissue_type": m.get("tissue_type", ""),
                    }
                )
            title = (
                "MDS on pairwise distance = 1 − %s, life stage: %s\n"
                "face = tissue, marker = histone, edge = reference; alpha=%.2f"
                % (sim_col, ls, EMBED_SCATTER_ALPHA)
            )
            fig = plot_embedding_scatter_figure(
                XY,
                labels,
                meta,
                title=title,
                xlabel="MDS 1",
                ylabel="MDS 2",
            )
            pdf.savefig(fig, dpi=150, bbox_inches="tight")
            plt.close(fig)
            pages += 1

    if coord_rows:
        pd.DataFrame(coord_rows).to_csv(out_pdf.parent / "mds_coordinates_by_lifestage.tsv", sep="\t", index=False)
    if pages:
        print(f"wrote {out_pdf} ({pages} page(s)), mds_coordinates_by_lifestage.tsv")
    else:
        print(f"warning: no MDS pages written: {out_pdf}", file=sys.stderr)


def plot_pca_lifestage_multipage_pdf(
    ann: pd.DataFrame,
    sim_col: str,
    out_pdf: Path,
    ls_vals: List[str],
) -> None:
    """
    One PDF page per life stage: PCA (SVD) on **rows** of the pairwise similarity matrix
    (columns centered), same similarity metric as the rest of the pipeline.
    """
    coord_rows: List[Dict[str, Any]] = []
    pages = 0
    with PdfPages(out_pdf) as pdf:
        for ls in ls_vals:
            mat, meta = build_similarity_matrix_lifestage(ann, sim_col, ls)
            if mat.empty or len(mat) < 2:
                print(
                    f"warning: skip PCA page for life stage {ls!r}: need >= 2 tracks",
                    file=sys.stderr,
                )
                continue
            try:
                XY = pca_two_components_similarity_rows(mat.values)
            except Exception as ex:  # pragma: no cover
                print(f"warning: PCA failed for life stage {ls!r}: {ex}", file=sys.stderr)
                continue
            labels = list(mat.index)
            for i, lab in enumerate(labels):
                m = meta.get(lab, {}) or {}
                coord_rows.append(
                    {
                        "life_stage": ls,
                        "track_id": lab,
                        "pc1": XY[i, 0],
                        "pc2": XY[i, 1],
                        "reference": m.get("reference", ""),
                        "histone_target": m.get("histone_target", ""),
                        "tissue_type": m.get("tissue_type", ""),
                    }
                )
            title = (
                "PCA (rows of similarity matrix), life stage: %s\n"
                "Columns of S centered; PC1/PC2 from SVD. "
                "face = tissue, marker = histone, edge = reference; alpha=%.2f"
                % (ls, EMBED_SCATTER_ALPHA)
            )
            fig = plot_embedding_scatter_figure(
                XY,
                labels,
                meta,
                title=title,
                xlabel="PC1",
                ylabel="PC2",
            )
            pdf.savefig(fig, dpi=150, bbox_inches="tight")
            plt.close(fig)
            pages += 1

    if coord_rows:
        pd.DataFrame(coord_rows).to_csv(out_pdf.parent / "pca_coordinates_by_lifestage.tsv", sep="\t", index=False)
    if pages:
        print(f"wrote {out_pdf} ({pages} page(s)), pca_coordinates_by_lifestage.tsv")
    else:
        print(f"warning: no PCA pages written (empty life-stage matrices?): {out_pdf}", file=sys.stderr)


def build_similarity_matrix_lifestage(
    df: pd.DataFrame,
    sim_col: str,
    life_stage: str,
    histone: Optional[str] = None,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, str]]]:
    """
    Pairwise similarity for one life stage, including multiple reference builds in the same matrix.
    Returns (symmetric matrix, label -> parsed uid fields).
    """
    ls = str(life_stage).strip()
    sub = df[
        (df["1_life_stage"].astype(str).str.strip() == ls)
        & (df["2_life_stage"].astype(str).str.strip() == ls)
    ].copy()
    if histone is not None:
        h = str(histone).strip()
        sub = sub[
            (sub["1_histone_target"].astype(str).str.strip() == h)
            & (sub["2_histone_target"].astype(str).str.strip() == h)
        ]
    if sub.empty:
        return pd.DataFrame(), {}

    rows_agg = []
    for _, row in sub.iterrows():
        if row["reference1"] == "unknown" or row["reference2"] == "unknown":
            continue
        u1 = sample_uid_with_reference(row, "1")
        u2 = sample_uid_with_reference(row, "2")
        if not u1.replace("|", "").strip() or not u2.replace("|", "").strip():
            continue
        a, b = sorted([u1, u2])
        rows_agg.append((a, b, float(row[sim_col])))

    if not rows_agg:
        return pd.DataFrame(), {}

    g = pd.DataFrame(rows_agg, columns=["a", "b", sim_col]).groupby(["a", "b"], as_index=False)[sim_col].mean()
    labels = sorted(set(g["a"]) | set(g["b"]))
    mat = pd.DataFrame(np.eye(len(labels)), index=labels, columns=labels)
    for _, r in g.iterrows():
        i, j = r["a"], r["b"]
        v = float(r[sim_col])
        mat.loc[i, j] = v
        mat.loc[j, i] = v

    meta = {lab: parse_sample_uid(lab) for lab in labels}
    return mat, meta


def _lifestage_slug(life_stage: str) -> str:
    s = re.sub(r"[^\w\-\.]+", "_", (life_stage or "").strip().lower())
    return s[:80] if s else "unknown"


def histone_slug(mark: str) -> str:
    s = re.sub(r"[^\w\-\.]+", "_", (mark or "").strip())
    return s[:80] if s else "unknown"


def iter_histone_marks(ann: pd.DataFrame) -> List[str]:
    if ann.empty or "1_histone_target" not in ann.columns:
        return []
    s1 = ann["1_histone_target"].astype(str).str.strip()
    s2 = ann["2_histone_target"].astype(str).str.strip()
    u = sorted(set(s1.unique()) | set(s2.unique()))
    return [x for x in u if x and str(x).lower() != "nan"]


def linkage_average(dist: pd.DataFrame):
    """Average linkage on condensed distances (matches seaborn clustermap default)."""
    condensed = squareform(dist.values, checks=False)
    return linkage(condensed, method="average")


def plot_dendrogram_standalone(
    dist: pd.DataFrame,
    meta: Dict[str, Dict[str, str]],
    title: str,
    out: Path,
    tick_label_fn,
    sim_label: str,
) -> None:
    """Standalone hierarchical dendrogram; branches by tissue, leaf label text by experiment group."""
    if dist.empty or len(dist) < 2:
        print(f"warning: skip dendrogram ({title}): matrix size {dist.shape}", file=sys.stderr)
        return

    labels = list(dist.index)
    tick_labels = [tick_label_fn(lab) for lab in labels]

    tissues = sorted({(meta.get(lab, {}) or {}).get("tissue_type", "") or "unknown" for lab in labels})
    n_t = max(len(tissues), 1)
    tissue_palette = ["#7b2cbf", "#f77f00", "#2a9d8f", "#e63946", "#6a994e", "#8f2d56"]
    t_to_color = {t: tissue_palette[i % len(tissue_palette)] for i, t in enumerate(tissues)}

    # Experiment identity for leaf text: same experiment_accession + histone_target share one color.
    exp_groups = []
    for lab in labels:
        m = meta.get(lab, {}) or {}
        exp = (m.get("experiment_accession") or "").strip() or "unknown_exp"
        hist = (m.get("histone_target") or "").strip() or "unknown_histone"
        exp_groups.append(f"{exp}|{hist}")
    exp_groups_unique = sorted(set(exp_groups))
    exp_palette = sns.color_palette("husl", n_colors=max(len(exp_groups_unique), 1))
    exp_to_color = {g: mcolors.to_hex(exp_palette[i % len(exp_palette)]) for i, g in enumerate(exp_groups_unique)}
    leaf_group = {lab: exp_groups[i] for i, lab in enumerate(labels)}

    Z = linkage_average(dist)
    color_thresh = float(0.70 * np.max(Z[:, 2])) if Z.shape[0] else 0.0
    fig_w = max(14.0, len(labels) * 0.28)
    fig, ax = plt.subplots(figsize=(fig_w, 7.0))

    # Build node -> dominant tissue mapping so each branch color follows tissue composition.
    n_leaves = len(labels)
    node_tissues: Dict[int, List[str]] = {}
    for i, lab in enumerate(labels):
        t = (meta.get(lab, {}) or {}).get("tissue_type", "") or "unknown"
        node_tissues[i] = [t]
    for i, row in enumerate(Z):
        left = int(row[0])
        right = int(row[1])
        node_id = n_leaves + i
        node_tissues[node_id] = node_tissues.get(left, []) + node_tissues.get(right, [])

    def dominant_tissue(node_id: int) -> str:
        vals = node_tissues.get(node_id, [])
        if not vals:
            return "unknown"
        return Counter(vals).most_common(1)[0][0]

    def link_color_func(node_id: int) -> str:
        t = dominant_tissue(int(node_id))
        return t_to_color.get(t, "#333333")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ddata = dendrogram(
            Z,
            labels=tick_labels,
            leaf_rotation=90.0,
            ax=ax,
            color_threshold=color_thresh,
            above_threshold_color="#222222",
            link_color_func=link_color_func,
        )

    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Height (average linkage; pairwise input = 1 − %s)" % sim_label)

    # Leaf order left-to-right: ddata["leaves"] are indices into original ``labels``
    for xpos, tick in enumerate(ax.get_xticklabels()):
        li = ddata["leaves"][xpos]
        lab = labels[li]
        g = leaf_group.get(lab, "unknown_exp|unknown_histone")
        tick.set_color(exp_to_color.get(g, "#333333"))
        fs = max(5, min(8, 100 // max(len(labels), 1)))
        tick.set_fontsize(fs)

    tissue_handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="None",
            color="none",
            markerfacecolor=t_to_color[t],
            markeredgecolor=t_to_color[t],
            markeredgewidth=0.6,
            markersize=9,
            label=t,
        )
        for t in tissues
    ]
    leg_tissue = ax.legend(
        handles=tissue_handles,
        title="Branch line color = dominant tissue",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
        fontsize=8,
        title_fontsize=9,
    )
    ax.add_artist(leg_tissue)

    # Keep legend compact: show experiment|histone groups in leaf label colors.
    exp_handles = []
    for g in exp_groups_unique[:12]:
        exp, hist = g.split("|", 1)
        label = f"{exp[-8:]}|{hist}"
        exp_handles.append(
            Line2D([0], [0], marker="o", linestyle="None", color=exp_to_color[g], markersize=6, label=label)
        )
    if len(exp_groups_unique) > 12:
        exp_handles.append(Line2D([0], [0], linestyle="None", color="none", label=f"... +{len(exp_groups_unique)-12} more"))
    ax.legend(
        handles=exp_handles,
        title="Leaf label color = experiment|histone",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.62),
        frameon=True,
        fontsize=7,
        title_fontsize=9,
    )

    fig.subplots_adjust(bottom=0.30, right=0.74, top=0.92)
    fig.text(
        0.5,
        0.06,
        "Branch colors reflect dominant tissue in each merged branch. "
        "Leaf label text color reflects experiment_accession|histone_target.",
        ha="center",
        fontsize=7,
        color="#333333",
    )
    fig.savefig(out, dpi=150, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)


def plot_clustermap_faceted_lifestage(
    mat: pd.DataFrame,
    meta: Dict[str, Dict[str, str]],
    life_stage: str,
    sim_label: str,
    out: Path,
) -> None:
    """
    Hierarchical clustering on 1 − Jaccard distance; tick labels emphasize reference build;
    row/column color strips encode tissue type.
    """
    if mat.empty or len(mat) < 2:
        print(
            f"warning: skip lifestage clustermap ({life_stage!r}): matrix size {mat.shape}",
            file=sys.stderr,
        )
        return

    dist = 1.0 - mat.clip(0, 1)
    np.fill_diagonal(dist.values, 0.0)
    vmin, vmax = heatmap_vmin_data_vmax_one(dist.values, ignore_diagonal=True)

    tissues = sorted({meta.get(i, {}).get("tissue_type", "") or "unknown" for i in mat.index})
    n_t = max(len(tissues), 1)
    pal = sns.color_palette("tab20", n_colors=n_t)
    t_to_color = {t: mcolors.to_hex(pal[i % len(pal)]) for i, t in enumerate(tissues)}

    def tissue_hex(lab: str) -> str:
        t = (meta.get(lab, {}) or {}).get("tissue_type", "") or "unknown"
        return t_to_color.get(t, "#808080")

    row_colors = pd.DataFrame({"tissue": [tissue_hex(i) for i in mat.index]}, index=mat.index)
    col_colors = pd.DataFrame({"tissue": [tissue_hex(i) for i in mat.columns]}, index=mat.columns)

    def tick_label(lab: str) -> str:
        m = meta.get(lab, {}) or {}
        ref = (m.get("reference") or "?").strip()
        hist = (m.get("histone_target") or "").strip()
        exp = (m.get("experiment_accession") or "").strip()
        exp_short = exp[-10:] if len(exp) > 10 else exp
        return f"{ref}\n{hist}\n{exp_short}"

    tl = [tick_label(i) for i in mat.index]
    fs = max(5, min(9, 110 // max(len(mat), 1)))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cg = sns.clustermap(
            dist,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            row_colors=row_colors,
            col_colors=col_colors,
            xticklabels=tl,
            yticklabels=tl,
            figsize=(max(11, len(mat) * 0.42), max(10, len(mat) * 0.42)),
            dendrogram_ratio=0.14,
            cbar_kws={
                "label": f"1 − {sim_label} (mean pairwise); vmin = data min (off-diag), vmax = 1"
            },
            linewidths=0.0,
        )

    for ax in (cg.ax_heatmap,):
        ax.tick_params(axis="both", labelsize=fs)

    cg.fig.suptitle(
        f"Life stage: {life_stage}\n(distance = 1 − Jaccard; strip color = tissue; label lines = reference, histone, experiment)",
        fontsize=11,
        y=1.02,
    )

    handles = [Patch(facecolor=t_to_color[t], edgecolor="none", label=t) for t in tissues]
    cg.fig.legend(
        handles=handles,
        title="Tissue",
        loc="center left",
        bbox_to_anchor=(1.02, 0.55),
        frameon=True,
        fontsize=8,
        title_fontsize=9,
    )

    cg.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(cg.fig)

    dend_out = out.parent / out.name.replace("clustermap_", "dendrogram_", 1)
    plot_dendrogram_standalone(
        dist,
        meta,
        "Dendrogram — life stage: %s" % life_stage,
        dend_out,
        tick_label,
        sim_label,
    )


def tissue_mean_block(
    df: pd.DataFrame,
    sim_col: str,
    reference: str,
    histone: Optional[str] = None,
) -> pd.DataFrame:
    both = df[(df["reference1"] == reference) & (df["reference2"] == reference)].copy()
    if histone is not None:
        h = str(histone).strip()
        both = both[
            (both["1_histone_target"].astype(str).str.strip() == h)
            & (both["2_histone_target"].astype(str).str.strip() == h)
        ]
    both = both[(both["1_tissue_type"] != "") & (both["2_tissue_type"] != "")]
    if both.empty:
        return pd.DataFrame()
    both["t1"] = both["1_tissue_type"].astype(str)
    both["t2"] = both["2_tissue_type"].astype(str)
    g = both.groupby(["t1", "t2"], as_index=False)[sim_col].mean()
    pivot = g.pivot(index="t1", columns="t2", values=sim_col)
    return pivot


def plot_tissue_heatmap(block: pd.DataFrame, title: str, out: Path) -> None:
    if block.empty:
        print(f"warning: skip tissue heatmap ({title}): no data", file=sys.stderr)
        return
    vmin, vmax = heatmap_vmin_data_vmax_one(block.values.astype(float), ignore_diagonal=False)
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(
        block,
        ax=ax,
        cmap="RdYlBu_r",
        vmin=vmin,
        vmax=vmax,
        square=True,
        cbar_kws={"label": "Mean similarity; vmin = data min, vmax = 1"},
    )
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("hammock_csv", type=Path, help="Hammock output CSV")
    p.add_argument("--manifest", type=Path, required=True, help="histone_fastq_manifest.tsv")
    p.add_argument(
        "--samplesheet",
        type=Path,
        default=None,
        help="nf-core CSV with sample + encode_ip_file_accession (recommended)",
    )
    p.add_argument("--outdir", type=Path, default=Path("hammock_histone_plots"))
    p.add_argument("--organism", type=str, default="", help='e.g. "Homo sapiens" (filter annotated rows)')
    p.add_argument(
        "--path-filter",
        type=str,
        default="",
        help="If set, keep only rows where file1 and file2 contain this substring (e.g. broadPeak)",
    )
    p.add_argument(
        "--similarity-col",
        type=str,
        default=None,
        help="Override CSV column (default: jaccard_similarity_with_ends only)",
    )
    p.add_argument(
        "--biological-replicate",
        type=str,
        default="all",
        metavar="N",
        help=(
            "If set to a number (e.g. 1), keep only pairs where both files map to that ENCODE "
            "biological replicate. Default: all (include all replicates; matrices and violins "
            "average replicate combinations)."
        ),
    )
    p.add_argument(
        "--exclude-experiment-accession",
        type=str,
        default="ENCSR330MAM,ENCSR158WBG",
        help=(
            "Comma-separated ENCODE experiment_accession values to drop (any pairwise row "
            "where side 1 or 2 matches is removed). Default: ENCSR330MAM,ENCSR158WBG "
            "(PCA/MDS outliers in typical histone reference runs). Pass '-' to keep all experiments."
        ),
    )
    args = p.parse_args(argv)

    raw = pd.read_csv(args.hammock_csv)
    if "file1" not in raw.columns or "file2" not in raw.columns:
        raise SystemExit("CSV must contain file1 and file2 columns")

    sim_col = pick_similarity_column(raw, args.similarity_col)
    raw[sim_col] = pd.to_numeric(raw[sim_col], errors="coerce")
    raw = raw.dropna(subset=[sim_col])

    if args.path_filter:
        pf = args.path_filter
        raw = raw[raw["file1"].astype(str).str.contains(pf, regex=False) & raw["file2"].astype(str).str.contains(pf, regex=False)]

    manifest = load_manifest(args.manifest)
    samplesheet = load_samplesheet_map(args.samplesheet)

    ann = annotate_pairs(raw, manifest, samplesheet, sim_col)

    if args.organism.strip() and not ann.empty and "1_organism" in ann.columns:
        org = args.organism.strip()
        m1 = ann["1_organism"].astype(str).str.contains(org, case=False, na=False)
        m2 = ann["2_organism"].astype(str).str.contains(org, case=False, na=False)
        ann = ann[m1 & m2]

    br_arg = (args.biological_replicate or "").strip().lower()
    if br_arg not in ("", "all", "none") and not ann.empty and "1_biological_replicate" in ann.columns:
        br = args.biological_replicate.strip()
        b1 = ann["1_biological_replicate"].astype(str).str.strip() == br
        b2 = ann["2_biological_replicate"].astype(str).str.strip() == br
        ann = ann[b1 & b2]

    excl_arg = (args.exclude_experiment_accession or "").strip()
    if excl_arg not in ("", "-", "none", "false") and not ann.empty:
        exclude_set = {x.strip() for x in excl_arg.split(",") if x.strip()}
        if exclude_set:
            e1 = ann["1_experiment_accession"].astype(str).str.strip()
            e2 = ann["2_experiment_accession"].astype(str).str.strip()
            touch = e1.isin(exclude_set) | e2.isin(exclude_set)
            n_drop = int(touch.sum())
            if n_drop:
                ann = ann.loc[~touch].copy()
                print(
                    f"excluded {n_drop} pairwise rows involving experiment_accession "
                    f"{sorted(exclude_set)} (--exclude-experiment-accession)",
                    file=sys.stderr,
                )

    args.outdir.mkdir(parents=True, exist_ok=True)
    tsv_out = args.outdir / "annotated_pairs.tsv"
    ann.to_csv(tsv_out, sep="\t", index=False)
    print(f"wrote {tsv_out}")

    plot_cross_reference_violin_by_lifestage_pdf(
        ann,
        sim_col,
        args.outdir / "cross_reference_violin_by_lifestage.pdf",
    )

    mat_all, _meta_all = build_similarity_matrix_all_tracks(ann, sim_col)
    if not mat_all.empty:
        mat_all.to_csv(args.outdir / "similarity_matrix_all_tracks.tsv", sep="\t")

    ls_vals = sorted(
        {x.strip() for x in ann["1_life_stage"].astype(str).unique()}
        | {x.strip() for x in ann["2_life_stage"].astype(str).unique()}
    )
    ls_vals = [x for x in ls_vals if x and str(x).lower() != "nan"]

    plot_mds_lifestage_multipage_pdf(ann, sim_col, args.outdir / "mds_embedding_by_lifestage.pdf", ls_vals)
    plot_pca_lifestage_multipage_pdf(ann, sim_col, args.outdir / "pca_by_lifestage.pdf", ls_vals)

    refs = sorted({r for r in ann["reference1"].tolist() + ann["reference2"].tolist() if r and r != "unknown"})
    for ref in refs:
        mat, _ = build_similarity_matrix(ann, sim_col, ref)
        if not mat.empty:
            mat.to_csv(args.outdir / f"similarity_matrix_{ref}.tsv", sep="\t")
        plot_clustermap(
            mat,
            f"Clustered distance (1 − Jaccard), {ref}",
            args.outdir / f"clustermap_{ref}.png",
            sim_col,
        )
        block = tissue_mean_block(ann, sim_col, ref)
        if not block.empty:
            block.to_csv(args.outdir / f"tissue_mean_block_{ref}.tsv", sep="\t")
        plot_tissue_heatmap(
            block,
            f"Mean {sim_col} by tissue × tissue ({ref})",
            args.outdir / f"tissue_block_heatmap_{ref}.png",
        )

    for ls in ls_vals:
        mat_ls, meta_ls = build_similarity_matrix_lifestage(ann, sim_col, ls)
        if not mat_ls.empty:
            mat_ls.to_csv(args.outdir / f"similarity_matrix_lifestage_{_lifestage_slug(ls)}.tsv", sep="\t")
        slug = _lifestage_slug(ls)
        plot_clustermap_faceted_lifestage(
            mat_ls,
            meta_ls,
            ls,
            sim_col,
            args.outdir / f"clustermap_lifestage_{slug}.png",
        )

    print(f"plots under {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
