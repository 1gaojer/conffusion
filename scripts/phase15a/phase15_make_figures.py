#!/usr/bin/env python3
"""Make Phase 1.5 retrieval/readout figures.

The inputs are Jerry-owned Phase 1.5 summary tables plus the saved Phase 1.4
MCA tensor shards for the PCA readout panel. Gaeun/shared files are read-only.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


DEFAULT_PHASE15B_DIR = Path("/project/liulab/jg1920/conffusion/phase15b_20260627_failure_analysis")
DEFAULT_PHASE15C_DIR = Path("/project/liulab/jg1920/conffusion/phase15c_20260627_region_readout")
DEFAULT_OUT_DIR = Path("/project/liulab/jg1920/conffusion/figures/phase15_20260627")

READOUT_ORDER = [
    "global_hl",
    "h_cdr3",
    "all_cdrs_hl",
    "global_plus_h_cdr3",
    "global_plus_all_cdrs",
    "all_cdrs_mean_std",
]

HEATMAP_CONDITIONS = [
    "full_128",
    "single_first",
    "h_cdr3_greedy_kcenter_k32",
    "h_cdr3_greedy_kcenter_k64",
    "all_cdrs_greedy_kcenter_k32",
    "all_cdrs_greedy_kcenter_k64",
    "first_k32",
    "first_k64",
    "random_k32_mean",
    "random_k64_mean",
]

GLOBAL_CONDITIONS = [
    "full_128",
    "single_first",
    "h_cdr3_greedy_kcenter_k32",
    "h_cdr3_greedy_kcenter_k64",
    "all_cdrs_greedy_kcenter_k32",
    "all_cdrs_greedy_kcenter_k64",
    "first_k32",
    "first_k64",
    "random_k32_mean",
    "random_k64_mean",
]

READOUT_LABELS = {
    "global_hl": "Global H/L",
    "h_cdr3": "H-CDR3",
    "all_cdrs_hl": "All CDRs",
    "global_plus_h_cdr3": "Global + H-CDR3",
    "global_plus_all_cdrs": "Global + CDRs",
    "all_cdrs_mean_std": "CDR mean+std",
}

CONDITION_LABELS = {
    "full_128": "Full 128",
    "single_first": "Single",
    "h_cdr3_greedy_kcenter_k32": "H-CDR3 K32",
    "h_cdr3_greedy_kcenter_k64": "H-CDR3 K64",
    "all_cdrs_greedy_kcenter_k32": "All-CDR K32",
    "all_cdrs_greedy_kcenter_k64": "All-CDR K64",
    "first_k32": "First K32",
    "first_k64": "First K64",
    "random_k32_mean": "Random K32",
    "random_k64_mean": "Random K64",
    "random_k32": "Random K32",
    "random_k64": "Random K64",
}


@dataclass
class FigureRecord:
    number: int
    stem: str
    title: str
    sources: list[str]
    note: str


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def label_readout(value: str) -> str:
    return READOUT_LABELS.get(value, value)


def label_condition(value: str) -> str:
    return CONDITION_LABELS.get(value, value.replace("_", " "))


def ensure_numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def set_clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", color="#e8e8e8", linewidth=0.8)
    ax.set_axisbelow(True)


def condensed_region_summary(summary: pd.DataFrame, random_rollup: pd.DataFrame) -> pd.DataFrame:
    nonrandom = summary[summary["family"] != "random"].copy()
    rows = []
    for _, row in nonrandom.iterrows():
        rows.append(row.to_dict())
    for _, row in random_rollup.iterrows():
        rows.append(
            {
                "readout": row["readout"],
                "condition": f"random_k{int(row['k'])}_mean",
                "family": "random_mean",
                "k": int(row["k"]),
                "replicate": -1,
                "n_queries": np.nan,
                "n_eval": np.nan,
                "recall_at_10": row["recall_at_10_mean"],
                "mrr": row["mrr_mean"],
                "median_first_positive_rank": row["median_first_positive_rank_mean"],
                "mean_cosine_to_readout_full": row["mean_cosine_to_readout_full_mean"],
            }
        )
    out = pd.DataFrame(rows)
    out["readout"] = pd.Categorical(out["readout"], READOUT_ORDER, ordered=True)
    out["condition"] = pd.Categorical(out["condition"], HEATMAP_CONDITIONS, ordered=True)
    return out.sort_values(["readout", "condition"])


def condensed_global_summary(summary: pd.DataFrame, random_rollup: pd.DataFrame) -> pd.DataFrame:
    nonrandom = summary[summary["family"] != "random"].copy()
    rows = []
    for _, row in nonrandom.iterrows():
        item = row.to_dict()
        item["is_random_mean"] = False
        item["recall_at_10_sd"] = 0.0
        item["mrr_sd"] = 0.0
        rows.append(item)
    for _, row in random_rollup.iterrows():
        rows.append(
            {
                "condition": f"random_k{int(row['k'])}_mean",
                "family": "random_mean",
                "k": int(row["k"]),
                "replicate": -1,
                "recall_at_10": row["recall_at_10_mean"],
                "recall_at_10_sd": row["recall_at_10_sd"],
                "mrr": row["mrr_mean"],
                "mrr_sd": row["mrr_sd"],
                "median_first_positive_rank": row["median_first_positive_rank_mean"],
                "is_random_mean": True,
            }
        )
    out = pd.DataFrame(rows)
    out["condition"] = pd.Categorical(out["condition"], GLOBAL_CONDITIONS, ordered=True)
    return out.sort_values("condition")


def figure_01_heatmap(summary: pd.DataFrame, random_rollup: pd.DataFrame, out_dir: Path) -> FigureRecord:
    df = condensed_region_summary(summary, random_rollup)
    mat = df.pivot_table(index="readout", columns="condition", values="mrr", observed=False)
    mat = mat.reindex(index=READOUT_ORDER, columns=HEATMAP_CONDITIONS)
    values = mat.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(12.5, 5.2))
    im = ax.imshow(values, vmin=0.0, vmax=np.nanmax(values), cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(HEATMAP_CONDITIONS)))
    ax.set_xticklabels([label_condition(c) for c in HEATMAP_CONDITIONS], rotation=35, ha="right")
    ax.set_yticks(range(len(READOUT_ORDER)))
    ax.set_yticklabels([label_readout(r) for r in READOUT_ORDER])
    ax.set_title("MRR by readout and conformer subset", pad=12)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            val = values[i, j]
            if np.isfinite(val):
                color = "white" if val > 0.45 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color=color)
    cbar = fig.colorbar(im, ax=ax, shrink=0.88)
    cbar.set_label("MRR")
    save_figure(fig, out_dir, "01_readout_condition_heatmap_mrr")
    return FigureRecord(
        1,
        "01_readout_condition_heatmap_mrr",
        "Readout x condition heatmap",
        ["region_readout_condition_summary.tsv", "region_readout_random_rollup.tsv"],
        "MRR summarized across readouts and conformer subset conditions.",
    )


def figure_02_global_dotplot(summary: pd.DataFrame, random_rollup: pd.DataFrame, out_dir: Path) -> FigureRecord:
    df = condensed_global_summary(summary, random_rollup)
    x = np.arange(len(df))
    labels = [label_condition(str(c)) for c in df["condition"]]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharex=True)
    for ax, metric, sd_col, ylabel in [
        (axes[0], "mrr", "mrr_sd", "MRR"),
        (axes[1], "recall_at_10", "recall_at_10_sd", "Recall@10"),
    ]:
        y = pd.to_numeric(df[metric], errors="coerce")
        yerr = pd.to_numeric(df.get(sd_col, 0.0), errors="coerce").fillna(0)
        colors = ["#4C78A8" if not bool(v) else "#F58518" for v in df["is_random_mean"]]
        ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor="#555555", elinewidth=1, capsize=3, zorder=1)
        ax.scatter(x, y, s=64, c=colors, edgecolor="white", linewidth=0.8, zorder=2)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        set_clean_axis(ax)
    axes[0].set_title("Global endpoint: rank quality")
    axes[1].set_title("Global endpoint: hit rate")
    fig.tight_layout()
    save_figure(fig, out_dir, "02_global_endpoint_condition_dotplot")
    return FigureRecord(
        2,
        "02_global_endpoint_condition_dotplot",
        "Global endpoint dot plot",
        ["condition_summary.tsv", "random_control_rollup.tsv"],
        "Full ensemble, subset selectors, first-K controls, and random-control means under global H/L pooling.",
    )


def figure_03_delta_lollipop(delta: pd.DataFrame, out_dir: Path) -> FigureRecord:
    df = delta[(delta["condition"] == "full_128") & (delta["readout"] != "global_hl")].copy()
    df = df.sort_values("mrr_minus_global")
    labels = [label_readout(r) for r in df["readout"]]
    y = np.arange(len(df))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, metric, title in [
        (axes[0], "mrr_minus_global", "Delta MRR"),
        (axes[1], "recall10_minus_global", "Delta Recall@10"),
    ]:
        vals = pd.to_numeric(df[metric], errors="coerce").to_numpy()
        ax.axvline(0, color="#777777", linewidth=1)
        for yi, val in zip(y, vals):
            ax.plot([0, val], [yi, yi], color="#999999", linewidth=1.4)
        ax.scatter(vals, y, s=70, color="#54A24B", edgecolor="white", linewidth=0.8, zorder=3)
        ax.set_xlabel(title)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        set_clean_axis(ax)
    axes[0].set_title("Full ensemble readout gain")
    axes[1].set_title("Full ensemble hit-rate gain")
    fig.tight_layout()
    save_figure(fig, out_dir, "03_cdr_readout_delta_vs_global")
    return FigureRecord(
        3,
        "03_cdr_readout_delta_vs_global",
        "CDR-aware readout improvement",
        ["region_readout_delta_vs_global.tsv"],
        "Full-ensemble readout gains relative to global H/L mean pooling.",
    )


def figure_04_rank_distribution(results: pd.DataFrame, out_dir: Path) -> FigureRecord:
    readouts = ["global_hl", "h_cdr3", "all_cdrs_hl", "all_cdrs_mean_std"]
    df = results[(results["condition"] == "full_128") & (results["readout"].isin(readouts))].copy()
    df = df[pd.to_numeric(df["n_positive_after_self_exclusion"], errors="coerce") > 0]
    groups = [
        pd.to_numeric(df.loc[df["readout"] == r, "first_positive_rank"], errors="coerce").dropna().to_numpy()
        for r in readouts
    ]

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    parts = ax.violinplot(groups, showmeans=False, showmedians=True, widths=0.72)
    for body in parts["bodies"]:
        body.set_facecolor("#72B7B2")
        body.set_edgecolor("#3C5F5B")
        body.set_alpha(0.75)
    parts["cmedians"].set_color("#222222")
    rng = np.random.default_rng(42)
    for i, vals in enumerate(groups, start=1):
        jitter = rng.uniform(-0.08, 0.08, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=10, color="#1f1f1f", alpha=0.35, linewidth=0)
    ax.set_xticks(range(1, len(readouts) + 1))
    ax.set_xticklabels([label_readout(r) for r in readouts], rotation=20, ha="right")
    ax.set_yscale("log")
    ax.set_yticks([1, 2, 5, 10, 20, 50, 100, 150])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_ylabel("First positive rank, log scale")
    ax.set_title("True antigen neighbors move upward with CDR-aware pooling")
    set_clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, out_dir, "04_first_positive_rank_distribution")
    return FigureRecord(
        4,
        "04_first_positive_rank_distribution",
        "First positive rank distribution",
        ["region_readout_retrieval_results.tsv"],
        "Distribution of first same-antigen neighbor ranks for full-ensemble readouts.",
    )


def figure_05_target_condition_heatmap(results: pd.DataFrame, out_dir: Path) -> FigureRecord:
    combos = [
        ("global_hl", "full_128", "Global full"),
        ("global_hl", "single_first", "Global single"),
        ("global_hl", "h_cdr3_greedy_kcenter_k64", "Global H-CDR3 K64"),
        ("global_hl", "all_cdrs_greedy_kcenter_k64", "Global all-CDR K64"),
        ("h_cdr3", "full_128", "H-CDR3 full"),
        ("h_cdr3", "h_cdr3_greedy_kcenter_k64", "H-CDR3 K64"),
        ("all_cdrs_mean_std", "full_128", "CDR mean+std full"),
        ("all_cdrs_mean_std", "h_cdr3_greedy_kcenter_k64", "CDR mean+std K64"),
    ]
    eval_targets = results[
        (results["readout"] == "global_hl")
        & (results["condition"] == "full_128")
        & (pd.to_numeric(results["n_positive_after_self_exclusion"], errors="coerce") > 0)
    ]["target_id"].tolist()
    mat = pd.DataFrame(index=eval_targets)
    for readout, condition, label in combos:
        sub = results[(results["readout"] == readout) & (results["condition"] == condition)]
        sub = sub.set_index("target_id")
        mat[label] = pd.to_numeric(sub["reciprocal_rank"], errors="coerce")
    mat["sort_key"] = mat.mean(axis=1, skipna=True)
    mat = mat.sort_values("sort_key", ascending=False).drop(columns=["sort_key"])

    fig, ax = plt.subplots(figsize=(9.6, 8.4))
    masked = np.ma.masked_invalid(mat.to_numpy(dtype=float))
    im = ax.imshow(masked, cmap="magma", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right")
    ax.set_yticks([])
    ax.set_ylabel(f"Targets, n={mat.shape[0]}")
    ax.set_title("Target-level retrieval success patterns")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Reciprocal rank")
    fig.tight_layout()
    save_figure(fig, out_dir, "05_target_condition_failure_heatmap")
    return FigureRecord(
        5,
        "05_target_condition_failure_heatmap",
        "Target x condition failure heatmap",
        ["region_readout_retrieval_results.tsv"],
        "Target-level reciprocal-rank patterns across global and CDR-aware readouts.",
    )


def figure_06_selector_vs_random(selector: pd.DataFrame, out_dir: Path) -> FigureRecord:
    conditions = [
        "h_cdr3_greedy_kcenter_k32",
        "h_cdr3_greedy_kcenter_k64",
        "all_cdrs_greedy_kcenter_k32",
        "all_cdrs_greedy_kcenter_k64",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.5), sharex=True, sharey=True)
    for ax, condition in zip(axes.ravel(), conditions):
        sub = selector[selector["condition"] == condition].copy()
        x = pd.to_numeric(sub["random_rr_mean"], errors="coerce")
        y = pd.to_numeric(sub["reciprocal_rank"], errors="coerce")
        ax.scatter(x, y, s=24, color="#4C78A8", alpha=0.68, edgecolor="white", linewidth=0.3)
        ax.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=1)
        mean_delta = pd.to_numeric(sub["rr_minus_random_mean"], errors="coerce").mean()
        ax.set_title(f"{label_condition(condition)}  mean diff {mean_delta:+.3f}", fontsize=10)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        set_clean_axis(ax)
    for ax in axes[-1, :]:
        ax.set_xlabel("Random-control mean reciprocal rank")
    for ax in axes[:, 0]:
        ax.set_ylabel("Selector reciprocal rank")
    fig.suptitle("Structural selectors versus random subsets", y=0.98)
    fig.tight_layout()
    save_figure(fig, out_dir, "06_selector_vs_random_paired_scatter")
    return FigureRecord(
        6,
        "06_selector_vs_random_paired_scatter",
        "Selector versus random paired plot",
        ["selector_vs_random_by_target.tsv"],
        "Per-target selector reciprocal rank compared with same-budget random-control mean.",
    )


def figure_07_cosine_vs_delta(features: pd.DataFrame, out_dir: Path) -> FigureRecord:
    keep = [
        "single_first",
        "h_cdr3_greedy_kcenter_k32",
        "h_cdr3_greedy_kcenter_k64",
        "all_cdrs_greedy_kcenter_k32",
        "all_cdrs_greedy_kcenter_k64",
        "first_k64",
    ]
    df = features[features["condition"].isin(keep)].copy()
    df["one_minus_cosine"] = 1.0 - pd.to_numeric(df["cosine_to_full"], errors="coerce")
    df["one_minus_cosine"] = df["one_minus_cosine"].clip(lower=1e-7)
    df["rr_delta_vs_full"] = pd.to_numeric(df["rr_delta_vs_full"], errors="coerce")
    color_map = {
        "single_first": "#E45756",
        "h_cdr3_greedy_kcenter_k32": "#4C78A8",
        "h_cdr3_greedy_kcenter_k64": "#72B7B2",
        "all_cdrs_greedy_kcenter_k32": "#F58518",
        "all_cdrs_greedy_kcenter_k64": "#54A24B",
        "first_k64": "#B279A2",
    }

    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    for condition in keep:
        sub = df[df["condition"] == condition]
        ax.scatter(
            sub["one_minus_cosine"],
            sub["rr_delta_vs_full"],
            s=18,
            alpha=0.45,
            label=label_condition(condition),
            color=color_map[condition],
            linewidth=0,
        )
    ax.axhline(0, color="#333333", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel("1 - cosine to full ensemble vector")
    ax.set_ylabel("Reciprocal-rank change vs full")
    ax.set_title("Tiny vector changes can move retrieval ranks")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    set_clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, out_dir, "07_cosine_preservation_vs_rr_change")
    return FigureRecord(
        7,
        "07_cosine_preservation_vs_rr_change",
        "Cosine preservation vs retrieval change",
        ["target_failure_features.tsv"],
        "Near-identical global vectors can still produce different nearest-neighbor ranks.",
    )


def top_label_series(labels: pd.Series, n: int = 6) -> pd.Series:
    clean = labels.fillna("").astype(str)
    clean = clean.where(clean != "", "unlabeled")
    short = clean.str.split(r"\|\|", regex=True).str[0]
    top = set(short.value_counts().head(n).index)
    return short.where(short.isin(top), "other")


def load_phase15c_module():
    import phase15c_region_readout as p15c

    return p15c


def compute_readout_vectors(results: pd.DataFrame, readouts: list[str]) -> tuple[list[str], pd.Series, dict[str, np.ndarray], dict[str, float]]:
    from sklearn.decomposition import PCA

    p15c = load_phase15c_module()
    labels_df = results[
        (results["readout"] == "global_hl")
        & (results["condition"] == "full_128")
        & (pd.to_numeric(results["n_positive_after_self_exclusion"], errors="coerce") > 0)
    ][["target_id", "positive_key"]].drop_duplicates()
    target_ids = labels_df["target_id"].tolist()
    labels = top_label_series(labels_df["positive_key"])
    region_map = p15c.load_region_map(p15c.DEFAULT_CDR_ASSIGNMENTS)
    vectors: dict[str, list[np.ndarray]] = {r: [] for r in readouts}

    for target_id in target_ids:
        shard = p15c.load_shard(p15c.DEFAULT_EMB_DIR / f"{target_id}.pt")
        reps = p15c.chain_reps(shard)
        n_conf = int(reps["H"]["mca_repr"].shape[0])
        selected = list(range(n_conf))
        per_readout = p15c.build_readout_vectors(shard, region_map.get(target_id, {}), selected)
        for readout in readouts:
            vectors[readout].append(per_readout[readout])

    matrix_by_readout = {readout: np.vstack(vals) for readout, vals in vectors.items()}
    explained: dict[str, float] = {}
    pca_by_readout: dict[str, np.ndarray] = {}
    for readout, matrix in matrix_by_readout.items():
        x = p15c.normalize_rows(matrix.astype(np.float32))
        pca = PCA(n_components=2, random_state=42)
        xy = pca.fit_transform(x)
        pca_by_readout[readout] = xy
        explained[readout] = float(pca.explained_variance_ratio_[:2].sum())
    return target_ids, labels.reset_index(drop=True), pca_by_readout, explained


def write_pca_points_table(
    out_dir: Path,
    target_ids: list[str],
    labels: pd.Series,
    xy_by_readout: dict[str, np.ndarray],
    explained: dict[str, float],
) -> Path:
    rows = []
    for readout, xy in xy_by_readout.items():
        for idx, target_id in enumerate(target_ids):
            rows.append(
                {
                    "target_id": target_id,
                    "antigen_label": labels.iloc[idx],
                    "readout": readout,
                    "pc1": float(xy[idx, 0]),
                    "pc2": float(xy[idx, 1]),
                    "explained_var_pc1_pc2": explained[readout],
                }
            )
    path = out_dir / "pca_readout_points.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def load_or_compute_pca_points(
    results: pd.DataFrame,
    out_dir: Path,
    pca_points_table: Path | None,
    readouts: list[str],
) -> tuple[pd.DataFrame, dict[str, float]]:
    if pca_points_table is not None and pca_points_table.exists():
        points = read_table(pca_points_table)
        points = ensure_numeric(points, ["pc1", "pc2", "explained_var_pc1_pc2"])
        explained = points.groupby("readout")["explained_var_pc1_pc2"].first().to_dict()
        return points, {str(k): float(v) for k, v in explained.items()}

    target_ids, labels, xy_by_readout, explained = compute_readout_vectors(results, readouts)
    table_path = write_pca_points_table(out_dir, target_ids, labels, xy_by_readout, explained)
    points = read_table(table_path)
    points = ensure_numeric(points, ["pc1", "pc2", "explained_var_pc1_pc2"])
    return points, explained


def figure_08_pca(results: pd.DataFrame, out_dir: Path, pca_points_table: Path | None) -> FigureRecord:
    readouts = ["global_hl", "h_cdr3", "all_cdrs_mean_std"]
    points, explained = load_or_compute_pca_points(results, out_dir, pca_points_table, readouts)
    labels = points.loc[points["readout"] == readouts[0], "antigen_label"].reset_index(drop=True)
    categories = labels.value_counts().index.tolist()
    palette = {
        cat: color
        for cat, color in zip(
            categories,
            ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#B279A2", "#9D755D", "#BAB0AC"],
        )
    }
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=False, sharey=False)
    for ax, readout in zip(axes, readouts):
        sub_points = points[points["readout"] == readout].reset_index(drop=True)
        xy = sub_points[["pc1", "pc2"]].to_numpy(dtype=float)
        sub_labels = sub_points["antigen_label"]
        for cat in categories:
            mask = sub_labels == cat
            ax.scatter(
                xy[mask, 0],
                xy[mask, 1],
                s=28 if cat != "other" else 16,
                color=palette[cat],
                alpha=0.78 if cat != "other" else 0.35,
                linewidth=0,
                label=cat,
            )
        ax.set_title(f"{label_readout(readout)}  PC1+2 {explained[readout]:.0%}")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles, legend_labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="center right", frameon=False, fontsize=7, title="Antigen label")
    fig.suptitle("PCA of full-ensemble readout vectors", y=1.02)
    fig.tight_layout(rect=(0, 0, 0.88, 1))
    save_figure(fig, out_dir, "08_pca_readout_embeddings")
    return FigureRecord(
        8,
        "08_pca_readout_embeddings",
        "PCA of readout embeddings",
        ["region_readout_retrieval_results.tsv", "Phase 1.4 MCA tensor shards"],
        "PCA projections of full-ensemble global, H-CDR3, and CDR mean+std vectors.",
    )


def figure_09_rank_transition(results: pd.DataFrame, out_dir: Path) -> FigureRecord:
    base = results[(results["readout"] == "global_hl") & (results["condition"] == "full_128")].copy()
    cdr = results[(results["readout"] == "h_cdr3") & (results["condition"] == "full_128")].copy()
    base = base[pd.to_numeric(base["n_positive_after_self_exclusion"], errors="coerce") > 0]
    cdr = cdr[pd.to_numeric(cdr["n_positive_after_self_exclusion"], errors="coerce") > 0]
    merged = base[["target_id", "first_positive_rank"]].merge(
        cdr[["target_id", "first_positive_rank"]],
        on="target_id",
        suffixes=("_global", "_hcdr3"),
    )
    merged["first_positive_rank_global"] = pd.to_numeric(merged["first_positive_rank_global"], errors="coerce")
    merged["first_positive_rank_hcdr3"] = pd.to_numeric(merged["first_positive_rank_hcdr3"], errors="coerce")
    merged = merged.dropna()
    improved = merged["first_positive_rank_hcdr3"] < merged["first_positive_rank_global"]

    fig, ax = plt.subplots(figsize=(6.4, 5.8))
    ax.scatter(
        merged.loc[~improved, "first_positive_rank_global"],
        merged.loc[~improved, "first_positive_rank_hcdr3"],
        s=34,
        color="#E45756",
        alpha=0.7,
        label="same/worse",
        linewidth=0,
    )
    ax.scatter(
        merged.loc[improved, "first_positive_rank_global"],
        merged.loc[improved, "first_positive_rank_hcdr3"],
        s=34,
        color="#54A24B",
        alpha=0.75,
        label="improved",
        linewidth=0,
    )
    max_rank = max(float(merged["first_positive_rank_global"].max()), float(merged["first_positive_rank_hcdr3"].max()))
    ax.plot([1, max_rank], [1, max_rank], color="#666666", linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ticks = [1, 2, 5, 10, 20, 50, 100, 150]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("Global H/L first positive rank")
    ax.set_ylabel("H-CDR3 first positive rank")
    ax.set_title("Rank transition: global to H-CDR3 readout")
    ax.legend(frameon=False)
    set_clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, out_dir, "09_rank_transition_global_to_hcdr3")
    return FigureRecord(
        9,
        "09_rank_transition_global_to_hcdr3",
        "Rank transition plot",
        ["region_readout_retrieval_results.tsv"],
        "Target-level rank changes when switching from global H/L pooling to H-CDR3 pooling.",
    )


def rounded_box(ax: plt.Axes, xy: tuple[float, float], wh: tuple[float, float], text: str, color: str) -> None:
    box = FancyBboxPatch(
        xy,
        wh[0],
        wh[1],
        boxstyle="round,pad=0.02,rounding_size=0.03",
        linewidth=1.2,
        edgecolor="#333333",
        facecolor=color,
        alpha=0.92,
    )
    ax.add_patch(box)
    ax.text(xy[0] + wh[0] / 2, xy[1] + wh[1] / 2, text, ha="center", va="center", fontsize=10)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, linewidth=1.2, color="#333333"))


def figure_10_schematic(global_summary: pd.DataFrame, region_summary: pd.DataFrame, out_dir: Path) -> FigureRecord:
    global_full = float(global_summary.loc[global_summary["condition"] == "full_128", "mrr"].iloc[0])
    global_single = float(global_summary.loc[global_summary["condition"] == "single_first", "mrr"].iloc[0])
    region_full = region_summary[region_summary["condition"] == "full_128"].set_index("readout")
    vals = {
        "Global": float(region_full.loc["global_hl", "mrr"]),
        "H-CDR3": float(region_full.loc["h_cdr3", "mrr"]),
        "CDR mean+std": float(region_full.loc["all_cdrs_mean_std", "mrr"]),
    }

    fig = plt.figure(figsize=(11.5, 4.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.35, 1.0], wspace=0.25)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rounded_box(ax, (0.03, 0.55), (0.2, 0.2), "128\nconformers", "#DCEBFF")
    rounded_box(ax, (0.33, 0.55), (0.2, 0.2), "MCA\ntensors", "#E6F4EA")
    rounded_box(ax, (0.64, 0.70), (0.26, 0.16), "global\npooling", "#FDEBD0")
    rounded_box(ax, (0.64, 0.42), (0.26, 0.16), "CDR-aware\npooling", "#EADCF8")
    rounded_box(ax, (0.64, 0.14), (0.26, 0.16), "retrieval\nranking", "#F5F5F5")
    arrow(ax, (0.23, 0.65), (0.33, 0.65))
    arrow(ax, (0.53, 0.65), (0.64, 0.78))
    arrow(ax, (0.53, 0.65), (0.64, 0.50))
    arrow(ax, (0.77, 0.42), (0.77, 0.30))
    arrow(ax, (0.77, 0.70), (0.77, 0.30))
    ax.text(0.03, 0.17, f"Global MRR: single {global_single:.2f}  full {global_full:.2f}", fontsize=10)

    ax2 = fig.add_subplot(gs[0, 1])
    names = list(vals.keys())
    y = list(vals.values())
    ax2.bar(names, y, color=["#4C78A8", "#54A24B", "#B279A2"])
    ax2.set_ylim(0, max(y) * 1.25)
    ax2.set_ylabel("MRR")
    ax2.set_title("Full-ensemble readout MRR")
    for i, val in enumerate(y):
        ax2.text(i, val + 0.015, f"{val:.2f}", ha="center", va="bottom", fontsize=10)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(True, axis="y", color="#e8e8e8")
    fig.suptitle("Phase 1.5 readout bottleneck diagnostic", y=1.02)
    fig.tight_layout()
    save_figure(fig, out_dir, "10_pipeline_result_schematic")
    return FigureRecord(
        10,
        "10_pipeline_result_schematic",
        "Compact pipeline/result schematic",
        ["condition_summary.tsv", "region_readout_condition_summary.tsv"],
        "Minimal schematic tying conformer ensembles, MCA tensors, pooling, and retrieval MRR.",
    )


def write_manifest(out_dir: Path, records: list[FigureRecord], args: argparse.Namespace) -> None:
    rows = []
    for record in records:
        rows.append(
            {
                "number": record.number,
                "stem": record.stem,
                "png": f"{record.stem}.png",
                "svg": f"{record.stem}.svg",
                "title": record.title,
                "sources": "; ".join(record.sources),
                "note": record.note,
            }
        )
    pd.DataFrame(rows).to_csv(out_dir / "figure_manifest.tsv", sep="\t", index=False)
    with (out_dir / "README.md").open("w") as handle:
        handle.write("# Phase 1.5 Figure Pack\n\n")
        handle.write("Generated from Jerry-owned Phase 1.5 outputs.\n\n")
        handle.write("## Inputs\n\n")
        handle.write(f"- Phase 1.5b: `{args.phase15b_dir}`\n")
        handle.write(f"- Phase 1.5c: `{args.phase15c_dir}`\n")
        if args.pca_points_table:
            handle.write(f"- PCA points: `{args.pca_points_table}`\n\n")
        else:
            p15c = load_phase15c_module()
            handle.write(f"- MCA shards for PCA: `{p15c.DEFAULT_EMB_DIR}`\n\n")
        handle.write("## Figures\n\n")
        for record in records:
            handle.write(f"{record.number}. `{record.stem}.png` / `.svg` - {record.note}\n")
        handle.write("\n## Caveat\n\n")
        handle.write(
            "These are internal diagnostic figures from the 149-target Phase 1.5 set. "
            "They are not independent external validation figures.\n"
        )
    with (out_dir / "figure_inputs.json").open("w") as handle:
        figure_inputs = {
            "phase15b_dir": str(args.phase15b_dir),
            "phase15c_dir": str(args.phase15c_dir),
            "pca_points_table": str(args.pca_points_table) if args.pca_points_table else "",
        }
        if not args.pca_points_table:
            p15c = load_phase15c_module()
            figure_inputs["embedding_dir"] = str(p15c.DEFAULT_EMB_DIR)
            figure_inputs["cdr_assignments"] = str(p15c.DEFAULT_CDR_ASSIGNMENTS)
        json.dump(figure_inputs, handle, indent=2, sort_keys=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase15b-dir", type=Path, default=DEFAULT_PHASE15B_DIR)
    parser.add_argument("--phase15c-dir", type=Path, default=DEFAULT_PHASE15C_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--pca-points-table",
        type=Path,
        default=None,
        help="Optional precomputed PCA points table for figure 08.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8,
            "figure.dpi": 120,
        }
    )

    condition_summary = ensure_numeric(
        read_table(args.phase15b_dir / "condition_summary.tsv"),
        ["recall_at_10", "mrr", "median_first_positive_rank"],
    )
    random_rollup = ensure_numeric(
        read_table(args.phase15b_dir / "random_control_rollup.tsv"),
        ["recall_at_10_mean", "recall_at_10_sd", "mrr_mean", "mrr_sd", "median_first_positive_rank_mean"],
    )
    selector = ensure_numeric(
        read_table(args.phase15b_dir / "selector_vs_random_by_target.tsv"),
        ["random_rr_mean", "reciprocal_rank", "rr_minus_random_mean", "recall10_minus_random_mean"],
    )
    features = ensure_numeric(
        read_table(args.phase15b_dir / "target_failure_features.tsv"),
        ["cosine_to_full", "rr_delta_vs_full"],
    )
    region_summary = ensure_numeric(
        read_table(args.phase15c_dir / "region_readout_condition_summary.tsv"),
        ["recall_at_10", "mrr", "median_first_positive_rank", "mean_cosine_to_readout_full"],
    )
    region_delta = ensure_numeric(
        read_table(args.phase15c_dir / "region_readout_delta_vs_global.tsv"),
        ["recall10_minus_global", "mrr_minus_global", "median_rank_minus_global"],
    )
    region_random = ensure_numeric(
        read_table(args.phase15c_dir / "region_readout_random_rollup.tsv"),
        ["recall_at_10_mean", "recall_at_10_sd", "mrr_mean", "mrr_sd", "median_first_positive_rank_mean"],
    )
    retrieval = ensure_numeric(
        read_table(args.phase15c_dir / "region_readout_retrieval_results.tsv"),
        ["n_positive_after_self_exclusion", "first_positive_rank", "reciprocal_rank"],
    )

    records = [
        figure_01_heatmap(region_summary, region_random, args.out_dir),
        figure_02_global_dotplot(condition_summary, random_rollup, args.out_dir),
        figure_03_delta_lollipop(region_delta, args.out_dir),
        figure_04_rank_distribution(retrieval, args.out_dir),
        figure_05_target_condition_heatmap(retrieval, args.out_dir),
        figure_06_selector_vs_random(selector, args.out_dir),
        figure_07_cosine_vs_delta(features, args.out_dir),
        figure_08_pca(retrieval, args.out_dir, args.pca_points_table),
        figure_09_rank_transition(retrieval, args.out_dir),
        figure_10_schematic(condition_summary, region_summary, args.out_dir),
    ]
    write_manifest(args.out_dir, records, args)
    print(f"[phase15-figures] wrote {len(records)} figures to {args.out_dir}")
    for record in records:
        print(f"{record.number:02d} {record.stem}")


if __name__ == "__main__":
    main()
