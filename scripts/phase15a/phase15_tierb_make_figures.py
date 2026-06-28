#!/usr/bin/env python3
"""Make Tier B Step 4 retrieval/readout figures.

Inputs are copied Jerry-owned Step 4 tables under figure_inputs/. The optional
PCA panel uses a precomputed PCA points TSV exported from MCA readout vectors.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


DEFAULT_INPUT_DIR = Path("figure_inputs/phase15_tierb_step4_20260627")
DEFAULT_INTERNAL_DIR = Path("figure_inputs/phase15d_20260627")
DEFAULT_OUT_DIR = Path("figures/phase15_tierb_step4_20260627")

READOUT_LABELS = {
    "global_hl": "Global H/L",
    "framework_hl": "Framework",
    "h_cdr1": "H-CDR1",
    "h_cdr2": "H-CDR2",
    "h_cdr3": "H-CDR3",
    "l_cdr1": "L-CDR1",
    "l_cdr2": "L-CDR2",
    "l_cdr3": "L-CDR3",
    "all_cdrs_hl": "All CDRs H/L",
    "all_cdrs_mean_std": "All CDRs mean+std",
    "global_plus_h_cdr3": "Global + H-CDR3",
    "global_plus_all_cdrs": "Global + CDRs",
}

FIG_COLORS = {
    "global": "#4C78A8",
    "full": "#4C78A8",
    "cdr": "#54A24B",
    "light": "#72B7B2",
    "control": "#9D755D",
    "random": "#F58518",
    "shuffle": "#BAB0AC",
    "bad": "#E45756",
    "purple": "#B279A2",
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


def read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def short_label(text: str, max_len: int = 38) -> str:
    value = str(text).replace("_", " ").replace("||", "|").strip()
    value = value.split("|")[0].strip() if value else "unlabeled"
    if len(value) > max_len:
        return value[: max_len - 1] + "."
    return value


def label_readout(value: str) -> str:
    return READOUT_LABELS.get(str(value), str(value).replace("_", " "))


def set_clean_axis(ax: plt.Axes, *, xgrid: bool = False) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="x" if xgrid else "y", color="#e7e7e7", linewidth=0.8)
    ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def full_condition(summary: pd.DataFrame) -> str:
    vals = [str(v) for v in summary["condition"].dropna().unique()]
    for value in vals:
        if value.startswith("full_"):
            return value
    raise ValueError("No full_* condition found")


def bool_to_float(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": 1.0, "false": 0.0}).astype(float)


def load_inputs(input_dir: Path, internal_dir: Path) -> dict[str, object]:
    data: dict[str, object] = {}
    data["endpoint"] = read_json(input_dir / "endpoint" / "tierb_endpoint_prep_summary.json")
    data["target_overlap"] = read_table(input_dir / "endpoint" / "target_overlap.tsv")
    data["region_summary"] = numeric(
        read_table(input_dir / "phase15c" / "region_readout_condition_summary.tsv"),
        [
            "k",
            "replicate",
            "n_queries",
            "n_eval",
            "recall_at_1",
            "recall_at_5",
            "recall_at_10",
            "mrr",
            "median_first_positive_rank",
            "mean_cosine_to_readout_full",
        ],
    )
    data["region_results"] = numeric(
        read_table(input_dir / "phase15c" / "region_readout_retrieval_results.tsv"),
        [
            "n_positive_after_self_exclusion",
            "first_positive_rank",
            "reciprocal_rank",
            "top1_similarity",
        ],
    )
    data["region_random"] = numeric(
        read_table(input_dir / "phase15c" / "region_readout_random_rollup.tsv"),
        [
            "k",
            "recall_at_10_mean",
            "recall_at_10_sd",
            "mrr_mean",
            "mrr_sd",
            "median_first_positive_rank_mean",
        ],
    )
    data["control_summary"] = numeric(
        read_table(input_dir / "phase15d" / "region_control_summary.tsv"),
        [
            "n_queries",
            "n_eval",
            "recall_at_1",
            "recall_at_5",
            "recall_at_10",
            "mrr",
            "median_first_positive_rank",
        ],
    )
    data["control_results"] = numeric(
        read_table(input_dir / "phase15d" / "region_control_retrieval_results.tsv"),
        [
            "n_positive_after_self_exclusion",
            "first_positive_rank",
            "reciprocal_rank",
            "top1_similarity",
        ],
    )
    data["random_window"] = numeric(
        read_table(input_dir / "phase15d" / "region_control_random_window_rollup.tsv"),
        [
            "recall_at_10_mean",
            "recall_at_10_sd",
            "mrr_mean",
            "mrr_sd",
            "median_first_positive_rank_mean",
        ],
    )
    data["shuffled"] = numeric(
        read_table(input_dir / "phase15d" / "region_control_shuffled_label_summary.tsv"),
        [
            "recall_at_10_mean",
            "recall_at_10_sd",
            "mrr_mean",
            "mrr_sd",
            "median_first_positive_rank_mean",
        ],
    )
    data["internal_control_summary"] = numeric(
        read_table(internal_dir / "region_control_summary.tsv"),
        [
            "n_queries",
            "n_eval",
            "recall_at_10",
            "mrr",
            "median_first_positive_rank",
        ],
    )
    data["internal_random_window"] = numeric(
        read_table(internal_dir / "region_control_random_window_rollup.tsv"),
        [
            "recall_at_10_mean",
            "mrr_mean",
            "median_first_positive_rank_mean",
        ],
    )
    return data


def metric_row(summary: pd.DataFrame, readout: str) -> pd.Series:
    sub = summary[summary["readout"] == readout]
    if sub.empty:
        raise KeyError(readout)
    return sub.iloc[0]


def figure_01_ladder(data: dict[str, object], out_dir: Path) -> FigureRecord:
    summary = data["control_summary"]  # type: ignore[assignment]
    random_window = data["random_window"]  # type: ignore[assignment]
    shuffled = data["shuffled"]  # type: ignore[assignment]
    rows = []
    order = [
        ("Shuffled global", "shuffle", float(shuffled.loc[shuffled["readout"] == "global_hl", "mrr_mean"].iloc[0]), float(shuffled.loc[shuffled["readout"] == "global_hl", "recall_at_10_mean"].iloc[0])),
        ("Random framework", "random", float(random_window["mrr_mean"].iloc[0]), float(random_window["recall_at_10_mean"].iloc[0])),
        ("Framework", "control", float(metric_row(summary, "framework_hl")["mrr"]), float(metric_row(summary, "framework_hl")["recall_at_10"])),
        ("Global H/L", "global", float(metric_row(summary, "global_hl")["mrr"]), float(metric_row(summary, "global_hl")["recall_at_10"])),
        ("H-CDR3", "cdr", float(metric_row(summary, "h_cdr3")["mrr"]), float(metric_row(summary, "h_cdr3")["recall_at_10"])),
        ("L-CDR3", "light", float(metric_row(summary, "l_cdr3")["mrr"]), float(metric_row(summary, "l_cdr3")["recall_at_10"])),
        ("All CDRs mean+std", "purple", float(metric_row(summary, "all_cdrs_mean_std")["mrr"]), float(metric_row(summary, "all_cdrs_mean_std")["recall_at_10"])),
    ]
    for label, kind, mrr, recall in order:
        rows.append({"label": label, "kind": kind, "mrr": mrr, "recall": recall})
    df = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
    y = np.arange(len(df))
    colors = [FIG_COLORS[k] for k in df["kind"]]
    for ax, metric, title in [(axes[0], "mrr", "MRR"), (axes[1], "recall", "Recall@10")]:
        vals = df[metric].to_numpy(float)
        ax.barh(y, vals, color=colors, alpha=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels(df["label"])
        ax.set_xlabel(title)
        ax.set_xlim(0, max(vals) * 1.25)
        for yi, val in zip(y, vals):
            ax.text(val + max(vals) * 0.025, yi, f"{val:.3f}", va="center", fontsize=8)
        set_clean_axis(ax, xgrid=True)
    fig.suptitle("Tier B Step 4: CDR-aware readouts retrieve same-antigen neighbors better", y=1.02)
    fig.tight_layout()
    save_figure(fig, out_dir, "01_readout_performance_ladder")
    return FigureRecord(
        1,
        "01_readout_performance_ladder",
        "Readout performance ladder",
        ["phase15d/region_control_summary.tsv", "phase15d/region_control_random_window_rollup.tsv", "phase15d/region_control_shuffled_label_summary.tsv"],
        "Global, framework, CDR-aware, random-window, and shuffled-label readouts ranked by MRR and Recall@10.",
    )


def paired_full_results(region_results: pd.DataFrame, readout_a: str, readout_b: str, full: str) -> pd.DataFrame:
    cols = ["target_id", "positive_key", "n_positive_after_self_exclusion", "first_positive_rank", "reciprocal_rank", "recall_at_10"]
    a = region_results[(region_results["readout"] == readout_a) & (region_results["condition"] == full)][cols].copy()
    b = region_results[(region_results["readout"] == readout_b) & (region_results["condition"] == full)][cols].copy()
    a = a[a["n_positive_after_self_exclusion"] > 0]
    b = b[b["n_positive_after_self_exclusion"] > 0]
    merged = a.merge(b, on=["target_id", "positive_key"], suffixes=("_a", "_b"))
    return numeric(merged, ["first_positive_rank_a", "first_positive_rank_b", "reciprocal_rank_a", "reciprocal_rank_b"])


def figure_02_rank_shift(data: dict[str, object], out_dir: Path) -> FigureRecord:
    summary = data["region_summary"]  # type: ignore[assignment]
    results = data["region_results"]  # type: ignore[assignment]
    full = full_condition(summary)
    merged = paired_full_results(results, "global_hl", "all_cdrs_mean_std", full)
    merged = merged.sort_values("first_positive_rank_a", ascending=False).reset_index(drop=True)
    n = len(merged)
    y = np.arange(n)
    improved = merged["first_positive_rank_b"] < merged["first_positive_rank_a"]

    fig, ax = plt.subplots(figsize=(7.8, 8.8))
    for idx, row in merged.iterrows():
        color = FIG_COLORS["cdr"] if improved.iloc[idx] else FIG_COLORS["bad"]
        alpha = 0.55 if improved.iloc[idx] else 0.35
        ax.plot([0, 1], [row["first_positive_rank_a"], row["first_positive_rank_b"]], color=color, alpha=alpha, linewidth=1)
    ax.scatter(np.zeros(n), merged["first_positive_rank_a"], s=18, color=FIG_COLORS["global"], label="Global H/L", zorder=3)
    ax.scatter(np.ones(n), merged["first_positive_rank_b"], s=18, color=FIG_COLORS["purple"], label="All CDRs mean+std", zorder=3)
    ax.axhline(10, color="#333333", linestyle="--", linewidth=1, alpha=0.75)
    ax.text(1.03, 10, "rank 10", va="center", fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Global H/L", "All CDRs\nmean+std"])
    ax.set_yscale("log")
    max_rank = max(float(merged["first_positive_rank_a"].max()), float(merged["first_positive_rank_b"].max()))
    ticks = [1, 2, 5, 10, 20, 50, 100, 160]
    ax.set_yticks([t for t in ticks if t <= max_rank * 1.1])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.invert_yaxis()
    ax.set_ylabel("First same-antigen neighbor rank, lower is better")
    ax.set_title(f"Rank shift per evaluable target, n={n}")
    ax.legend(frameon=False, loc="lower center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", color="#e8e8e8")
    fig.tight_layout()
    save_figure(fig, out_dir, "02_global_to_cdr_rank_shift")
    return FigureRecord(
        2,
        "02_global_to_cdr_rank_shift",
        "Global to CDR rank shift",
        ["phase15c/region_readout_retrieval_results.tsv"],
        "Target-level first-positive ranks when switching from global pooling to all-CDR mean+std pooling.",
    )


def ecdf(values: pd.Series, max_rank: int) -> tuple[np.ndarray, np.ndarray]:
    vals = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    xs = np.arange(1, max_rank + 1)
    ys = np.array([(vals <= x).mean() for x in xs])
    return xs, ys


def figure_03_first_neighbor_cdf(data: dict[str, object], out_dir: Path) -> FigureRecord:
    region_summary = data["region_summary"]  # type: ignore[assignment]
    region_results = data["region_results"]  # type: ignore[assignment]
    control_results = data["control_results"]  # type: ignore[assignment]
    full = full_condition(region_summary)
    specs = [
        ("Global H/L", region_results[(region_results["condition"] == full) & (region_results["readout"] == "global_hl")], FIG_COLORS["global"]),
        ("H-CDR3", region_results[(region_results["condition"] == full) & (region_results["readout"] == "h_cdr3")], FIG_COLORS["cdr"]),
        ("L-CDR3", control_results[control_results["readout"] == "l_cdr3"], FIG_COLORS["light"]),
        ("All CDRs mean+std", region_results[(region_results["condition"] == full) & (region_results["readout"] == "all_cdrs_mean_std")], FIG_COLORS["purple"]),
    ]
    max_rank = 160
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for label, df, color in specs:
        valid = df[df["n_positive_after_self_exclusion"] > 0]
        xs, ys = ecdf(valid["first_positive_rank"], max_rank)
        ax.step(xs, ys, where="post", label=label, color=color, linewidth=2)
    ax.axvline(10, color="#333333", linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100, 160])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Rank cutoff")
    ax.set_ylabel("Fraction with a same-antigen neighbor by this rank")
    ax.set_title("First correct neighbor moves earlier under CDR-aware readouts")
    ax.legend(frameon=False)
    set_clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, out_dir, "03_first_correct_neighbor_cdf")
    return FigureRecord(
        3,
        "03_first_correct_neighbor_cdf",
        "First correct neighbor CDF",
        ["phase15c/region_readout_retrieval_results.tsv", "phase15d/region_control_retrieval_results.tsv"],
        "Cumulative distribution of first same-antigen neighbor rank for global and CDR-aware readouts.",
    )


def lookup_metric(summary: pd.DataFrame, random_window: pd.DataFrame, readout: str, metric: str) -> float:
    if readout == "random_framework_window":
        return float(random_window[f"{metric}_mean"].iloc[0])
    return float(metric_row(summary, readout)[metric])


def figure_04_replication_heatmap(data: dict[str, object], out_dir: Path) -> FigureRecord:
    internal = data["internal_control_summary"]  # type: ignore[assignment]
    internal_random = data["internal_random_window"]  # type: ignore[assignment]
    tierb = data["control_summary"]  # type: ignore[assignment]
    tierb_random = data["random_window"]  # type: ignore[assignment]
    readouts = ["global_hl", "framework_hl", "random_framework_window", "h_cdr3", "l_cdr3", "all_cdrs_mean_std"]
    labels = ["Global", "Framework", "Random\nframework", "H-CDR3", "L-CDR3", "All CDRs\nmean+std"]
    mat = np.array(
        [
            [lookup_metric(internal, internal_random, r, "mrr") for r in readouts],
            [lookup_metric(tierb, tierb_random, r, "mrr") for r in readouts],
        ],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(8.8, 3.3))
    im = ax.imshow(mat, cmap="viridis", vmin=0, vmax=np.nanmax(mat), aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Internal 149", "Tier B 160"])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            color = "white" if mat[i, j] > np.nanmax(mat) * 0.55 else "black"
            ax.text(j, i, f"{mat[i, j]:.3f}", ha="center", va="center", color=color, fontsize=9)
    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("MRR")
    ax.set_title("CDR-aware readout pattern replicates on the stricter endpoint")
    fig.tight_layout()
    save_figure(fig, out_dir, "04_internal_vs_tierb_mrr_heatmap")
    return FigureRecord(
        4,
        "04_internal_vs_tierb_mrr_heatmap",
        "Internal vs Tier B replication heatmap",
        ["figure_inputs/phase15d_20260627/region_control_summary.tsv", "phase15d/region_control_summary.tsv"],
        "MRR pattern across readouts in the original internal diagnostic and stricter Tier B endpoint.",
    )


def figure_05_controls(data: dict[str, object], out_dir: Path) -> FigureRecord:
    summary = data["control_summary"]  # type: ignore[assignment]
    random_window = data["random_window"]  # type: ignore[assignment]
    shuffled = data["shuffled"]  # type: ignore[assignment]
    bars = [
        ("Global H/L", "global", float(metric_row(summary, "global_hl")["mrr"]), None),
        ("Framework", "control", float(metric_row(summary, "framework_hl")["mrr"]), None),
        ("Random framework", "random", float(random_window["mrr_mean"].iloc[0]), float(random_window["mrr_sd"].iloc[0])),
        ("Shuffled global", "shuffle", float(shuffled.loc[shuffled["readout"] == "global_hl", "mrr_mean"].iloc[0]), float(shuffled.loc[shuffled["readout"] == "global_hl", "mrr_sd"].iloc[0])),
        ("H-CDR3", "cdr", float(metric_row(summary, "h_cdr3")["mrr"]), None),
        ("L-CDR3", "light", float(metric_row(summary, "l_cdr3")["mrr"]), None),
        ("All CDRs mean+std", "purple", float(metric_row(summary, "all_cdrs_mean_std")["mrr"]), None),
    ]
    df = pd.DataFrame(bars, columns=["label", "kind", "mrr", "sd"])
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(9.4, 4.5))
    ax.bar(x, df["mrr"], color=[FIG_COLORS[k] for k in df["kind"]], alpha=0.88)
    yerr = df["sd"].fillna(0).to_numpy(float)
    ax.errorbar(x, df["mrr"], yerr=yerr, fmt="none", color="#333333", capsize=3, linewidth=1)
    for xi, val in zip(x, df["mrr"]):
        ax.text(xi, val + 0.012, f"{val:.3f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=28, ha="right")
    ax.set_ylabel("MRR")
    ax.set_title("Controls: true CDR readouts beat framework and shuffled baselines")
    set_clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, out_dir, "05_controls_cdr_vs_framework_shuffle")
    return FigureRecord(
        5,
        "05_controls_cdr_vs_framework_shuffle",
        "Controls: CDR vs framework and shuffled labels",
        ["phase15d/region_control_summary.tsv", "phase15d/region_control_random_window_rollup.tsv", "phase15d/region_control_shuffled_label_summary.tsv"],
        "Focused control comparison showing CDR readouts against framework, random-window, and shuffled-label baselines.",
    )


def figure_06_compression(data: dict[str, object], out_dir: Path) -> FigureRecord:
    summary = data["region_summary"]  # type: ignore[assignment]
    random = data["region_random"]  # type: ignore[assignment]
    full = full_condition(summary)
    readout = "all_cdrs_mean_std"
    selectors = [
        ("Full 100", "full", summary[(summary["readout"] == readout) & (summary["condition"] == full)].iloc[0], None),
        ("H-CDR3 K64", "cdr", summary[(summary["readout"] == readout) & (summary["condition"] == "h_cdr3_greedy_kcenter_k64")].iloc[0], None),
        ("All-CDR K64", "purple", summary[(summary["readout"] == readout) & (summary["condition"] == "all_cdrs_greedy_kcenter_k64")].iloc[0], None),
    ]
    random64 = random[(random["readout"] == readout) & (random["k"] == 64)].iloc[0]
    rows = []
    for label, kind, row, err in selectors:
        rows.append({"label": label, "kind": kind, "mrr": float(row["mrr"]), "mrr_sd": 0.0, "recall": float(row["recall_at_10"]), "recall_sd": 0.0})
    rows.append({"label": "Random K64 mean", "kind": "random", "mrr": float(random64["mrr_mean"]), "mrr_sd": float(random64["mrr_sd"]), "recall": float(random64["recall_at_10_mean"]), "recall_sd": float(random64["recall_at_10_sd"])})
    df = pd.DataFrame(rows)
    x = np.arange(len(df))

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3), sharex=True)
    for ax, metric, sd_col, title in [(axes[0], "mrr", "mrr_sd", "MRR"), (axes[1], "recall", "recall_sd", "Recall@10")]:
        ax.bar(x, df[metric], color=[FIG_COLORS[k] for k in df["kind"]], alpha=0.9)
        ax.errorbar(x, df[metric], yerr=df[sd_col], fmt="none", color="#333333", capsize=3, linewidth=1)
        for xi, val in zip(x, df[metric]):
            ax.text(xi, val + 0.012, f"{val:.3f}", ha="center", va="bottom", fontsize=8)
        ax.set_ylabel(title)
        ax.set_xticks(x)
        ax.set_xticklabels(df["label"], rotation=25, ha="right")
        set_clean_axis(ax)
    fig.suptitle("Compression honesty check: K64 selectors are close, but random is competitive", y=1.02)
    fig.tight_layout()
    save_figure(fig, out_dir, "06_compression_honesty_k64")
    return FigureRecord(
        6,
        "06_compression_honesty_k64",
        "Compression honesty plot",
        ["phase15c/region_readout_condition_summary.tsv", "phase15c/region_readout_random_rollup.tsv"],
        "Full 100 conformers versus K64 structural selectors and same-budget random K64 under the strongest readout.",
    )


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def figure_07_endpoint_funnel(data: dict[str, object], out_dir: Path, repo_root: Path) -> FigureRecord:
    endpoint = data["endpoint"]  # type: ignore[assignment]
    manifest_dir = repo_root / "manifests"
    generated_test = count_rows(manifest_dir / "gaeun_conformer_ensembles_generated_non_1000_test_20260627.tsv")
    generated_all = count_rows(manifest_dir / "gaeun_conformer_ensembles_generated_non_1000_all_20260627.tsv")
    tierb_pass = count_rows(manifest_dir / "tierb_sequence_deoverlap_20260627" / "gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_pass_20260627.tsv")
    rows = [
        ("Generated, not 1000 test IDs", generated_test),
        ("Exact-PDB unseen vs full checkpoint", generated_all),
        ("Tier B VH+VL <0.85 pass", tierb_pass),
        ("Usable Step 4 endpoint", int(endpoint["target_count"])),
        ("Evaluable retrieval queries", 74),
    ]
    df = pd.DataFrame(rows, columns=["stage", "count"])
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    y = np.arange(len(df))
    colors = ["#4C78A8", "#72B7B2", "#54A24B", "#F58518", "#B279A2"]
    ax.barh(y, df["count"], color=colors, alpha=0.88)
    ax.set_yticks(y)
    ax.set_yticklabels(df["stage"])
    ax.invert_yaxis()
    ax.set_xlabel("Targets")
    for yi, val in zip(y, df["count"]):
        ax.text(val + max(df["count"]) * 0.015, yi, f"{val:,}", va="center", fontsize=9)
    ax.set_title("Endpoint funnel for Tier B Step 4")
    set_clean_axis(ax, xgrid=True)
    fig.tight_layout()
    save_figure(fig, out_dir, "07_endpoint_funnel")
    return FigureRecord(
        7,
        "07_endpoint_funnel",
        "Endpoint funnel",
        ["manifests/*.tsv", "endpoint/tierb_endpoint_prep_summary.json"],
        "How the Gaeun-generated pool narrows into the Tier B Step 4 evaluable retrieval endpoint.",
    )


def per_label_metrics(region_results: pd.DataFrame, control_results: pd.DataFrame, full: str) -> pd.DataFrame:
    sources = [
        ("Global H/L", region_results[(region_results["condition"] == full) & (region_results["readout"] == "global_hl")]),
        ("All CDRs mean+std", region_results[(region_results["condition"] == full) & (region_results["readout"] == "all_cdrs_mean_std")]),
        ("L-CDR3", control_results[control_results["readout"] == "l_cdr3"]),
    ]
    rows = []
    for readout, df in sources:
        valid = df[df["n_positive_after_self_exclusion"] > 0].copy()
        valid["recall_at_10_num"] = bool_to_float(valid["recall_at_10"])
        for label, group in valid.groupby("positive_key"):
            rows.append(
                {
                    "positive_key": label,
                    "label": short_label(label, 42),
                    "readout": readout,
                    "n": int(len(group)),
                    "mrr": float(group["reciprocal_rank"].mean()),
                    "recall_at_10": float(group["recall_at_10_num"].mean()),
                }
            )
    return pd.DataFrame(rows)


def figure_08_per_antigen_heatmap(data: dict[str, object], out_dir: Path) -> FigureRecord:
    summary = data["region_summary"]  # type: ignore[assignment]
    region_results = data["region_results"]  # type: ignore[assignment]
    control_results = data["control_results"]  # type: ignore[assignment]
    full = full_condition(summary)
    df = per_label_metrics(region_results, control_results, full)
    counts = df.groupby(["positive_key", "label"])["n"].max().reset_index().sort_values("n", ascending=False)
    keep = counts[counts["n"] >= 2].head(14)
    sub = df.merge(keep[["positive_key", "label"]], on=["positive_key", "label"])
    pivot = sub.pivot_table(index="label", columns="readout", values="mrr", aggfunc="first")
    pivot = pivot.reindex(index=keep["label"].tolist(), columns=["Global H/L", "All CDRs mean+std", "L-CDR3"])

    fig, ax = plt.subplots(figsize=(8.8, 6.7))
    values = pivot.to_numpy(float)
    im = ax.imshow(values, cmap="viridis", vmin=0, vmax=max(0.55, np.nanmax(values)), aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=20, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ylabels = [f"{label}  n={int(keep.loc[keep['label'] == label, 'n'].iloc[0])}" for label in pivot.index]
    ax.set_yticklabels(ylabels)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if np.isfinite(values[i, j]):
                color = "white" if values[i, j] > np.nanmax(values) * 0.55 else "black"
                ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center", fontsize=8, color=color)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Mean reciprocal rank")
    ax.set_title("Per-antigen readout improvements, labels with at least two queries")
    fig.tight_layout()
    save_figure(fig, out_dir, "08_per_antigen_readout_heatmap")
    return FigureRecord(
        8,
        "08_per_antigen_readout_heatmap",
        "Per-antigen improvement heatmap",
        ["phase15c/region_readout_retrieval_results.tsv", "phase15d/region_control_retrieval_results.tsv"],
        "MRR by antigen label for global pooling, all-CDR mean+std, and L-CDR3.",
    )


def figure_09_pca(data: dict[str, object], out_dir: Path, pca_points_table: Path) -> FigureRecord:
    points = read_table(pca_points_table)
    points = numeric(points, ["pc1", "pc2", "explained_var_pc1_pc2"])
    readouts = ["global_hl", "all_cdrs_mean_std"]
    labels = points[points["readout"] == readouts[0]]["label_short"].value_counts()
    top = labels.head(6).index.tolist()
    points["plot_label"] = points["label_short"].where(points["label_short"].isin(top), "other")
    categories = top + (["other"] if "other" in set(points["plot_label"]) else [])
    palette = dict(
        zip(
            categories,
            ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#B279A2", "#BAB0AC"],
        )
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
    for ax, readout in zip(axes, readouts):
        sub = points[points["readout"] == readout]
        for category in categories:
            part = sub[sub["plot_label"] == category]
            ax.scatter(
                part["pc1"],
                part["pc2"],
                s=30 if category != "other" else 16,
                color=palette[category],
                alpha=0.78 if category != "other" else 0.35,
                linewidth=0,
                label=category,
            )
        explained = float(sub["explained_var_pc1_pc2"].iloc[0]) if not sub.empty else float("nan")
        ax.set_title(f"{label_readout(readout)}  PC1+2 {explained:.0%}")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles, labels_ = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="center right", frameon=False, fontsize=7, title="Antigen label")
    fig.suptitle("PCA of Tier B full-ensemble readout vectors", y=1.02)
    fig.tight_layout(rect=(0, 0, 0.83, 1))
    save_figure(fig, out_dir, "09_pca_global_vs_cdr_readouts")
    return FigureRecord(
        9,
        "09_pca_global_vs_cdr_readouts",
        "PCA side-by-side",
        [str(pca_points_table)],
        "PCA of actual full-ensemble global and all-CDR mean+std readout vectors, colored by antigen label.",
    )


def rounded_box(ax: plt.Axes, xy: tuple[float, float], wh: tuple[float, float], text: str, color: str) -> None:
    box = FancyBboxPatch(
        xy,
        wh[0],
        wh[1],
        boxstyle="round,pad=0.02,rounding_size=0.03",
        linewidth=1.0,
        edgecolor="#333333",
        facecolor=color,
        alpha=0.95,
    )
    ax.add_patch(box)
    ax.text(xy[0] + wh[0] / 2, xy[1] + wh[1] / 2, text, ha="center", va="center", fontsize=9)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12, linewidth=1.1, color="#333333"))


def figure_10_schematic(data: dict[str, object], out_dir: Path) -> FigureRecord:
    summary = data["control_summary"]  # type: ignore[assignment]
    global_mrr = float(metric_row(summary, "global_hl")["mrr"])
    cdr_mrr = float(metric_row(summary, "all_cdrs_mean_std")["mrr"])
    l3_mrr = float(metric_row(summary, "l_cdr3")["mrr"])
    fig = plt.figure(figsize=(11, 4.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.35, 1.0], wspace=0.25)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rounded_box(ax, (0.04, 0.58), (0.19, 0.18), "100\nconformers", "#DCEBFF")
    rounded_box(ax, (0.33, 0.58), (0.18, 0.18), "MCA\nresidue tensors", "#E6F4EA")
    rounded_box(ax, (0.62, 0.72), (0.26, 0.14), "global H/L\npooling", "#FDEBD0")
    rounded_box(ax, (0.62, 0.45), (0.26, 0.14), "CDR/paratope\npooling", "#EADCF8")
    rounded_box(ax, (0.62, 0.18), (0.26, 0.14), "nearest-neighbor\nretrieval", "#F5F5F5")
    arrow(ax, (0.23, 0.67), (0.33, 0.67))
    arrow(ax, (0.51, 0.67), (0.62, 0.79))
    arrow(ax, (0.51, 0.67), (0.62, 0.52))
    arrow(ax, (0.75, 0.72), (0.75, 0.32))
    arrow(ax, (0.75, 0.45), (0.75, 0.32))
    ax.text(0.04, 0.14, "Same tensors, different readout.", fontsize=10)

    ax2 = fig.add_subplot(gs[0, 1])
    labels = ["Global", "L-CDR3", "All CDRs\nmean+std"]
    vals = [global_mrr, l3_mrr, cdr_mrr]
    ax2.bar(labels, vals, color=[FIG_COLORS["global"], FIG_COLORS["light"], FIG_COLORS["purple"]])
    ax2.set_ylim(0, max(vals) * 1.3)
    ax2.set_ylabel("MRR")
    ax2.set_title("Tier B Step 4 full-ensemble MRR")
    for i, val in enumerate(vals):
        ax2.text(i, val + 0.012, f"{val:.3f}", ha="center", va="bottom", fontsize=9)
    set_clean_axis(ax2)
    fig.suptitle("Readout bottleneck: CDR-local signal is easier to retrieve than global average signal", y=1.02)
    fig.tight_layout()
    save_figure(fig, out_dir, "10_readout_bottleneck_schematic")
    return FigureRecord(
        10,
        "10_readout_bottleneck_schematic",
        "Readout bottleneck schematic",
        ["phase15d/region_control_summary.tsv"],
        "Minimal schematic linking conformer ensembles, MCA tensors, pooling choice, and Tier B retrieval MRR.",
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
        handle.write("# Tier B Step 4 Figure Pack\n\n")
        handle.write("Generated from Jerry-owned Tier B Step 4 outputs.\n\n")
        handle.write("## Inputs\n\n")
        handle.write(f"- Tier B Step 4 tables: `{args.input_dir}`\n")
        handle.write(f"- Internal comparison controls: `{args.internal_dir}`\n")
        handle.write(f"- PCA points: `{args.pca_points_table}`\n\n")
        handle.write("## Figures\n\n")
        for record in records:
            handle.write(f"{record.number}. `{record.stem}.png` / `.svg` - {record.note}\n")
        handle.write("\n## Caveat\n\n")
        handle.write(
            "Tier B is exact-PDB-unseen and antibody-sequence-de-overlapped against the checkpoint universe, "
            "but antigen/source de-overlap is not yet enforced. Treat this as cleaner readout evidence, not final external validation.\n"
        )
    with (out_dir / "figure_inputs.json").open("w") as handle:
        json.dump(
            {
                "input_dir": str(args.input_dir),
                "internal_dir": str(args.internal_dir),
                "out_dir": str(args.out_dir),
                "pca_points_table": str(args.pca_points_table),
            },
            handle,
            indent=2,
            sort_keys=True,
        )


def make_contact_sheet(out_dir: Path, records: list[FigureRecord]) -> None:
    images = []
    for record in records:
        path = out_dir / f"{record.stem}.png"
        images.append((record.number, record.title, mpimg.imread(path)))
    ncols = 2
    nrows = math.ceil(len(images) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4.8 * nrows))
    axes_arr = np.asarray(axes).reshape(-1)
    for ax, (number, title, img) in zip(axes_arr, images):
        ax.imshow(img)
        ax.set_title(f"{number:02d}. {title}", fontsize=11)
        ax.axis("off")
    for ax in axes_arr[len(images) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "contact_sheet.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--internal-dir", type=Path, default=DEFAULT_INTERNAL_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--pca-points-table", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    if not args.input_dir.is_absolute():
        args.input_dir = repo_root / args.input_dir
    if not args.internal_dir.is_absolute():
        args.internal_dir = repo_root / args.internal_dir
    if not args.out_dir.is_absolute():
        args.out_dir = repo_root / args.out_dir
    if not args.pca_points_table.is_absolute():
        args.pca_points_table = repo_root / args.pca_points_table
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

    data = load_inputs(args.input_dir, args.internal_dir)
    records = [
        figure_01_ladder(data, args.out_dir),
        figure_02_rank_shift(data, args.out_dir),
        figure_03_first_neighbor_cdf(data, args.out_dir),
        figure_04_replication_heatmap(data, args.out_dir),
        figure_05_controls(data, args.out_dir),
        figure_06_compression(data, args.out_dir),
        figure_07_endpoint_funnel(data, args.out_dir, repo_root),
        figure_08_per_antigen_heatmap(data, args.out_dir),
        figure_09_pca(data, args.out_dir, args.pca_points_table),
        figure_10_schematic(data, args.out_dir),
    ]
    write_manifest(args.out_dir, records, args)
    make_contact_sheet(args.out_dir, records)
    print(f"[tierb-figures] wrote {len(records)} figures to {args.out_dir}")
    for record in records:
        print(f"{record.number:02d} {record.stem}")


if __name__ == "__main__":
    main()
