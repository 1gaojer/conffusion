#!/usr/bin/env python3
"""Phase 1.5b failure/random-control analysis for the Phase 1.5a run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_PHASE15A_RUN = Path(
    "/project/liulab/jg1920/conffusion/phase15a_20260627_175558_endpoint_audit"
)
DEFAULT_CDR_DIR = Path(
    "/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627"
)
DEFAULT_OUT = Path("/project/liulab/jg1920/conffusion/phase15b_20260627_failure_analysis")


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    if df.empty:
        return "No rows."
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in df.iterrows():
        vals = [str(row.get(c, "")).replace("\n", " ") for c in df.columns]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def pearson(x: pd.Series, y: pd.Series) -> float:
    data = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if len(data) < 3 or data["x"].nunique() < 2 or data["y"].nunique() < 2:
        return float("nan")
    return float(data["x"].corr(data["y"], method="pearson"))


def spearman(x: pd.Series, y: pd.Series) -> float:
    data = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if len(data) < 3 or data["x"].nunique() < 2 or data["y"].nunique() < 2:
        return float("nan")
    return float(data["x"].rank().corr(data["y"].rank(), method="pearson"))


def summarize_conditions(results: pd.DataFrame, preservation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for condition, group in results.groupby("condition", sort=False):
        valid = group[group["n_positive_after_self_exclusion"] > 0]
        pres = preservation[preservation["condition"] == condition]
        row: dict[str, object] = {
            "condition": condition,
            "family": group["family"].iloc[0],
            "k": group["k"].iloc[0],
            "replicate": group["replicate"].iloc[0],
            "n_queries": len(group),
            "n_eval": len(valid),
            "recall_at_1": bool_series(valid["recall_at_1"]).mean() if not valid.empty else np.nan,
            "recall_at_5": bool_series(valid["recall_at_5"]).mean() if not valid.empty else np.nan,
            "recall_at_10": bool_series(valid["recall_at_10"]).mean() if not valid.empty else np.nan,
            "mrr": pd.to_numeric(valid["reciprocal_rank"], errors="coerce").mean(),
            "median_first_positive_rank": pd.to_numeric(valid["first_positive_rank"], errors="coerce").median(),
            "mean_cosine_to_full": pd.to_numeric(pres["cosine_to_full"], errors="coerce").mean(),
            "median_cosine_to_full": pd.to_numeric(pres["cosine_to_full"], errors="coerce").median(),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def random_rollup(condition_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    random_rows = condition_summary[condition_summary["family"] == "random"].copy()
    for k, group in random_rows.groupby("k", sort=True):
        row: dict[str, object] = {"random_family": f"random_k{int(k)}", "k": int(k), "n_replicates": len(group)}
        for metric in ["recall_at_1", "recall_at_5", "recall_at_10", "mrr", "median_first_positive_rank", "mean_cosine_to_full"]:
            vals = pd.to_numeric(group[metric], errors="coerce").dropna()
            row[f"{metric}_mean"] = float(vals.mean()) if len(vals) else np.nan
            row[f"{metric}_sd"] = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
            row[f"{metric}_min"] = float(vals.min()) if len(vals) else np.nan
            row[f"{metric}_max"] = float(vals.max()) if len(vals) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def condition_deltas(results: pd.DataFrame, preservation: pd.DataFrame) -> pd.DataFrame:
    full = results[results["condition"] == "full_128"][
        ["target_id", "first_positive_rank", "reciprocal_rank", "recall_at_1", "recall_at_5", "recall_at_10"]
    ].copy()
    full.columns = [
        "target_id",
        "full_first_positive_rank",
        "full_reciprocal_rank",
        "full_recall_at_1",
        "full_recall_at_5",
        "full_recall_at_10",
    ]
    merged = results.merge(full, on="target_id", how="left")
    merged = merged.merge(preservation, on=["target_id", "condition"], how="left")
    for col in ["first_positive_rank", "reciprocal_rank", "full_first_positive_rank", "full_reciprocal_rank"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")
    for col in ["recall_at_1", "recall_at_5", "recall_at_10", "full_recall_at_1", "full_recall_at_5", "full_recall_at_10"]:
        merged[col] = bool_series(merged[col])
    merged["rank_delta_vs_full"] = merged["first_positive_rank"] - merged["full_first_positive_rank"]
    merged["rr_delta_vs_full"] = merged["reciprocal_rank"] - merged["full_reciprocal_rank"]
    merged["lost_recall10_vs_full"] = merged["full_recall_at_10"] & ~merged["recall_at_10"]
    merged["gained_recall10_vs_full"] = ~merged["full_recall_at_10"] & merged["recall_at_10"]
    return merged


def add_structural_features(deltas: pd.DataFrame, cdr_dir: Path) -> pd.DataFrame:
    whole = read_tsv(cdr_dir / "whole_vs_cdr_summary.tsv")
    keep = [
        "target_id",
        "h_cdr3_pairwise_mean",
        "h_cdr3_pairwise_p90",
        "h_cdr3_pairwise_max",
        "all_cdr_pairwise_mean",
        "all_cdr_pairwise_p90",
        "all_cdr_pairwise_max",
    ]
    return deltas.merge(whole[keep], on="target_id", how="left")


def selector_vs_random(deltas: pd.DataFrame) -> pd.DataFrame:
    rows = []
    valid = deltas[deltas["n_positive_after_self_exclusion"] > 0].copy()
    random = valid[valid["family"] == "random"].copy()
    random["k"] = pd.to_numeric(random["k"], errors="coerce")
    random_by_target = (
        random.groupby(["target_id", "k"], as_index=False)
        .agg(
            random_rr_mean=("reciprocal_rank", "mean"),
            random_rr_sd=("reciprocal_rank", "std"),
            random_recall10_mean=("recall_at_10", "mean"),
            random_rank_median=("first_positive_rank", "median"),
            random_cosine_mean=("cosine_to_full", "mean"),
        )
    )
    selector_names = [
        "h_cdr3_greedy_kcenter_k32",
        "h_cdr3_greedy_kcenter_k64",
        "all_cdrs_greedy_kcenter_k32",
        "all_cdrs_greedy_kcenter_k64",
        "first_k32",
        "first_k64",
    ]
    selectors = valid[valid["condition"].isin(selector_names)].copy()
    selectors["k"] = pd.to_numeric(selectors["k"], errors="coerce")
    merged = selectors.merge(random_by_target, on=["target_id", "k"], how="left")
    merged["rr_minus_random_mean"] = merged["reciprocal_rank"] - merged["random_rr_mean"]
    merged["recall10_minus_random_mean"] = merged["recall_at_10"].astype(float) - merged["random_recall10_mean"]
    for condition, group in merged.groupby("condition", sort=False):
        rows.append(
            {
                "condition": condition,
                "k": int(group["k"].iloc[0]),
                "n_targets": len(group),
                "mean_rr_minus_random": group["rr_minus_random_mean"].mean(),
                "median_rr_minus_random": group["rr_minus_random_mean"].median(),
                "frac_rr_above_random_mean": (group["rr_minus_random_mean"] > 0).mean(),
                "mean_recall10_minus_random": group["recall10_minus_random_mean"].mean(),
                "frac_recall10_above_random_mean": (group["recall10_minus_random_mean"] > 0).mean(),
            }
        )
    return pd.DataFrame(rows), merged


def correlation_summary(feature_deltas: pd.DataFrame) -> pd.DataFrame:
    features = [
        "h_cdr3_pairwise_mean",
        "h_cdr3_pairwise_p90",
        "h_cdr3_pairwise_max",
        "all_cdr_pairwise_mean",
        "all_cdr_pairwise_p90",
        "all_cdr_pairwise_max",
        "cosine_to_full",
    ]
    outcomes = [
        "rr_delta_vs_full",
        "rank_delta_vs_full",
        "lost_recall10_vs_full",
    ]
    conditions = [
        "single_first",
        "h_cdr3_greedy_kcenter_k32",
        "h_cdr3_greedy_kcenter_k64",
        "all_cdrs_greedy_kcenter_k32",
        "all_cdrs_greedy_kcenter_k64",
        "first_k32",
        "first_k64",
    ]
    rows = []
    valid = feature_deltas[feature_deltas["n_positive_after_self_exclusion"] > 0].copy()
    valid["lost_recall10_vs_full"] = valid["lost_recall10_vs_full"].astype(float)
    for condition in conditions:
        group = valid[valid["condition"] == condition]
        for feature in features:
            for outcome in outcomes:
                rows.append(
                    {
                        "condition": condition,
                        "feature": feature,
                        "outcome": outcome,
                        "n": group[[feature, outcome]].dropna().shape[0],
                        "pearson": pearson(group[feature], group[outcome]),
                        "spearman": spearman(group[feature], group[outcome]),
                    }
                )
    return pd.DataFrame(rows)


def write_report(
    path: Path,
    condition_summary: pd.DataFrame,
    random_summary: pd.DataFrame,
    selector_summary: pd.DataFrame,
    corr: pd.DataFrame,
    audit: dict[str, object],
) -> None:
    full = condition_summary[condition_summary["condition"] == "full_128"].iloc[0]
    single = condition_summary[condition_summary["condition"] == "single_first"].iloc[0]
    rows = [
        "# Phase 1.5b Failure And Selector Analysis",
        "",
        "## Inputs",
        "",
        f"- Phase 1.5a run: `{audit['phase15a_run']}`",
        f"- CDR structural directory: `{audit['cdr_dir']}`",
        "",
        "## Main Readout",
        "",
        f"- Full 128 conformers: Recall@10 `{full['recall_at_10']:.3f}`, MRR `{full['mrr']:.3f}`.",
        f"- Single first conformer: Recall@10 `{single['recall_at_10']:.3f}`, MRR `{single['mrr']:.3f}`.",
        "- Random controls were condensed across five replicates before interpretation.",
        "",
        "## Condition Summary",
        "",
        markdown_table(
            condition_summary[
                [
                    "condition",
                    "n_eval",
                    "recall_at_10",
                    "mrr",
                    "median_first_positive_rank",
                    "mean_cosine_to_full",
                ]
            ]
        ),
        "",
        "## Random-Control Rollup",
        "",
        markdown_table(random_summary),
        "",
        "## Selector Versus Random Mean",
        "",
        markdown_table(selector_summary),
        "",
        "## Strongest Absolute Spearman Correlations",
        "",
        markdown_table(
            corr.assign(abs_spearman=lambda d: d["spearman"].abs())
            .sort_values("abs_spearman", ascending=False)
            .drop(columns=["abs_spearman"])
            .head(12)
        ),
        "",
        "## Interpretation",
        "",
        "- Full ensemble retrieval is meaningfully better than a single conformer.",
        "- K=64 generally outperforms K=32, but CDR structural k-center does not cleanly beat random K=64 in the global mean-pooled endpoint.",
        "- The near-1.0 cosine preservation values show the global pooled vector barely changes, yet retrieval ranks still move; small vector changes can matter near decision boundaries.",
        "- Treat this as evidence for a readout/selector problem, not yet as proof that conformer generation itself should be optimized by CDR k-center alone.",
        "",
    ]
    path.write_text("\n".join(rows))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase15a-run", type=Path, default=DEFAULT_PHASE15A_RUN)
    parser.add_argument("--cdr-dir", type=Path, default=DEFAULT_CDR_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    results = read_tsv(args.phase15a_run / "retrieval_smoke_results.tsv")
    preservation = read_tsv(args.phase15a_run / "subset_vector_preservation.tsv")
    for col in ["n_positive_after_self_exclusion", "k", "replicate", "selected_n"]:
        results[col] = pd.to_numeric(results[col], errors="coerce")

    condition_summary = summarize_conditions(results, preservation)
    random_summary = random_rollup(condition_summary)
    deltas = condition_deltas(results, preservation)
    feature_deltas = add_structural_features(deltas, args.cdr_dir)
    selector_summary, selector_target = selector_vs_random(deltas)
    corr = correlation_summary(feature_deltas)

    audit = {
        "phase15a_run": str(args.phase15a_run),
        "cdr_dir": str(args.cdr_dir),
        "out_dir": str(args.out_dir),
    }
    (args.out_dir / "phase15b_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True))
    condition_summary.to_csv(args.out_dir / "condition_summary.tsv", sep="\t", index=False)
    random_summary.to_csv(args.out_dir / "random_control_rollup.tsv", sep="\t", index=False)
    deltas.to_csv(args.out_dir / "target_condition_deltas.tsv", sep="\t", index=False)
    feature_deltas.to_csv(args.out_dir / "target_failure_features.tsv", sep="\t", index=False)
    selector_summary.to_csv(args.out_dir / "selector_vs_random_summary.tsv", sep="\t", index=False)
    selector_target.to_csv(args.out_dir / "selector_vs_random_by_target.tsv", sep="\t", index=False)
    corr.to_csv(args.out_dir / "loss_correlation_summary.tsv", sep="\t", index=False)
    write_report(
        args.out_dir / "phase15b_failure_analysis_report.md",
        condition_summary,
        random_summary,
        selector_summary,
        corr,
        audit,
    )
    print(f"[phase15b] wrote {args.out_dir}")
    print(condition_summary[["condition", "n_eval", "recall_at_10", "mrr", "mean_cosine_to_full"]].to_string(index=False))
    print("\n[random rollup]")
    print(random_summary.to_string(index=False))
    print("\n[selector vs random]")
    print(selector_summary.to_string(index=False))


if __name__ == "__main__":
    main()
