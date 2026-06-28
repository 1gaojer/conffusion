#!/usr/bin/env python3
"""Analyze Phase 1.4 CDR geometry to MCA embedding sensitivity.

This script reads Jerry-owned Phase 1.4 embedding outputs and staged copied
PH/AF3 CIFs. It reuses the Phase 1.3 CDR structural definitions so CDR-H3 and
all-CDR geometry are measured consistently with the previous result.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import random
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

import numpy as np
import torch
from scipy import stats


K_VALUES = (1, 2, 4, 8, 16, 32, 64, 128)
RANDOM_SEEDS = (0, 1, 2, 3, 4)
EMBEDDING_MODES = ("hl_concat", "h_global", "l_global", "h_cdr3", "all_cdr_concat")
PRIMARY_EMBEDDING_MODES = ("hl_concat", "h_cdr3", "all_cdr_concat")
CDR_MODULE_CACHE = {}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def cached_cdr_module(path: Path):
    key = str(path.resolve())
    if key not in CDR_MODULE_CACHE:
        CDR_MODULE_CACHE[key] = load_module(path.resolve(), f"phase14_cdr_structural_{os.getpid()}")
    return CDR_MODULE_CACHE[key]


def upper_triangle_values(matrix: np.ndarray) -> np.ndarray:
    if matrix.shape[0] < 2:
        return np.array([], dtype=np.float64)
    return matrix[np.triu_indices(matrix.shape[0], k=1)]


def safe_float(value: float) -> float | str:
    if value is None or not np.isfinite(value):
        return ""
    return float(value)


def cosine_distance_between(a: np.ndarray, b: np.ndarray) -> float:
    an = float(np.linalg.norm(a))
    bn = float(np.linalg.norm(b))
    if an == 0.0 or bn == 0.0:
        return float("nan")
    return float(1.0 - np.dot(a, b) / (an * bn))


def pairwise_cosine_distance(vectors: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=np.float64)
    norms = np.linalg.norm(vectors, axis=1)
    denom = np.outer(norms, norms)
    sim = np.zeros((vectors.shape[0], vectors.shape[0]), dtype=np.float64)
    valid = denom > 0
    sim[valid] = (vectors @ vectors.T)[valid] / denom[valid]
    sim = np.clip(sim, -1.0, 1.0)
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    return dist


def pairwise_euclidean_distance(vectors: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=np.float64)
    sq = np.sum(vectors * vectors, axis=1)
    dist_sq = np.maximum(sq[:, None] + sq[None, :] - 2.0 * (vectors @ vectors.T), 0.0)
    dist = np.sqrt(dist_sq)
    np.fill_diagonal(dist, 0.0)
    return dist


def evenly_spaced_indices(n: int, k: int) -> list[int]:
    if k >= n:
        return list(range(n))
    if k == 1:
        return [0]
    return sorted({round(i * (n - 1) / (k - 1)) for i in range(k)})[:k]


def greedy_kcenter_indices(dist: np.ndarray, k: int) -> list[int]:
    n = dist.shape[0]
    if k >= n:
        return list(range(n))
    medoid = int(np.argmin(dist.mean(axis=1)))
    selected = [medoid]
    nearest = dist[:, medoid].copy()
    while len(selected) < k:
        next_idx = int(np.argmax(nearest))
        if next_idx in selected:
            break
        selected.append(next_idx)
        nearest = np.minimum(nearest, dist[:, next_idx])
    return selected


def coverage_stats(dist: np.ndarray | None, selected: list[int], prefix: str) -> dict[str, object]:
    if dist is None or dist.shape[0] == 0 or not selected:
        return {
            f"{prefix}_coverage_mean_nearest": "",
            f"{prefix}_coverage_p90_nearest": "",
            f"{prefix}_coverage_max_nearest": "",
        }
    nearest = dist[:, selected].min(axis=1)
    return {
        f"{prefix}_coverage_mean_nearest": float(nearest.mean()),
        f"{prefix}_coverage_p90_nearest": float(np.quantile(nearest, 0.90)),
        f"{prefix}_coverage_max_nearest": float(nearest.max()),
    }


def select_indices(
    *,
    strategy: str,
    k: int,
    n: int,
    target_index: int,
    seed: int,
    structural: dict[str, np.ndarray],
    embedding_dist: dict[tuple[str, str], np.ndarray],
) -> list[int]:
    if strategy == "first":
        return list(range(k))
    if strategy == "evenly_spaced":
        return evenly_spaced_indices(n, k)
    if strategy == "random":
        rng = random.Random(seed + target_index * 1000 + k)
        return sorted(rng.sample(range(n), k))
    if strategy == "h_cdr3_kcenter":
        return greedy_kcenter_indices(structural["h_cdr3_frame_ca"], k)
    if strategy == "all_cdr_kcenter":
        return greedy_kcenter_indices(structural["all_cdr_frame_ca"], k)
    if strategy == "embedding_hl_kcenter":
        return greedy_kcenter_indices(embedding_dist[("hl_concat", "cosine")], k)
    raise ValueError(f"Unknown strategy: {strategy}")


def region_indices(regions: dict[str, object], names: list[str], max_len: int) -> list[int]:
    out: list[int] = []
    for name in names:
        region = regions.get(name)
        if region is None:
            continue
        for idx in getattr(region, "seq_indices"):
            idx = int(idx)
            if 0 <= idx < max_len:
                out.append(idx)
    return sorted(set(out))


def mean_region(arr: np.ndarray, indices: list[int]) -> np.ndarray:
    if not indices:
        return arr.mean(axis=1)
    return arr[:, indices, :].mean(axis=1)


def load_embedding_vectors(
    *,
    shard_path: Path,
    target_id: str,
    target_meta: dict[str, str],
    cdr_mod,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    shard = torch.load(shard_path, map_location="cpu")
    chains = shard["chain_representations"]
    h_repr = chains["H"]["mca_repr"].float().numpy()
    l_repr = chains["L"]["mca_repr"].float().numpy()
    regions = cdr_mod.target_region_defs(target_id, target_meta, "imgt")

    h_mean = h_repr.mean(axis=1)
    l_mean = l_repr.mean(axis=1)
    h_cdr3_idx = region_indices(regions, ["H-CDR3"], h_repr.shape[1])
    h_cdr_idx = region_indices(regions, ["H-CDR1", "H-CDR2", "H-CDR3"], h_repr.shape[1])
    l_cdr_idx = region_indices(regions, ["L-CDR1", "L-CDR2", "L-CDR3"], l_repr.shape[1])
    h_cdr_mean = mean_region(h_repr, h_cdr_idx)
    l_cdr_mean = mean_region(l_repr, l_cdr_idx)
    vectors = {
        "h_global": h_mean,
        "l_global": l_mean,
        "hl_concat": np.concatenate([h_mean, l_mean], axis=1),
        "h_cdr3": mean_region(h_repr, h_cdr3_idx),
        "all_cdr_concat": np.concatenate([h_cdr_mean, l_cdr_mean], axis=1),
    }
    meta = {
        "source_checkpoint": shard.get("source_checkpoint", ""),
        "n_embedding_conformers": int(h_repr.shape[0]),
        "h_len": int(h_repr.shape[1]),
        "l_len": int(l_repr.shape[1]),
        "h_cdr3_embedding_residues": len(h_cdr3_idx),
        "h_cdr_embedding_residues": len(h_cdr_idx),
        "l_cdr_embedding_residues": len(l_cdr_idx),
    }
    return vectors, meta


def make_region_groups(cdr_mod, target_id: str, regions: dict[str, object]) -> dict[str, list[str]]:
    groups = cdr_mod.make_region_groups(target_id, regions)
    groups["paired-VHVL"] = [
        "H-FR1",
        "H-CDR1",
        "H-FR2",
        "H-CDR2",
        "H-FR3",
        "H-CDR3",
        "H-FR4",
        "L-FR1",
        "L-CDR1",
        "L-FR2",
        "L-CDR2",
        "L-FR3",
        "L-CDR3",
        "L-FR4",
    ]
    return groups


def build_structural_matrices(
    *,
    dataset_root: Path,
    target_id: str,
    target_index: int,
    rows: list[dict[str, str]],
    target_meta: dict[str, str],
    cdr_mod,
) -> tuple[dict[str, np.ndarray], list[int], list[dict[str, object]], dict[str, object]]:
    regions = cdr_mod.target_region_defs(target_id, target_meta, "imgt")
    groups = make_region_groups(cdr_mod, target_id, regions)
    assigned = [cdr_mod.parse_and_assign(dataset_root, row, target_meta) for row in rows]
    mapped = [
        cdr_mod.map_conformer_regions(conf, target_meta, regions, groups)
        for conf in assigned
        if conf.assignment_ok
    ]
    valid_pairs = [
        (idx, conf)
        for idx, conf in enumerate(mapped)
        if conf.role_map_score.get("heavy", 0.0) >= cdr_mod.CHAIN_MAP_THRESHOLD
        and conf.role_map_score.get("light", 0.0) >= cdr_mod.CHAIN_MAP_THRESHOLD
    ]
    valid_positions = [idx for idx, _conf in valid_pairs]
    mapped_ok = [conf for _idx, conf in valid_pairs]
    mapping_rows = []
    for idx, conf in valid_pairs:
        mapping_rows.append(
            {
                "target_id": target_id,
                "target_index": target_index,
                "conformer_index": rows[idx].get("conformer_index", idx),
                "pick_rank": conf.assigned.pick_rank,
                "heavy_chain": conf.assigned.heavy_chain,
                "light_chain": conf.assigned.light_chain,
                "heavy_score": conf.assigned.heavy_score,
                "light_score": conf.assigned.light_score,
                "heavy_map_score": conf.role_map_score.get("heavy", 0.0),
                "light_map_score": conf.role_map_score.get("light", 0.0),
                "mapping_error_count": len(conf.mapping_errors),
                "mapping_errors": "; ".join(conf.mapping_errors[:4]),
            }
        )

    specs = [
        ("h_cdr3_frame_ca", "H-CDR3", "H-FR"),
        ("h_cdrs_frame_ca", "H-CDRs", "H-FR"),
        ("l_cdrs_frame_ca", "L-CDRs", "L-FR"),
        ("all_cdr_frame_ca", "all-CDRs", "all-FR"),
        ("paired_vhvl_frame_ca", "paired-VHVL", "all-FR"),
    ]
    matrices: dict[str, np.ndarray] = {}
    for name, region, align_region in specs:
        region_arrays, region_keys = cdr_mod.common_region_arrays(mapped_ok, region)
        frame_arrays, frame_keys = cdr_mod.common_region_arrays(mapped_ok, align_region)
        if len(region_arrays) < 2 or len(region_keys) == 0 or len(frame_keys) < 3:
            continue
        matrices[name] = cdr_mod.frame_aligned_distance_matrix(region_arrays, frame_arrays, "ca")

    summary = {
        "target_id": target_id,
        "target_index": target_index,
        "n_manifest_conformers": len(rows),
        "n_assignment_ok": sum(1 for conf in assigned if conf.assignment_ok),
        "n_mapped_ok": len(mapped_ok),
        "available_structural_matrices": ",".join(sorted(matrices)),
    }
    return matrices, valid_positions, mapping_rows, summary


def target_strata(cdr_summary: Path, target_ids: list[str]) -> dict[str, dict[str, object]]:
    allowed = set(target_ids)
    rows = []
    for row in read_tsv(cdr_summary):
        if (
            row.get("target_id") in allowed
            and row.get("metric") == "frame_aligned_ca"
            and row.get("region") == "H-CDR3"
        ):
            rows.append((float(row["pairwise_rmsd_mean"]), row["target_id"], row))
    rows.sort()
    out: dict[str, dict[str, object]] = {}
    if not rows:
        return out
    for rank, (value, target_id, row) in enumerate(rows):
        frac = rank / max(len(rows) - 1, 1)
        if frac < 1 / 3:
            stratum = "low"
        elif frac < 2 / 3:
            stratum = "medium"
        else:
            stratum = "high"
        out[target_id] = {
            "target_id": target_id,
            "h_cdr3_diversity_rank": rank + 1,
            "h_cdr3_diversity_stratum": stratum,
            "h_cdr3_pairwise_rmsd_mean": value,
            "h_cdr3_pairwise_rmsd_p90": row.get("pairwise_rmsd_p90", ""),
            "h_cdr3_pairwise_rmsd_max": row.get("pairwise_rmsd_max", ""),
        }
    return out


def correlation_rows_for_target(
    *,
    target_id: str,
    stratum: str,
    structural: dict[str, np.ndarray],
    embedding_dist: dict[tuple[str, str], np.ndarray],
) -> list[dict[str, object]]:
    rows = []
    for structural_name, smat in structural.items():
        svals = upper_triangle_values(smat)
        if len(svals) < 3 or float(np.std(svals)) == 0.0:
            continue
        for (mode, distance_metric), emat in embedding_dist.items():
            evals = upper_triangle_values(emat)
            if len(evals) != len(svals) or float(np.std(evals)) == 0.0:
                continue
            spearman = stats.spearmanr(svals, evals)
            pearson = stats.pearsonr(svals, evals)
            rows.append(
                {
                    "target_id": target_id,
                    "h_cdr3_diversity_stratum": stratum,
                    "structural_metric": structural_name,
                    "embedding_mode": mode,
                    "embedding_distance": distance_metric,
                    "n_conformers": smat.shape[0],
                    "n_pairs": len(svals),
                    "structural_mean": float(svals.mean()),
                    "embedding_mean": float(evals.mean()),
                    "spearman_r": safe_float(float(spearman.statistic)),
                    "spearman_p": safe_float(float(spearman.pvalue)),
                    "pearson_r": safe_float(float(pearson.statistic)),
                    "pearson_p": safe_float(float(pearson.pvalue)),
                }
            )
    return rows


def subset_rows_for_target(
    *,
    target_id: str,
    target_index: int,
    stratum: str,
    structural: dict[str, np.ndarray],
    vectors: dict[str, np.ndarray],
    embedding_dist: dict[tuple[str, str], np.ndarray],
) -> list[dict[str, object]]:
    required = ["h_cdr3_frame_ca", "all_cdr_frame_ca"]
    if any(name not in structural for name in required):
        return []
    n = next(iter(structural.values())).shape[0]
    rows = []
    strategies = ["first", "evenly_spaced", "h_cdr3_kcenter", "all_cdr_kcenter", "embedding_hl_kcenter"]
    for k in [k for k in K_VALUES if k <= n]:
        strategy_seed_pairs = [(strategy, 0) for strategy in strategies]
        strategy_seed_pairs.extend(("random", seed) for seed in RANDOM_SEEDS)
        for strategy, seed in strategy_seed_pairs:
            selected = select_indices(
                strategy=strategy,
                k=k,
                n=n,
                target_index=target_index,
                seed=seed,
                structural=structural,
                embedding_dist=embedding_dist,
            )
            row: dict[str, object] = {
                "target_id": target_id,
                "target_index": target_index,
                "h_cdr3_diversity_stratum": stratum,
                "strategy": strategy,
                "k": k,
                "replicate": seed,
                "selected_n": len(selected),
                "selected_indices": ",".join(str(i) for i in selected),
            }
            row.update(coverage_stats(structural.get("h_cdr3_frame_ca"), selected, "h_cdr3"))
            row.update(coverage_stats(structural.get("all_cdr_frame_ca"), selected, "all_cdr"))
            row.update(coverage_stats(embedding_dist.get(("hl_concat", "cosine")), selected, "hl_embedding_cosine"))
            for mode in PRIMARY_EMBEDDING_MODES:
                full_mean = vectors[mode].mean(axis=0)
                subset_mean = vectors[mode][selected].mean(axis=0)
                row[f"{mode}_mean_cosine_distance"] = safe_float(
                    cosine_distance_between(full_mean, subset_mean)
                )
                row[f"{mode}_mean_euclidean_distance"] = float(
                    np.linalg.norm(full_mean - subset_mean)
                )
            rows.append(row)
    return rows


def control_rows_for_target(
    *,
    target_id: str,
    stratum: str,
    structural: dict[str, np.ndarray],
    vectors: dict[str, np.ndarray],
) -> list[dict[str, object]]:
    rows = []
    n = vectors["hl_concat"].shape[0]
    h_medoid = 0
    if "h_cdr3_frame_ca" in structural:
        h_medoid = int(np.argmin(structural["h_cdr3_frame_ca"].mean(axis=1)))
    cases = [("permutation_mean", None), ("repeat_first_128", 0), ("repeat_hcdr3_medoid_128", h_medoid)]
    for mode in PRIMARY_EMBEDDING_MODES:
        full_mean = vectors[mode].mean(axis=0)
        for control, idx in cases:
            if idx is None:
                control_mean = full_mean.copy()
            else:
                control_mean = vectors[mode][int(idx)]
            rows.append(
                {
                    "target_id": target_id,
                    "h_cdr3_diversity_stratum": stratum,
                    "control": control,
                    "embedding_mode": mode,
                    "base_conformer_index": "" if idx is None else int(idx),
                    "mean_cosine_distance": safe_float(cosine_distance_between(full_mean, control_mean)),
                    "mean_euclidean_distance": float(np.linalg.norm(full_mean - control_mean)),
                }
            )
        for factor in (2, 10, 100):
            weighted = (vectors[mode].sum(axis=0) + (factor - 1) * n * vectors[mode][0]) / (factor * n)
            rows.append(
                {
                    "target_id": target_id,
                    "h_cdr3_diversity_stratum": stratum,
                    "control": f"duplicate_first_{factor}x_weight",
                    "embedding_mode": mode,
                    "base_conformer_index": 0,
                    "mean_cosine_distance": safe_float(cosine_distance_between(full_mean, weighted)),
                    "mean_euclidean_distance": float(np.linalg.norm(full_mean - weighted)),
                }
            )
    return rows


def summarize_correlations(correlation_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    for row in correlation_rows:
        value = row.get("spearman_r", "")
        if value == "":
            continue
        grouped[
            (
                str(row["structural_metric"]),
                str(row["embedding_mode"]),
                str(row["embedding_distance"]),
                str(row["h_cdr3_diversity_stratum"]),
            )
        ].append(float(value))
    out = []
    for (structural_metric, mode, distance_metric, stratum), vals in sorted(grouped.items()):
        arr = np.asarray(vals, dtype=np.float64)
        out.append(
            {
                "structural_metric": structural_metric,
                "embedding_mode": mode,
                "embedding_distance": distance_metric,
                "h_cdr3_diversity_stratum": stratum,
                "n_targets": len(vals),
                "spearman_r_mean": float(arr.mean()),
                "spearman_r_median": float(np.median(arr)),
                "spearman_r_p10": float(np.quantile(arr, 0.10)),
                "spearman_r_p90": float(np.quantile(arr, 0.90)),
                "fraction_positive": float((arr > 0).mean()),
            }
        )
    return out


def summarize_subsets(subset_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = [
        "h_cdr3_coverage_mean_nearest",
        "all_cdr_coverage_mean_nearest",
        "hl_embedding_cosine_coverage_mean_nearest",
        "hl_concat_mean_cosine_distance",
        "h_cdr3_mean_cosine_distance",
        "all_cdr_concat_mean_cosine_distance",
    ]
    grouped: dict[tuple[str, int, str], list[dict[str, object]]] = defaultdict(list)
    for row in subset_rows:
        grouped[(str(row["strategy"]), int(row["k"]), str(row["h_cdr3_diversity_stratum"]))].append(row)
    out = []
    for (strategy, k, stratum), rows in sorted(grouped.items()):
        payload: dict[str, object] = {
            "strategy": strategy,
            "k": k,
            "h_cdr3_diversity_stratum": stratum,
            "n_rows": len(rows),
            "n_targets": len({str(row["target_id"]) for row in rows}),
        }
        for metric in metrics:
            vals = [float(row[metric]) for row in rows if row.get(metric, "") != ""]
            if vals:
                payload[f"{metric}_median"] = float(np.median(vals))
                payload[f"{metric}_mean"] = float(np.mean(vals))
        out.append(payload)
    return out


def make_plots(out_dir: Path, correlation_rows: list[dict[str, object]], subset_summary: list[dict[str, object]]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    primary = [
        row
        for row in correlation_rows
        if row["structural_metric"] in {"h_cdr3_frame_ca", "all_cdr_frame_ca"}
        and row["embedding_mode"] in {"hl_concat", "h_cdr3", "all_cdr_concat"}
        and row["embedding_distance"] == "cosine"
        and row.get("spearman_r", "") != ""
    ]
    labels = sorted({(row["structural_metric"], row["embedding_mode"]) for row in primary})
    if labels:
        plt.figure(figsize=(10, 4))
        data = [
            [float(row["spearman_r"]) for row in primary if (row["structural_metric"], row["embedding_mode"]) == label]
            for label in labels
        ]
        plt.boxplot(data, labels=[f"{a}\n{b}" for a, b in labels], showfliers=False)
        plt.ylabel("Per-target Spearman r")
        plt.xticks(rotation=30, ha="right")
        plt.tight_layout()
        plt.savefig(fig_dir / "correlation_boxplot_primary.png", dpi=180)
        plt.close()

    rows = [
        row
        for row in subset_summary
        if row["h_cdr3_diversity_stratum"] == "high"
        and row["strategy"] in {"first", "random", "h_cdr3_kcenter", "all_cdr_kcenter", "embedding_hl_kcenter"}
    ]
    if rows:
        plt.figure(figsize=(7, 4))
        for strategy in sorted({row["strategy"] for row in rows}):
            xs, ys = [], []
            for k in K_VALUES:
                vals = [
                    row.get("h_cdr3_coverage_mean_nearest_median", "")
                    for row in rows
                    if row["strategy"] == strategy and int(row["k"]) == k
                ]
                vals = [float(v) for v in vals if v != ""]
                if vals:
                    xs.append(k)
                    ys.append(float(np.median(vals)))
            if xs:
                plt.plot(xs, ys, marker="o", label=strategy)
        plt.xscale("log", base=2)
        plt.xlabel("Selected conformers per target (K)")
        plt.ylabel("Median H-CDR3 mean-nearest RMSD")
        plt.title("High H-CDR3 targets")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(fig_dir / "subset_hcdr3_coverage_high_targets.png", dpi=180)
        plt.close()


def write_report(
    *,
    out_dir: Path,
    args: argparse.Namespace,
    run_summary: dict[str, object],
    correlation_summary: list[dict[str, object]],
    subset_summary: list[dict[str, object]],
) -> None:
    primary_corr = [
        row
        for row in correlation_summary
        if row["structural_metric"] == "h_cdr3_frame_ca"
        and row["embedding_distance"] == "cosine"
        and row["h_cdr3_diversity_stratum"] in {"low", "medium", "high"}
        and row["embedding_mode"] in {"hl_concat", "h_cdr3", "all_cdr_concat"}
    ]
    primary_subset = [
        row
        for row in subset_summary
        if row["h_cdr3_diversity_stratum"] == "high"
        and int(row["k"]) in {8, 16, 32, 64}
        and row["strategy"] in {"first", "random", "h_cdr3_kcenter", "all_cdr_kcenter", "embedding_hl_kcenter"}
    ]
    lines = [
        "# Aim 1 Phase 1.4 Embedding Sensitivity Report",
        "",
        f"Created: {run_summary['created_at_utc']}",
        "",
        "## Inputs",
        "",
        f"- Run root: `{args.run_root}`",
        f"- Dataset root: `{args.dataset_root}`",
        f"- CDR summary: `{args.cdr_summary}`",
        f"- CDR script: `{args.cdr_script}`",
        "",
        "## Run Summary",
        "",
        f"- Scope label: `{args.scope_label}`",
        f"- Targets requested: {run_summary['n_targets_requested']}",
        f"- Targets analyzed: {run_summary['n_targets_analyzed']}",
        f"- Targets with H-CDR3 structural matrix: {run_summary['targets_with_h_cdr3_matrix']}",
        f"- Targets with all-CDR structural matrix: {run_summary['targets_with_all_cdr_matrix']}",
        f"- Elapsed seconds: {run_summary['elapsed_seconds']}",
        "",
        "## Primary H-CDR3 To Embedding Correlations",
        "",
        "| Embedding mode | Stratum | Targets | Median Spearman r | Mean Spearman r | Fraction positive |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in primary_corr:
        lines.append(
            f"| {row['embedding_mode']} | {row['h_cdr3_diversity_stratum']} | {row['n_targets']} | "
            f"{float(row['spearman_r_median']):.3f} | {float(row['spearman_r_mean']):.3f} | "
            f"{float(row['fraction_positive']):.2f} |"
        )
    lines.extend(
        [
            "",
            "## High-Stratum Subset Snapshot",
            "",
            "| Strategy | K | Rows | Median H-CDR3 coverage | Median all-CDR coverage | Median HL mean cosine error |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in primary_subset:
        lines.append(
            f"| {row['strategy']} | {row['k']} | {row['n_rows']} | "
            f"{float(row.get('h_cdr3_coverage_mean_nearest_median', math.nan)):.3f} | "
            f"{float(row.get('all_cdr_coverage_mean_nearest_median', math.nan)):.3f} | "
            f"{float(row.get('hl_concat_mean_cosine_distance_median', math.nan)):.6f} |"
        )
    lines.extend(
        [
            "",
            "## Output Files",
            "",
            "- `input_validation.json`",
            "- `target_embedding_summary.tsv`",
            "- `target_strata.tsv`",
            "- `structural_mapping_summary.tsv`",
            "- `cdr_to_embedding_correlations.tsv`",
            "- `correlation_summary.tsv`",
            "- `subset_embedding_preservation.tsv`",
            "- `subset_strategy_summary.tsv`",
            "- `order_duplication_controls.tsv`",
            "- `figures/`",
            "",
            "Interpretation should stay conservative: this tests whether current MCA embeddings are sensitive to generated CDR geometry. It does not by itself prove antigen-retrieval benefit.",
            "",
        ]
    )
    (out_dir / "aim1_phase1_4_embedding_sensitivity_report.md").write_text("\n".join(lines))


def process_target_task(payload: dict[str, object]) -> dict[str, object]:
    run_root = Path(str(payload["run_root"]))
    dataset_root = Path(str(payload["dataset_root"]))
    cdr_script = Path(str(payload["cdr_script"]))
    target_id = str(payload["target_id"])
    target_index = int(payload["target_index"])
    n_targets = int(payload["n_targets"])
    target_meta = payload["target_meta"]
    rows = payload["rows"]
    stratum_row = payload["stratum_row"]
    if not isinstance(target_meta, dict) or not isinstance(rows, list) or not isinstance(stratum_row, dict):
        raise TypeError("Invalid target payload")

    cdr_mod = cached_cdr_module(cdr_script)
    shard_path = run_root / "embeddings" / "mca1000_selected" / f"{target_id}.pt"
    if not shard_path.exists():
        return {
            "target_id": target_id,
            "missing_embedding_shard": target_id,
            "manifest_count_mismatch": None,
            "target_summary_rows": [],
            "mapping_rows": [],
            "correlation_rows": [],
            "subset_rows": [],
            "control_rows": [],
            "targets_analyzed": 0,
            "targets_with_h": 0,
            "targets_with_all": 0,
            "log": f"[{target_index}/{n_targets}] {target_id}: missing embedding shard",
        }

    vectors_all, emb_meta = load_embedding_vectors(
        shard_path=shard_path,
        target_id=target_id,
        target_meta=target_meta,
        cdr_mod=cdr_mod,
    )
    mismatch = None
    if vectors_all["hl_concat"].shape[0] != len(rows):
        mismatch = {
            "target_id": target_id,
            "embedding_n": int(vectors_all["hl_concat"].shape[0]),
            "manifest_n": len(rows),
        }
    structural, valid_positions, target_mapping_rows, structural_summary = build_structural_matrices(
        dataset_root=dataset_root,
        target_id=target_id,
        target_index=target_index,
        rows=rows,
        target_meta=target_meta,
        cdr_mod=cdr_mod,
    )
    if not valid_positions:
        return {
            "target_id": target_id,
            "missing_embedding_shard": None,
            "manifest_count_mismatch": mismatch,
            "target_summary_rows": [],
            "mapping_rows": target_mapping_rows,
            "correlation_rows": [],
            "subset_rows": [],
            "control_rows": [],
            "targets_analyzed": 0,
            "targets_with_h": 0,
            "targets_with_all": 0,
            "log": f"[{target_index}/{n_targets}] {target_id}: no valid mapped conformers",
        }

    vectors = {mode: value[valid_positions] for mode, value in vectors_all.items()}
    embedding_dist = {}
    for mode in EMBEDDING_MODES:
        embedding_dist[(mode, "cosine")] = pairwise_cosine_distance(vectors[mode])
        embedding_dist[(mode, "euclidean")] = pairwise_euclidean_distance(vectors[mode])

    stratum = str(stratum_row.get("h_cdr3_diversity_stratum", "unknown"))
    target_summary = {
        **structural_summary,
        **emb_meta,
        "h_cdr3_diversity_stratum": stratum,
        "h_cdr3_pairwise_rmsd_mean_phase13": stratum_row.get("h_cdr3_pairwise_rmsd_mean", ""),
    }
    return {
        "target_id": target_id,
        "missing_embedding_shard": None,
        "manifest_count_mismatch": mismatch,
        "target_summary_rows": [target_summary],
        "mapping_rows": target_mapping_rows,
        "correlation_rows": correlation_rows_for_target(
            target_id=target_id,
            stratum=stratum,
            structural=structural,
            embedding_dist=embedding_dist,
        ),
        "subset_rows": subset_rows_for_target(
            target_id=target_id,
            target_index=target_index,
            stratum=stratum,
            structural=structural,
            vectors=vectors,
            embedding_dist=embedding_dist,
        ),
        "control_rows": control_rows_for_target(
            target_id=target_id,
            stratum=stratum,
            structural=structural,
            vectors=vectors,
        ),
        "targets_analyzed": 1,
        "targets_with_h": 1 if "h_cdr3_frame_ca" in structural else 0,
        "targets_with_all": 1 if "all_cdr_frame_ca" in structural else 0,
        "log": (
            f"[{target_index}/{n_targets}] {target_id}: "
            f"mapped={structural_summary['n_mapped_ok']} "
            f"matrices={structural_summary['available_structural_matrices']}"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--dataset-root", required=True, type=Path)
    parser.add_argument("--cdr-summary", required=True, type=Path)
    parser.add_argument("--cdr-script", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--scope-label", default="full")
    parser.add_argument("--target-ids-file", type=Path)
    parser.add_argument("--max-targets", type=int)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    start = time.time()
    run_root = args.run_root.resolve()
    dataset_root = args.dataset_root.resolve()
    out_dir = args.out_dir.resolve()
    if out_dir.exists() and any(out_dir.iterdir()) and not args.force:
        raise SystemExit(f"Output dir exists and is not empty: {out_dir}; use --force")
    out_dir.mkdir(parents=True, exist_ok=True)

    target_rows = read_tsv(dataset_root / "manifests" / "medium_targets.tsv")
    target_meta = {row["target_id"]: row for row in target_rows}
    selected_rows = read_tsv(run_root / "manifests" / "selected_conformers.tsv")
    by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in selected_rows:
        by_target[row["target_id"]].append(row)
    for rows in by_target.values():
        rows.sort(key=lambda row: int(row.get("conformer_index") or row.get("pick_rank") or 0))

    if args.target_ids_file:
        target_ids = [
            line.strip()
            for line in args.target_ids_file.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
    else:
        target_ids = [row["target_id"] for row in target_rows if row["target_id"] in by_target]
    if args.max_targets:
        target_ids = target_ids[: args.max_targets]

    strata = target_strata(args.cdr_summary, target_ids)
    write_tsv(
        out_dir / "target_strata.tsv",
        [strata[target_id] for target_id in target_ids if target_id in strata],
        [
            "target_id",
            "h_cdr3_diversity_rank",
            "h_cdr3_diversity_stratum",
            "h_cdr3_pairwise_rmsd_mean",
            "h_cdr3_pairwise_rmsd_p90",
            "h_cdr3_pairwise_rmsd_max",
        ],
    )

    validation = {
        "run_root": str(run_root),
        "dataset_root": str(dataset_root),
        "embedding_dir": str(run_root / "embeddings" / "mca1000_selected"),
        "manifest": str(run_root / "manifests" / "selected_conformers.tsv"),
        "n_targets_requested": len(target_ids),
        "missing_embedding_shards": [],
        "manifest_count_mismatches": [],
    }

    target_summary_rows = []
    mapping_rows = []
    correlation_rows = []
    subset_rows = []
    control_rows = []
    targets_with_h = 0
    targets_with_all = 0
    targets_analyzed = 0

    tasks = [
        {
            "run_root": str(run_root),
            "dataset_root": str(dataset_root),
            "cdr_script": str(args.cdr_script.resolve()),
            "target_id": target_id,
            "target_index": target_index,
            "n_targets": len(target_ids),
            "target_meta": target_meta[target_id],
            "rows": by_target.get(target_id, []),
            "stratum_row": strata.get(target_id, {}),
        }
        for target_index, target_id in enumerate(target_ids, start=1)
    ]

    worker_count = max(1, min(int(args.workers), len(tasks) or 1))
    if worker_count == 1:
        results = []
        for task in tasks:
            result = process_target_task(task)
            results.append(result)
            print(result["log"], flush=True)
    else:
        results = []
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            futures = [executor.submit(process_target_task, task) for task in tasks]
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                print(result["log"], flush=True)

    results.sort(key=lambda result: target_ids.index(str(result["target_id"])))
    for result in results:
        if result["missing_embedding_shard"]:
            validation["missing_embedding_shards"].append(result["missing_embedding_shard"])
        if result["manifest_count_mismatch"]:
            validation["manifest_count_mismatches"].append(result["manifest_count_mismatch"])
        target_summary_rows.extend(result["target_summary_rows"])
        mapping_rows.extend(result["mapping_rows"])
        correlation_rows.extend(result["correlation_rows"])
        subset_rows.extend(result["subset_rows"])
        control_rows.extend(result["control_rows"])
        targets_analyzed += int(result["targets_analyzed"])
        targets_with_h += int(result["targets_with_h"])
        targets_with_all += int(result["targets_with_all"])

    validation["ok"] = not validation["missing_embedding_shards"] and not validation["manifest_count_mismatches"]
    (out_dir / "input_validation.json").write_text(json.dumps(validation, indent=2))

    correlation_summary = summarize_correlations(correlation_rows)
    subset_summary = summarize_subsets(subset_rows)
    run_summary = {
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "scope_label": args.scope_label,
        "n_targets_requested": len(target_ids),
        "n_targets_analyzed": targets_analyzed,
        "targets_with_h_cdr3_matrix": targets_with_h,
        "targets_with_all_cdr_matrix": targets_with_all,
        "elapsed_seconds": round(time.time() - start, 3),
    }
    (out_dir / "run_summary.json").write_text(json.dumps(run_summary, indent=2))

    write_tsv(
        out_dir / "target_embedding_summary.tsv",
        target_summary_rows,
        [
            "target_id",
            "target_index",
            "h_cdr3_diversity_stratum",
            "h_cdr3_pairwise_rmsd_mean_phase13",
            "n_manifest_conformers",
            "n_assignment_ok",
            "n_mapped_ok",
            "n_embedding_conformers",
            "h_len",
            "l_len",
            "h_cdr3_embedding_residues",
            "h_cdr_embedding_residues",
            "l_cdr_embedding_residues",
            "source_checkpoint",
            "available_structural_matrices",
        ],
    )
    write_tsv(
        out_dir / "structural_mapping_summary.tsv",
        mapping_rows,
        [
            "target_id",
            "target_index",
            "conformer_index",
            "pick_rank",
            "heavy_chain",
            "light_chain",
            "heavy_score",
            "light_score",
            "heavy_map_score",
            "light_map_score",
            "mapping_error_count",
            "mapping_errors",
        ],
    )
    write_tsv(
        out_dir / "cdr_to_embedding_correlations.tsv",
        correlation_rows,
        [
            "target_id",
            "h_cdr3_diversity_stratum",
            "structural_metric",
            "embedding_mode",
            "embedding_distance",
            "n_conformers",
            "n_pairs",
            "structural_mean",
            "embedding_mean",
            "spearman_r",
            "spearman_p",
            "pearson_r",
            "pearson_p",
        ],
    )
    write_tsv(
        out_dir / "correlation_summary.tsv",
        correlation_summary,
        [
            "structural_metric",
            "embedding_mode",
            "embedding_distance",
            "h_cdr3_diversity_stratum",
            "n_targets",
            "spearman_r_mean",
            "spearman_r_median",
            "spearman_r_p10",
            "spearman_r_p90",
            "fraction_positive",
        ],
    )
    subset_fieldnames = [
        "target_id",
        "target_index",
        "h_cdr3_diversity_stratum",
        "strategy",
        "k",
        "replicate",
        "selected_n",
        "h_cdr3_coverage_mean_nearest",
        "h_cdr3_coverage_p90_nearest",
        "h_cdr3_coverage_max_nearest",
        "all_cdr_coverage_mean_nearest",
        "all_cdr_coverage_p90_nearest",
        "all_cdr_coverage_max_nearest",
        "hl_embedding_cosine_coverage_mean_nearest",
        "hl_embedding_cosine_coverage_p90_nearest",
        "hl_embedding_cosine_coverage_max_nearest",
        "hl_concat_mean_cosine_distance",
        "hl_concat_mean_euclidean_distance",
        "h_cdr3_mean_cosine_distance",
        "h_cdr3_mean_euclidean_distance",
        "all_cdr_concat_mean_cosine_distance",
        "all_cdr_concat_mean_euclidean_distance",
        "selected_indices",
    ]
    write_tsv(out_dir / "subset_embedding_preservation.tsv", subset_rows, subset_fieldnames)
    write_tsv(
        out_dir / "subset_strategy_summary.tsv",
        subset_summary,
        sorted({key for row in subset_summary for key in row}),
    )
    write_tsv(
        out_dir / "order_duplication_controls.tsv",
        control_rows,
        [
            "target_id",
            "h_cdr3_diversity_stratum",
            "control",
            "embedding_mode",
            "base_conformer_index",
            "mean_cosine_distance",
            "mean_euclidean_distance",
        ],
    )
    make_plots(out_dir, correlation_rows, subset_summary)
    write_report(
        out_dir=out_dir,
        args=args,
        run_summary=run_summary,
        correlation_summary=correlation_summary,
        subset_summary=subset_summary,
    )
    print("[done] " + json.dumps(run_summary, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
