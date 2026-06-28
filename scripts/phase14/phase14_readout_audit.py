#!/usr/bin/env python3
"""Audit richer MCA readouts against CDR structural diversity.

This is a CPU-only follow-up to the Phase 1.4 embedding-sensitivity analysis.
It reads Jerry-owned exported MCA shards and copied PH/AF3 conformer manifests,
then tests whether CDR/paratope signal is more visible in richer readouts than
in the simple full-H/L average.

The current exported shards contain final per-conformer ``mca_repr`` and final
chain-local ``pair_repr``. They do not contain layer-wise conformer weights or
subset-specific pair representations; those require a later GPU re-encode.
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


def configure_torch_threads() -> None:
    threads = int(os.environ.get("PHASE14_TORCH_THREADS", "1"))
    try:
        torch.set_num_threads(max(1, threads))
    except RuntimeError:
        pass
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass


configure_torch_threads()


K_VALUES = (8, 16, 32, 64, 128)
RANDOM_SEEDS = (0, 1, 2, 3, 4)
PRIMARY_READOUTS = (
    "hl_global_mean",
    "h_cdr3_mean",
    "h_cdr3_flat",
    "h_cdr3_mean_std",
    "all_cdr_mean_concat",
    "all_cdr_flat_concat",
    "all_cdr_mean_std",
)
STRUCTURAL_METRICS = (
    "h_cdr3_frame_ca",
    "h_cdrs_frame_ca",
    "l_cdrs_frame_ca",
    "all_cdr_frame_ca",
    "paired_vhvl_frame_ca",
)

BASE_MODULE = None


def load_base_module():
    global BASE_MODULE
    if BASE_MODULE is not None:
        return BASE_MODULE
    path = Path(__file__).resolve().with_name("phase14_analyze_embeddings.py")
    spec = importlib.util.spec_from_file_location("phase14_analyze_embeddings_base", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import base analyzer from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["phase14_analyze_embeddings_base"] = module
    spec.loader.exec_module(module)
    BASE_MODULE = module
    return module


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


def safe_float(value: object) -> float | str:
    if value is None:
        return ""
    try:
        out = float(value)
    except (TypeError, ValueError):
        return ""
    if not np.isfinite(out):
        return ""
    return out


def upper_triangle_values(matrix: np.ndarray) -> np.ndarray:
    if matrix.shape[0] < 2:
        return np.array([], dtype=np.float64)
    return matrix[np.triu_indices(matrix.shape[0], k=1)]


def cosine_distance_between(a: np.ndarray, b: np.ndarray) -> float:
    an = float(np.linalg.norm(a))
    bn = float(np.linalg.norm(b))
    if an == 0.0 or bn == 0.0:
        return float("nan")
    return float(1.0 - np.dot(a, b) / (an * bn))


def vector_region(arr: np.ndarray, indices: list[int]) -> np.ndarray:
    if not indices:
        return arr
    return arr[:, indices, :]


def mean_region(arr: np.ndarray, indices: list[int]) -> np.ndarray:
    return vector_region(arr, indices).mean(axis=1)


def flat_region(arr: np.ndarray, indices: list[int]) -> np.ndarray:
    selected = vector_region(arr, indices)
    return selected.reshape(selected.shape[0], -1)


def mean_std_region(arr: np.ndarray, indices: list[int]) -> np.ndarray:
    selected = vector_region(arr, indices)
    return np.concatenate([selected.mean(axis=1), selected.std(axis=1)], axis=1)


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
    readout_dist: dict[str, np.ndarray],
) -> list[int]:
    if strategy == "first":
        return list(range(k))
    if strategy == "random":
        rng = random.Random(seed + target_index * 1000 + k)
        return sorted(rng.sample(range(n), k))
    if strategy == "h_cdr3_kcenter":
        return greedy_kcenter_indices(structural["h_cdr3_frame_ca"], k)
    if strategy == "all_cdr_kcenter":
        return greedy_kcenter_indices(structural["all_cdr_frame_ca"], k)
    if strategy == "readout_h_cdr3_flat_kcenter":
        return greedy_kcenter_indices(readout_dist["h_cdr3_flat"], k)
    if strategy == "readout_all_cdr_flat_kcenter":
        return greedy_kcenter_indices(readout_dist["all_cdr_flat_concat"], k)
    raise ValueError(f"Unknown strategy: {strategy}")


def nearest_weighted_mean(vectors: np.ndarray, selected: list[int], dist: np.ndarray) -> np.ndarray:
    nearest = np.argmin(dist[:, selected], axis=1)
    counts = np.bincount(nearest, minlength=len(selected)).astype(np.float64)
    weights = counts / counts.sum()
    return np.sum(vectors[selected] * weights[:, None], axis=0)


def pair_block_summary(
    *,
    target_id: str,
    chain: str,
    block: str,
    pair_repr: np.ndarray,
    idx_a: list[int],
    idx_b: list[int],
) -> dict[str, object]:
    row: dict[str, object] = {
        "target_id": target_id,
        "chain": chain,
        "block": block,
        "n_residues_a": len(idx_a),
        "n_residues_b": len(idx_b),
        "n_pairs": len(idx_a) * len(idx_b),
    }
    if not idx_a or not idx_b:
        row.update(
            {
                "block_scalar_mean": "",
                "block_scalar_std": "",
                "block_mean_abs": "",
                "block_mean_l2_norm": "",
                "block_dim": pair_repr.shape[-1],
            }
        )
        return row
    block_arr = pair_repr[np.ix_(idx_a, idx_b)]
    flat = block_arr.reshape(-1, pair_repr.shape[-1]).astype(np.float64)
    mean_vec = flat.mean(axis=0)
    row.update(
        {
            "block_scalar_mean": float(flat.mean()),
            "block_scalar_std": float(flat.std()),
            "block_mean_abs": float(np.abs(flat).mean()),
            "block_mean_l2_norm": float(np.linalg.norm(mean_vec)),
            "block_dim": int(pair_repr.shape[-1]),
        }
    )
    return row


def load_readouts(
    *,
    shard_path: Path,
    target_id: str,
    target_meta: dict[str, str],
    cdr_mod,
) -> tuple[dict[str, np.ndarray], list[dict[str, object]], dict[str, object], dict[str, object]]:
    base = load_base_module()
    shard = torch.load(shard_path, map_location="cpu")
    chains = shard["chain_representations"]
    h_repr = chains["H"]["mca_repr"].float().numpy().astype(np.float64, copy=False)
    l_repr = chains["L"]["mca_repr"].float().numpy().astype(np.float64, copy=False)
    regions = cdr_mod.target_region_defs(target_id, target_meta, "imgt")

    h_cdr3 = base.region_indices(regions, ["H-CDR3"], h_repr.shape[1])
    h_cdrs = base.region_indices(regions, ["H-CDR1", "H-CDR2", "H-CDR3"], h_repr.shape[1])
    l_cdrs = base.region_indices(regions, ["L-CDR1", "L-CDR2", "L-CDR3"], l_repr.shape[1])
    h_fr = base.region_indices(regions, ["H-FR1", "H-FR2", "H-FR3", "H-FR4"], h_repr.shape[1])
    l_fr = base.region_indices(regions, ["L-FR1", "L-FR2", "L-FR3", "L-FR4"], l_repr.shape[1])

    readouts = {
        "h_global_mean": h_repr.mean(axis=1),
        "l_global_mean": l_repr.mean(axis=1),
        "hl_global_mean": np.concatenate([h_repr.mean(axis=1), l_repr.mean(axis=1)], axis=1),
        "h_cdr3_mean": mean_region(h_repr, h_cdr3),
        "h_cdr3_flat": flat_region(h_repr, h_cdr3),
        "h_cdr3_mean_std": mean_std_region(h_repr, h_cdr3),
        "h_cdrs_mean": mean_region(h_repr, h_cdrs),
        "h_cdrs_flat": flat_region(h_repr, h_cdrs),
        "l_cdrs_mean": mean_region(l_repr, l_cdrs),
        "l_cdrs_flat": flat_region(l_repr, l_cdrs),
        "all_cdr_mean_concat": np.concatenate(
            [mean_region(h_repr, h_cdrs), mean_region(l_repr, l_cdrs)], axis=1
        ),
        "all_cdr_flat_concat": np.concatenate(
            [flat_region(h_repr, h_cdrs), flat_region(l_repr, l_cdrs)], axis=1
        ),
        "all_cdr_mean_std": np.concatenate(
            [mean_std_region(h_repr, h_cdrs), mean_std_region(l_repr, l_cdrs)], axis=1
        ),
    }

    pair_rows = []
    for chain, chain_len, cdr3_idx, cdr_idx, fr_idx in (
        ("H", h_repr.shape[1], h_cdr3, h_cdrs, h_fr),
        ("L", l_repr.shape[1], [], l_cdrs, l_fr),
    ):
        pair_repr = chains[chain]["pair_repr"].float().numpy().astype(np.float64, copy=False)
        all_idx = list(range(chain_len))
        block_defs = [
            ("all_all", all_idx, all_idx),
            ("cdrs_cdrs", cdr_idx, cdr_idx),
            ("cdrs_framework", cdr_idx, fr_idx),
            ("framework_framework", fr_idx, fr_idx),
        ]
        if chain == "H":
            block_defs.extend(
                [
                    ("cdr3_cdr3", cdr3_idx, cdr3_idx),
                    ("cdr3_framework", cdr3_idx, fr_idx),
                    ("cdr3_cdrs", cdr3_idx, cdr_idx),
                ]
            )
        for block, idx_a, idx_b in block_defs:
            pair_rows.append(
                pair_block_summary(
                    target_id=target_id,
                    chain=chain,
                    block=block,
                    pair_repr=pair_repr,
                    idx_a=idx_a,
                    idx_b=idx_b,
                )
            )

    sanity = {
        "target_id": target_id,
        "h_len": int(h_repr.shape[1]),
        "l_len": int(l_repr.shape[1]),
        "n_conformers": int(h_repr.shape[0]),
        "h_cdr3_n": len(h_cdr3),
        "h_cdr3_first_index": h_cdr3[0] if h_cdr3 else "",
        "h_cdr3_last_index": h_cdr3[-1] if h_cdr3 else "",
        "h_cdr3_sequence": getattr(regions.get("H-CDR3"), "sequence", ""),
        "h_cdrs_n": len(h_cdrs),
        "l_cdrs_n": len(l_cdrs),
        "source_checkpoint": shard.get("source_checkpoint", ""),
        "schema_version": shard.get("schema_version", ""),
    }
    meta = {
        "source_checkpoint": shard.get("source_checkpoint", ""),
        "schema_version": shard.get("schema_version", ""),
        "h_len": int(h_repr.shape[1]),
        "l_len": int(l_repr.shape[1]),
        "n_embedding_conformers": int(h_repr.shape[0]),
        "readout_dims_json": json.dumps({key: int(value.shape[1]) for key, value in readouts.items()}, sort_keys=True),
    }
    return readouts, pair_rows, sanity, meta


def correlation_rows_for_target(
    *,
    target_id: str,
    stratum: str,
    structural: dict[str, np.ndarray],
    readout_distances: dict[tuple[str, str], np.ndarray],
    readout_dims: dict[str, int],
) -> list[dict[str, object]]:
    rows = []
    for structural_metric in STRUCTURAL_METRICS:
        smat = structural.get(structural_metric)
        if smat is None:
            continue
        svals = upper_triangle_values(smat)
        if len(svals) < 3 or float(np.std(svals)) == 0.0:
            continue
        for (readout_mode, distance_metric), emat in readout_distances.items():
            evals = upper_triangle_values(emat)
            if len(evals) != len(svals) or float(np.std(evals)) == 0.0:
                continue
            spearman = stats.spearmanr(svals, evals)
            pearson = stats.pearsonr(svals, evals)
            rows.append(
                {
                    "target_id": target_id,
                    "h_cdr3_diversity_stratum": stratum,
                    "structural_metric": structural_metric,
                    "readout_mode": readout_mode,
                    "readout_distance": distance_metric,
                    "readout_dim": readout_dims[readout_mode],
                    "n_conformers": smat.shape[0],
                    "n_pairs": len(svals),
                    "structural_mean": float(svals.mean()),
                    "readout_distance_mean": float(evals.mean()),
                    "spearman_r": safe_float(spearman.statistic),
                    "spearman_p": safe_float(spearman.pvalue),
                    "pearson_r": safe_float(pearson.statistic),
                    "pearson_p": safe_float(pearson.pvalue),
                }
            )
    return rows


def subset_rows_for_target(
    *,
    target_id: str,
    target_index: int,
    stratum: str,
    structural: dict[str, np.ndarray],
    readouts: dict[str, np.ndarray],
    readout_dist: dict[str, np.ndarray],
) -> list[dict[str, object]]:
    if "h_cdr3_frame_ca" not in structural or "all_cdr_frame_ca" not in structural:
        return []
    n = structural["h_cdr3_frame_ca"].shape[0]
    strategies = [
        "first",
        "h_cdr3_kcenter",
        "all_cdr_kcenter",
        "readout_h_cdr3_flat_kcenter",
        "readout_all_cdr_flat_kcenter",
    ]
    strategy_seed_pairs = [(strategy, 0) for strategy in strategies]
    for seed in RANDOM_SEEDS:
        strategy_seed_pairs.append(("random", seed))

    rows = []
    for k in [value for value in K_VALUES if value <= n]:
        for strategy, seed in strategy_seed_pairs:
            selected = select_indices(
                strategy=strategy,
                k=k,
                n=n,
                target_index=target_index,
                seed=seed,
                structural=structural,
                readout_dist=readout_dist,
            )
            coverage = {}
            coverage.update(coverage_stats(structural.get("h_cdr3_frame_ca"), selected, "h_cdr3"))
            coverage.update(coverage_stats(structural.get("all_cdr_frame_ca"), selected, "all_cdr"))
            for readout_mode in PRIMARY_READOUTS:
                vectors = readouts[readout_mode]
                full_mean = vectors.mean(axis=0)
                equal_mean = vectors[selected].mean(axis=0)
                weighted_h = nearest_weighted_mean(vectors, selected, structural["h_cdr3_frame_ca"])
                weighted_all = nearest_weighted_mean(vectors, selected, structural["all_cdr_frame_ca"])
                full_var = vectors.var(axis=0)
                subset_var = vectors[selected].var(axis=0)
                row: dict[str, object] = {
                    "target_id": target_id,
                    "target_index": target_index,
                    "h_cdr3_diversity_stratum": stratum,
                    "strategy": strategy,
                    "k": k,
                    "replicate": seed,
                    "selected_n": len(selected),
                    "readout_mode": readout_mode,
                    "readout_dim": int(vectors.shape[1]),
                    "equal_mean_cosine_distance": safe_float(cosine_distance_between(full_mean, equal_mean)),
                    "equal_mean_l2_distance": float(np.linalg.norm(full_mean - equal_mean)),
                    "h_cdr3_weighted_mean_cosine_distance": safe_float(
                        cosine_distance_between(full_mean, weighted_h)
                    ),
                    "h_cdr3_weighted_mean_l2_distance": float(np.linalg.norm(full_mean - weighted_h)),
                    "all_cdr_weighted_mean_cosine_distance": safe_float(
                        cosine_distance_between(full_mean, weighted_all)
                    ),
                    "all_cdr_weighted_mean_l2_distance": float(np.linalg.norm(full_mean - weighted_all)),
                    "variance_cosine_distance": safe_float(cosine_distance_between(full_var, subset_var)),
                    "variance_l2_distance": float(np.linalg.norm(full_var - subset_var)),
                    "selected_indices": ",".join(str(idx) for idx in selected),
                }
                row.update(coverage)
                rows.append(row)
    return rows


def summarize_correlations(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    dims: dict[tuple[str, str, str, str], int] = {}
    for row in rows:
        value = row.get("spearman_r", "")
        if value == "":
            continue
        key = (
            str(row["structural_metric"]),
            str(row["readout_mode"]),
            str(row["readout_distance"]),
            str(row["h_cdr3_diversity_stratum"]),
        )
        grouped[key].append(float(value))
        dims[key] = int(row["readout_dim"])
    out = []
    for key, vals in sorted(grouped.items()):
        structural_metric, readout_mode, distance_metric, stratum = key
        arr = np.asarray(vals, dtype=np.float64)
        out.append(
            {
                "structural_metric": structural_metric,
                "readout_mode": readout_mode,
                "readout_distance": distance_metric,
                "h_cdr3_diversity_stratum": stratum,
                "readout_dim": dims[key],
                "n_targets": len(vals),
                "spearman_r_mean": float(arr.mean()),
                "spearman_r_median": float(np.median(arr)),
                "spearman_r_p10": float(np.quantile(arr, 0.10)),
                "spearman_r_p90": float(np.quantile(arr, 0.90)),
                "fraction_positive": float((arr > 0).mean()),
            }
        )
    return out


def summarize_subsets(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = [
        "h_cdr3_coverage_mean_nearest",
        "all_cdr_coverage_mean_nearest",
        "equal_mean_cosine_distance",
        "h_cdr3_weighted_mean_cosine_distance",
        "all_cdr_weighted_mean_cosine_distance",
        "variance_cosine_distance",
        "variance_l2_distance",
    ]
    grouped: dict[tuple[str, int, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[
            (
                str(row["strategy"]),
                int(row["k"]),
                str(row["readout_mode"]),
                str(row["h_cdr3_diversity_stratum"]),
            )
        ].append(row)
    out = []
    for (strategy, k, readout_mode, stratum), group_rows in sorted(grouped.items()):
        payload: dict[str, object] = {
            "strategy": strategy,
            "k": k,
            "readout_mode": readout_mode,
            "h_cdr3_diversity_stratum": stratum,
            "n_rows": len(group_rows),
            "n_targets": len({str(row["target_id"]) for row in group_rows}),
        }
        for metric in metrics:
            vals = [float(row[metric]) for row in group_rows if row.get(metric, "") != ""]
            if vals:
                payload[f"{metric}_mean"] = float(np.mean(vals))
                payload[f"{metric}_median"] = float(np.median(vals))
                payload[f"{metric}_p90"] = float(np.quantile(vals, 0.90))
        out.append(payload)
    return out


def pair_repr_cross_target_correlations(
    pair_rows: list[dict[str, object]],
    strata: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    metrics = ("block_scalar_std", "block_mean_abs", "block_mean_l2_norm")
    grouped: dict[tuple[str, str, str], list[tuple[float, float]]] = defaultdict(list)
    for row in pair_rows:
        target_id = str(row["target_id"])
        stratum = strata.get(target_id, {})
        diversity = stratum.get("h_cdr3_pairwise_rmsd_mean", "")
        if diversity == "":
            continue
        for metric in metrics:
            value = row.get(metric, "")
            if value == "":
                continue
            grouped[(str(row["chain"]), str(row["block"]), metric)].append((float(diversity), float(value)))
    out = []
    for (chain, block, metric), pairs in sorted(grouped.items()):
        if len(pairs) < 5:
            continue
        xs = np.asarray([x for x, _y in pairs], dtype=np.float64)
        ys = np.asarray([y for _x, y in pairs], dtype=np.float64)
        if float(np.std(xs)) == 0.0 or float(np.std(ys)) == 0.0:
            continue
        spearman = stats.spearmanr(xs, ys)
        pearson = stats.pearsonr(xs, ys)
        out.append(
            {
                "chain": chain,
                "block": block,
                "pair_repr_metric": metric,
                "target_diversity_metric": "h_cdr3_pairwise_rmsd_mean",
                "n_targets": len(pairs),
                "spearman_r": safe_float(spearman.statistic),
                "spearman_p": safe_float(spearman.pvalue),
                "pearson_r": safe_float(pearson.statistic),
                "pearson_p": safe_float(pearson.pvalue),
            }
        )
    return out


def process_target(payload: dict[str, object]) -> dict[str, object]:
    base = load_base_module()
    run_root = Path(str(payload["run_root"]))
    dataset_root = Path(str(payload["dataset_root"]))
    cdr_script = Path(str(payload["cdr_script"]))
    target_id = str(payload["target_id"])
    target_index = int(payload["target_index"])
    n_targets = int(payload["n_targets"])
    target_meta = payload["target_meta"]
    conformer_rows = payload["rows"]
    stratum_row = payload["stratum_row"]
    if not isinstance(target_meta, dict) or not isinstance(conformer_rows, list) or not isinstance(stratum_row, dict):
        raise TypeError("Invalid target payload")

    shard_path = run_root / "embeddings" / "mca1000_selected" / f"{target_id}.pt"
    if not shard_path.exists():
        return {
            "target_id": target_id,
            "missing_embedding_shard": target_id,
            "manifest_count_mismatch": None,
            "region_sanity_rows": [],
            "target_summary_rows": [],
            "pair_repr_rows": [],
            "correlation_rows": [],
            "subset_rows": [],
            "targets_analyzed": 0,
            "targets_with_h": 0,
            "targets_with_all": 0,
            "log": f"[{target_index}/{n_targets}] {target_id}: missing embedding shard",
        }

    cdr_mod = base.cached_cdr_module(cdr_script)
    readouts_all, pair_repr_rows, sanity, readout_meta = load_readouts(
        shard_path=shard_path,
        target_id=target_id,
        target_meta=target_meta,
        cdr_mod=cdr_mod,
    )
    mismatch = None
    if readouts_all["hl_global_mean"].shape[0] != len(conformer_rows):
        mismatch = {
            "target_id": target_id,
            "embedding_n": int(readouts_all["hl_global_mean"].shape[0]),
            "manifest_n": len(conformer_rows),
        }

    structural, valid_positions, mapping_rows, structural_summary = base.build_structural_matrices(
        dataset_root=dataset_root,
        target_id=target_id,
        target_index=target_index,
        rows=conformer_rows,
        target_meta=target_meta,
        cdr_mod=cdr_mod,
    )
    if not valid_positions:
        return {
            "target_id": target_id,
            "missing_embedding_shard": None,
            "manifest_count_mismatch": mismatch,
            "region_sanity_rows": [sanity],
            "target_summary_rows": [],
            "pair_repr_rows": pair_repr_rows,
            "correlation_rows": [],
            "subset_rows": [],
            "targets_analyzed": 0,
            "targets_with_h": 0,
            "targets_with_all": 0,
            "log": f"[{target_index}/{n_targets}] {target_id}: no valid mapped conformers",
        }

    readouts = {mode: value[valid_positions] for mode, value in readouts_all.items()}
    readout_dims = {mode: int(value.shape[1]) for mode, value in readouts.items()}
    readout_distances: dict[tuple[str, str], np.ndarray] = {}
    readout_cosine_for_selection: dict[str, np.ndarray] = {}
    for mode, vectors in readouts.items():
        cosine = pairwise_cosine_distance(vectors)
        readout_distances[(mode, "cosine")] = cosine
        readout_distances[(mode, "euclidean")] = pairwise_euclidean_distance(vectors)
        readout_cosine_for_selection[mode] = cosine

    stratum = str(stratum_row.get("h_cdr3_diversity_stratum", "unknown"))
    target_summary = {
        **structural_summary,
        **readout_meta,
        "h_cdr3_diversity_stratum": stratum,
        "h_cdr3_pairwise_rmsd_mean_phase13": stratum_row.get("h_cdr3_pairwise_rmsd_mean", ""),
        "n_valid_positions": len(valid_positions),
    }
    correlation_rows = correlation_rows_for_target(
        target_id=target_id,
        stratum=stratum,
        structural=structural,
        readout_distances=readout_distances,
        readout_dims=readout_dims,
    )
    subset_rows = subset_rows_for_target(
        target_id=target_id,
        target_index=target_index,
        stratum=stratum,
        structural=structural,
        readouts=readouts,
        readout_dist=readout_cosine_for_selection,
    )

    return {
        "target_id": target_id,
        "missing_embedding_shard": None,
        "manifest_count_mismatch": mismatch,
        "region_sanity_rows": [sanity],
        "target_summary_rows": [target_summary],
        "pair_repr_rows": pair_repr_rows,
        "correlation_rows": correlation_rows,
        "subset_rows": subset_rows,
        "targets_analyzed": 1,
        "targets_with_h": 1 if "h_cdr3_frame_ca" in structural else 0,
        "targets_with_all": 1 if "all_cdr_frame_ca" in structural else 0,
        "log": (
            f"[{target_index}/{n_targets}] {target_id}: "
            f"mapped={structural_summary['n_mapped_ok']} "
            f"h_cdr3_n={sanity['h_cdr3_n']} "
            f"readouts={len(readouts)}"
        ),
    }


def write_report(
    *,
    out_dir: Path,
    args: argparse.Namespace,
    run_summary: dict[str, object],
    validation: dict[str, object],
    correlation_summary: list[dict[str, object]],
    subset_summary: list[dict[str, object]],
    pair_cross_rows: list[dict[str, object]],
) -> None:
    primary_corr = [
        row
        for row in correlation_summary
        if row["structural_metric"] == "h_cdr3_frame_ca"
        and row["readout_distance"] == "cosine"
        and row["readout_mode"] in PRIMARY_READOUTS
    ]
    primary_subset = [
        row
        for row in subset_summary
        if row["k"] in {32, 64}
        and row["readout_mode"] in {"h_cdr3_flat", "all_cdr_flat_concat", "all_cdr_mean_std"}
        and row["strategy"] in {"random", "h_cdr3_kcenter", "all_cdr_kcenter"}
    ]
    top_pair = [
        row
        for row in sorted(
            pair_cross_rows,
            key=lambda r: abs(float(r.get("spearman_r") or 0.0)),
            reverse=True,
        )[:12]
    ]

    lines = [
        "# Aim 1 Phase 1.4 Rich MCA Readout Audit",
        "",
        f"Created: {run_summary['created_at_utc']}",
        "",
        "## Inputs",
        "",
        f"- Run root: `{args.run_root}`",
        f"- Dataset root: `{args.dataset_root}`",
        f"- CDR summary: `{args.cdr_summary}`",
        f"- CDR script: `{args.cdr_script}`",
        f"- Scope: `{args.scope_label}`",
        "",
        "## Run Summary",
        "",
        f"- Targets requested: {run_summary['n_targets_requested']}",
        f"- Targets analyzed: {run_summary['n_targets_analyzed']}",
        f"- Targets with H-CDR3 matrix: {run_summary['targets_with_h_cdr3_matrix']}",
        f"- Targets with all-CDR matrix: {run_summary['targets_with_all_cdr_matrix']}",
        f"- Elapsed seconds: {run_summary['elapsed_seconds']}",
        f"- Validation OK: {validation['ok']}",
        "",
        "## What This Tests",
        "",
        "This tests whether MCA1000 final per-conformer representations encode generated CDR geometry more strongly when read out through CDR-specific, residue-flattened, mean+variance, duplicate-normalized, or cluster-aware summaries.",
        "",
        "It does not test layer-wise conformer weights or subset-specific final pair representations, because those tensors were not saved in the current embedding shards. That requires a later GPU re-encode.",
        "",
        "## Primary H-CDR3 Correlations",
        "",
        "| Readout | Stratum | Targets | Median Spearman r | Mean Spearman r | Fraction positive |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in primary_corr:
        lines.append(
            f"| {row['readout_mode']} | {row['h_cdr3_diversity_stratum']} | {row['n_targets']} | "
            f"{float(row['spearman_r_median']):.3f} | {float(row['spearman_r_mean']):.3f} | "
            f"{float(row['fraction_positive']):.2f} |"
        )
    lines.extend(
        [
            "",
            "## K=32/K=64 Subset Preservation Snapshot",
            "",
            "| Strategy | K | Readout | Stratum | Targets | Median equal-mean cosine error | Median variance cosine error |",
            "|---|---:|---|---|---:|---:|---:|",
        ]
    )
    for row in primary_subset:
        lines.append(
            f"| {row['strategy']} | {row['k']} | {row['readout_mode']} | "
            f"{row['h_cdr3_diversity_stratum']} | {row['n_targets']} | "
            f"{float(row.get('equal_mean_cosine_distance_median', math.nan)):.6f} | "
            f"{float(row.get('variance_cosine_distance_median', math.nan)):.6f} |"
        )
    lines.extend(
        [
            "",
            "## Pair Representation Cross-Target Checks",
            "",
            "| Chain | Block | Metric | Targets | Spearman r vs H-CDR3 diversity |",
            "|---|---|---|---:|---:|",
        ]
    )
    for row in top_pair:
        lines.append(
            f"| {row['chain']} | {row['block']} | {row['pair_repr_metric']} | "
            f"{row['n_targets']} | {float(row['spearman_r']):.3f} |"
        )
    lines.extend(
        [
            "",
            "## Output Files",
            "",
            "- `input_validation.json`",
            "- `run_summary.json`",
            "- `region_sanity.tsv`",
            "- `target_readout_summary.tsv`",
            "- `readout_correlations.tsv`",
            "- `readout_correlation_summary.tsv`",
            "- `subset_readout_preservation.tsv`",
            "- `subset_readout_summary.tsv`",
            "- `pair_repr_region_summary.tsv`",
            "- `pair_repr_cross_target_correlations.tsv`",
            "",
        ]
    )
    (out_dir / "aim1_phase1_4_readout_audit_report.md").write_text("\n".join(lines))


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

    base = load_base_module()
    strata = base.target_strata(args.cdr_summary, target_ids)
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

    validation: dict[str, object] = {
        "run_root": str(run_root),
        "dataset_root": str(dataset_root),
        "embedding_dir": str(run_root / "embeddings" / "mca1000_selected"),
        "manifest": str(run_root / "manifests" / "selected_conformers.tsv"),
        "current_shard_limitation": "mca_shard_v2 contains final mca_repr and pair_repr only; no layer-wise conformer weights or subset-specific pair_repr",
        "n_targets_requested": len(target_ids),
        "missing_embedding_shards": [],
        "manifest_count_mismatches": [],
    }

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

    results = []
    worker_count = max(1, min(int(args.workers), len(tasks) or 1))
    if worker_count == 1:
        for task in tasks:
            result = process_target(task)
            results.append(result)
            print(result["log"], flush=True)
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            futures = [executor.submit(process_target, task) for task in tasks]
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                print(result["log"], flush=True)

    order = {target_id: idx for idx, target_id in enumerate(target_ids)}
    results.sort(key=lambda row: order.get(str(row["target_id"]), 999999))

    region_sanity_rows = []
    target_summary_rows = []
    pair_repr_rows = []
    correlation_rows = []
    subset_rows = []
    targets_analyzed = 0
    targets_with_h = 0
    targets_with_all = 0
    for result in results:
        if result["missing_embedding_shard"]:
            validation["missing_embedding_shards"].append(result["missing_embedding_shard"])
        if result["manifest_count_mismatch"]:
            validation["manifest_count_mismatches"].append(result["manifest_count_mismatch"])
        region_sanity_rows.extend(result["region_sanity_rows"])
        target_summary_rows.extend(result["target_summary_rows"])
        pair_repr_rows.extend(result["pair_repr_rows"])
        correlation_rows.extend(result["correlation_rows"])
        subset_rows.extend(result["subset_rows"])
        targets_analyzed += int(result["targets_analyzed"])
        targets_with_h += int(result["targets_with_h"])
        targets_with_all += int(result["targets_with_all"])

    validation["ok"] = not validation["missing_embedding_shards"] and not validation["manifest_count_mismatches"]
    (out_dir / "input_validation.json").write_text(json.dumps(validation, indent=2, sort_keys=True))

    correlation_summary = summarize_correlations(correlation_rows)
    subset_summary = summarize_subsets(subset_rows)
    pair_cross_rows = pair_repr_cross_target_correlations(pair_repr_rows, strata)

    run_summary = {
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "scope_label": args.scope_label,
        "n_targets_requested": len(target_ids),
        "n_targets_analyzed": targets_analyzed,
        "targets_with_h_cdr3_matrix": targets_with_h,
        "targets_with_all_cdr_matrix": targets_with_all,
        "elapsed_seconds": round(time.time() - start, 3),
    }
    (out_dir / "run_summary.json").write_text(json.dumps(run_summary, indent=2, sort_keys=True))

    write_tsv(
        out_dir / "region_sanity.tsv",
        region_sanity_rows,
        [
            "target_id",
            "h_len",
            "l_len",
            "n_conformers",
            "h_cdr3_n",
            "h_cdr3_first_index",
            "h_cdr3_last_index",
            "h_cdr3_sequence",
            "h_cdrs_n",
            "l_cdrs_n",
            "schema_version",
            "source_checkpoint",
        ],
    )
    write_tsv(
        out_dir / "target_readout_summary.tsv",
        target_summary_rows,
        [
            "target_id",
            "target_index",
            "h_cdr3_diversity_stratum",
            "h_cdr3_pairwise_rmsd_mean_phase13",
            "n_manifest_conformers",
            "n_assignment_ok",
            "n_mapped_ok",
            "n_valid_positions",
            "n_embedding_conformers",
            "h_len",
            "l_len",
            "schema_version",
            "source_checkpoint",
            "available_structural_matrices",
            "readout_dims_json",
        ],
    )
    write_tsv(
        out_dir / "readout_correlations.tsv",
        correlation_rows,
        [
            "target_id",
            "h_cdr3_diversity_stratum",
            "structural_metric",
            "readout_mode",
            "readout_distance",
            "readout_dim",
            "n_conformers",
            "n_pairs",
            "structural_mean",
            "readout_distance_mean",
            "spearman_r",
            "spearman_p",
            "pearson_r",
            "pearson_p",
        ],
    )
    write_tsv(
        out_dir / "readout_correlation_summary.tsv",
        correlation_summary,
        [
            "structural_metric",
            "readout_mode",
            "readout_distance",
            "h_cdr3_diversity_stratum",
            "readout_dim",
            "n_targets",
            "spearman_r_mean",
            "spearman_r_median",
            "spearman_r_p10",
            "spearman_r_p90",
            "fraction_positive",
        ],
    )
    write_tsv(
        out_dir / "subset_readout_preservation.tsv",
        subset_rows,
        [
            "target_id",
            "target_index",
            "h_cdr3_diversity_stratum",
            "strategy",
            "k",
            "replicate",
            "selected_n",
            "readout_mode",
            "readout_dim",
            "h_cdr3_coverage_mean_nearest",
            "h_cdr3_coverage_p90_nearest",
            "h_cdr3_coverage_max_nearest",
            "all_cdr_coverage_mean_nearest",
            "all_cdr_coverage_p90_nearest",
            "all_cdr_coverage_max_nearest",
            "equal_mean_cosine_distance",
            "equal_mean_l2_distance",
            "h_cdr3_weighted_mean_cosine_distance",
            "h_cdr3_weighted_mean_l2_distance",
            "all_cdr_weighted_mean_cosine_distance",
            "all_cdr_weighted_mean_l2_distance",
            "variance_cosine_distance",
            "variance_l2_distance",
            "selected_indices",
        ],
    )
    write_tsv(
        out_dir / "subset_readout_summary.tsv",
        subset_summary,
        sorted({key for row in subset_summary for key in row}),
    )
    write_tsv(
        out_dir / "pair_repr_region_summary.tsv",
        pair_repr_rows,
        [
            "target_id",
            "chain",
            "block",
            "n_residues_a",
            "n_residues_b",
            "n_pairs",
            "block_scalar_mean",
            "block_scalar_std",
            "block_mean_abs",
            "block_mean_l2_norm",
            "block_dim",
        ],
    )
    write_tsv(
        out_dir / "pair_repr_cross_target_correlations.tsv",
        pair_cross_rows,
        [
            "chain",
            "block",
            "pair_repr_metric",
            "target_diversity_metric",
            "n_targets",
            "spearman_r",
            "spearman_p",
            "pearson_r",
            "pearson_p",
        ],
    )
    write_report(
        out_dir=out_dir,
        args=args,
        run_summary=run_summary,
        validation=validation,
        correlation_summary=correlation_summary,
        subset_summary=subset_summary,
        pair_cross_rows=pair_cross_rows,
    )
    print("[done] " + json.dumps(run_summary, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
