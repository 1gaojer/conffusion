#!/usr/bin/env python3
"""Antibody-focused structural analysis for copied PH/AF3 ensembles.

Phase 1.2 validates whether the global diversity observed in Phase 1.1 remains
when measuring only the antibody heavy/light chains.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import json
import math
import random
import shlex
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np


K_VALUES = (1, 2, 4, 8, 16, 32, 64, 128)
RANDOM_SEEDS = (0, 1, 2, 3, 4)
ROLE_SCORE_THRESHOLD = 0.95

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass
class ChainParsed:
    ca_by_key: dict[tuple[str, int], np.ndarray]
    seq_by_chain: dict[str, str]
    chain_counts: dict[str, int]
    chain_centroids: dict[str, np.ndarray]


@dataclass
class AssignedConformer:
    target_id: str
    pick_rank: int
    copied_relative_path: str
    parse_ok: bool
    assignment_ok: bool
    error: str
    heavy_chain: str
    light_chain: str
    heavy_score: float
    light_score: float
    parsed: ChainParsed | None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def find_manifest(dataset_root: Path, suffix: str) -> Path:
    matches = sorted((dataset_root / "manifests").glob(f"*{suffix}.tsv"))
    if len(matches) != 1:
        raise SystemExit(f"Expected exactly one *{suffix}.tsv under {dataset_root / 'manifests'}, found {matches}")
    return matches[0]


def parse_atom_site_chains(path: Path) -> ChainParsed:
    lines = path.read_text(errors="replace").splitlines()
    headers: list[str] = []
    data_start = None
    for i, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        j = i + 1
        loop_headers = []
        while j < len(lines) and lines[j].startswith("_"):
            loop_headers.append(lines[j].strip())
            j += 1
        if loop_headers and loop_headers[0].startswith("_atom_site."):
            headers = loop_headers
            data_start = j
            break
    if not headers or data_start is None:
        raise ValueError("missing _atom_site loop")
    col = {h.replace("_atom_site.", ""): idx for idx, h in enumerate(headers)}
    required = ["label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"]
    missing = [name for name in required if name not in col]
    if missing:
        raise ValueError(f"missing atom_site columns: {','.join(missing)}")

    ca_by_key: dict[tuple[str, int], np.ndarray] = {}
    residues: dict[tuple[str, int], str] = {}
    chain_points: dict[str, list[np.ndarray]] = defaultdict(list)
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped or stripped == "#" or stripped == "loop_" or stripped.startswith("_"):
            break
        parts = shlex.split(stripped)
        if len(parts) < len(headers):
            continue
        if parts[col["label_atom_id"]] != "CA":
            continue
        chain = parts[col["label_asym_id"]]
        seq_raw = parts[col["label_seq_id"]]
        if seq_raw in {".", "?"}:
            continue
        try:
            seq_id = int(seq_raw)
            coord = np.array(
                [
                    float(parts[col["Cartn_x"]]),
                    float(parts[col["Cartn_y"]]),
                    float(parts[col["Cartn_z"]]),
                ],
                dtype=np.float64,
            )
        except ValueError:
            continue
        key = (chain, seq_id)
        if key in ca_by_key:
            continue
        ca_by_key[key] = coord
        residues[key] = AA3_TO_1.get(parts[col["label_comp_id"]], "X")
        chain_points[chain].append(coord)

    if not ca_by_key:
        raise ValueError("no CA coordinates parsed")
    seq_by_chain: dict[str, str] = defaultdict(str)
    for chain, seq_id in sorted(residues):
        seq_by_chain[chain] += residues[(chain, seq_id)]
    chain_counts = {chain: len(points) for chain, points in chain_points.items()}
    chain_centroids = {chain: np.vstack(points).mean(axis=0) for chain, points in chain_points.items()}
    return ChainParsed(
        ca_by_key=ca_by_key,
        seq_by_chain=dict(seq_by_chain),
        chain_counts=chain_counts,
        chain_centroids=chain_centroids,
    )


def seq_score(observed: str, expected: str) -> float:
    if not observed or not expected:
        return 0.0
    if observed == expected:
        return 1.0
    if observed in expected or expected in observed:
        return min(len(observed), len(expected)) / max(len(observed), len(expected))
    return difflib.SequenceMatcher(None, observed, expected).ratio()


def assign_heavy_light(parsed: ChainParsed, vh_seq: str, vl_seq: str) -> tuple[str, str, float, float]:
    chains = sorted(parsed.seq_by_chain)
    if len(chains) < 2:
        raise ValueError("fewer than two chains parsed")
    best: tuple[float, str, str, float, float] | None = None
    for heavy_chain in chains:
        heavy_score = seq_score(parsed.seq_by_chain[heavy_chain], vh_seq)
        for light_chain in chains:
            if light_chain == heavy_chain:
                continue
            light_score = seq_score(parsed.seq_by_chain[light_chain], vl_seq)
            total = heavy_score + light_score
            candidate = (total, heavy_chain, light_chain, heavy_score, light_score)
            if best is None or candidate > best:
                best = candidate
    if best is None:
        raise ValueError("no distinct heavy/light chain pair")
    _total, heavy_chain, light_chain, heavy_score, light_score = best
    if heavy_score < ROLE_SCORE_THRESHOLD or light_score < ROLE_SCORE_THRESHOLD:
        raise ValueError(f"low chain assignment score heavy={heavy_score:.3f} light={light_score:.3f}")
    return heavy_chain, light_chain, heavy_score, light_score


def parse_and_assign(dataset_root: Path, row: dict[str, str], target_meta: dict[str, str]) -> AssignedConformer:
    target_id = row["target_id"]
    pick_rank = int(row.get("pick_rank") or 0)
    copied_relative_path = row["copied_relative_path"]
    path = dataset_root / copied_relative_path
    try:
        parsed = parse_atom_site_chains(path)
        heavy_chain, light_chain, heavy_score, light_score = assign_heavy_light(
            parsed,
            target_meta.get("vh_seq", ""),
            target_meta.get("vl_seq", ""),
        )
        return AssignedConformer(
            target_id=target_id,
            pick_rank=pick_rank,
            copied_relative_path=copied_relative_path,
            parse_ok=True,
            assignment_ok=True,
            error="",
            heavy_chain=heavy_chain,
            light_chain=light_chain,
            heavy_score=heavy_score,
            light_score=light_score,
            parsed=parsed,
        )
    except Exception as exc:  # noqa: BLE001 - keep failure rows for audit.
        return AssignedConformer(
            target_id=target_id,
            pick_rank=pick_rank,
            copied_relative_path=copied_relative_path,
            parse_ok=False,
            assignment_ok=False,
            error=f"{type(exc).__name__}: {exc}",
            heavy_chain="",
            light_chain="",
            heavy_score=0.0,
            light_score=0.0,
            parsed=None,
        )


def kabsch_rmsd(points_a: np.ndarray, points_b: np.ndarray) -> float:
    a = points_a - points_a.mean(axis=0)
    b = points_b - points_b.mean(axis=0)
    cov = a.T @ b
    u, _s, vt = np.linalg.svd(cov)
    correction = np.eye(3)
    correction[2, 2] = np.sign(np.linalg.det(u @ vt))
    rot = u @ correction @ vt
    aligned = a @ rot
    diff = aligned - b
    return float(np.sqrt((diff * diff).sum(axis=1).mean()))


def pairwise_distance_matrix(arrays: list[np.ndarray]) -> np.ndarray:
    n = len(arrays)
    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            d = kabsch_rmsd(arrays[i], arrays[j])
            dist[i, j] = d
            dist[j, i] = d
    return dist


def summarize_distances(dist: np.ndarray) -> dict[str, float]:
    if dist.shape[0] < 2:
        return {
            "pairwise_rmsd_mean": float("nan"),
            "pairwise_rmsd_median": float("nan"),
            "pairwise_rmsd_p90": float("nan"),
            "pairwise_rmsd_max": float("nan"),
        }
    vals = dist[np.triu_indices(dist.shape[0], k=1)]
    return {
        "pairwise_rmsd_mean": float(vals.mean()),
        "pairwise_rmsd_median": float(np.median(vals)),
        "pairwise_rmsd_p90": float(np.quantile(vals, 0.90)),
        "pairwise_rmsd_max": float(vals.max()),
    }


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


def coverage_row(
    target_id: str,
    metric: str,
    strategy: str,
    k: int,
    replicate: str,
    selected: list[int],
    dist: np.ndarray,
) -> dict[str, object]:
    selected = sorted(set(selected))
    nearest = dist[:, selected].min(axis=1)
    return {
        "target_id": target_id,
        "metric": metric,
        "strategy": strategy,
        "k": k,
        "replicate": replicate,
        "selected_n": len(selected),
        "coverage_mean_nearest_rmsd": float(nearest.mean()),
        "coverage_median_nearest_rmsd": float(np.median(nearest)),
        "coverage_p90_nearest_rmsd": float(np.quantile(nearest, 0.90)),
        "coverage_max_nearest_rmsd": float(nearest.max()),
        "selected_indices": ",".join(str(i) for i in selected),
    }


def common_role_arrays(confs: list[AssignedConformer], roles: tuple[str, ...]) -> tuple[list[np.ndarray], int]:
    role_key_sets = []
    for conf in confs:
        assert conf.parsed is not None
        role_to_chain = {"heavy": conf.heavy_chain, "light": conf.light_chain}
        keys = set()
        for role in roles:
            chain = role_to_chain[role]
            keys.update((role, seq_id) for c, seq_id in conf.parsed.ca_by_key if c == chain)
        role_key_sets.append(keys)
    common_keys = set.intersection(*role_key_sets)
    common_keys_sorted = sorted(common_keys)
    arrays = []
    for conf in confs:
        assert conf.parsed is not None
        role_to_chain = {"heavy": conf.heavy_chain, "light": conf.light_chain}
        coords = []
        for role, seq_id in common_keys_sorted:
            coords.append(conf.parsed.ca_by_key[(role_to_chain[role], seq_id)])
        arrays.append(np.vstack(coords))
    return arrays, len(common_keys_sorted)


def orientation_summary(confs: list[AssignedConformer]) -> dict[str, float]:
    distances = []
    unit_vectors = []
    for conf in confs:
        assert conf.parsed is not None
        hc = conf.parsed.chain_centroids[conf.heavy_chain]
        lc = conf.parsed.chain_centroids[conf.light_chain]
        vec = lc - hc
        norm = float(np.linalg.norm(vec))
        if norm == 0:
            continue
        distances.append(norm)
        unit_vectors.append(vec / norm)
    if not distances:
        return {
            "vh_vl_centroid_distance_mean": float("nan"),
            "vh_vl_centroid_distance_std": float("nan"),
            "vh_vl_axis_angle_mean_deg": float("nan"),
            "vh_vl_axis_angle_p90_deg": float("nan"),
            "vh_vl_axis_angle_max_deg": float("nan"),
        }
    angles = []
    for i in range(len(unit_vectors)):
        for j in range(i + 1, len(unit_vectors)):
            cos = float(np.clip(np.dot(unit_vectors[i], unit_vectors[j]), -1.0, 1.0))
            angles.append(math.degrees(math.acos(cos)))
    return {
        "vh_vl_centroid_distance_mean": float(np.mean(distances)),
        "vh_vl_centroid_distance_std": float(np.std(distances)),
        "vh_vl_axis_angle_mean_deg": float(np.mean(angles)) if angles else float("nan"),
        "vh_vl_axis_angle_p90_deg": float(np.quantile(angles, 0.90)) if angles else float("nan"),
        "vh_vl_axis_angle_max_deg": float(np.max(angles)) if angles else float("nan"),
    }


def make_plots(out_dir: Path, saturation_rows: list[dict[str, object]]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    grouped: dict[tuple[str, str, int], list[float]] = defaultdict(list)
    for row in saturation_rows:
        grouped[(str(row["metric"]), str(row["strategy"]), int(row["k"]))].append(float(row["coverage_mean_nearest_rmsd"]))
    for metric in sorted({key[0] for key in grouped}):
        plt.figure(figsize=(7, 4))
        for strategy in sorted({key[1] for key in grouped if key[0] == metric}):
            xs = []
            ys = []
            for k in K_VALUES:
                vals = grouped.get((metric, strategy, k), [])
                if vals:
                    xs.append(k)
                    ys.append(float(np.median(vals)))
            if xs:
                plt.plot(xs, ys, marker="o", label=strategy)
        plt.xscale("log", base=2)
        plt.xlabel("Selected conformers per target (K)")
        plt.ylabel("Median mean-nearest RMSD to selected set (A)")
        plt.title(metric)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(fig_dir / f"saturation_{metric}_median_mean_nearest_rmsd.png", dpi=180)
        plt.close()


def read_global_diversity(global_results_dir: Path | None) -> dict[str, dict[str, str]]:
    if global_results_dir is None:
        return {}
    path = global_results_dir / "target_diversity_summary.tsv"
    if not path.exists():
        return {}
    rows = read_tsv(path)
    return {row["target_id"]: row for row in rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--global-results-dir", type=Path)
    parser.add_argument("--max-targets", type=int)
    parser.add_argument("--max-conformers-per-target", type=int)
    args = parser.parse_args()

    dataset_root = args.dataset_root.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    targets_path = find_manifest(dataset_root, "targets")
    conformers_path = find_manifest(dataset_root, "conformers")
    target_rows = read_tsv(targets_path)
    conformer_rows = read_tsv(conformers_path)
    target_meta = {row["target_id"]: row for row in target_rows}
    target_order = [row["target_id"] for row in target_rows]
    if args.max_targets is not None:
        target_order = target_order[: args.max_targets]
    target_order_set = set(target_order)
    by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in conformer_rows:
        if row["target_id"] in target_order_set:
            by_target[row["target_id"]].append(row)
    for rows in by_target.values():
        rows.sort(key=lambda r: int(r.get("pick_rank") or 0))
        if args.max_conformers_per_target is not None:
            del rows[args.max_conformers_per_target :]

    global_div = read_global_diversity(args.global_results_dir.resolve() if args.global_results_dir else None)

    assignment_rows = []
    target_assignment_rows = []
    diversity_rows = []
    saturation_rows = []
    compare_rows = []
    random.seed(0)

    for target_index, target_id in enumerate(target_order, start=1):
        rows = by_target[target_id]
        assigned = [parse_and_assign(dataset_root, row, target_meta[target_id]) for row in rows]
        ok = [conf for conf in assigned if conf.assignment_ok]
        heavy_counts = Counter(conf.heavy_chain for conf in ok)
        light_counts = Counter(conf.light_chain for conf in ok)
        for conf in assigned:
            parsed = conf.parsed
            chain_lengths = ""
            chain_scores = ""
            if parsed is not None:
                chain_lengths = ";".join(f"{chain}:{len(seq)}" for chain, seq in sorted(parsed.seq_by_chain.items()))
                vh = target_meta[target_id].get("vh_seq", "")
                vl = target_meta[target_id].get("vl_seq", "")
                chain_scores = ";".join(
                    f"{chain}:H={seq_score(seq, vh):.3f},L={seq_score(seq, vl):.3f}"
                    for chain, seq in sorted(parsed.seq_by_chain.items())
                )
            assignment_rows.append(
                {
                    "target_id": target_id,
                    "pick_rank": conf.pick_rank,
                    "copied_relative_path": conf.copied_relative_path,
                    "parse_ok": str(conf.parse_ok).lower(),
                    "assignment_ok": str(conf.assignment_ok).lower(),
                    "assignment_error": conf.error,
                    "heavy_chain": conf.heavy_chain,
                    "light_chain": conf.light_chain,
                    "heavy_score": conf.heavy_score,
                    "light_score": conf.light_score,
                    "chain_lengths": chain_lengths,
                    "chain_scores": chain_scores,
                }
            )
        target_assignment_rows.append(
            {
                "target_id": target_id,
                "manifest_conformers": len(rows),
                "assignment_ok_conformers": len(ok),
                "assignment_failed_conformers": len(assigned) - len(ok),
                "heavy_chain_counts": ";".join(f"{k}:{v}" for k, v in sorted(heavy_counts.items())),
                "light_chain_counts": ";".join(f"{k}:{v}" for k, v in sorted(light_counts.items())),
                "stable_heavy_chain": heavy_counts.most_common(1)[0][0] if heavy_counts else "",
                "stable_light_chain": light_counts.most_common(1)[0][0] if light_counts else "",
                "vh_seq_len": len(target_meta[target_id].get("vh_seq", "")),
                "vl_seq_len": len(target_meta[target_id].get("vl_seq", "")),
            }
        )
        if len(ok) < 2:
            continue

        orient = orientation_summary(ok)
        metric_distances: dict[str, np.ndarray] = {}
        metric_common_counts: dict[str, int] = {}
        for metric, roles in [
            ("heavy_ca", ("heavy",)),
            ("light_ca", ("light",)),
            ("antibody_ca", ("heavy", "light")),
        ]:
            arrays, n_common = common_role_arrays(ok, roles)
            if n_common < 3:
                continue
            dist = pairwise_distance_matrix(arrays)
            metric_distances[metric] = dist
            metric_common_counts[metric] = n_common
            medoid_idx = int(np.argmin(dist.mean(axis=1)))
            diversity_rows.append(
                {
                    "target_id": target_id,
                    "target_index": target_index,
                    "metric": metric,
                    "n_conformers": len(ok),
                    "common_ca_count": n_common,
                    "medoid_pick_rank": ok[medoid_idx].pick_rank,
                    "medoid_mean_distance": float(dist[medoid_idx].mean()),
                    **summarize_distances(dist),
                    **orient,
                }
            )

            for k in [k for k in K_VALUES if k <= len(ok)]:
                saturation_rows.append(coverage_row(target_id, metric, "first", k, "0", list(range(k)), dist))
                saturation_rows.append(coverage_row(target_id, metric, "evenly_spaced", k, "0", evenly_spaced_indices(len(ok), k), dist))
                saturation_rows.append(coverage_row(target_id, metric, "greedy_kcenter", k, "0", greedy_kcenter_indices(dist, k), dist))
                for seed in RANDOM_SEEDS:
                    rng = random.Random(seed + target_index * 1000 + k)
                    selected = sorted(rng.sample(range(len(ok)), k))
                    saturation_rows.append(coverage_row(target_id, metric, "random", k, str(seed), selected, dist))

        if target_id in global_div and "antibody_ca" in metric_distances:
            antibody_summary = summarize_distances(metric_distances["antibody_ca"])
            global_mean = float(global_div[target_id]["pairwise_rmsd_mean"])
            antibody_mean = antibody_summary["pairwise_rmsd_mean"]
            compare_rows.append(
                {
                    "target_id": target_id,
                    "global_pairwise_rmsd_mean": global_mean,
                    "antibody_pairwise_rmsd_mean": antibody_mean,
                    "antibody_minus_global": antibody_mean - global_mean,
                    "antibody_div_global": antibody_mean / global_mean if global_mean else float("nan"),
                    "global_common_ca_count": global_div[target_id].get("common_ca_count", ""),
                    "antibody_common_ca_count": metric_common_counts.get("antibody_ca", ""),
                }
            )

        print(f"processed target {target_index}/{len(target_order)} {target_id} assignment_ok={len(ok)}", flush=True)

    aggregate_rows = []
    grouped: dict[tuple[str, str, int], list[dict[str, object]]] = defaultdict(list)
    for row in saturation_rows:
        grouped[(str(row["metric"]), str(row["strategy"]), int(row["k"]))].append(row)
    for (metric, strategy, k), rows in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1], item[0][2])):
        vals = [float(row["coverage_mean_nearest_rmsd"]) for row in rows]
        targets = {str(row["target_id"]) for row in rows}
        aggregate_rows.append(
            {
                "metric": metric,
                "strategy": strategy,
                "k": k,
                "n_rows": len(rows),
                "n_targets": len(targets),
                "mean_coverage_mean_nearest_rmsd": float(np.mean(vals)),
                "median_coverage_mean_nearest_rmsd": float(np.median(vals)),
                "p90_coverage_mean_nearest_rmsd": float(np.quantile(vals, 0.90)),
            }
        )

    write_tsv(
        out_dir / "chain_role_assignments.tsv",
        assignment_rows,
        [
            "target_id",
            "pick_rank",
            "copied_relative_path",
            "parse_ok",
            "assignment_ok",
            "assignment_error",
            "heavy_chain",
            "light_chain",
            "heavy_score",
            "light_score",
            "chain_lengths",
            "chain_scores",
        ],
    )
    write_tsv(
        out_dir / "target_chain_role_summary.tsv",
        target_assignment_rows,
        [
            "target_id",
            "manifest_conformers",
            "assignment_ok_conformers",
            "assignment_failed_conformers",
            "heavy_chain_counts",
            "light_chain_counts",
            "stable_heavy_chain",
            "stable_light_chain",
            "vh_seq_len",
            "vl_seq_len",
        ],
    )
    write_tsv(
        out_dir / "antibody_diversity_summary.tsv",
        diversity_rows,
        [
            "target_id",
            "target_index",
            "metric",
            "n_conformers",
            "common_ca_count",
            "medoid_pick_rank",
            "medoid_mean_distance",
            "pairwise_rmsd_mean",
            "pairwise_rmsd_median",
            "pairwise_rmsd_p90",
            "pairwise_rmsd_max",
            "vh_vl_centroid_distance_mean",
            "vh_vl_centroid_distance_std",
            "vh_vl_axis_angle_mean_deg",
            "vh_vl_axis_angle_p90_deg",
            "vh_vl_axis_angle_max_deg",
        ],
    )
    write_tsv(
        out_dir / "antibody_saturation_curves.tsv",
        saturation_rows,
        [
            "target_id",
            "metric",
            "strategy",
            "k",
            "replicate",
            "selected_n",
            "coverage_mean_nearest_rmsd",
            "coverage_median_nearest_rmsd",
            "coverage_p90_nearest_rmsd",
            "coverage_max_nearest_rmsd",
            "selected_indices",
        ],
    )
    write_tsv(
        out_dir / "antibody_selection_strategy_summary.tsv",
        aggregate_rows,
        [
            "metric",
            "strategy",
            "k",
            "n_rows",
            "n_targets",
            "mean_coverage_mean_nearest_rmsd",
            "median_coverage_mean_nearest_rmsd",
            "p90_coverage_mean_nearest_rmsd",
        ],
    )
    write_tsv(
        out_dir / "global_vs_antibody_rmsd.tsv",
        compare_rows,
        [
            "target_id",
            "global_pairwise_rmsd_mean",
            "antibody_pairwise_rmsd_mean",
            "antibody_minus_global",
            "antibody_div_global",
            "global_common_ca_count",
            "antibody_common_ca_count",
        ],
    )

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_root": str(dataset_root),
        "out_dir": str(out_dir),
        "targets_in_manifest": len(target_rows),
        "targets_analyzed": len(target_order),
        "conformers_in_manifest": sum(len(rows) for rows in by_target.values()),
        "assignment_ok_conformers": sum(1 for row in assignment_rows if row["assignment_ok"] == "true"),
        "assignment_failed_conformers": sum(1 for row in assignment_rows if row["assignment_ok"] == "false"),
        "targets_with_antibody_distances": len({row["target_id"] for row in diversity_rows if row["metric"] == "antibody_ca"}),
        "role_score_threshold": ROLE_SCORE_THRESHOLD,
        "k_values": list(K_VALUES),
        "random_seeds": list(RANDOM_SEEDS),
    }
    (out_dir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    make_plots(out_dir, saturation_rows)
    write_report(out_dir, dataset_root, summary, aggregate_rows, compare_rows)


def write_report(
    out_dir: Path,
    dataset_root: Path,
    summary: dict[str, object],
    aggregate_rows: list[dict[str, object]],
    compare_rows: list[dict[str, object]],
) -> None:
    report = out_dir / "aim1_phase1_2_antibody_structural_report.md"
    antibody_rows = [
        row
        for row in aggregate_rows
        if row["metric"] == "antibody_ca" and row["strategy"] in {"evenly_spaced", "random", "greedy_kcenter", "first"}
    ]
    compare_ratios = [float(row["antibody_div_global"]) for row in compare_rows]
    lines = [
        "# Aim 1 Phase 1.2 Antibody Structural Analysis",
        "",
        f"Created: {summary['created_at_utc']}",
        "",
        f"Dataset root: `{dataset_root}`",
        "",
        "## Chain Assignment Summary",
        "",
        f"- Targets analyzed: {summary['targets_analyzed']}",
        f"- Conformers in manifest: {summary['conformers_in_manifest']}",
        f"- Assignment-ok conformers: {summary['assignment_ok_conformers']}",
        f"- Assignment-failed conformers: {summary['assignment_failed_conformers']}",
        f"- Targets with antibody distance matrices: {summary['targets_with_antibody_distances']}",
        "",
    ]
    if compare_ratios:
        lines.extend(
            [
                "## Global vs Antibody-Only",
                "",
                f"- Median antibody/global mean pairwise RMSD ratio: {float(np.median(compare_ratios)):.3f}",
                f"- Mean antibody/global mean pairwise RMSD ratio: {float(np.mean(compare_ratios)):.3f}",
                "",
            ]
        )
    lines.extend(
        [
            "## Antibody-CA Saturation Summary",
            "",
            "Median coverage mean-nearest RMSD across targets. Lower is better.",
            "",
            "| Strategy | K | Median mean-nearest RMSD (A) | Mean mean-nearest RMSD (A) | Targets |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in antibody_rows:
        lines.append(
            f"| {row['strategy']} | {row['k']} | {float(row['median_coverage_mean_nearest_rmsd']):.3f} | "
            f"{float(row['mean_coverage_mean_nearest_rmsd']):.3f} | {row['n_targets']} |"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `chain_role_assignments.tsv`",
            "- `target_chain_role_summary.tsv`",
            "- `antibody_diversity_summary.tsv`",
            "- `antibody_saturation_curves.tsv`",
            "- `antibody_selection_strategy_summary.tsv`",
            "- `global_vs_antibody_rmsd.tsv`",
            "- `figures/`",
            "",
            "Caveat: this remains structural-only and does not include CDR-H3 numbering, embeddings, or strict-300 retrieval.",
            "",
        ]
    )
    report.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
