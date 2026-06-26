#!/usr/bin/env python3
"""CPU-only structural saturation analysis for copied PH/AF3 ensembles.

This script intentionally reads from a Jerry-owned copied dataset root and writes
to an explicit output directory. It does not need GPU access.
"""

from __future__ import annotations

import argparse
import csv
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


@dataclass
class ParsedConformer:
    target_id: str
    pick_rank: int
    copied_relative_path: str
    path: Path
    parse_ok: bool
    error: str
    ca_by_key: dict[tuple[str, int], np.ndarray]
    chain_counts: dict[str, int]
    chain_centroids: dict[str, np.ndarray]


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


def percentile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    arr = sorted(values)
    if len(arr) == 1:
        return arr[0]
    pos = (len(arr) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return arr[lo]
    frac = pos - lo
    return arr[lo] * (1.0 - frac) + arr[hi] * frac


def find_manifest(dataset_root: Path, suffix: str) -> Path:
    matches = sorted((dataset_root / "manifests").glob(f"*{suffix}.tsv"))
    if len(matches) != 1:
        raise SystemExit(f"Expected exactly one *{suffix}.tsv under {dataset_root / 'manifests'}, found {matches}")
    return matches[0]


def parse_atom_site_ca(path: Path) -> tuple[dict[tuple[str, int], np.ndarray], dict[str, int], dict[str, np.ndarray]]:
    """Parse CA coordinates from an AlphaFold/ModelCIF-style mmCIF file."""

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
    required = ["label_atom_id", "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"]
    missing = [name for name in required if name not in col]
    if missing:
        raise ValueError(f"missing atom_site columns: {','.join(missing)}")

    ca_by_key: dict[tuple[str, int], np.ndarray] = {}
    chain_points: dict[str, list[np.ndarray]] = defaultdict(list)
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped or stripped == "#":
            break
        if stripped == "loop_" or stripped.startswith("_"):
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
        chain_points[chain].append(coord)

    if not ca_by_key:
        raise ValueError("no CA coordinates parsed")
    chain_counts = {chain: len(points) for chain, points in chain_points.items()}
    chain_centroids = {chain: np.vstack(points).mean(axis=0) for chain, points in chain_points.items()}
    return ca_by_key, chain_counts, chain_centroids


def parse_conformer(dataset_root: Path, row: dict[str, str]) -> ParsedConformer:
    target_id = row["target_id"]
    pick_rank = int(row.get("pick_rank") or 0)
    copied_relative_path = row["copied_relative_path"]
    path = dataset_root / copied_relative_path
    try:
        ca_by_key, chain_counts, chain_centroids = parse_atom_site_ca(path)
        return ParsedConformer(
            target_id=target_id,
            pick_rank=pick_rank,
            copied_relative_path=copied_relative_path,
            path=path,
            parse_ok=True,
            error="",
            ca_by_key=ca_by_key,
            chain_counts=chain_counts,
            chain_centroids=chain_centroids,
        )
    except Exception as exc:  # noqa: BLE001 - report parse failures rather than aborting the whole dataset.
        return ParsedConformer(
            target_id=target_id,
            pick_rank=pick_rank,
            copied_relative_path=copied_relative_path,
            path=path,
            parse_ok=False,
            error=f"{type(exc).__name__}: {exc}",
            ca_by_key={},
            chain_counts={},
            chain_centroids={},
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


def coverage_row(target_id: str, strategy: str, k: int, replicate: str, selected: list[int], dist: np.ndarray) -> dict[str, object]:
    selected = sorted(set(selected))
    nearest = dist[:, selected].min(axis=1)
    return {
        "target_id": target_id,
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


def top_two_chain_ids(parsed: list[ParsedConformer]) -> list[str]:
    counts: Counter[str] = Counter()
    seen: Counter[str] = Counter()
    for conf in parsed:
        for chain, n_ca in conf.chain_counts.items():
            counts[chain] += n_ca
            seen[chain] += 1
    common = [chain for chain, n_seen in seen.items() if n_seen == len(parsed)]
    return sorted(common, key=lambda c: (-counts[c], c))[:2]


def target_arrays(parsed: list[ParsedConformer], keys: list[tuple[str, int]]) -> list[np.ndarray]:
    return [np.vstack([conf.ca_by_key[key] for key in keys]) for conf in parsed]


def make_plots(out_dir: Path, target_rows: list[dict[str, object]], saturation_rows: list[dict[str, object]]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return

    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    means = [float(row["pairwise_rmsd_mean"]) for row in target_rows if row.get("pairwise_rmsd_mean") not in {"", "nan"}]
    if means:
        plt.figure(figsize=(7, 4))
        plt.hist(means, bins=30, color="#315c72", edgecolor="white")
        plt.xlabel("Mean pairwise CA RMSD per target (A)")
        plt.ylabel("Target count")
        plt.tight_layout()
        plt.savefig(fig_dir / "target_mean_pairwise_rmsd_hist.png", dpi=180)
        plt.close()

    grouped: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in saturation_rows:
        grouped[(str(row["strategy"]), int(row["k"]))].append(float(row["coverage_mean_nearest_rmsd"]))
    if grouped:
        plt.figure(figsize=(7, 4))
        for strategy in sorted({key[0] for key in grouped}):
            xs = []
            ys = []
            for k in K_VALUES:
                vals = grouped.get((strategy, k), [])
                if vals:
                    xs.append(k)
                    ys.append(float(np.median(vals)))
            if xs:
                plt.plot(xs, ys, marker="o", label=strategy)
        plt.xscale("log", base=2)
        plt.xlabel("Selected conformers per target (K)")
        plt.ylabel("Median mean-nearest RMSD to selected set (A)")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(fig_dir / "saturation_median_mean_nearest_rmsd.png", dpi=180)
        plt.close()


def write_report(out_dir: Path, dataset_root: Path, summary: dict[str, object], aggregate_rows: list[dict[str, object]]) -> None:
    report = out_dir / "aim1_phase1_structural_report.md"
    top_lines = [
        "# Aim 1 Phase 1 Structural Analysis",
        "",
        f"Created: {summary['created_at_utc']}",
        "",
        f"Dataset root: `{dataset_root}`",
        "",
        "## Parse Summary",
        "",
        f"- Targets analyzed: {summary['targets_analyzed']}",
        f"- Conformers in manifest: {summary['conformers_in_manifest']}",
        f"- Parse-ok conformers: {summary['parse_ok_conformers']}",
        f"- Parse-failed conformers: {summary['parse_failed_conformers']}",
        f"- Targets with distance matrices: {summary['targets_with_distance_matrices']}",
        "",
        "## Saturation Summary",
        "",
        "Median coverage mean-nearest RMSD across targets. Lower is better.",
        "",
        "| Strategy | K | Median mean-nearest RMSD (A) | Mean mean-nearest RMSD (A) | Targets |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        top_lines.append(
            f"| {row['strategy']} | {row['k']} | {float(row['median_coverage_mean_nearest_rmsd']):.3f} | "
            f"{float(row['mean_coverage_mean_nearest_rmsd']):.3f} | {row['n_targets']} |"
        )
    top_lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `analysis_catalog.tsv`",
            "- `target_parse_summary.tsv`",
            "- `target_diversity_summary.tsv`",
            "- `saturation_curves.tsv`",
            "- `selection_strategy_summary.tsv`",
            "- `figures/`",
            "",
            "Caveat: this is a structural-only readout on copied SAbDab/PDB-like PH/AF3 ensembles, not the strict-300 disease/retrieval benchmark.",
            "",
        ]
    )
    report.write_text("\n".join(top_lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--max-targets", type=int, default=None)
    parser.add_argument("--max-conformers-per-target", type=int, default=None)
    parser.add_argument("--random-seed", type=int, default=0)
    args = parser.parse_args()

    dataset_root = args.dataset_root.resolve()
    out_dir = args.out_dir.resolve()
    if not dataset_root.exists():
        raise SystemExit(f"dataset root not found: {dataset_root}")
    out_dir.mkdir(parents=True, exist_ok=True)

    conformers_path = find_manifest(dataset_root, "conformers")
    targets_path = find_manifest(dataset_root, "targets")
    conformer_rows = read_tsv(conformers_path)
    target_rows = read_tsv(targets_path)
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

    catalog_rows: list[dict[str, object]] = []
    parse_summary_rows: list[dict[str, object]] = []
    diversity_rows: list[dict[str, object]] = []
    saturation_rows: list[dict[str, object]] = []
    parse_failures: list[dict[str, object]] = []
    random.seed(args.random_seed)

    for target_index, target_id in enumerate(target_order, start=1):
        rows = by_target.get(target_id, [])
        parsed = [parse_conformer(dataset_root, row) for row in rows]
        ok = [conf for conf in parsed if conf.parse_ok]
        for conf in parsed:
            chain_counts_text = ";".join(f"{chain}:{n}" for chain, n in sorted(conf.chain_counts.items()))
            catalog_rows.append(
                {
                    "target_id": conf.target_id,
                    "pick_rank": conf.pick_rank,
                    "copied_relative_path": conf.copied_relative_path,
                    "parse_ok": str(conf.parse_ok).lower(),
                    "parse_error": conf.error,
                    "n_ca": len(conf.ca_by_key),
                    "n_chains": len(conf.chain_counts),
                    "chain_ids": ",".join(sorted(conf.chain_counts)),
                    "chain_ca_counts": chain_counts_text,
                }
            )
            if not conf.parse_ok:
                parse_failures.append(
                    {
                        "target_id": conf.target_id,
                        "pick_rank": conf.pick_rank,
                        "copied_relative_path": conf.copied_relative_path,
                        "parse_error": conf.error,
                    }
                )
        parse_summary_rows.append(
            {
                "target_id": target_id,
                "manifest_conformers": len(rows),
                "parse_ok_conformers": len(ok),
                "parse_failed_conformers": len(parsed) - len(ok),
                "chain_sets": "|".join(sorted({",".join(sorted(conf.chain_counts)) for conf in ok})),
            }
        )
        if len(ok) < 2:
            continue

        common_keys = set(ok[0].ca_by_key)
        for conf in ok[1:]:
            common_keys &= set(conf.ca_by_key)
        common_keys_sorted = sorted(common_keys)
        if len(common_keys_sorted) < 3:
            continue
        arrays = target_arrays(ok, common_keys_sorted)
        dist = pairwise_distance_matrix(arrays)
        dist_summary = summarize_distances(dist)

        common_chain_ids = sorted({chain for chain, _seq in common_keys_sorted})
        top2 = top_two_chain_ids(ok)
        top2_distance_mean = float("nan")
        top2_distance_std = float("nan")
        if len(top2) == 2:
            distances = [
                float(np.linalg.norm(conf.chain_centroids[top2[0]] - conf.chain_centroids[top2[1]]))
                for conf in ok
                if top2[0] in conf.chain_centroids and top2[1] in conf.chain_centroids
            ]
            if distances:
                top2_distance_mean = float(np.mean(distances))
                top2_distance_std = float(np.std(distances))

        medoid_idx = int(np.argmin(dist.mean(axis=1)))
        diversity_rows.append(
            {
                "target_id": target_id,
                "target_index": target_index,
                "n_conformers_parse_ok": len(ok),
                "common_ca_count": len(common_keys_sorted),
                "common_chain_ids": ",".join(common_chain_ids),
                "top2_chain_ids": ",".join(top2),
                "top2_centroid_distance_mean": top2_distance_mean,
                "top2_centroid_distance_std": top2_distance_std,
                "medoid_pick_rank": ok[medoid_idx].pick_rank,
                "medoid_mean_distance": float(dist[medoid_idx].mean()),
                **dist_summary,
            }
        )

        n = len(ok)
        effective_k_values = [k for k in K_VALUES if k <= n]
        for k in effective_k_values:
            saturation_rows.append(coverage_row(target_id, "first", k, "0", list(range(k)), dist))
            saturation_rows.append(coverage_row(target_id, "evenly_spaced", k, "0", evenly_spaced_indices(n, k), dist))
            saturation_rows.append(coverage_row(target_id, "greedy_kcenter", k, "0", greedy_kcenter_indices(dist, k), dist))
            for seed in RANDOM_SEEDS:
                rng = random.Random(seed + target_index * 1000 + k)
                selected = sorted(rng.sample(range(n), k))
                saturation_rows.append(coverage_row(target_id, "random", k, str(seed), selected, dist))

        print(f"processed target {target_index}/{len(target_order)} {target_id} parse_ok={len(ok)}", flush=True)

    aggregate_rows = []
    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for row in saturation_rows:
        grouped[(str(row["strategy"]), int(row["k"]))].append(row)
    for (strategy, k), rows in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1])):
        vals = [float(row["coverage_mean_nearest_rmsd"]) for row in rows]
        targets = {str(row["target_id"]) for row in rows}
        aggregate_rows.append(
            {
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
        out_dir / "analysis_catalog.tsv",
        catalog_rows,
        ["target_id", "pick_rank", "copied_relative_path", "parse_ok", "parse_error", "n_ca", "n_chains", "chain_ids", "chain_ca_counts"],
    )
    write_tsv(
        out_dir / "target_parse_summary.tsv",
        parse_summary_rows,
        ["target_id", "manifest_conformers", "parse_ok_conformers", "parse_failed_conformers", "chain_sets"],
    )
    write_tsv(
        out_dir / "parse_failures.tsv",
        parse_failures,
        ["target_id", "pick_rank", "copied_relative_path", "parse_error"],
    )
    write_tsv(
        out_dir / "target_diversity_summary.tsv",
        diversity_rows,
        [
            "target_id",
            "target_index",
            "n_conformers_parse_ok",
            "common_ca_count",
            "common_chain_ids",
            "top2_chain_ids",
            "top2_centroid_distance_mean",
            "top2_centroid_distance_std",
            "medoid_pick_rank",
            "medoid_mean_distance",
            "pairwise_rmsd_mean",
            "pairwise_rmsd_median",
            "pairwise_rmsd_p90",
            "pairwise_rmsd_max",
        ],
    )
    write_tsv(
        out_dir / "saturation_curves.tsv",
        saturation_rows,
        [
            "target_id",
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
        out_dir / "selection_strategy_summary.tsv",
        aggregate_rows,
        [
            "strategy",
            "k",
            "n_rows",
            "n_targets",
            "mean_coverage_mean_nearest_rmsd",
            "median_coverage_mean_nearest_rmsd",
            "p90_coverage_mean_nearest_rmsd",
        ],
    )

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_root": str(dataset_root),
        "out_dir": str(out_dir),
        "targets_in_manifest": len(target_rows),
        "targets_analyzed": len(target_order),
        "conformers_in_manifest": sum(len(rows) for rows in by_target.values()),
        "parse_ok_conformers": sum(1 for row in catalog_rows if row["parse_ok"] == "true"),
        "parse_failed_conformers": sum(1 for row in catalog_rows if row["parse_ok"] == "false"),
        "targets_with_distance_matrices": len(diversity_rows),
        "k_values": list(K_VALUES),
        "random_seeds": list(RANDOM_SEEDS),
    }
    (out_dir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    make_plots(out_dir, diversity_rows, saturation_rows)
    write_report(out_dir, dataset_root, summary, aggregate_rows)


if __name__ == "__main__":
    main()
