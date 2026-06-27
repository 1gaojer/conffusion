#!/usr/bin/env python3
"""CDR-focused structural analysis for copied PH/AF3 antibody ensembles.

Phase 1.3 tests whether conformer diversity is concentrated in antibody CDR
loops after the Phase 1.2 whole-heavy/light analysis showed low global VH/VL
motion.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import json
import random
import shlex
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
from abnumber import Chain


K_VALUES = (1, 2, 4, 8, 16, 32, 64, 128)
RANDOM_SEEDS = (0, 1, 2, 3, 4)
ROLE_SCORE_THRESHOLD = 0.90
CHAIN_MAP_THRESHOLD = 0.90

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
    seq_ids_by_chain: dict[str, list[int]]
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


@dataclass(frozen=True)
class RegionDef:
    target_id: str
    chain_role: str
    chain_type: str
    region: str
    region_type: str
    seq_indices: tuple[int, ...]
    imgt_positions: tuple[str, ...]
    sequence: str


@dataclass
class MappedConformer:
    assigned: AssignedConformer
    role_map_score: dict[str, float]
    coords_by_region: dict[str, dict[tuple[str, int], np.ndarray]]
    sequence_by_region: dict[str, str]
    seq_ids_by_region: dict[str, list[int]]
    mapping_errors: list[str]


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

    seq_by_chain: dict[str, str] = {}
    seq_ids_by_chain: dict[str, list[int]] = {}
    by_chain: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for chain, seq_id in sorted(residues):
        by_chain[chain].append((seq_id, residues[(chain, seq_id)]))
    for chain, pairs in by_chain.items():
        seq_ids_by_chain[chain] = [seq_id for seq_id, _aa in pairs]
        seq_by_chain[chain] = "".join(aa for _seq_id, aa in pairs)

    chain_counts = {chain: len(points) for chain, points in chain_points.items()}
    chain_centroids = {chain: np.vstack(points).mean(axis=0) for chain, points in chain_points.items()}
    return ChainParsed(
        ca_by_key=ca_by_key,
        seq_by_chain=seq_by_chain,
        seq_ids_by_chain=seq_ids_by_chain,
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
    return difflib.SequenceMatcher(None, observed, expected, autojunk=False).ratio()


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
    except Exception as exc:  # noqa: BLE001 - preserve audit rows.
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


def chain_region_defs(target_id: str, chain_role: str, sequence: str, scheme: str) -> list[RegionDef]:
    chain = Chain(sequence, scheme=scheme, cdr_definition=scheme, use_anarcii=True)
    pos_to_seq_index = {pos: i for i, (pos, _aa) in enumerate(chain)}
    prefix = "H" if chain_role == "heavy" else "L"
    expected_type = "H" if chain_role == "heavy" else None
    if expected_type is not None and chain.chain_type != expected_type:
        raise ValueError(f"expected heavy chain but AbNumber assigned {chain.chain_type}")
    if chain_role == "light" and chain.chain_type not in {"K", "L"}:
        raise ValueError(f"expected light chain but AbNumber assigned {chain.chain_type}")

    region_attrs = [
        ("FR1", "framework", chain.fr1_dict),
        ("CDR1", "cdr", chain.cdr1_dict),
        ("FR2", "framework", chain.fr2_dict),
        ("CDR2", "cdr", chain.cdr2_dict),
        ("FR3", "framework", chain.fr3_dict),
        ("CDR3", "cdr", chain.cdr3_dict),
        ("FR4", "framework", chain.fr4_dict),
    ]
    regions = []
    for short_name, region_type, region_dict in region_attrs:
        seq_indices = []
        imgt_positions = []
        letters = []
        for pos, aa in region_dict.items():
            if pos not in pos_to_seq_index:
                continue
            seq_indices.append(pos_to_seq_index[pos])
            imgt_positions.append(str(pos))
            letters.append(aa)
        regions.append(
            RegionDef(
                target_id=target_id,
                chain_role=chain_role,
                chain_type=chain.chain_type,
                region=f"{prefix}-{short_name}",
                region_type=region_type,
                seq_indices=tuple(seq_indices),
                imgt_positions=tuple(imgt_positions),
                sequence="".join(letters),
            )
        )
    return regions


def target_region_defs(target_id: str, target_meta: dict[str, str], scheme: str) -> dict[str, RegionDef]:
    regions = {}
    for region in chain_region_defs(target_id, "heavy", target_meta.get("vh_seq", ""), scheme):
        regions[region.region] = region
    for region in chain_region_defs(target_id, "light", target_meta.get("vl_seq", ""), scheme):
        regions[region.region] = region
    return regions


def expected_to_observed_map(expected: str, observed: str) -> dict[int, int]:
    if expected == observed:
        return {i: i for i in range(len(expected))}
    matcher = difflib.SequenceMatcher(None, expected, observed, autojunk=False)
    mapping: dict[int, int] = {}
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag != "equal":
            continue
        for offset in range(min(i2 - i1, j2 - j1)):
            mapping[i1 + offset] = j1 + offset
    return mapping


def merge_region_defs(target_id: str, region: str, defs: list[RegionDef]) -> RegionDef:
    seq_indices = []
    imgt_positions = []
    sequence = []
    chain_roles = []
    chain_types = []
    for region_def in defs:
        seq_indices.extend(region_def.seq_indices)
        imgt_positions.extend(region_def.imgt_positions)
        sequence.append(region_def.sequence)
        chain_roles.append(region_def.chain_role)
        chain_types.append(region_def.chain_type)
    return RegionDef(
        target_id=target_id,
        chain_role="+".join(chain_roles),
        chain_type="+".join(chain_types),
        region=region,
        region_type="combined",
        seq_indices=tuple(seq_indices),
        imgt_positions=tuple(imgt_positions),
        sequence="|".join(sequence),
    )


def make_region_groups(target_id: str, regions: dict[str, RegionDef]) -> dict[str, list[str]]:
    groups = {
        "H-CDR1": ["H-CDR1"],
        "H-CDR2": ["H-CDR2"],
        "H-CDR3": ["H-CDR3"],
        "L-CDR1": ["L-CDR1"],
        "L-CDR2": ["L-CDR2"],
        "L-CDR3": ["L-CDR3"],
        "H-CDRs": ["H-CDR1", "H-CDR2", "H-CDR3"],
        "L-CDRs": ["L-CDR1", "L-CDR2", "L-CDR3"],
        "all-CDRs": ["H-CDR1", "H-CDR2", "H-CDR3", "L-CDR1", "L-CDR2", "L-CDR3"],
        "H-FR": ["H-FR1", "H-FR2", "H-FR3", "H-FR4"],
        "L-FR": ["L-FR1", "L-FR2", "L-FR3", "L-FR4"],
        "all-FR": ["H-FR1", "H-FR2", "H-FR3", "H-FR4", "L-FR1", "L-FR2", "L-FR3", "L-FR4"],
    }
    for name, parts in list(groups.items()):
        missing = [part for part in parts if part not in regions]
        if missing:
            raise ValueError(f"{target_id} missing region definitions for {name}: {','.join(missing)}")
    return groups


def map_conformer_regions(
    conf: AssignedConformer,
    target_meta: dict[str, str],
    regions: dict[str, RegionDef],
    region_groups: dict[str, list[str]],
) -> MappedConformer:
    if conf.parsed is None:
        return MappedConformer(conf, {}, {}, {}, {}, [conf.error or "parse/assignment failed"])

    role_to_chain = {"heavy": conf.heavy_chain, "light": conf.light_chain}
    role_to_expected = {"heavy": target_meta.get("vh_seq", ""), "light": target_meta.get("vl_seq", "")}
    expected_to_cif: dict[str, dict[int, tuple[int, np.ndarray, str]]] = {}
    role_map_score = {}
    errors = []
    for role, chain_id in role_to_chain.items():
        expected = role_to_expected[role]
        observed = conf.parsed.seq_by_chain.get(chain_id, "")
        seq_ids = conf.parsed.seq_ids_by_chain.get(chain_id, [])
        seq_map = expected_to_observed_map(expected, observed)
        role_map_score[role] = len(seq_map) / len(expected) if expected else 0.0
        if role_map_score[role] < CHAIN_MAP_THRESHOLD:
            errors.append(f"{role} sequence map score {role_map_score[role]:.3f} below {CHAIN_MAP_THRESHOLD:.3f}")
        mapped = {}
        for expected_index, observed_index in seq_map.items():
            if observed_index >= len(seq_ids):
                continue
            seq_id = seq_ids[observed_index]
            coord = conf.parsed.ca_by_key.get((chain_id, seq_id))
            if coord is None:
                continue
            mapped[expected_index] = (seq_id, coord, observed[observed_index])
        expected_to_cif[role] = mapped

    coords_by_atomic_region: dict[str, dict[tuple[str, int], np.ndarray]] = {}
    sequence_by_atomic_region: dict[str, str] = {}
    seq_ids_by_atomic_region: dict[str, list[int]] = {}
    for name, region_def in regions.items():
        coords = {}
        letters = []
        seq_ids = []
        role_map = expected_to_cif.get(region_def.chain_role, {})
        for expected_index in region_def.seq_indices:
            mapped = role_map.get(expected_index)
            if mapped is None:
                continue
            seq_id, coord, aa = mapped
            coords[(region_def.chain_role, expected_index)] = coord
            seq_ids.append(seq_id)
            letters.append(aa)
        coords_by_atomic_region[name] = coords
        sequence_by_atomic_region[name] = "".join(letters)
        seq_ids_by_atomic_region[name] = seq_ids
        if len(coords) != len(region_def.seq_indices):
            errors.append(
                f"{name} mapped {len(coords)}/{len(region_def.seq_indices)} residues"
            )

    coords_by_region = dict(coords_by_atomic_region)
    sequence_by_region = dict(sequence_by_atomic_region)
    seq_ids_by_region = dict(seq_ids_by_atomic_region)
    for group_name, parts in region_groups.items():
        coords = {}
        seq = []
        seq_ids = []
        for part in parts:
            coords.update(coords_by_atomic_region[part])
            seq.append(sequence_by_atomic_region[part])
            seq_ids.extend(seq_ids_by_atomic_region[part])
        coords_by_region[group_name] = coords
        sequence_by_region[group_name] = "|".join(seq)
        seq_ids_by_region[group_name] = seq_ids

    return MappedConformer(conf, role_map_score, coords_by_region, sequence_by_region, seq_ids_by_region, errors)


def kabsch_fit_transform(mobile_frame: np.ndarray, reference_frame: np.ndarray, mobile_points: np.ndarray) -> np.ndarray:
    mobile_center = mobile_frame.mean(axis=0)
    reference_center = reference_frame.mean(axis=0)
    mobile_zero = mobile_frame - mobile_center
    reference_zero = reference_frame - reference_center
    cov = mobile_zero.T @ reference_zero
    u, _s, vt = np.linalg.svd(cov)
    correction = np.eye(3)
    correction[2, 2] = np.sign(np.linalg.det(u @ vt))
    rot = u @ correction @ vt
    return (mobile_points - mobile_center) @ rot + reference_center


def ca_rmsd(points_a: np.ndarray, points_b: np.ndarray) -> float:
    diff = points_a - points_b
    return float(np.sqrt((diff * diff).sum(axis=1).mean()))


def kabsch_rmsd(points_a: np.ndarray, points_b: np.ndarray) -> float:
    aligned = kabsch_fit_transform(points_a, points_b, points_a)
    return ca_rmsd(aligned, points_b)


def common_region_arrays(
    mapped: list[MappedConformer],
    region: str,
) -> tuple[list[np.ndarray], list[tuple[str, int]]]:
    key_sets = [set(conf.coords_by_region.get(region, {})) for conf in mapped]
    common_keys = sorted(set.intersection(*key_sets)) if key_sets else []
    if not common_keys:
        return [], []
    arrays = []
    for conf in mapped:
        coords = conf.coords_by_region.get(region, {})
        arrays.append(np.vstack([coords[key] for key in common_keys]))
    return arrays, common_keys


def frame_aligned_distance_matrix(
    region_arrays: list[np.ndarray],
    frame_arrays: list[np.ndarray],
    reducer: str,
) -> np.ndarray:
    n = len(region_arrays)
    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            moved = kabsch_fit_transform(frame_arrays[i], frame_arrays[j], region_arrays[i])
            if reducer == "centroid":
                d = float(np.linalg.norm(moved.mean(axis=0) - region_arrays[j].mean(axis=0)))
            else:
                d = ca_rmsd(moved, region_arrays[j])
            dist[i, j] = d
            dist[j, i] = d
    return dist


def self_aligned_distance_matrix(region_arrays: list[np.ndarray]) -> np.ndarray:
    n = len(region_arrays)
    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            d = kabsch_rmsd(region_arrays[i], region_arrays[j])
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
    region: str,
    align_region: str,
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
        "region": region,
        "align_region": align_region,
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


def add_saturation_rows(
    target_index: int,
    target_id: str,
    metric: str,
    region: str,
    align_region: str,
    dist: np.ndarray,
    saturation_rows: list[dict[str, object]],
) -> None:
    for k in [k for k in K_VALUES if k <= dist.shape[0]]:
        saturation_rows.append(coverage_row(target_id, metric, region, align_region, "first", k, "0", list(range(k)), dist))
        saturation_rows.append(
            coverage_row(target_id, metric, region, align_region, "evenly_spaced", k, "0", evenly_spaced_indices(dist.shape[0], k), dist)
        )
        saturation_rows.append(
            coverage_row(target_id, metric, region, align_region, "greedy_kcenter", k, "0", greedy_kcenter_indices(dist, k), dist)
        )
        for seed in RANDOM_SEEDS:
            rng = random.Random(seed + target_index * 1000 + k)
            selected = sorted(rng.sample(range(dist.shape[0]), k))
            saturation_rows.append(coverage_row(target_id, metric, region, align_region, "random", k, str(seed), selected, dist))


def region_alignments() -> list[tuple[str, str]]:
    return [
        ("H-CDR1", "H-FR"),
        ("H-CDR2", "H-FR"),
        ("H-CDR3", "H-FR"),
        ("L-CDR1", "L-FR"),
        ("L-CDR2", "L-FR"),
        ("L-CDR3", "L-FR"),
        ("H-CDRs", "H-FR"),
        ("L-CDRs", "L-FR"),
        ("all-CDRs", "all-FR"),
    ]


def numbering_rows_for_target(target_id: str, target_meta: dict[str, str], scheme: str) -> list[dict[str, object]]:
    regions = target_region_defs(target_id, target_meta, scheme)
    rows = []
    for region in regions.values():
        rows.append(
            {
                "target_id": target_id,
                "chain_role": region.chain_role,
                "chain_type": region.chain_type,
                "region": region.region,
                "region_type": region.region_type,
                "numbering_scheme": scheme,
                "residue_count": len(region.seq_indices),
                "seq_index_start_1based": min(region.seq_indices) + 1 if region.seq_indices else "",
                "seq_index_end_1based": max(region.seq_indices) + 1 if region.seq_indices else "",
                "imgt_start": region.imgt_positions[0] if region.imgt_positions else "",
                "imgt_end": region.imgt_positions[-1] if region.imgt_positions else "",
                "imgt_positions": ",".join(region.imgt_positions),
                "sequence": region.sequence,
            }
        )
    return rows


def aggregate_saturation_rows(saturation_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    aggregate_rows = []
    grouped: dict[tuple[str, str, str, str, int], list[dict[str, object]]] = defaultdict(list)
    for row in saturation_rows:
        grouped[
            (
                str(row["metric"]),
                str(row["region"]),
                str(row["align_region"]),
                str(row["strategy"]),
                int(row["k"]),
            )
        ].append(row)
    for (metric, region, align_region, strategy, k), rows in sorted(grouped.items()):
        vals = [float(row["coverage_mean_nearest_rmsd"]) for row in rows]
        targets = {str(row["target_id"]) for row in rows}
        aggregate_rows.append(
            {
                "metric": metric,
                "region": region,
                "align_region": align_region,
                "strategy": strategy,
                "k": k,
                "n_rows": len(rows),
                "n_targets": len(targets),
                "mean_coverage_mean_nearest_rmsd": float(np.mean(vals)),
                "median_coverage_mean_nearest_rmsd": float(np.median(vals)),
                "p90_coverage_mean_nearest_rmsd": float(np.quantile(vals, 0.90)),
            }
        )
    return aggregate_rows


def per_residue_variability(
    target_id: str,
    target_index: int,
    mapped: list[MappedConformer],
    region: str,
    align_region: str,
) -> list[dict[str, object]]:
    region_arrays, region_keys = common_region_arrays(mapped, region)
    frame_arrays, frame_keys = common_region_arrays(mapped, align_region)
    if len(region_arrays) < 2 or len(region_keys) == 0 or len(frame_keys) < 3:
        return []
    frame_dist = self_aligned_distance_matrix(frame_arrays)
    medoid_idx = int(np.argmin(frame_dist.mean(axis=1)))
    reference_frame = frame_arrays[medoid_idx]
    aligned_region_arrays = []
    for frame, arr in zip(frame_arrays, region_arrays, strict=True):
        aligned_region_arrays.append(kabsch_fit_transform(frame, reference_frame, arr))
    stack = np.stack(aligned_region_arrays, axis=0)
    mean_coords = stack.mean(axis=0)
    deviations = np.sqrt(((stack - mean_coords) ** 2).sum(axis=2))
    rows = []
    for idx, key in enumerate(region_keys):
        role, expected_index = key
        rows.append(
            {
                "target_id": target_id,
                "target_index": target_index,
                "region": region,
                "align_region": align_region,
                "chain_role": role,
                "seq_index_1based": expected_index + 1,
                "n_conformers": len(mapped),
                "rmsf_like_ca": float(np.sqrt(np.mean(deviations[:, idx] ** 2))),
                "mean_distance_to_mean_ca": float(np.mean(deviations[:, idx])),
                "p90_distance_to_mean_ca": float(np.quantile(deviations[:, idx], 0.90)),
                "max_distance_to_mean_ca": float(np.max(deviations[:, idx])),
            }
        )
    return rows


def make_plots(out_dir: Path, aggregate_rows: list[dict[str, object]]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return

    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    key_regions = ["H-CDR3", "all-CDRs", "H-CDRs", "L-CDRs"]
    for metric in ["frame_aligned_ca", "self_aligned_ca", "frame_aligned_centroid"]:
        rows = [row for row in aggregate_rows if row["metric"] == metric and row["region"] in key_regions]
        if not rows:
            continue
        for region in key_regions:
            region_rows = [row for row in rows if row["region"] == region]
            if not region_rows:
                continue
            plt.figure(figsize=(7, 4))
            for strategy in sorted({str(row["strategy"]) for row in region_rows}):
                xs = []
                ys = []
                for k in K_VALUES:
                    vals = [
                        float(row["median_coverage_mean_nearest_rmsd"])
                        for row in region_rows
                        if row["strategy"] == strategy and int(row["k"]) == k
                    ]
                    if vals:
                        xs.append(k)
                        ys.append(float(np.median(vals)))
                if xs:
                    plt.plot(xs, ys, marker="o", label=strategy)
            plt.xscale("log", base=2)
            plt.xlabel("Selected conformers per target (K)")
            plt.ylabel("Median mean-nearest distance (A)")
            plt.title(f"{region} {metric}")
            plt.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig(fig_dir / f"saturation_{region}_{metric}.png", dpi=180)
            plt.close()


def write_report(
    out_dir: Path,
    dataset_root: Path,
    summary: dict[str, object],
    aggregate_rows: list[dict[str, object]],
    diversity_rows: list[dict[str, object]],
) -> None:
    report = out_dir / "aim1_phase1_3_cdr_structural_report.md"
    headline = [
        row
        for row in diversity_rows
        if row["metric"] == "frame_aligned_ca" and row["region"] in {"H-CDR3", "all-CDRs", "H-CDRs", "L-CDRs"}
    ]
    headline_grouped: dict[str, list[float]] = defaultdict(list)
    for row in headline:
        headline_grouped[str(row["region"])].append(float(row["pairwise_rmsd_mean"]))

    sat_rows = [
        row
        for row in aggregate_rows
        if row["metric"] == "frame_aligned_ca"
        and row["region"] in {"H-CDR3", "all-CDRs"}
        and row["strategy"] in {"first", "evenly_spaced", "greedy_kcenter", "random"}
        and int(row["k"]) in {8, 16, 32, 64, 128}
    ]

    lines = [
        "# Aim 1 Phase 1.3 CDR Structural Analysis",
        "",
        f"Created: {summary['created_at_utc']}",
        "",
        f"Dataset root: `{dataset_root}`",
        f"Numbering scheme: `{summary['numbering_scheme']}`",
        "",
        "## Summary",
        "",
        f"- Targets analyzed: {summary['targets_analyzed']}",
        f"- Conformers in manifest subset: {summary['conformers_in_manifest']}",
        f"- Assignment-ok conformers: {summary['assignment_ok_conformers']}",
        f"- Mapped conformers used: {summary['mapped_ok_conformers']}",
        f"- Targets with CDR distance matrices: {summary['targets_with_cdr_distances']}",
        "",
        "## Headline Pairwise Diversity",
        "",
        "Target-level mean pairwise frame-aligned CA RMSD, summarized across targets.",
        "",
        "| Region | Median target mean RMSD (A) | Mean target mean RMSD (A) | Targets |",
        "|---|---:|---:|---:|",
    ]
    for region in ["H-CDR3", "H-CDRs", "L-CDRs", "all-CDRs"]:
        vals = headline_grouped.get(region, [])
        if vals:
            lines.append(f"| {region} | {float(np.median(vals)):.3f} | {float(np.mean(vals)):.3f} | {len(vals)} |")

    lines.extend(
        [
            "",
            "## Saturation Snapshot",
            "",
            "Median across targets of mean-nearest distance to the selected subset. Lower is better.",
            "",
            "| Metric | Region | Strategy | K | Median mean-nearest distance (A) | Targets |",
            "|---|---|---|---:|---:|---:|",
        ]
    )
    for row in sat_rows:
        lines.append(
            f"| {row['metric']} | {row['region']} | {row['strategy']} | {row['k']} | "
            f"{float(row['median_coverage_mean_nearest_rmsd']):.3f} | {row['n_targets']} |"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `cdr_numbering_assignments.tsv`",
            "- `cdr_mapping_summary.tsv`",
            "- `cdr_diversity_summary.tsv`",
            "- `cdr_saturation_curves.tsv`",
            "- `cdr_selection_strategy_summary.tsv`",
            "- `cdr_per_residue_variability.tsv`",
            "- `whole_vs_cdr_summary.tsv`",
            "- `figures/`",
            "",
            "Caveat: this is a structural CDR readout. It does not yet test MCA embedding or antigen-retrieval sensitivity.",
            "",
        ]
    )
    report.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--antibody-results-dir", type=Path)
    parser.add_argument("--scheme", default="imgt")
    parser.add_argument("--max-targets", type=int)
    parser.add_argument("--max-conformers-per-target", type=int)
    parser.add_argument("--validate-only", action="store_true")
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

    numbering_rows = []
    mapping_rows = []
    diversity_rows = []
    saturation_rows = []
    per_residue_rows = []
    whole_vs_cdr_rows = []
    assignment_ok_count = 0
    mapped_ok_count = 0

    for target_index, target_id in enumerate(target_order, start=1):
        target_regions = target_region_defs(target_id, target_meta[target_id], args.scheme)
        region_groups = make_region_groups(target_id, target_regions)
        numbering_rows.extend(numbering_rows_for_target(target_id, target_meta[target_id], args.scheme))

        assigned = [parse_and_assign(dataset_root, row, target_meta[target_id]) for row in by_target[target_id]]
        ok_assigned = [conf for conf in assigned if conf.assignment_ok]
        assignment_ok_count += len(ok_assigned)
        mapped = [map_conformer_regions(conf, target_meta[target_id], target_regions, region_groups) for conf in ok_assigned]
        mapped_ok = [
            conf
            for conf in mapped
            if conf.role_map_score.get("heavy", 0.0) >= CHAIN_MAP_THRESHOLD
            and conf.role_map_score.get("light", 0.0) >= CHAIN_MAP_THRESHOLD
        ]
        mapped_ok_count += len(mapped_ok)

        for conf in mapped:
            for region in ["H-CDR1", "H-CDR2", "H-CDR3", "L-CDR1", "L-CDR2", "L-CDR3", "H-CDRs", "L-CDRs", "all-CDRs"]:
                target_def = None
                if region in target_regions:
                    target_def = target_regions[region]
                elif region in region_groups:
                    target_def = merge_region_defs(target_id, region, [target_regions[name] for name in region_groups[region]])
                expected_count = len(target_def.seq_indices) if target_def is not None else ""
                expected_int = int(expected_count) if expected_count != "" else 0
                mapped_count = len(conf.coords_by_region.get(region, {}))
                chain_mapping_ok = (
                    conf.role_map_score.get("heavy", 0.0) >= CHAIN_MAP_THRESHOLD
                    and conf.role_map_score.get("light", 0.0) >= CHAIN_MAP_THRESHOLD
                )
                mapping_rows.append(
                    {
                        "target_id": target_id,
                        "pick_rank": conf.assigned.pick_rank,
                        "copied_relative_path": conf.assigned.copied_relative_path,
                        "region": region,
                        "expected_residues": expected_count,
                        "mapped_residues": mapped_count,
                        "mapped_sequence": conf.sequence_by_region.get(region, ""),
                        "cif_seq_ids": ",".join(str(x) for x in conf.seq_ids_by_region.get(region, [])),
                        "heavy_map_score": conf.role_map_score.get("heavy", 0.0),
                        "light_map_score": conf.role_map_score.get("light", 0.0),
                        "mapping_ok": str(chain_mapping_ok and mapped_count == expected_int).lower(),
                        "mapping_errors": "; ".join(conf.mapping_errors),
                    }
                )

        if args.validate_only:
            print(
                f"validated target {target_index}/{len(target_order)} {target_id} "
                f"assignment_ok={len(ok_assigned)} mapped_ok={len(mapped_ok)}",
                flush=True,
            )
            continue

        if len(mapped_ok) < 2:
            print(
                f"processed target {target_index}/{len(target_order)} {target_id} "
                f"assignment_ok={len(ok_assigned)} mapped_ok={len(mapped_ok)} skipped_distances",
                flush=True,
            )
            continue

        target_has_cdr_distance = False
        for region, align_region in region_alignments():
            region_arrays, region_keys = common_region_arrays(mapped_ok, region)
            frame_arrays, frame_keys = common_region_arrays(mapped_ok, align_region)
            if len(region_arrays) < 2 or len(region_keys) == 0 or len(frame_keys) < 3:
                continue
            target_has_cdr_distance = True

            frame_dist = frame_aligned_distance_matrix(region_arrays, frame_arrays, "ca")
            medoid_idx = int(np.argmin(frame_dist.mean(axis=1)))
            diversity_rows.append(
                {
                    "target_id": target_id,
                    "target_index": target_index,
                    "metric": "frame_aligned_ca",
                    "region": region,
                    "align_region": align_region,
                    "n_conformers": len(mapped_ok),
                    "region_ca_count": len(region_keys),
                    "align_ca_count": len(frame_keys),
                    "medoid_pick_rank": mapped_ok[medoid_idx].assigned.pick_rank,
                    "medoid_mean_distance": float(frame_dist[medoid_idx].mean()),
                    **summarize_distances(frame_dist),
                }
            )
            add_saturation_rows(target_index, target_id, "frame_aligned_ca", region, align_region, frame_dist, saturation_rows)

            centroid_dist = frame_aligned_distance_matrix(region_arrays, frame_arrays, "centroid")
            diversity_rows.append(
                {
                    "target_id": target_id,
                    "target_index": target_index,
                    "metric": "frame_aligned_centroid",
                    "region": region,
                    "align_region": align_region,
                    "n_conformers": len(mapped_ok),
                    "region_ca_count": len(region_keys),
                    "align_ca_count": len(frame_keys),
                    "medoid_pick_rank": mapped_ok[int(np.argmin(centroid_dist.mean(axis=1)))].assigned.pick_rank,
                    "medoid_mean_distance": float(centroid_dist.mean(axis=1).min()),
                    **summarize_distances(centroid_dist),
                }
            )
            add_saturation_rows(
                target_index,
                target_id,
                "frame_aligned_centroid",
                region,
                align_region,
                centroid_dist,
                saturation_rows,
            )

            if len(region_keys) >= 3:
                self_dist = self_aligned_distance_matrix(region_arrays)
                diversity_rows.append(
                    {
                        "target_id": target_id,
                        "target_index": target_index,
                        "metric": "self_aligned_ca",
                        "region": region,
                        "align_region": region,
                        "n_conformers": len(mapped_ok),
                        "region_ca_count": len(region_keys),
                        "align_ca_count": len(region_keys),
                        "medoid_pick_rank": mapped_ok[int(np.argmin(self_dist.mean(axis=1)))].assigned.pick_rank,
                        "medoid_mean_distance": float(self_dist.mean(axis=1).min()),
                        **summarize_distances(self_dist),
                    }
                )
                add_saturation_rows(target_index, target_id, "self_aligned_ca", region, region, self_dist, saturation_rows)

            if region in {"H-CDR3", "all-CDRs"}:
                per_residue_rows.extend(per_residue_variability(target_id, target_index, mapped_ok, region, align_region))

        whole_rows = [row for row in diversity_rows if row["target_id"] == target_id and row["metric"] == "frame_aligned_ca"]
        whole_lookup = {row["region"]: row for row in whole_rows}
        if "H-CDR3" in whole_lookup:
            h3 = whole_lookup["H-CDR3"]
            all_cdr = whole_lookup.get("all-CDRs")
            whole_vs_cdr_rows.append(
                {
                    "target_id": target_id,
                    "target_index": target_index,
                    "n_conformers": h3["n_conformers"],
                    "h_cdr3_region_ca_count": h3["region_ca_count"],
                    "h_cdr3_pairwise_mean": h3["pairwise_rmsd_mean"],
                    "h_cdr3_pairwise_median": h3["pairwise_rmsd_median"],
                    "h_cdr3_pairwise_p90": h3["pairwise_rmsd_p90"],
                    "h_cdr3_pairwise_max": h3["pairwise_rmsd_max"],
                    "all_cdr_pairwise_mean": all_cdr["pairwise_rmsd_mean"] if all_cdr else "",
                    "all_cdr_pairwise_median": all_cdr["pairwise_rmsd_median"] if all_cdr else "",
                    "all_cdr_pairwise_p90": all_cdr["pairwise_rmsd_p90"] if all_cdr else "",
                    "all_cdr_pairwise_max": all_cdr["pairwise_rmsd_max"] if all_cdr else "",
                }
            )

        print(
            f"processed target {target_index}/{len(target_order)} {target_id} "
            f"assignment_ok={len(ok_assigned)} mapped_ok={len(mapped_ok)} cdr_distances={target_has_cdr_distance}",
            flush=True,
        )

    aggregate_rows = aggregate_saturation_rows(saturation_rows)

    write_tsv(
        out_dir / "cdr_numbering_assignments.tsv",
        numbering_rows,
        [
            "target_id",
            "chain_role",
            "chain_type",
            "region",
            "region_type",
            "numbering_scheme",
            "residue_count",
            "seq_index_start_1based",
            "seq_index_end_1based",
            "imgt_start",
            "imgt_end",
            "imgt_positions",
            "sequence",
        ],
    )
    write_tsv(
        out_dir / "cdr_mapping_summary.tsv",
        mapping_rows,
        [
            "target_id",
            "pick_rank",
            "copied_relative_path",
            "region",
            "expected_residues",
            "mapped_residues",
            "mapped_sequence",
            "cif_seq_ids",
            "heavy_map_score",
            "light_map_score",
            "mapping_ok",
            "mapping_errors",
        ],
    )
    write_tsv(
        out_dir / "cdr_diversity_summary.tsv",
        diversity_rows,
        [
            "target_id",
            "target_index",
            "metric",
            "region",
            "align_region",
            "n_conformers",
            "region_ca_count",
            "align_ca_count",
            "medoid_pick_rank",
            "medoid_mean_distance",
            "pairwise_rmsd_mean",
            "pairwise_rmsd_median",
            "pairwise_rmsd_p90",
            "pairwise_rmsd_max",
        ],
    )
    write_tsv(
        out_dir / "cdr_saturation_curves.tsv",
        saturation_rows,
        [
            "target_id",
            "metric",
            "region",
            "align_region",
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
        out_dir / "cdr_selection_strategy_summary.tsv",
        aggregate_rows,
        [
            "metric",
            "region",
            "align_region",
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
        out_dir / "cdr_per_residue_variability.tsv",
        per_residue_rows,
        [
            "target_id",
            "target_index",
            "region",
            "align_region",
            "chain_role",
            "seq_index_1based",
            "n_conformers",
            "rmsf_like_ca",
            "mean_distance_to_mean_ca",
            "p90_distance_to_mean_ca",
            "max_distance_to_mean_ca",
        ],
    )
    write_tsv(
        out_dir / "whole_vs_cdr_summary.tsv",
        whole_vs_cdr_rows,
        [
            "target_id",
            "target_index",
            "n_conformers",
            "h_cdr3_region_ca_count",
            "h_cdr3_pairwise_mean",
            "h_cdr3_pairwise_median",
            "h_cdr3_pairwise_p90",
            "h_cdr3_pairwise_max",
            "all_cdr_pairwise_mean",
            "all_cdr_pairwise_median",
            "all_cdr_pairwise_p90",
            "all_cdr_pairwise_max",
        ],
    )

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_root": str(dataset_root),
        "out_dir": str(out_dir),
        "numbering_scheme": args.scheme,
        "numbering_backend": "abnumber_use_anarcii",
        "targets_in_manifest": len(target_rows),
        "targets_analyzed": len(target_order),
        "conformers_in_manifest": sum(len(rows) for rows in by_target.values()),
        "assignment_ok_conformers": assignment_ok_count,
        "mapped_ok_conformers": mapped_ok_count,
        "targets_with_cdr_distances": len({row["target_id"] for row in diversity_rows if row["metric"] == "frame_aligned_ca"}),
        "role_score_threshold": ROLE_SCORE_THRESHOLD,
        "chain_map_threshold": CHAIN_MAP_THRESHOLD,
        "k_values": list(K_VALUES),
        "random_seeds": list(RANDOM_SEEDS),
        "validate_only": args.validate_only,
    }
    (out_dir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if not args.validate_only:
        make_plots(out_dir, aggregate_rows)
        write_report(out_dir, dataset_root, summary, aggregate_rows, diversity_rows)


if __name__ == "__main__":
    main()
