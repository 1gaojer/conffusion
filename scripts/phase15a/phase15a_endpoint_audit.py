#!/usr/bin/env python3
"""Phase 1.5a endpoint audit and retrieval smoke for Conffusion.

This is a Jerry-owned wrapper. It reads Gaeun/shared inputs as immutable
artifacts and writes all outputs into a caller-provided Jerry-owned directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch


DEFAULT_CONFUSION_ROOT = Path("/project/liulab/jg1920/conffusion")
DEFAULT_PHASE14_RUN = DEFAULT_CONFUSION_ROOT / "phase14_20260626_2301" / "runs" / "full"
DEFAULT_PHASE14_EMB = DEFAULT_PHASE14_RUN / "embeddings" / "mca1000_selected"
DEFAULT_PHASE14_CKPT = (
    DEFAULT_CONFUSION_ROOT
    / "phase14_20260626_2301"
    / "shared"
    / "checkpoints"
    / "mca1000_checkpoint.pt"
)
DEFAULT_CONFIRMED_CKPT = Path(
    "/external/liulab/gkim/antigen_prediction/"
    "2026.06.19_retrain_mca_1000_confs/train_v1/checkpoint.pt"
)
DEFAULT_REF_NPZ = Path(
    "/external/liulab/gkim/antigen_prediction/"
    "2026.06.24_hiv_flu_covid_nn_retrieval/outputs_1000/reference_embeddings.npz"
)
DEFAULT_JUNE24_QUERY_INDEX = Path(
    "/external/liulab/gkim/antigen_prediction/"
    "2026.06.24_hiv_flu_covid_nn_retrieval/query_mca_embeddings/index.json"
)
DEFAULT_CDR_CURVES = Path(
    "/external/liulab/jg1920/conffusion/"
    "aim1_phase1_4_cdr_coresets_20260627/cdr_saturation_curves.tsv"
)
DEFAULT_MEDIUM_TARGETS = (
    DEFAULT_CONFUSION_ROOT
    / "phase14_20260626_2301"
    / "shared"
    / "gaeun_ph_af3_medium_20260626"
    / "manifests"
    / "medium_targets.tsv"
)
DEFAULT_STRICT300 = Path(
    "/project/liulab/jg1920/bcr-conformer-runs/"
    "model3_iedb_latest_pipeline_inputs_20260624/"
    "recommended_strict_300/input/targets.tsv"
)
DEFAULT_OUT = DEFAULT_CONFUSION_ROOT / "phase15a_20260627_endpoint_audit"

CHAINS = ("H", "L")
MCA_DIM = 256


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def file_info(path: Path, with_hash: bool = True) -> dict[str, object]:
    info: dict[str, object] = {
        "path": str(path),
        "exists": path.exists(),
    }
    if not path.exists():
        return info
    stat = path.stat()
    info.update(
        {
            "size_bytes": stat.st_size,
            "mtime_epoch": stat.st_mtime,
        }
    )
    if with_hash:
        info["sha256"] = sha256_file(path)
    return info


def load_index_groups(path: Path) -> list[str]:
    data = json.loads(path.read_text())
    groups = data.get("group_names") or data.get("groups") or []
    if groups and isinstance(groups[0], dict):
        return [str(row.get("group_name") or row.get("target_id")) for row in groups]
    return [str(x) for x in groups]


def parse_indices(text: str) -> list[int]:
    if text is None:
        return []
    text = str(text).strip()
    if not text:
        return []
    return [int(x) for x in text.split(",") if x != ""]


def selected_indices_from_curves(path: Path) -> dict[tuple[str, str, str, int, int], list[int]]:
    """Map (target, region, strategy, k, replicate) to selected indices."""
    out: dict[tuple[str, str, str, int, int], list[int]] = {}
    for row in read_tsv(path):
        if row.get("metric") != "frame_aligned_ca":
            continue
        try:
            key = (
                row["target_id"],
                row["region"],
                row["strategy"],
                int(row["k"]),
                int(row.get("replicate") or 0),
            )
        except (KeyError, ValueError):
            continue
        out[key] = parse_indices(row.get("selected_indices", ""))
    return out


def load_ref_bank(path: Path) -> pd.DataFrame:
    arr = np.load(path, allow_pickle=False)
    required = ["X", "group_names", "pdb", "antigen_name", "antigen_species"]
    missing = [key for key in required if key not in arr.files]
    if missing:
        raise RuntimeError(f"Reference NPZ missing required arrays: {missing}")
    df = pd.DataFrame(
        {
            "group_name": arr["group_names"].astype(str),
            "pdb": arr["pdb"].astype(str),
            "antigen_name": arr["antigen_name"].astype(str),
            "antigen_species": arr["antigen_species"].astype(str),
        }
    )
    df["positive_key"] = (
        df["antigen_name"].fillna("").astype(str).str.lower()
        + "||"
        + df["antigen_species"].fillna("").astype(str).str.lower()
    )
    df.attrs["X"] = arr["X"].astype(np.float32)
    return df


def normalize_rows(x: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(x, axis=1, keepdims=True)
    denom[denom == 0] = 1.0
    return x / denom


def load_shard(path: Path) -> dict[str, object]:
    try:
        return torch.load(path, map_location="cpu")
    except TypeError:
        return torch.load(path, map_location="cpu")


def chain_representations(shard: dict[str, object]) -> dict[str, dict[str, torch.Tensor]]:
    for key in ("chain_representations", "chain_repr", "chains"):
        value = shard.get(key)
        if isinstance(value, dict):
            return value  # type: ignore[return-value]
    return {k: v for k, v in shard.items() if k in CHAINS and isinstance(v, dict)}  # type: ignore[return-value]


def n_conformers_in_shard(shard: dict[str, object]) -> int:
    reps = chain_representations(shard)
    for chain in CHAINS:
        if chain in reps and "mca_repr" in reps[chain]:
            return int(reps[chain]["mca_repr"].shape[0])
    raise RuntimeError("No H/L mca_repr found in shard")


def pool_shard(shard: dict[str, object], selected: list[int]) -> np.ndarray:
    reps = chain_representations(shard)
    vecs: list[np.ndarray] = []
    for chain in CHAINS:
        if chain not in reps or "mca_repr" not in reps[chain]:
            vecs.append(np.zeros(MCA_DIM, dtype=np.float32))
            continue
        tensor = reps[chain]["mca_repr"]
        n_conf = int(tensor.shape[0])
        valid = [idx for idx in selected if 0 <= idx < n_conf]
        if not valid:
            raise RuntimeError(f"No valid selected indices for chain {chain}; n_conf={n_conf}")
        pooled = tensor[valid].to(torch.float32).mean(dim=(0, 1)).cpu().numpy()
        vecs.append(pooled.astype(np.float32))
    return np.concatenate(vecs, axis=0).astype(np.float32)


def condition_indices(
    target_id: str,
    n_conf: int,
    cdr_index: dict[tuple[str, str, str, int, int], list[int]],
    random_replicates: int,
) -> list[dict[str, object]]:
    conditions: list[dict[str, object]] = [
        {
            "condition": "full_128",
            "family": "full",
            "region": "all",
            "strategy": "all",
            "k": n_conf,
            "replicate": 0,
            "indices": list(range(n_conf)),
        },
        {
            "condition": "single_first",
            "family": "single",
            "region": "none",
            "strategy": "first",
            "k": 1,
            "replicate": 0,
            "indices": [0],
        },
    ]
    for region in ("H-CDR3", "all-CDRs"):
        for k in (32, 64):
            strategy = "greedy_kcenter"
            indices = cdr_index.get((target_id, region, strategy, k, 0), [])
            if indices:
                conditions.append(
                    {
                        "condition": f"{region.lower().replace('-', '_')}_{strategy}_k{k}",
                        "family": "cdr_kcenter",
                        "region": region,
                        "strategy": strategy,
                        "k": k,
                        "replicate": 0,
                        "indices": indices,
                    }
                )
    for k in (32, 64):
        indices = cdr_index.get((target_id, "H-CDR3", "first", k, 0), [])
        if indices:
            conditions.append(
                {
                    "condition": f"first_k{k}",
                    "family": "first",
                    "region": "none",
                    "strategy": "first",
                    "k": k,
                    "replicate": 0,
                    "indices": indices,
                }
            )
        for rep in range(random_replicates):
            indices = cdr_index.get((target_id, "H-CDR3", "random", k, rep), [])
            if indices:
                conditions.append(
                    {
                        "condition": f"random_k{k}_rep{rep}",
                        "family": "random",
                        "region": "none",
                        "strategy": "random",
                        "k": k,
                        "replicate": rep,
                        "indices": indices,
                    }
                )
    return conditions


def rank_query(
    query_vec: np.ndarray,
    target_id: str,
    ref_df: pd.DataFrame,
    ref_norm: np.ndarray,
) -> dict[str, object]:
    target_l = target_id.lower()
    query_rows = ref_df[ref_df["group_name"].str.lower() == target_l]
    if query_rows.empty:
        query_rows = ref_df[ref_df["pdb"].str.lower() == target_l]
    has_ref_label = not query_rows.empty
    positive_key = str(query_rows.iloc[0]["positive_key"]) if has_ref_label else ""

    q = query_vec.astype(np.float32)
    denom = np.linalg.norm(q)
    if denom == 0:
        denom = 1.0
    sims = ref_norm @ (q / denom)

    self_mask = (ref_df["group_name"].str.lower().to_numpy() == target_l) | (
        ref_df["pdb"].str.lower().to_numpy() == target_l
    )
    candidate_mask = ~self_mask
    positive_mask = (
        (ref_df["positive_key"].to_numpy() == positive_key) & candidate_mask
        if positive_key
        else np.zeros(len(ref_df), dtype=bool)
    )
    candidate_indices = np.where(candidate_mask)[0]
    ordered = candidate_indices[np.argsort(-sims[candidate_indices])]
    top = int(ordered[0]) if ordered.size else -1

    positive_positions = np.where(positive_mask[ordered])[0]
    if positive_positions.size:
        first_positive_rank = int(positive_positions[0] + 1)
        reciprocal_rank = 1.0 / first_positive_rank
    else:
        first_positive_rank = ""
        reciprocal_rank = ""

    def recall_at(k: int) -> object:
        if not positive_mask.any():
            return ""
        return bool(positive_mask[ordered[:k]].any())

    result: dict[str, object] = {
        "has_ref_label": has_ref_label,
        "positive_key": positive_key,
        "n_ref_candidates": len(ref_df),
        "n_candidates_after_self_exclusion": int(candidate_mask.sum()),
        "n_self_excluded": int(self_mask.sum()),
        "n_positive_after_self_exclusion": int(positive_mask.sum()),
        "first_positive_rank": first_positive_rank,
        "reciprocal_rank": reciprocal_rank,
        "recall_at_1": recall_at(1),
        "recall_at_5": recall_at(5),
        "recall_at_10": recall_at(10),
    }
    if top >= 0:
        top_row = ref_df.iloc[top]
        result.update(
            {
                "top1_group": top_row["group_name"],
                "top1_pdb": top_row["pdb"],
                "top1_antigen_name": top_row["antigen_name"],
                "top1_antigen_species": top_row["antigen_species"],
                "top1_similarity": float(sims[top]),
            }
        )
    return result


def summarize_results(results: pd.DataFrame, preservation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for condition, group in results.groupby("condition", sort=False):
        valid = group[group["n_positive_after_self_exclusion"].astype(int) > 0]
        pres = preservation[preservation["condition"] == condition]
        row: dict[str, object] = {
            "condition": condition,
            "n_queries": len(group),
            "n_with_ref_label": int(group["has_ref_label"].astype(bool).sum()),
            "n_with_positive_after_self_exclusion": len(valid),
            "median_selected_n": float(pd.to_numeric(group["selected_n"]).median()),
            "mean_full_vector_cosine": "",
            "median_full_vector_cosine": "",
        }
        if not valid.empty:
            for col in ("recall_at_1", "recall_at_5", "recall_at_10"):
                row[f"mean_{col}"] = float(valid[col].astype(bool).mean())
            rr = pd.to_numeric(valid["reciprocal_rank"], errors="coerce")
            row["mrr"] = float(rr.mean())
            ranks = pd.to_numeric(valid["first_positive_rank"], errors="coerce")
            row["median_first_positive_rank"] = float(ranks.median())
        else:
            for col in ("recall_at_1", "recall_at_5", "recall_at_10"):
                row[f"mean_{col}"] = ""
            row["mrr"] = ""
            row["median_first_positive_rank"] = ""
        if not pres.empty:
            values = pd.to_numeric(pres["cosine_to_full"], errors="coerce").dropna()
            if not values.empty:
                row["mean_full_vector_cosine"] = float(values.mean())
                row["median_full_vector_cosine"] = float(values.median())
        rows.append(row)
    return pd.DataFrame(rows)


def write_summary_md(
    path: Path,
    audit: dict[str, object],
    summary: pd.DataFrame,
    overlap_counts: dict[str, object],
) -> None:
    def markdown_table(df: pd.DataFrame) -> list[str]:
        if df.empty:
            return ["No summary rows were produced."]
        columns = [str(c) for c in df.columns]
        rows = ["| " + " | ".join(columns) + " |"]
        rows.append("| " + " | ".join("---" for _ in columns) + " |")
        for _, row in df.iterrows():
            values = [str(row.get(col, "")).replace("\n", " ") for col in df.columns]
            rows.append("| " + " | ".join(values) + " |")
        return rows

    lines = [
        "# Phase 1.5a Endpoint Audit And Retrieval Smoke",
        "",
        "## Scope",
        "",
        "CPU-only Jerry-owned wrapper. Gaeun/shared files were read as immutable inputs; outputs were written only under the run directory.",
        "",
        "## Checkpoint",
        "",
        f"- Confirmed checkpoint: `{audit['confirmed_checkpoint']['path']}`",
        f"- Staged Phase 1.4 checkpoint: `{audit['phase14_checkpoint']['path']}`",
        f"- Hash match: `{audit['checkpoint_hash_match']}`",
        f"- Existing Phase 1.4 embeddings reused: `{audit['reuse_phase14_embeddings']}`",
        "",
        "## Endpoint",
        "",
        f"- Reference bank: `{audit['reference_npz']}`",
        f"- Positive definition: exact lowercased `antigen_name || antigen_species`, excluding the query PDB/group itself.",
        f"- Retrieval metric: cosine nearest neighbor over pooled H|L MCA vectors.",
        "",
        "## Overlap",
        "",
        f"- Phase 1.4 targets: `{overlap_counts['phase14_target_count']}`",
        f"- Reference-bank overlap: `{overlap_counts['phase14_ref_overlap_count']}`",
        f"- Strict-300 query-id overlap with June24 query embeddings: `{overlap_counts['strict300_june24_query_overlap_count']}`",
        f"- Strict-300 PDB overlap with reference bank: `{overlap_counts['strict300_ref_pdb_overlap_count']}`",
        "",
        "## Condition Summary",
        "",
    ]
    lines.extend(markdown_table(summary))
    lines.extend(
        [
            "",
            "## Interpretation Guardrail",
            "",
            "This is an endpoint-feasibility and plumbing result, not an independent validation claim. The 149-target PH/AF3 set largely overlaps the reference bank, so self-matches are excluded and metrics depend on whether other antibodies share the same antigen label.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase14-emb-dir", type=Path, default=DEFAULT_PHASE14_EMB)
    parser.add_argument("--phase14-checkpoint", type=Path, default=DEFAULT_PHASE14_CKPT)
    parser.add_argument("--confirmed-checkpoint", type=Path, default=DEFAULT_CONFIRMED_CKPT)
    parser.add_argument("--reference-npz", type=Path, default=DEFAULT_REF_NPZ)
    parser.add_argument("--june24-query-index", type=Path, default=DEFAULT_JUNE24_QUERY_INDEX)
    parser.add_argument("--cdr-curves", type=Path, default=DEFAULT_CDR_CURVES)
    parser.add_argument("--medium-targets", type=Path, default=DEFAULT_MEDIUM_TARGETS)
    parser.add_argument("--strict300", type=Path, default=DEFAULT_STRICT300)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--limit-targets", type=int, default=0)
    parser.add_argument("--random-replicates", type=int, default=5)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    phase14_index = args.phase14_emb_dir / "index.json"
    target_ids = load_index_groups(phase14_index)
    if args.limit_targets:
        target_ids = target_ids[: args.limit_targets]

    ref_df = load_ref_bank(args.reference_npz)
    ref_x = ref_df.attrs["X"]
    ref_norm = normalize_rows(ref_x)
    ref_ids = set(ref_df["group_name"].str.lower()) | set(ref_df["pdb"].str.lower())

    june24_query_ids = set()
    if args.june24_query_index.exists():
        june24_query_ids = {x.lower() for x in load_index_groups(args.june24_query_index)}

    strict_rows = read_tsv(args.strict300) if args.strict300.exists() else []
    strict_query_ids = {row.get("query_id", "").lower() for row in strict_rows if row.get("query_id")}
    strict_pdb_ids = {row.get("pdb_id", "").lower() for row in strict_rows if row.get("pdb_id")}

    cdr_index = selected_indices_from_curves(args.cdr_curves)

    confirmed_info = file_info(args.confirmed_checkpoint)
    phase14_ckpt_info = file_info(args.phase14_checkpoint)
    checkpoint_hash_match = (
        confirmed_info.get("sha256") == phase14_ckpt_info.get("sha256")
        and confirmed_info.get("sha256") is not None
    )

    overlap_rows: list[dict[str, object]] = []
    result_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    preservation_rows: list[dict[str, object]] = []

    for idx, target in enumerate(target_ids, start=1):
        shard_path = args.phase14_emb_dir / f"{target}.pt"
        if not shard_path.exists():
            raise FileNotFoundError(f"Missing Phase 1.4 shard: {shard_path}")
        shard = load_shard(shard_path)
        n_conf = n_conformers_in_shard(shard)

        ref_matches = ref_df[
            (ref_df["group_name"].str.lower() == target.lower())
            | (ref_df["pdb"].str.lower() == target.lower())
        ]
        if not ref_matches.empty:
            ref_row = ref_matches.iloc[0]
            positive_key = ref_row["positive_key"]
            n_same_label = int((ref_df["positive_key"] == positive_key).sum())
        else:
            ref_row = None
            positive_key = ""
            n_same_label = 0
        overlap_rows.append(
            {
                "target_id": target,
                "in_reference_bank": bool(target.lower() in ref_ids),
                "reference_antigen_name": "" if ref_row is None else ref_row["antigen_name"],
                "reference_antigen_species": "" if ref_row is None else ref_row["antigen_species"],
                "same_label_reference_count": n_same_label,
                "same_label_after_self_exclusion": max(0, n_same_label - len(ref_matches)),
                "in_june24_query_index": bool(target.lower() in june24_query_ids),
                "in_strict300_query_id": bool(target.lower() in strict_query_ids),
                "in_strict300_pdb_id": bool(target.lower() in strict_pdb_ids),
            }
        )

        condition_defs = condition_indices(target, n_conf, cdr_index, args.random_replicates)
        vectors: dict[str, np.ndarray] = {}
        for condition in condition_defs:
            selected = condition["indices"]
            assert isinstance(selected, list)
            vec = pool_shard(shard, selected)
            condition_name = str(condition["condition"])
            vectors[condition_name] = vec

            manifest_rows.append(
                {
                    "target_id": target,
                    "condition": condition_name,
                    "family": condition["family"],
                    "region": condition["region"],
                    "strategy": condition["strategy"],
                    "k": condition["k"],
                    "replicate": condition["replicate"],
                    "selected_n": len(selected),
                    "selected_indices": ",".join(str(i) for i in selected),
                }
            )

            ranked = rank_query(vec, target, ref_df, ref_norm)
            result_rows.append(
                {
                    "target_id": target,
                    "condition": condition_name,
                    "family": condition["family"],
                    "region": condition["region"],
                    "strategy": condition["strategy"],
                    "k": condition["k"],
                    "replicate": condition["replicate"],
                    "selected_n": len(selected),
                    **ranked,
                }
            )

        full_vec = vectors.get("full_128")
        if full_vec is not None:
            full_norm = np.linalg.norm(full_vec) or 1.0
            for condition_name, vec in vectors.items():
                denom = (np.linalg.norm(vec) or 1.0) * full_norm
                preservation_rows.append(
                    {
                        "target_id": target,
                        "condition": condition_name,
                        "cosine_to_full": float(np.dot(vec, full_vec) / denom),
                        "l2_to_full": float(np.linalg.norm(vec - full_vec)),
                    }
                )

        print(f"[phase15a] processed {idx}/{len(target_ids)} {target}", flush=True)

    overlap_counts = {
        "phase14_target_count": len(target_ids),
        "phase14_ref_overlap_count": sum(bool(r["in_reference_bank"]) for r in overlap_rows),
        "phase14_june24_query_overlap_count": sum(bool(r["in_june24_query_index"]) for r in overlap_rows),
        "strict300_count": len(strict_rows),
        "strict300_june24_query_overlap_count": len(strict_query_ids & june24_query_ids),
        "strict300_ref_pdb_overlap_count": len(strict_pdb_ids & ref_ids),
    }

    audit = {
        "run_dir": str(args.out_dir),
        "phase14_embedding_dir": str(args.phase14_emb_dir),
        "phase14_index": str(phase14_index),
        "phase14_checkpoint": phase14_ckpt_info,
        "confirmed_checkpoint": confirmed_info,
        "checkpoint_hash_match": checkpoint_hash_match,
        "reuse_phase14_embeddings": checkpoint_hash_match,
        "reference_npz": str(args.reference_npz),
        "reference_rows": len(ref_df),
        "june24_query_index": str(args.june24_query_index),
        "june24_query_count": len(june24_query_ids),
        "cdr_curves": str(args.cdr_curves),
        "strict300": str(args.strict300),
        "overlap_counts": overlap_counts,
        "positive_definition": "exact lowercased antigen_name || antigen_species, self PDB/group excluded",
        "retrieval_metric": "cosine similarity over pooled H|L MCA mca_repr vectors",
        "read_only_inputs": [
            str(args.confirmed_checkpoint),
            str(args.reference_npz),
            str(args.june24_query_index),
            str(args.cdr_curves),
        ],
    }

    results = pd.DataFrame(result_rows)
    preservation = pd.DataFrame(preservation_rows)
    summary = summarize_results(results, preservation)

    (args.out_dir / "endpoint_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True))
    write_tsv(
        args.out_dir / "target_overlap.tsv",
        overlap_rows,
        [
            "target_id",
            "in_reference_bank",
            "reference_antigen_name",
            "reference_antigen_species",
            "same_label_reference_count",
            "same_label_after_self_exclusion",
            "in_june24_query_index",
            "in_strict300_query_id",
            "in_strict300_pdb_id",
        ],
    )
    write_tsv(
        args.out_dir / "condition_manifest.tsv",
        manifest_rows,
        [
            "target_id",
            "condition",
            "family",
            "region",
            "strategy",
            "k",
            "replicate",
            "selected_n",
            "selected_indices",
        ],
    )
    results.to_csv(args.out_dir / "retrieval_smoke_results.tsv", sep="\t", index=False)
    preservation.to_csv(args.out_dir / "subset_vector_preservation.tsv", sep="\t", index=False)
    summary.to_csv(args.out_dir / "retrieval_smoke_summary.tsv", sep="\t", index=False)
    write_summary_md(args.out_dir / "retrieval_smoke_summary.md", audit, summary, overlap_counts)

    print("[phase15a] wrote", args.out_dir, flush=True)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
