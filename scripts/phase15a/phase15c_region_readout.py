#!/usr/bin/env python3
"""Phase 1.5c CDR/paratope-aware retrieval readout diagnostic.

This runner uses already-exported MCA shards and CDR numbering assignments.
Candidate vectors are full-ensemble vectors from the same target bank; query
vectors vary by subset condition.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch


DEFAULT_PHASE15A_RUN = Path(
    "/project/liulab/jg1920/conffusion/phase15a_20260627_175558_endpoint_audit"
)
DEFAULT_EMB_DIR = Path(
    "/project/liulab/jg1920/conffusion/phase14_20260626_2301/runs/full/embeddings/mca1000_selected"
)
DEFAULT_CDR_ASSIGNMENTS = Path(
    "/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627/cdr_numbering_assignments.tsv"
)
DEFAULT_OUT = Path("/project/liulab/jg1920/conffusion/phase15c_20260627_region_readout")

CHAINS = ("H", "L")
READOUTS = (
    "global_hl",
    "h_cdr3",
    "all_cdrs_hl",
    "global_plus_h_cdr3",
    "global_plus_all_cdrs",
    "all_cdrs_mean_std",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def parse_indices(text: str) -> list[int]:
    if text is None:
        return []
    text = str(text).strip()
    if not text:
        return []
    return [int(x) for x in text.split(",") if x != ""]


def load_region_map(path: Path) -> dict[str, dict[str, dict[str, list[int]]]]:
    """Return target -> H/L -> region -> zero-based sequence indices."""
    out: dict[str, dict[str, dict[str, list[int]]]] = {}
    for row in read_tsv(path):
        target = row["target_id"]
        chain_role = row.get("chain_role", "").lower()
        if chain_role == "heavy":
            chain = "H"
        elif chain_role == "light":
            chain = "L"
        else:
            chain = row["chain_type"]
            if chain == "K":
                chain = "L"
        region = row["region"]
        try:
            start = int(row["seq_index_start_1based"]) - 1
            end = int(row["seq_index_end_1based"])
        except ValueError:
            continue
        out.setdefault(target, {}).setdefault(chain, {})[region] = list(range(start, end))
    return out


def load_shard(path: Path) -> dict[str, object]:
    return torch.load(path, map_location="cpu")


def chain_reps(shard: dict[str, object]) -> dict[str, dict[str, torch.Tensor]]:
    reps = shard.get("chain_representations")
    if isinstance(reps, dict):
        return reps  # type: ignore[return-value]
    return {k: v for k, v in shard.items() if k in CHAINS and isinstance(v, dict)}  # type: ignore[return-value]


def valid_selected(selected: list[int], n_conf: int) -> list[int]:
    vals = [idx for idx in selected if 0 <= idx < n_conf]
    if not vals:
        raise RuntimeError(f"No valid selected conformer indices; n_conf={n_conf}")
    return vals


def region_positions(region_map: dict[str, dict[str, list[int]]], chain: str, regions: list[str], seq_len: int) -> list[int]:
    positions: list[int] = []
    for region in regions:
        positions.extend(region_map.get(chain, {}).get(region, []))
    uniq = sorted({idx for idx in positions if 0 <= idx < seq_len})
    return uniq


def pool_chain_global(tensor: torch.Tensor, selected: list[int]) -> np.ndarray:
    selected = valid_selected(selected, int(tensor.shape[0]))
    return tensor[selected].to(torch.float32).mean(dim=(0, 1)).cpu().numpy().astype(np.float32)


def pool_chain_region(
    tensor: torch.Tensor,
    selected: list[int],
    positions: list[int],
    *,
    stat: str = "mean",
) -> np.ndarray:
    selected = valid_selected(selected, int(tensor.shape[0]))
    if not positions:
        return np.zeros(int(tensor.shape[-1]), dtype=np.float32)
    data = tensor[selected][:, positions, :].to(torch.float32).reshape(-1, int(tensor.shape[-1]))
    if stat == "std":
        return data.std(dim=0, unbiased=False).cpu().numpy().astype(np.float32)
    return data.mean(dim=0).cpu().numpy().astype(np.float32)


def build_readout_vectors(
    shard: dict[str, object],
    region_map: dict[str, dict[str, list[int]]],
    selected: list[int],
) -> dict[str, np.ndarray]:
    reps = chain_reps(shard)
    global_parts = []
    all_cdr_parts = []
    all_cdr_std_parts = []
    h_cdr3 = np.zeros(256, dtype=np.float32)
    for chain in CHAINS:
        tensor = reps[chain]["mca_repr"]
        global_parts.append(pool_chain_global(tensor, selected))
        cdr_regions = [f"{chain}-CDR1", f"{chain}-CDR2", f"{chain}-CDR3"]
        positions = region_positions(region_map, chain, cdr_regions, int(tensor.shape[1]))
        all_cdr_parts.append(pool_chain_region(tensor, selected, positions, stat="mean"))
        all_cdr_std_parts.append(pool_chain_region(tensor, selected, positions, stat="std"))
        if chain == "H":
            h_positions = region_positions(region_map, "H", ["H-CDR3"], int(tensor.shape[1]))
            h_cdr3 = pool_chain_region(tensor, selected, h_positions, stat="mean")
    global_hl = np.concatenate(global_parts).astype(np.float32)
    all_cdrs_hl = np.concatenate(all_cdr_parts).astype(np.float32)
    all_cdrs_mean_std = np.concatenate(all_cdr_parts + all_cdr_std_parts).astype(np.float32)
    return {
        "global_hl": global_hl,
        "h_cdr3": h_cdr3,
        "all_cdrs_hl": all_cdrs_hl,
        "global_plus_h_cdr3": np.concatenate([global_hl, h_cdr3]).astype(np.float32),
        "global_plus_all_cdrs": np.concatenate([global_hl, all_cdrs_hl]).astype(np.float32),
        "all_cdrs_mean_std": all_cdrs_mean_std,
    }


def normalize_rows(x: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(x, axis=1, keepdims=True)
    denom[denom == 0] = 1.0
    return x / denom


def rank_against_bank(
    query_vec: np.ndarray,
    target_id: str,
    bank_ids: list[str],
    bank_labels: list[str],
    bank_matrix_norm: np.ndarray,
) -> dict[str, object]:
    label_by_target = dict(zip(bank_ids, bank_labels, strict=False))
    positive_key = label_by_target.get(target_id, "")
    q = query_vec.astype(np.float32)
    denom = np.linalg.norm(q)
    if denom == 0:
        denom = 1.0
    sims = bank_matrix_norm @ (q / denom)
    bank_ids_array = np.array(bank_ids)
    self_mask = bank_ids_array == target_id
    candidate_mask = ~self_mask
    positive_mask = np.array([(lab == positive_key and lab != "") for lab in bank_labels]) & candidate_mask
    candidate_indices = np.where(candidate_mask)[0]
    ordered = candidate_indices[np.argsort(-sims[candidate_indices])]
    pos = np.where(positive_mask[ordered])[0]
    if len(pos):
        first_positive_rank: object = int(pos[0] + 1)
        rr: object = 1.0 / int(first_positive_rank)
    else:
        first_positive_rank = ""
        rr = ""
    top = int(ordered[0]) if len(ordered) else -1
    result: dict[str, object] = {
        "positive_key": positive_key,
        "n_candidates_after_self_exclusion": int(candidate_mask.sum()),
        "n_positive_after_self_exclusion": int(positive_mask.sum()),
        "first_positive_rank": first_positive_rank,
        "reciprocal_rank": rr,
        "recall_at_1": bool(len(pos) and pos[0] < 1) if positive_mask.any() else "",
        "recall_at_5": bool(len(pos) and pos[0] < 5) if positive_mask.any() else "",
        "recall_at_10": bool(len(pos) and pos[0] < 10) if positive_mask.any() else "",
    }
    if top >= 0:
        result.update(
            {
                "top1_target_id": bank_ids[top],
                "top1_positive_key": bank_labels[top],
                "top1_similarity": float(sims[top]),
            }
        )
    return result


def summarize(results: pd.DataFrame, preservation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (readout, condition), group in results.groupby(["readout", "condition"], sort=False):
        valid = group[group["n_positive_after_self_exclusion"] > 0]
        pres = preservation[(preservation["readout"] == readout) & (preservation["condition"] == condition)]
        rows.append(
            {
                "readout": readout,
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
                "mean_cosine_to_readout_full": pd.to_numeric(pres["cosine_to_readout_full"], errors="coerce").mean(),
            }
        )
    return pd.DataFrame(rows)


def random_rollup(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    random = summary[summary["family"] == "random"]
    for (readout, k), group in random.groupby(["readout", "k"], sort=False):
        row: dict[str, object] = {"readout": readout, "k": k, "n_replicates": len(group)}
        for metric in ["recall_at_10", "mrr", "median_first_positive_rank", "mean_cosine_to_readout_full"]:
            vals = pd.to_numeric(group[metric], errors="coerce").dropna()
            row[f"{metric}_mean"] = vals.mean()
            row[f"{metric}_sd"] = vals.std(ddof=1) if len(vals) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def delta_vs_global(summary: pd.DataFrame) -> pd.DataFrame:
    base = summary[summary["readout"] == "global_hl"][
        ["condition", "recall_at_10", "mrr", "median_first_positive_rank"]
    ].copy()
    base = base.rename(
        columns={
            "recall_at_10": "global_recall_at_10",
            "mrr": "global_mrr",
            "median_first_positive_rank": "global_median_first_positive_rank",
        }
    )
    out = summary.merge(base, on="condition", how="left")
    out["recall10_minus_global"] = out["recall_at_10"] - out["global_recall_at_10"]
    out["mrr_minus_global"] = out["mrr"] - out["global_mrr"]
    out["median_rank_minus_global"] = (
        out["median_first_positive_rank"] - out["global_median_first_positive_rank"]
    )
    return out


def labels_from_target_overlap(target_overlap: pd.DataFrame) -> dict[str, str]:
    if "positive_key" in target_overlap.columns:
        positive = target_overlap["positive_key"].fillna("").astype(str).str.lower()
    else:
        positive = (
            target_overlap["reference_antigen_name"].fillna("").astype(str).str.lower()
            + "||"
            + target_overlap["reference_antigen_species"].fillna("").astype(str).str.lower()
        )
        positive.loc[target_overlap["reference_antigen_name"].isna()] = ""
    return dict(zip(target_overlap["target_id"].astype(str), positive, strict=False))


def write_report(
    path: Path,
    summary: pd.DataFrame,
    rand: pd.DataFrame,
    deltas: pd.DataFrame,
    audit: dict[str, object],
) -> None:
    full = summary[summary["family"] == "full"].sort_values("mrr", ascending=False)
    k64 = summary[summary["condition"].isin(["h_cdr3_greedy_kcenter_k64", "all_cdrs_greedy_kcenter_k64"])]
    rows = [
        "# Phase 1.5c CDR/Paratope-Aware Readout Diagnostic",
        "",
        "## Scope",
        "",
        str(audit.get("candidate_bank", "Internal retrieval diagnostic using full-ensemble candidate vectors; query vectors vary by subset condition.")),
        "",
        "## Inputs",
        "",
        f"- Phase 1.5a run: `{audit['phase15a_run']}`",
        f"- Embedding shards: `{audit['embedding_dir']}`",
        f"- CDR numbering assignments: `{audit['cdr_assignments']}`",
        "",
        "## Full-Ensemble Readout Ranking",
        "",
        markdown_table(full[["readout", "n_eval", "recall_at_10", "mrr", "median_first_positive_rank"]]),
        "",
        "## K=64 CDR Conditions",
        "",
        markdown_table(k64[["readout", "condition", "n_eval", "recall_at_10", "mrr", "median_first_positive_rank", "mean_cosine_to_readout_full"]]),
        "",
        "## Random Rollup",
        "",
        markdown_table(rand),
        "",
        "## Readout Delta Versus Global H/L",
        "",
        markdown_table(
            deltas[
                [
                    "readout",
                    "condition",
                    "recall10_minus_global",
                    "mrr_minus_global",
                    "median_rank_minus_global",
                ]
            ].sort_values("mrr_minus_global", ascending=False),
            max_rows=30,
        ),
        "",
        "## Interpretation Guardrail",
        "",
        "Use this as a readout diagnostic. It tests whether CDR/paratope-aware pooling changes retrieval behavior on the ready 149-target bank; it does not replace the stricter external or strict-300 endpoint.",
        "",
    ]
    path.write_text("\n".join(rows))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase15a-run", type=Path, default=DEFAULT_PHASE15A_RUN)
    parser.add_argument("--embedding-dir", type=Path, default=DEFAULT_EMB_DIR)
    parser.add_argument("--cdr-assignments", type=Path, default=DEFAULT_CDR_ASSIGNMENTS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--limit-targets", type=int, default=0)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    target_overlap = pd.read_csv(args.phase15a_run / "target_overlap.tsv", sep="\t")
    labels = labels_from_target_overlap(target_overlap)

    manifest = pd.read_csv(args.phase15a_run / "condition_manifest.tsv", sep="\t")
    target_ids = list(dict.fromkeys(manifest["target_id"].astype(str).tolist()))
    if args.limit_targets:
        target_ids = target_ids[: args.limit_targets]
        manifest = manifest[manifest["target_id"].isin(target_ids)].copy()

    region_map = load_region_map(args.cdr_assignments)
    condition_vectors: dict[tuple[str, str, str], np.ndarray] = {}
    full_vectors: dict[tuple[str, str], np.ndarray] = {}
    vector_rows = []

    for idx, target_id in enumerate(target_ids, start=1):
        shard = load_shard(args.embedding_dir / f"{target_id}.pt")
        per_target_regions = region_map.get(target_id, {})
        reps = chain_reps(shard)
        full_selected = list(range(int(reps["H"]["mca_repr"].shape[0])))
        full_readouts = build_readout_vectors(shard, per_target_regions, full_selected)
        for readout, vec in full_readouts.items():
            full_vectors[(readout, target_id)] = vec

        target_manifest = manifest[manifest["target_id"] == target_id]
        for _, row in target_manifest.iterrows():
            selected = parse_indices(row["selected_indices"])
            readouts = build_readout_vectors(shard, per_target_regions, selected)
            for readout, vec in readouts.items():
                condition_vectors[(readout, str(row["condition"]), target_id)] = vec
                full_vec = full_readouts[readout]
                denom = (np.linalg.norm(vec) or 1.0) * (np.linalg.norm(full_vec) or 1.0)
                vector_rows.append(
                    {
                        "target_id": target_id,
                        "readout": readout,
                        "condition": row["condition"],
                        "cosine_to_readout_full": float(np.dot(vec, full_vec) / denom),
                        "l2_to_readout_full": float(np.linalg.norm(vec - full_vec)),
                    }
                )
        print(f"[phase15c] pooled {idx}/{len(target_ids)} {target_id}", flush=True)

    result_rows = []
    bank_ids = target_ids
    bank_labels = [labels.get(target_id, "") for target_id in bank_ids]
    for readout in READOUTS:
        bank_matrix = np.stack([full_vectors[(readout, target_id)] for target_id in bank_ids]).astype(np.float32)
        bank_norm = normalize_rows(bank_matrix)
        for _, row in manifest.iterrows():
            target_id = str(row["target_id"])
            condition = str(row["condition"])
            vec = condition_vectors[(readout, condition, target_id)]
            ranked = rank_against_bank(vec, target_id, bank_ids, bank_labels, bank_norm)
            result_rows.append(
                {
                    "target_id": target_id,
                    "readout": readout,
                    "condition": condition,
                    "family": row["family"],
                    "region": row["region"],
                    "strategy": row["strategy"],
                    "k": row["k"],
                    "replicate": row["replicate"],
                    "selected_n": row["selected_n"],
                    **ranked,
                }
            )

    results = pd.DataFrame(result_rows)
    preservation = pd.DataFrame(vector_rows)
    summary = summarize(results, preservation)
    rand = random_rollup(summary)
    deltas = delta_vs_global(summary)

    audit = {
        "phase15a_run": str(args.phase15a_run),
        "embedding_dir": str(args.embedding_dir),
        "cdr_assignments": str(args.cdr_assignments),
        "out_dir": str(args.out_dir),
        "target_count": len(target_ids),
        "candidate_bank": "Internal target bank using full-ensemble candidate vectors; self-matches excluded.",
        "readouts": list(READOUTS),
    }
    (args.out_dir / "phase15c_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True))
    results.to_csv(args.out_dir / "region_readout_retrieval_results.tsv", sep="\t", index=False)
    preservation.to_csv(args.out_dir / "region_readout_vector_preservation.tsv", sep="\t", index=False)
    summary.to_csv(args.out_dir / "region_readout_condition_summary.tsv", sep="\t", index=False)
    rand.to_csv(args.out_dir / "region_readout_random_rollup.tsv", sep="\t", index=False)
    deltas.to_csv(args.out_dir / "region_readout_delta_vs_global.tsv", sep="\t", index=False)
    write_report(args.out_dir / "phase15c_region_readout_report.md", summary, rand, deltas, audit)
    print(f"[phase15c] wrote {args.out_dir}")
    print(summary[["readout", "condition", "n_eval", "recall_at_10", "mrr"]].to_string(index=False))


if __name__ == "__main__":
    main()
