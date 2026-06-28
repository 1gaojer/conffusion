#!/usr/bin/env python3
"""Phase 1.5d CDR readout sanity-control diagnostic.

This runner reuses the exported 149-target Phase 1.4 MCA tensors and CDR
numbering assignments. It tests whether the Phase 1.5c CDR-aware retrieval
gain is region-specific by comparing CDR readouts against framework-only,
random framework windows, and shuffled-label controls.
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
DEFAULT_OUT = Path("/project/liulab/jg1920/conffusion/phase15d_20260627_region_controls")

CHAINS = ("H", "L")
BASE_READOUTS = (
    "global_hl",
    "h_cdr1",
    "h_cdr2",
    "h_cdr3",
    "l_cdr1",
    "l_cdr2",
    "l_cdr3",
    "h_all_cdrs",
    "l_all_cdrs",
    "all_cdrs_hl",
    "all_cdrs_mean_std",
    "framework_hl",
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
    return sorted({idx for idx in positions if 0 <= idx < seq_len})


def framework_positions(region_map: dict[str, dict[str, list[int]]], chain: str, seq_len: int) -> list[int]:
    cdr_regions = [f"{chain}-CDR1", f"{chain}-CDR2", f"{chain}-CDR3"]
    cdr_positions = set(region_positions(region_map, chain, cdr_regions, seq_len))
    return [idx for idx in range(seq_len) if idx not in cdr_positions]


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


def choose_framework_window(
    seq_len: int,
    allowed_positions: list[int],
    window_len: int,
    rng: np.random.Generator,
) -> list[int]:
    if window_len <= 0:
        return []
    allowed = set(allowed_positions)
    starts = [
        start
        for start in range(0, max(seq_len - window_len + 1, 0))
        if all(pos in allowed for pos in range(start, start + window_len))
    ]
    if starts:
        start = int(rng.choice(starts))
        return list(range(start, start + window_len))
    if len(allowed_positions) >= window_len:
        return sorted(rng.choice(allowed_positions, size=window_len, replace=False).tolist())
    return sorted(allowed_positions)


def make_target_windows(
    reps: dict[str, dict[str, torch.Tensor]],
    per_target_regions: dict[str, dict[str, list[int]]],
    rng: np.random.Generator,
    n_random_reps: int,
) -> dict[str, list[int]]:
    h_tensor = reps["H"]["mca_repr"]
    h_len = int(h_tensor.shape[1])
    h_cdr3 = region_positions(per_target_regions, "H", ["H-CDR3"], h_len)
    h_framework = framework_positions(per_target_regions, "H", h_len)
    window_len = max(len(h_cdr3), 1)
    return {
        f"random_h_framework_window_hcdr3_len_r{rep}": choose_framework_window(
            h_len, h_framework, window_len, rng
        )
        for rep in range(1, n_random_reps + 1)
    }


def build_readout_vectors(
    shard: dict[str, object],
    region_map: dict[str, dict[str, list[int]]],
    selected: list[int],
    random_windows: dict[str, list[int]],
) -> dict[str, np.ndarray]:
    reps = chain_reps(shard)
    global_parts = []
    all_cdr_parts = []
    all_cdr_std_parts = []
    framework_parts = []
    single_regions: dict[str, np.ndarray] = {}

    for chain in CHAINS:
        tensor = reps[chain]["mca_repr"]
        seq_len = int(tensor.shape[1])
        global_parts.append(pool_chain_global(tensor, selected))
        cdr_regions = [f"{chain}-CDR1", f"{chain}-CDR2", f"{chain}-CDR3"]
        cdr_positions = region_positions(region_map, chain, cdr_regions, seq_len)
        fw_positions = framework_positions(region_map, chain, seq_len)
        all_cdr_parts.append(pool_chain_region(tensor, selected, cdr_positions, stat="mean"))
        all_cdr_std_parts.append(pool_chain_region(tensor, selected, cdr_positions, stat="std"))
        framework_parts.append(pool_chain_region(tensor, selected, fw_positions, stat="mean"))
        for cdr in cdr_regions:
            positions = region_positions(region_map, chain, [cdr], seq_len)
            single_regions[cdr.lower().replace("-", "_")] = pool_chain_region(tensor, selected, positions)

    global_hl = np.concatenate(global_parts).astype(np.float32)
    all_cdrs_hl = np.concatenate(all_cdr_parts).astype(np.float32)
    all_cdrs_mean_std = np.concatenate(all_cdr_parts + all_cdr_std_parts).astype(np.float32)
    framework_hl = np.concatenate(framework_parts).astype(np.float32)

    out = {
        "global_hl": global_hl,
        "h_cdr1": single_regions["h_cdr1"],
        "h_cdr2": single_regions["h_cdr2"],
        "h_cdr3": single_regions["h_cdr3"],
        "l_cdr1": single_regions["l_cdr1"],
        "l_cdr2": single_regions["l_cdr2"],
        "l_cdr3": single_regions["l_cdr3"],
        "h_all_cdrs": all_cdr_parts[0],
        "l_all_cdrs": all_cdr_parts[1],
        "all_cdrs_hl": all_cdrs_hl,
        "all_cdrs_mean_std": all_cdrs_mean_std,
        "framework_hl": framework_hl,
    }
    h_tensor = reps["H"]["mca_repr"]
    for name, positions in random_windows.items():
        out[name] = pool_chain_region(h_tensor, selected, positions)
    return out


def normalize_rows(x: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(x, axis=1, keepdims=True)
    denom[denom == 0] = 1.0
    return x / denom


def rank_against_bank(
    query_vec: np.ndarray,
    target_id: str,
    bank_ids: list[str],
    label_by_target: dict[str, str],
    bank_labels: list[str],
    bank_matrix_norm: np.ndarray,
) -> dict[str, object]:
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


def summarize(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (readout, family), group in results.groupby(["readout", "readout_family"], sort=False):
        valid = group[group["n_positive_after_self_exclusion"] > 0]
        rows.append(
            {
                "readout": readout,
                "readout_family": family,
                "n_queries": len(group),
                "n_eval": len(valid),
                "recall_at_1": bool_series(valid["recall_at_1"]).mean() if not valid.empty else np.nan,
                "recall_at_5": bool_series(valid["recall_at_5"]).mean() if not valid.empty else np.nan,
                "recall_at_10": bool_series(valid["recall_at_10"]).mean() if not valid.empty else np.nan,
                "mrr": pd.to_numeric(valid["reciprocal_rank"], errors="coerce").mean(),
                "median_first_positive_rank": pd.to_numeric(valid["first_positive_rank"], errors="coerce").median(),
            }
        )
    return pd.DataFrame(rows)


def shuffled_label_controls(
    readout_vectors: dict[tuple[str, str], np.ndarray],
    target_ids: list[str],
    label_by_target: dict[str, str],
    readouts: list[str],
    n_reps: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    labels = np.array([label_by_target.get(target_id, "") for target_id in target_ids], dtype=object)
    rows = []
    for rep in range(1, n_reps + 1):
        shuffled = labels.copy()
        rng.shuffle(shuffled)
        shuffled_by_target = dict(zip(target_ids, shuffled.tolist(), strict=False))
        for readout in readouts:
            bank_matrix = np.stack([readout_vectors[(readout, target_id)] for target_id in target_ids]).astype(np.float32)
            bank_norm = normalize_rows(bank_matrix)
            bank_labels = [shuffled_by_target.get(target_id, "") for target_id in target_ids]
            for target_id in target_ids:
                ranked = rank_against_bank(
                    readout_vectors[(readout, target_id)],
                    target_id,
                    target_ids,
                    shuffled_by_target,
                    bank_labels,
                    bank_norm,
                )
                rows.append(
                    {
                        "target_id": target_id,
                        "readout": readout,
                        "shuffle_replicate": rep,
                        **ranked,
                    }
                )
    return pd.DataFrame(rows)


def summarize_shuffled(results: pd.DataFrame) -> pd.DataFrame:
    per_rep = []
    for (readout, rep), group in results.groupby(["readout", "shuffle_replicate"], sort=False):
        valid = group[group["n_positive_after_self_exclusion"] > 0]
        per_rep.append(
            {
                "readout": readout,
                "shuffle_replicate": rep,
                "n_eval": len(valid),
                "recall_at_10": bool_series(valid["recall_at_10"]).mean() if not valid.empty else np.nan,
                "mrr": pd.to_numeric(valid["reciprocal_rank"], errors="coerce").mean(),
                "median_first_positive_rank": pd.to_numeric(valid["first_positive_rank"], errors="coerce").median(),
            }
        )
    per_rep_df = pd.DataFrame(per_rep)
    rows = []
    for readout, group in per_rep_df.groupby("readout", sort=False):
        row: dict[str, object] = {"readout": readout, "n_shuffle_replicates": len(group)}
        for metric in ("recall_at_10", "mrr", "median_first_positive_rank"):
            vals = pd.to_numeric(group[metric], errors="coerce").dropna()
            row[f"{metric}_mean"] = vals.mean()
            row[f"{metric}_sd"] = vals.std(ddof=1) if len(vals) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def readout_family(readout: str) -> str:
    if readout == "global_hl":
        return "global"
    if readout == "framework_hl":
        return "framework"
    if readout.startswith("random_h_framework_window"):
        return "random_framework_window"
    if readout.startswith("h_cdr"):
        return "h_cdr_single"
    if readout.startswith("l_cdr"):
        return "l_cdr_single"
    if readout in {"h_all_cdrs", "l_all_cdrs", "all_cdrs_hl", "all_cdrs_mean_std"}:
        return "cdr_combined"
    return "other"


def write_report(
    path: Path,
    summary: pd.DataFrame,
    random_rollup: pd.DataFrame,
    shuffled_summary: pd.DataFrame,
    region_audit: pd.DataFrame,
    audit: dict[str, object],
) -> None:
    sorted_summary = summary.sort_values("mrr", ascending=False)
    global_row = summary[summary["readout"] == "global_hl"][["readout", "recall_at_10", "mrr", "median_first_positive_rank"]]
    rows = [
        "# Phase 1.5d Region-Control Readout Diagnostic",
        "",
        "## Scope",
        "",
        str(audit.get("candidate_bank", "Internal diagnostic using saved MCA tensors. It compares CDR readouts against framework-only, random framework-window, and shuffled-label controls.")),
        "",
        "## Inputs",
        "",
        f"- Phase 1.5a run: `{audit['phase15a_run']}`",
        f"- Embedding shards: `{audit['embedding_dir']}`",
        f"- CDR numbering assignments: `{audit['cdr_assignments']}`",
        f"- Random framework-window replicates: {audit['n_random_window_reps']}",
        f"- Shuffled-label replicates: {audit['n_shuffle_reps']}",
        "",
        "## Global Baseline",
        "",
        markdown_table(global_row),
        "",
        "## Readout Summary",
        "",
        markdown_table(sorted_summary[["readout", "readout_family", "n_eval", "recall_at_10", "mrr", "median_first_positive_rank"]]),
        "",
        "## Random Framework-Window Rollup",
        "",
        markdown_table(random_rollup),
        "",
        "## Shuffled-Label Control Rollup",
        "",
        markdown_table(shuffled_summary),
        "",
        "## Region Extraction Audit Sample",
        "",
        markdown_table(region_audit.head(18)),
        "",
        "## Interpretation Guardrail",
        "",
        "The CDR-readout result is stronger if true CDR regions beat framework-only, random framework windows, and shuffled labels. Treat this as internal evidence that should be repeated on strict/de-overlapped targets before external validation claims.",
        "",
    ]
    path.write_text("\n".join(rows))


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase15a-run", type=Path, default=DEFAULT_PHASE15A_RUN)
    parser.add_argument("--embedding-dir", type=Path, default=DEFAULT_EMB_DIR)
    parser.add_argument("--cdr-assignments", type=Path, default=DEFAULT_CDR_ASSIGNMENTS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--limit-targets", type=int, default=0)
    parser.add_argument("--n-random-window-reps", type=int, default=10)
    parser.add_argument("--n-shuffle-reps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260627)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    target_overlap = pd.read_csv(args.phase15a_run / "target_overlap.tsv", sep="\t")
    label_by_target = labels_from_target_overlap(target_overlap)

    manifest = pd.read_csv(args.phase15a_run / "condition_manifest.tsv", sep="\t")
    target_ids = list(dict.fromkeys(manifest["target_id"].astype(str).tolist()))
    if args.limit_targets:
        target_ids = target_ids[: args.limit_targets]

    region_map = load_region_map(args.cdr_assignments)
    rng = np.random.default_rng(args.seed)
    readout_vectors: dict[tuple[str, str], np.ndarray] = {}
    region_rows = []
    readouts_seen: list[str] = []

    for idx, target_id in enumerate(target_ids, start=1):
        shard = load_shard(args.embedding_dir / f"{target_id}.pt")
        reps = chain_reps(shard)
        per_target_regions = region_map.get(target_id, {})
        n_conf = int(reps["H"]["mca_repr"].shape[0])
        full_selected = list(range(n_conf))
        random_windows = make_target_windows(reps, per_target_regions, rng, args.n_random_window_reps)
        readouts = build_readout_vectors(shard, per_target_regions, full_selected, random_windows)
        if not readouts_seen:
            readouts_seen = list(readouts)
        for readout, vec in readouts.items():
            readout_vectors[(readout, target_id)] = vec

        for chain in CHAINS:
            seq_len = int(reps[chain]["mca_repr"].shape[1])
            for region in [f"{chain}-CDR1", f"{chain}-CDR2", f"{chain}-CDR3"]:
                positions = region_positions(per_target_regions, chain, [region], seq_len)
                region_rows.append(
                    {
                        "target_id": target_id,
                        "chain": chain,
                        "region": region,
                        "seq_len": seq_len,
                        "n_positions": len(positions),
                        "start_1based": min(positions) + 1 if positions else "",
                        "end_1based": max(positions) + 1 if positions else "",
                        "positions_1based": ",".join(str(pos + 1) for pos in positions),
                    }
                )
            fw_positions = framework_positions(per_target_regions, chain, seq_len)
            region_rows.append(
                {
                    "target_id": target_id,
                    "chain": chain,
                    "region": f"{chain}-framework",
                    "seq_len": seq_len,
                    "n_positions": len(fw_positions),
                    "start_1based": min(fw_positions) + 1 if fw_positions else "",
                    "end_1based": max(fw_positions) + 1 if fw_positions else "",
                    "positions_1based": "",
                }
            )
        for name, positions in random_windows.items():
            region_rows.append(
                {
                    "target_id": target_id,
                    "chain": "H",
                    "region": name,
                    "seq_len": int(reps["H"]["mca_repr"].shape[1]),
                    "n_positions": len(positions),
                    "start_1based": min(positions) + 1 if positions else "",
                    "end_1based": max(positions) + 1 if positions else "",
                    "positions_1based": ",".join(str(pos + 1) for pos in positions),
                }
            )
        print(f"[phase15d] pooled {idx}/{len(target_ids)} {target_id}", flush=True)

    result_rows = []
    bank_ids = target_ids
    bank_labels = [label_by_target.get(target_id, "") for target_id in bank_ids]
    for readout in readouts_seen:
        bank_matrix = np.stack([readout_vectors[(readout, target_id)] for target_id in bank_ids]).astype(np.float32)
        bank_norm = normalize_rows(bank_matrix)
        for target_id in target_ids:
            ranked = rank_against_bank(
                readout_vectors[(readout, target_id)],
                target_id,
                bank_ids,
                label_by_target,
                bank_labels,
                bank_norm,
            )
            result_rows.append(
                {
                    "target_id": target_id,
                    "readout": readout,
                    "readout_family": readout_family(readout),
                    **ranked,
                }
            )

    results = pd.DataFrame(result_rows)
    summary = summarize(results)
    random_summary = summary[summary["readout_family"] == "random_framework_window"].copy()
    random_rollup_rows = []
    if not random_summary.empty:
        row: dict[str, object] = {
            "readout_family": "random_framework_window",
            "n_window_replicates": len(random_summary),
        }
        for metric in ("recall_at_10", "mrr", "median_first_positive_rank"):
            vals = pd.to_numeric(random_summary[metric], errors="coerce").dropna()
            row[f"{metric}_mean"] = vals.mean()
            row[f"{metric}_sd"] = vals.std(ddof=1) if len(vals) > 1 else 0.0
        random_rollup_rows.append(row)
    random_rollup = pd.DataFrame(random_rollup_rows)

    shuffle_readouts = [
        "global_hl",
        "h_cdr3",
        "all_cdrs_hl",
        "all_cdrs_mean_std",
        "framework_hl",
    ]
    shuffled = shuffled_label_controls(
        readout_vectors,
        target_ids,
        label_by_target,
        shuffle_readouts,
        args.n_shuffle_reps,
        args.seed + 1,
    )
    shuffled_summary = summarize_shuffled(shuffled)
    region_audit = pd.DataFrame(region_rows)

    audit = {
        "phase15a_run": str(args.phase15a_run),
        "embedding_dir": str(args.embedding_dir),
        "cdr_assignments": str(args.cdr_assignments),
        "out_dir": str(args.out_dir),
        "target_count": len(target_ids),
        "candidate_bank": "Internal target bank using full-ensemble vectors; self-matches excluded.",
        "readouts": readouts_seen,
        "n_random_window_reps": args.n_random_window_reps,
        "n_shuffle_reps": args.n_shuffle_reps,
        "seed": args.seed,
    }
    (args.out_dir / "phase15d_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True))
    results.to_csv(args.out_dir / "region_control_retrieval_results.tsv", sep="\t", index=False)
    summary.to_csv(args.out_dir / "region_control_summary.tsv", sep="\t", index=False)
    random_rollup.to_csv(args.out_dir / "region_control_random_window_rollup.tsv", sep="\t", index=False)
    shuffled.to_csv(args.out_dir / "region_control_shuffled_label_results.tsv", sep="\t", index=False)
    shuffled_summary.to_csv(args.out_dir / "region_control_shuffled_label_summary.tsv", sep="\t", index=False)
    region_audit.to_csv(args.out_dir / "region_extraction_audit.tsv", sep="\t", index=False)
    write_report(
        args.out_dir / "phase15d_region_control_report.md",
        summary,
        random_rollup,
        shuffled_summary,
        region_audit,
        audit,
    )
    print(f"[phase15d] wrote {args.out_dir}")
    print(summary[["readout", "readout_family", "n_eval", "recall_at_10", "mrr", "median_first_positive_rank"]].to_string(index=False))


if __name__ == "__main__":
    main()
