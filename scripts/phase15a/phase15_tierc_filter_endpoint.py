#!/usr/bin/env python3
"""Filter a Phase 1.5 endpoint to Tier C antigen/source-de-overlapped targets."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_ENDPOINT_DIR = Path(
    "/external/liulab/jg1920/conffusion/phase15_tierb_step4_full_20260627_233600/endpoint"
)
DEFAULT_CHECKPOINT_IDS = Path(
    "/external/liulab/gkim/antigen_prediction/2026.06.19_retrain_mca_1000_confs/split_confs/all_pdb_ids.txt"
)
DEFAULT_ANTIGEN_LABELS = Path(
    "/external/liulab/gkim/antigen_prediction/datasets/antigen_labels/full_sabdab_relabeled_antigens.csv"
)


def read_table(path: Path, delimiter: str = "\t") -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader), list(reader.fieldnames or [])


def write_table(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def norm_text(value: object) -> str:
    return str(value or "").strip().lower()


def label_pair(row: dict[str, str]) -> str:
    name = norm_text(row.get("antigen_name"))
    species = norm_text(row.get("antigen_species"))
    if not name and not species:
        return ""
    return f"{name} || {species}"


def load_label_pairs(path: Path) -> dict[str, set[str]]:
    pairs_by_pdb: dict[str, set[str]] = defaultdict(set)
    rows, _ = read_table(path, delimiter=",")
    for row in rows:
        pdb = norm_text(row.get("pdb"))
        pair = label_pair(row)
        if pdb and pair:
            pairs_by_pdb[pdb].add(pair)
    return pairs_by_pdb


def read_ids(path: Path) -> set[str]:
    return {norm_text(line) for line in path.read_text().splitlines() if line.strip()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint-dir", type=Path, default=DEFAULT_ENDPOINT_DIR)
    parser.add_argument("--checkpoint-ids", type=Path, default=DEFAULT_CHECKPOINT_IDS)
    parser.add_argument("--antigen-labels", type=Path, default=DEFAULT_ANTIGEN_LABELS)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    target_overlap, overlap_fields = read_table(args.endpoint_dir / "target_overlap.tsv")
    condition_manifest, condition_fields = read_table(args.endpoint_dir / "condition_manifest.tsv")
    label_pairs = load_label_pairs(args.antigen_labels)
    checkpoint_ids = read_ids(args.checkpoint_ids)
    checkpoint_pairs = set()
    for pdb_id in checkpoint_ids:
        checkpoint_pairs.update(label_pairs.get(pdb_id, set()))

    kept_ids: list[str] = []
    dropped_rows: list[dict[str, object]] = []
    missing_label_ids: list[str] = []
    overlap_by_target: dict[str, list[str]] = {}

    for row in target_overlap:
        target_id = norm_text(row["target_id"])
        target_pairs = sorted(label_pairs.get(target_id, set()))
        if not target_pairs:
            missing_label_ids.append(target_id)
        overlapping = [pair for pair in target_pairs if pair in checkpoint_pairs]
        overlap_by_target[target_id] = overlapping
        if overlapping:
            dropped_rows.append(
                {
                    "target_id": target_id,
                    "positive_key": row.get("positive_key", ""),
                    "n_target_pairs": len(target_pairs),
                    "n_overlapping_pairs": len(overlapping),
                    "overlapping_pairs": "; ".join(overlapping),
                }
            )
        else:
            kept_ids.append(target_id)

    kept_set = set(kept_ids)
    filtered_overlap = [row for row in target_overlap if norm_text(row["target_id"]) in kept_set]
    filtered_manifest = [row for row in condition_manifest if norm_text(row["target_id"]) in kept_set]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_table(args.out_dir / "target_overlap.tsv", filtered_overlap, overlap_fields)
    write_table(args.out_dir / "condition_manifest.tsv", filtered_manifest, condition_fields)
    (args.out_dir / "tierc_target_ids.txt").write_text("\n".join(kept_ids) + "\n")
    write_table(
        args.out_dir / "tierc_dropped_targets.tsv",
        dropped_rows,
        ["target_id", "positive_key", "n_target_pairs", "n_overlapping_pairs", "overlapping_pairs"],
    )

    label_counts = Counter(row.get("positive_key", "") for row in filtered_overlap)
    label_rows = [
        {"positive_key": label, "target_count": count}
        for label, count in sorted(label_counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    write_table(args.out_dir / "tierc_label_counts.tsv", label_rows, ["positive_key", "target_count"])
    evaluable_target_count = sum(count for count in label_counts.values() if count >= 2)

    summary = {
        "endpoint_dir": str(args.endpoint_dir),
        "checkpoint_ids": str(args.checkpoint_ids),
        "antigen_labels": str(args.antigen_labels),
        "out_dir": str(args.out_dir),
        "input_target_count": len(target_overlap),
        "input_condition_rows": len(condition_manifest),
        "checkpoint_id_count": len(checkpoint_ids),
        "checkpoint_antigen_source_pair_count": len(checkpoint_pairs),
        "kept_target_count": len(kept_ids),
        "dropped_target_count": len(dropped_rows),
        "condition_rows": len(filtered_manifest),
        "missing_label_target_count": len(missing_label_ids),
        "missing_label_ids": missing_label_ids,
        "labels_with_at_least_two_targets": sum(1 for count in label_counts.values() if count >= 2),
        "self_match_excluded_evaluable_target_count": evaluable_target_count,
        "top_retained_labels": label_rows[:12],
        "top_dropped_overlap_pairs": [
            {"overlapping_pair": pair, "target_count": count}
            for pair, count in Counter(
                pair for rows in overlap_by_target.values() for pair in rows
            ).most_common(12)
        ],
        "rule": "Keep endpoint targets whose exact normalized antigen_name || antigen_species pairs do not occur in the checkpoint universe.",
    }
    (args.out_dir / "tierc_endpoint_filter_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
