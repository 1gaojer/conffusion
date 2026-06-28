#!/usr/bin/env python3
"""Prepare a Tier B Phase 1.5 endpoint manifest.

This creates the minimal Phase 1.5-style files needed by the readout scripts:
target_overlap.tsv and condition_manifest.tsv. It uses Jerry-owned Tier B
manifests and saved CDR structural selections.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


DEFAULT_TARGET_IDS = Path(
    "/project/liulab/jg1920/conffusion/phase14_tierb_vhvl85_full100_20260627_215733/run/manifests/target_ids.txt"
)
DEFAULT_ANTIGEN_LABELS = Path(
    "/project/liulab/jg1920/conffusion/tierb_vhvl85_dataset_20260627/manifests/tierb_antigen_labels.tsv"
)
DEFAULT_CDR_CURVES = Path(
    "/project/liulab/jg1920/conffusion/tierb_vhvl85_cdr_structural_20260627/cdr_saturation_curves.tsv"
)
DEFAULT_OUT = Path("/project/liulab/jg1920/conffusion/phase15_tierb_step4_smoke")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def parse_indices(text: str) -> list[int]:
    text = str(text or "").strip()
    if not text:
        return []
    return [int(part) for part in text.split(",") if part != ""]


def selected_indices_from_curves(path: Path) -> dict[tuple[str, str, str, int, int], list[int]]:
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


def normalize_positive_key(row: dict[str, str]) -> str:
    value = row.get("positive_key", "")
    if value:
        return value.strip().lower()
    name = row.get("antigen_name", row.get("reference_antigen_name", "")).strip().lower()
    species = row.get("antigen_species", row.get("reference_antigen_species", "")).strip().lower()
    if not name and not species:
        return ""
    return f"{name}||{species}"


def choose_smoke_targets(target_ids: list[str], labels: dict[str, dict[str, str]], n: int) -> list[str]:
    if n <= 0 or n >= len(target_ids):
        return target_ids
    by_label: dict[str, list[str]] = defaultdict(list)
    for target_id in target_ids:
        key = normalize_positive_key(labels.get(target_id, {}))
        if key:
            by_label[key].append(target_id)
    grouped = sorted(
        (targets for targets in by_label.values() if len(targets) >= 2),
        key=lambda vals: (-len(vals), vals[0]),
    )
    selected: list[str] = []
    for group in grouped:
        for target_id in group:
            if target_id not in selected:
                selected.append(target_id)
            if len(selected) >= n:
                return selected
    for target_id in target_ids:
        if target_id not in selected:
            selected.append(target_id)
        if len(selected) >= n:
            return selected
    return selected


def condition_rows(
    target_id: str,
    n_conformers: int,
    cdr_index: dict[tuple[str, str, str, int, int], list[int]],
    random_replicates: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = [
        {
            "target_id": target_id,
            "condition": f"full_{n_conformers}",
            "family": "full",
            "region": "all",
            "strategy": "all",
            "k": n_conformers,
            "replicate": 0,
            "selected_n": n_conformers,
            "selected_indices": ",".join(str(i) for i in range(n_conformers)),
        },
        {
            "target_id": target_id,
            "condition": "single_first",
            "family": "single",
            "region": "none",
            "strategy": "first",
            "k": 1,
            "replicate": 0,
            "selected_n": 1,
            "selected_indices": "0",
        },
    ]
    for region in ("H-CDR3", "all-CDRs"):
        for k in (32, 64):
            selected = cdr_index.get((target_id, region, "greedy_kcenter", k, 0), [])
            if selected:
                rows.append(
                    {
                        "target_id": target_id,
                        "condition": f"{region.lower().replace('-', '_')}_greedy_kcenter_k{k}",
                        "family": "cdr_kcenter",
                        "region": region,
                        "strategy": "greedy_kcenter",
                        "k": k,
                        "replicate": 0,
                        "selected_n": len(selected),
                        "selected_indices": ",".join(str(i) for i in selected),
                    }
                )
    for k in (32, 64):
        selected = cdr_index.get((target_id, "H-CDR3", "first", k, 0), [])
        if selected:
            rows.append(
                {
                    "target_id": target_id,
                    "condition": f"first_k{k}",
                    "family": "first",
                    "region": "none",
                    "strategy": "first",
                    "k": k,
                    "replicate": 0,
                    "selected_n": len(selected),
                    "selected_indices": ",".join(str(i) for i in selected),
                }
            )
        for rep in range(random_replicates):
            selected = cdr_index.get((target_id, "H-CDR3", "random", k, rep), [])
            if selected:
                rows.append(
                    {
                        "target_id": target_id,
                        "condition": f"random_k{k}_rep{rep}",
                        "family": "random",
                        "region": "none",
                        "strategy": "random",
                        "k": k,
                        "replicate": rep,
                        "selected_n": len(selected),
                        "selected_indices": ",".join(str(i) for i in selected),
                    }
                )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-ids", type=Path, default=DEFAULT_TARGET_IDS)
    parser.add_argument("--antigen-labels", type=Path, default=DEFAULT_ANTIGEN_LABELS)
    parser.add_argument("--cdr-curves", type=Path, default=DEFAULT_CDR_CURVES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--n-conformers", type=int, default=100)
    parser.add_argument("--smoke-target-count", type=int, default=0)
    parser.add_argument("--random-replicates", type=int, default=5)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    target_ids = [line.strip() for line in args.target_ids.read_text().splitlines() if line.strip()]
    labels = {row["target_id"]: row for row in read_tsv(args.antigen_labels)}
    target_ids = [target_id for target_id in target_ids if target_id in labels]
    target_ids = choose_smoke_targets(target_ids, labels, args.smoke_target_count)
    cdr_index = selected_indices_from_curves(args.cdr_curves)

    overlap_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    for target_id in target_ids:
        label = labels[target_id]
        positive_key = normalize_positive_key(label)
        name, sep, species = positive_key.partition("||")
        overlap_rows.append(
            {
                "target_id": target_id,
                "positive_key": positive_key,
                "reference_antigen_name": name,
                "reference_antigen_species": species if sep else "",
                "in_reference_bank": "",
                "same_label_reference_count": "",
                "same_label_after_self_exclusion": "",
                "in_june24_query_index": "",
                "in_strict300_query_id": "",
                "in_strict300_pdb_id": "",
            }
        )
        manifest_rows.extend(condition_rows(target_id, args.n_conformers, cdr_index, args.random_replicates))

    write_tsv(
        args.out_dir / "target_overlap.tsv",
        overlap_rows,
        [
            "target_id",
            "positive_key",
            "reference_antigen_name",
            "reference_antigen_species",
            "in_reference_bank",
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
    label_counts: dict[str, int] = defaultdict(int)
    for row in overlap_rows:
        label_counts[str(row["positive_key"])] += 1
    summary = {
        "out_dir": str(args.out_dir),
        "target_count": len(target_ids),
        "condition_rows": len(manifest_rows),
        "n_conformers": args.n_conformers,
        "random_replicates": args.random_replicates,
        "labels_with_at_least_two_targets": sum(1 for count in label_counts.values() if count >= 2),
        "top_labels": sorted(label_counts.items(), key=lambda item: (-item[1], item[0]))[:12],
        "inputs": {
            "target_ids": str(args.target_ids),
            "antigen_labels": str(args.antigen_labels),
            "cdr_curves": str(args.cdr_curves),
        },
    }
    (args.out_dir / "tierb_endpoint_prep_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
