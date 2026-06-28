#!/usr/bin/env python3
"""Materialize Tier B pass/fail manifests from a scored identity table."""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import time
from pathlib import Path
from typing import Any


DEFAULT_SCORED = Path(
    "manifests/tierb_sequence_deoverlap_20260627/"
    "gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl80_scored_20260627.tsv"
)
DEFAULT_OUT_DIR = Path("manifests/tierb_sequence_deoverlap_20260627")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def threshold_label(threshold: float) -> str:
    return f"vhvl{int(round(threshold * 100)):02d}"


def strip_old_threshold_columns(columns: list[str]) -> list[str]:
    old_pass = re.compile(r"^tierb_pass_vhvl\d+$")
    old_fail = re.compile(r"^tierb_fail_reason(?:_vhvl\d+)?$")
    return [col for col in columns if not old_pass.match(col) and not old_fail.match(col)]


def quantiles(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"min": None, "p25": None, "median": None, "p75": None, "max": None}
    if len(values) == 1:
        val = values[0]
        return {"min": val, "p25": val, "median": val, "p75": val, "max": val}
    qs = statistics.quantiles(values, n=4)
    return {
        "min": min(values),
        "p25": qs[0],
        "median": statistics.median(values),
        "p75": qs[2],
        "max": max(values),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scored", type=Path, default=DEFAULT_SCORED)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--threshold", type=float, required=True)
    parser.add_argument("--date", default="20260627")
    args = parser.parse_args()

    label = threshold_label(args.threshold)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    base_columns, rows = read_tsv(args.scored)
    columns = strip_old_threshold_columns(base_columns)
    pass_col = f"tierb_pass_{label}"
    fail_col = f"tierb_fail_reason_{label}"
    for col in [pass_col, fail_col]:
        if col not in columns:
            columns.append(col)

    scored_rows: list[dict[str, Any]] = []
    missing_identity = 0
    for row in rows:
        out = {col: row.get(col, "") for col in columns if col not in {pass_col, fail_col}}
        raw_identity = row.get("max_vhvl_identity_to_1000_all", "").strip()
        if not raw_identity:
            missing_identity += 1
            passed = False
        else:
            passed = float(raw_identity) < args.threshold
        out[pass_col] = passed
        out[fail_col] = "" if passed else f"max_vhvl_identity_ge_{args.threshold:.2f}"
        scored_rows.append(out)

    passed_rows = [row for row in scored_rows if str(row[pass_col]).lower() == "true"]
    failed_rows = [row for row in scored_rows if str(row[pass_col]).lower() != "true"]

    stem = f"gaeun_conformer_ensembles_generated_non_1000_all_tierB_{label}"
    scored_path = args.out_dir / f"{stem}_scored_{args.date}.tsv"
    pass_path = args.out_dir / f"{stem}_pass_{args.date}.tsv"
    fail_path = args.out_dir / f"{stem}_fail_{args.date}.tsv"
    write_tsv(scored_path, scored_rows, columns)
    write_tsv(pass_path, passed_rows, columns)
    write_tsv(fail_path, failed_rows, columns)

    pass_vals = [float(row["max_vhvl_identity_to_1000_all"]) for row in passed_rows]
    fail_vals = [float(row["max_vhvl_identity_to_1000_all"]) for row in failed_rows]
    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime()),
        "source_scored_manifest": str(args.scored),
        "threshold": args.threshold,
        "threshold_label": label,
        "identity_rule": "global NW identical_aligned_columns / max(len), concatenated VH+VL",
        "candidates": len(rows),
        "tierb_pass_count": len(passed_rows),
        "tierb_fail_count": len(failed_rows),
        "missing_identity_values": missing_identity,
        "pass_identity_distribution": quantiles(pass_vals),
        "fail_identity_distribution": quantiles(fail_vals),
        "fail_count_ge_0_90": sum(val >= 0.90 for val in fail_vals),
        "fail_count_ge_0_95": sum(val >= 0.95 for val in fail_vals),
        "fail_count_eq_1_00": sum(abs(val - 1.0) < 1e-9 for val in fail_vals),
        "outputs": {
            "scored": str(scored_path),
            "pass": str(pass_path),
            "fail": str(fail_path),
        },
    }

    label_summary_json = args.out_dir / f"tierb_sequence_deoverlap_{label}_summary.json"
    label_summary_md = args.out_dir / f"tierb_sequence_deoverlap_{label}_summary.md"
    current_summary_json = args.out_dir / "tierb_sequence_deoverlap_summary.json"
    current_summary_md = args.out_dir / "tierb_sequence_deoverlap_summary.md"

    for path in [label_summary_json, current_summary_json]:
        path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    md = [
        f"# Tier B Sequence De-Overlap ({label})",
        "",
        f"Date: {time.strftime('%Y-%m-%d')}",
        "",
        "## Definition",
        "",
        "Tier B starts from the 359 generated-output targets that do not overlap",
        "the full June 19/21 `1000` checkpoint PDB universe. It keeps only",
        "targets whose concatenated VH+VL variable-region global identity is",
        f"below {args.threshold:.2f} to every checkpoint-universe target.",
        "",
        "Identity mirrors Gaeun's split script: global Needleman-Wunsch,",
        "`identical_aligned_columns / max(len_a, len_b)`, using concatenated",
        "VH+VL sequences.",
        "",
        "## Counts",
        "",
        f"- Tier A candidates scored: {len(rows)}",
        f"- Tier B pass: {len(passed_rows)}",
        f"- Tier B fail: {len(failed_rows)}",
        f"- Missing identity values: {missing_identity}",
        "",
        "## Outputs",
        "",
        f"- Scored all candidates: `{scored_path.name}`",
        f"- Tier B pass manifest: `{pass_path.name}`",
        f"- Tier B fail manifest: `{fail_path.name}`",
        f"- Machine-readable summary: `{label_summary_json.name}`",
        "",
        "## Caveat",
        "",
        "Tier B is antibody-sequence de-overlapped against the checkpoint universe.",
        "It does not yet enforce antigen-name, antigen-sequence, or same-antigen",
        "group de-overlap.",
        "",
    ]
    md_text = "\n".join(md)
    for path in [label_summary_md, current_summary_md]:
        path.write_text(md_text)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
