#!/usr/bin/env python3
"""Build Tier B sequence-deoverlapped conformer target manifest.

Tier A starts from exact-PDB-new generated outputs. Tier B additionally removes
targets whose concatenated VH+VL variable-region sequence has global
Needleman-Wunsch identity >= threshold to any target in the checkpoint universe.

The identity function mirrors Gaeun's split script:
identical aligned columns / max(sequence lengths), using Biopython's global
PairwiseAligner with match=1, mismatch=0, gap open=-1, gap extend=-0.5.
"""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import time
from pathlib import Path
from typing import Any

from Bio.Align import PairwiseAligner


DEFAULT_CANDIDATE_MANIFEST = Path(
    "/project/liulab/jg1920/conffusion/manifests/"
    "gaeun_conformer_ensembles_generated_non_1000_all_20260627.tsv"
)
DEFAULT_SEQUENCE_JSON = Path(
    "/external/liulab/gkim/antigen_prediction/datasets/seq_files/full_sabdab_seq_file.json"
)
DEFAULT_CHECKPOINT_IDS = Path(
    "/external/liulab/gkim/antigen_prediction/2026.06.19_retrain_mca_1000_confs/"
    "split_confs/all_pdb_ids.txt"
)
DEFAULT_OUT_DIR = Path("/project/liulab/jg1920/conffusion/tierb_sequence_deoverlap_20260627")

_TARGETS: list[dict[str, str]] = []
_ALIGNER: PairwiseAligner | None = None


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5
    return aligner


def identity(a: str, b: str) -> float:
    global _ALIGNER
    if not a or not b:
        return 0.0
    if a == b:
        return 1.0
    if min(len(a), len(b)) / max(len(a), len(b)) == 0:
        return 0.0
    if _ALIGNER is None:
        _ALIGNER = make_aligner()
    aln = _ALIGNER.align(a, b)[0]
    matches = 0
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1], strict=False):
        for k in range(a1 - a0):
            if a[a0 + k] == b[b0 + k]:
                matches += 1
    return matches / max(len(a), len(b))


def init_worker(targets: list[dict[str, str]]) -> None:
    global _TARGETS, _ALIGNER
    _TARGETS = targets
    _ALIGNER = make_aligner()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def read_ids(path: Path) -> list[str]:
    return sorted({line.strip().lower() for line in path.read_text().splitlines() if line.strip()})


def load_sequences(path: Path) -> dict[str, dict[str, str]]:
    raw = json.loads(path.read_text())
    if isinstance(raw, list):
        by_id = {str(rec["name"]).lower(): rec for rec in raw if rec.get("name")}
    else:
        by_id = {str(k).lower(): v for k, v in raw.items()}
    out: dict[str, dict[str, str]] = {}
    for target_id, rec in by_id.items():
        vh = str(rec.get("vh_seq", "")).strip().upper()
        vl = str(rec.get("vl_seq", "")).strip().upper()
        if vh and vl:
            out[target_id] = {"vh": vh, "vl": vl, "vhvl": vh + vl}
    return out


def best_identity(query: str, field: str) -> tuple[float, str]:
    best = -1.0
    best_id = ""
    q_len = len(query)
    for target in _TARGETS:
        target_seq = target[field]
        # Upper bound on identity with this denominator. If it cannot beat the
        # current best, skip the alignment.
        if q_len and target_seq and min(q_len, len(target_seq)) / max(q_len, len(target_seq)) < best:
            continue
        val = identity(query, target_seq)
        if val > best:
            best = val
            best_id = target["target_id"]
    return max(best, 0.0), best_id


def score_one(args: tuple[dict[str, str], float, str]) -> dict[str, Any]:
    row, threshold, threshold_label = args
    pass_col = f"tierb_pass_{threshold_label}"
    fail_col = f"tierb_fail_reason_{threshold_label}"
    target_id = row["target_id"].lower()
    seqs = row.get("_seqs")
    if not isinstance(seqs, dict):
        out: dict[str, Any] = dict(row)
        out.update(
            {
                "tierb_has_sequence": False,
                pass_col: False,
                fail_col: "missing_vh_or_vl_sequence",
                "max_vhvl_identity_to_1000_all": "",
                "nearest_vhvl_target_1000_all": "",
                "max_vh_identity_to_1000_all": "",
                "nearest_vh_target_1000_all": "",
                "max_vl_identity_to_1000_all": "",
                "nearest_vl_target_1000_all": "",
            }
        )
        return out

    max_vhvl, nearest_vhvl = best_identity(seqs["vhvl"], "vhvl")
    max_vh, nearest_vh = best_identity(seqs["vh"], "vh")
    max_vl, nearest_vl = best_identity(seqs["vl"], "vl")
    passed = max_vhvl < threshold
    out = {k: v for k, v in row.items() if k != "_seqs"}
    out.update(
        {
            "tierb_has_sequence": True,
            pass_col: passed,
            fail_col: "" if passed else f"max_vhvl_identity_ge_{threshold:.2f}",
            "max_vhvl_identity_to_1000_all": f"{max_vhvl:.6f}",
            "nearest_vhvl_target_1000_all": nearest_vhvl,
            "max_vh_identity_to_1000_all": f"{max_vh:.6f}",
            "nearest_vh_target_1000_all": nearest_vh,
            "max_vl_identity_to_1000_all": f"{max_vl:.6f}",
            "nearest_vl_target_1000_all": nearest_vl,
        }
    )
    return out


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-manifest", type=Path, default=DEFAULT_CANDIDATE_MANIFEST)
    parser.add_argument("--sequence-json", type=Path, default=DEFAULT_SEQUENCE_JSON)
    parser.add_argument("--checkpoint-ids", type=Path, default=DEFAULT_CHECKPOINT_IDS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--threshold", type=float, default=0.80)
    parser.add_argument("--workers", type=int, default=max(1, min(32, mp.cpu_count() or 1)))
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    threshold_label = f"vhvl{int(round(args.threshold * 100)):02d}"
    pass_col = f"tierb_pass_{threshold_label}"
    fail_col = f"tierb_fail_reason_{threshold_label}"

    base_columns, candidate_rows = read_tsv(args.candidate_manifest)
    seq_by_id = load_sequences(args.sequence_json)
    checkpoint_ids = read_ids(args.checkpoint_ids)
    checkpoint_targets = []
    missing_checkpoint = []
    for target_id in checkpoint_ids:
        seqs = seq_by_id.get(target_id)
        if seqs:
            checkpoint_targets.append({"target_id": target_id, **seqs})
        else:
            missing_checkpoint.append(target_id)

    scoring_rows = []
    missing_candidates = []
    for row in candidate_rows:
        target_id = row["target_id"].lower()
        seqs = seq_by_id.get(target_id)
        new_row = dict(row)
        if seqs:
            new_row["_seqs"] = seqs  # type: ignore[assignment]
        else:
            missing_candidates.append(target_id)
        scoring_rows.append(new_row)

    print(f"[tierB] candidates={len(candidate_rows)}", flush=True)
    print(f"[tierB] checkpoint ids={len(checkpoint_ids)} usable={len(checkpoint_targets)}", flush=True)
    print(f"[tierB] missing candidate sequences={len(missing_candidates)}", flush=True)
    print(f"[tierB] threshold={args.threshold:.2f} label={threshold_label}", flush=True)
    print(f"[tierB] workers={args.workers}", flush=True)
    t0 = time.time()

    with mp.Pool(processes=args.workers, initializer=init_worker, initargs=(checkpoint_targets,)) as pool:
        scored = list(
            pool.imap_unordered(
                score_one,
                [(row, args.threshold, threshold_label) for row in scoring_rows],
                chunksize=1,
            )
        )
    scored.sort(key=lambda row: str(row["target_id"]))
    passed = [row for row in scored if str(row.get(pass_col)).lower() == "true"]
    failed = [row for row in scored if str(row.get(pass_col)).lower() != "true"]

    extra_columns = [
        "tierb_has_sequence",
        pass_col,
        fail_col,
        "max_vhvl_identity_to_1000_all",
        "nearest_vhvl_target_1000_all",
        "max_vh_identity_to_1000_all",
        "nearest_vh_target_1000_all",
        "max_vl_identity_to_1000_all",
        "nearest_vl_target_1000_all",
    ]
    columns = base_columns + [col for col in extra_columns if col not in base_columns]
    stem = f"gaeun_conformer_ensembles_generated_non_1000_all_tierB_{threshold_label}"
    scored_path = args.out_dir / f"{stem}_scored_20260627.tsv"
    pass_path = args.out_dir / f"{stem}_pass_20260627.tsv"
    fail_path = args.out_dir / f"{stem}_fail_20260627.tsv"
    write_tsv(scored_path, scored, columns)
    write_tsv(pass_path, passed, columns)
    write_tsv(fail_path, failed, columns)

    def max_float(rows: list[dict[str, Any]], key: str) -> float | None:
        vals = [float(row[key]) for row in rows if str(row.get(key, "")).strip()]
        return max(vals) if vals else None

    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime()),
        "candidate_manifest": str(args.candidate_manifest),
        "sequence_json": str(args.sequence_json),
        "checkpoint_ids": str(args.checkpoint_ids),
        "identity_rule": "global NW identical_aligned_columns / max(len), concatenated VH+VL",
        "threshold": args.threshold,
        "threshold_label": threshold_label,
        "candidates": len(candidate_rows),
        "checkpoint_ids_count": len(checkpoint_ids),
        "checkpoint_ids_with_sequences": len(checkpoint_targets),
        "missing_checkpoint_sequences": len(missing_checkpoint),
        "missing_candidate_sequences": len(missing_candidates),
        "tierb_pass_count": len(passed),
        "tierb_fail_count": len(failed),
        "max_pass_vhvl_identity": max_float(passed, "max_vhvl_identity_to_1000_all"),
        "min_fail_vhvl_identity": min(
            [float(row["max_vhvl_identity_to_1000_all"]) for row in failed if str(row.get("max_vhvl_identity_to_1000_all", "")).strip()],
            default=None,
        ),
        "max_pass_vh_identity_annotation": max_float(passed, "max_vh_identity_to_1000_all"),
        "max_pass_vl_identity_annotation": max_float(passed, "max_vl_identity_to_1000_all"),
        "elapsed_seconds": round(time.time() - t0, 2),
        "outputs": {
            "scored": str(scored_path),
            "pass": str(pass_path),
            "fail": str(fail_path),
        },
    }
    summary_json = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    (args.out_dir / f"tierb_sequence_deoverlap_{threshold_label}_summary.json").write_text(summary_json)
    (args.out_dir / "tierb_sequence_deoverlap_summary.json").write_text(summary_json)
    md = [
        f"# Tier B Sequence De-Overlap ({threshold_label})",
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
        f"- Tier A candidates scored: {len(candidate_rows)}",
        f"- Checkpoint IDs compared against: {len(checkpoint_targets)}",
        f"- Tier B pass: {len(passed)}",
        f"- Tier B fail: {len(failed)}",
        f"- Missing candidate sequences: {len(missing_candidates)}",
        "",
        "## Outputs",
        "",
        f"- Scored all candidates: `{scored_path.name}`",
        f"- Tier B pass manifest: `{pass_path.name}`",
        f"- Tier B fail manifest: `{fail_path.name}`",
        f"- Machine-readable summary: `tierb_sequence_deoverlap_{threshold_label}_summary.json`",
        "",
        "## Caveat",
        "",
        "Tier B is antibody-sequence de-overlapped against the checkpoint universe.",
        "It does not yet enforce antigen-name, antigen-sequence, or same-antigen",
        "group de-overlap.",
        "",
    ]
    summary_md = "\n".join(md)
    (args.out_dir / f"tierb_sequence_deoverlap_{threshold_label}_summary.md").write_text(summary_md)
    (args.out_dir / "tierb_sequence_deoverlap_summary.md").write_text(summary_md)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
