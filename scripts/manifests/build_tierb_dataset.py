#!/usr/bin/env python3
"""Materialize a Jerry-owned Tier B conformer dataset manifest.

This script reads Gaeun-owned conformer outputs and labels as immutable inputs.
It writes a lightweight dataset root under Jerry-owned space. By default it
creates symlinks to source CIFs rather than copying large CIF files.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


DEFAULT_CONFFUSION_ROOT = Path("/project/liulab/jg1920/conffusion")
DEFAULT_TIERB_PASS = (
    DEFAULT_CONFFUSION_ROOT
    / "tierb_sequence_deoverlap_20260627"
    / "gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_pass_20260627.tsv"
)
DEFAULT_SEQ_JSON = Path(
    "/external/liulab/gkim/antigen_prediction/datasets/seq_files/full_sabdab_seq_file.json"
)
DEFAULT_LABELS = Path(
    "/external/liulab/gkim/antigen_prediction/datasets/antigen_labels/full_sabdab_relabeled_antigens.csv"
)
DEFAULT_CONF_GEN_ROOT = Path(
    "/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results"
)
DEFAULT_OUT_ROOT = DEFAULT_CONFFUSION_ROOT / "tierb_vhvl85_dataset_20260627"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def norm_text(value: object) -> str:
    return " ".join(str(value or "").strip().lower().split())


def load_sequence_meta(path: Path) -> dict[str, dict[str, object]]:
    data = json.loads(path.read_text())
    if not isinstance(data, list):
        raise TypeError(f"Expected list in sequence JSON: {path}")
    out = {}
    for row in data:
        name = norm_text(row.get("name"))
        if name:
            out[name] = row
    return out


def load_antigen_labels(path: Path) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pdb = norm_text(row.get("pdb"))
            if pdb:
                grouped[pdb].append(row)
    return grouped


def label_summary(target_id: str, labels_by_pdb: dict[str, list[dict[str, str]]]) -> dict[str, object]:
    rows = labels_by_pdb.get(norm_text(target_id), [])
    if not rows:
        return {
            "target_id": target_id,
            "has_antigen_label": False,
            "label_row_count": 0,
            "positive_key": "",
            "positive_key_count": 0,
        }

    keys = []
    names = []
    species = []
    antigen_types = []
    chains = []
    short_headers = []
    compounds = []
    organisms = []
    for row in rows:
        name = norm_text(row.get("antigen_name"))
        specie = norm_text(row.get("antigen_species"))
        key = f"{name} || {specie}".strip()
        keys.append(key)
        names.append(str(row.get("antigen_name") or "").strip())
        species.append(str(row.get("antigen_species") or "").strip())
        antigen_types.append(str(row.get("antigen_type") or "").strip())
        chains.append(str(row.get("antigen_chain") or "").strip())
        short_headers.append(str(row.get("short_header") or "").strip())
        compounds.append(str(row.get("compound") or "").strip())
        organisms.append(str(row.get("organism") or "").strip())

    uniq_keys = sorted({key for key in keys if key and key != "||"})
    return {
        "target_id": target_id,
        "has_antigen_label": bool(uniq_keys),
        "label_row_count": len(rows),
        "positive_key": " + ".join(uniq_keys),
        "positive_key_count": len(uniq_keys),
        "antigen_names": " | ".join(sorted({x for x in names if x})),
        "antigen_species": " | ".join(sorted({x for x in species if x})),
        "antigen_types": " | ".join(sorted({x for x in antigen_types if x})),
        "antigen_chains": " | ".join(sorted({x for x in chains if x})),
        "short_headers": " | ".join(sorted({x for x in short_headers if x})),
        "compounds": " | ".join(sorted({x for x in compounds if x})),
        "organisms": " | ".join(sorted({x for x in organisms if x})),
    }


def parse_run_cycle(target_id: str, path: Path) -> tuple[int, int]:
    pattern = re.compile(rf"^{re.escape(target_id)}_run(?P<run>\d+)_cycle(?P<cycle>\d+)$")
    match = pattern.match(path.name)
    if not match:
        return (10**9, 10**9)
    return int(match.group("run")), int(match.group("cycle"))


def choose_model_cif(result_dir: Path) -> tuple[Path | None, str]:
    direct = result_dir / f"{result_dir.name}_model.cif"
    if direct.exists():
        return direct, "af3_result_root_model"
    nested = sorted(result_dir.glob("**/*_model.cif"))
    if nested:
        return nested[0], "af3_result_nested_model"
    return None, "missing_model_cif"


def target_conformer_sources(
    *,
    target_id: str,
    conf_gen_root: Path,
    max_conformers: int,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    af3_root = conf_gen_root / "af3_results"
    result_dirs = [
        path
        for path in af3_root.glob(f"{target_id}_run*_cycle*")
        if path.is_dir()
    ]
    result_dirs.sort(key=lambda path: (*parse_run_cycle(target_id, path), path.name))

    rows = []
    missing_model_dirs = 0
    for result_dir in result_dirs:
        model_path, model_kind = choose_model_cif(result_dir)
        if model_path is None:
            missing_model_dirs += 1
            continue
        run_id, cycle_id = parse_run_cycle(target_id, result_dir)
        rows.append(
            {
                "target_id": target_id,
                "run_id": run_id if run_id < 10**9 else "",
                "cycle_id": cycle_id if cycle_id < 10**9 else "",
                "source_af3_result_dir": str(result_dir),
                "source_af3_model_path": str(model_path),
                "source_model_kind": model_kind,
                "source_size_bytes": model_path.stat().st_size,
                "source_mtime_epoch": round(model_path.stat().st_mtime, 3),
            }
        )
        if len(rows) >= max_conformers:
            break

    coverage = {
        "target_id": target_id,
        "af3_result_dir_count": len(result_dirs),
        "model_cif_count": len(rows),
        "missing_model_dirs": missing_model_dirs,
        "selected_conformer_count": len(rows),
        "selected_max_conformers": max_conformers,
    }
    return rows, coverage


def link_or_copy(src: Path, dst: Path, mode: str, force: bool) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        if not force:
            if dst.is_symlink() and os.readlink(dst) == str(src):
                return
            raise FileExistsError(f"Destination already exists: {dst}")
        if dst.is_dir() and not dst.is_symlink():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    if mode == "symlink":
        dst.symlink_to(src)
    elif mode == "copy":
        shutil.copy2(src, dst)
    elif mode == "none":
        return
    else:
        raise ValueError(f"Unknown materialize mode: {mode}")


def write_report(out_root: Path, summary: dict[str, object]) -> None:
    lines = [
        "# Tier B VH/VL85 Dataset Preflight",
        "",
        f"Created: {summary['created_at_utc']}",
        "",
        "## Inputs",
        "",
        f"- Tier B pass manifest: `{summary['tierb_pass_manifest']}`",
        f"- Sequence JSON: `{summary['sequence_json']}`",
        f"- Antigen labels: `{summary['antigen_labels_csv']}`",
        f"- Gaeun conformer root: `{summary['conf_gen_root']}`",
        "",
        "## Outputs",
        "",
        f"- Dataset root: `{summary['out_root']}`",
        "- `manifests/tierb_targets.tsv`",
        "- `manifests/tierb_conformers.tsv`",
        "- `manifests/tierb_antigen_labels.tsv`",
        "- `manifests/tierb_conformer_coverage.tsv`",
        "- `manifests/source_summary.json`",
        "",
        "## Coverage",
        "",
        f"- Tier B targets requested: {summary['tierb_targets_requested']}",
        f"- Targets with sequence metadata: {summary['targets_with_sequence']}",
        f"- Targets with antigen labels: {summary['targets_with_antigen_label']}",
        f"- Targets with at least one selected conformer: {summary['targets_with_conformers']}",
        f"- Total selected conformers: {summary['selected_conformer_count']}",
        f"- Materialization mode: `{summary['materialize_mode']}`",
        "",
        "## Conformer Count Distribution",
        "",
        "| Selected conformers | Targets |",
        "|---:|---:|",
    ]
    for n, count in summary["selected_conformer_count_distribution"]:
        lines.append(f"| {n} | {count} |")
    lines.extend(
        [
            "",
            "Read-only guarantee: source paths under Gaeun-owned directories were only read; all generated manifests and links/copies live under the Jerry-owned dataset root.",
            "",
        ]
    )
    (out_root / "tierb_dataset_preflight_summary.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tierb-pass-manifest", type=Path, default=DEFAULT_TIERB_PASS)
    parser.add_argument("--sequence-json", type=Path, default=DEFAULT_SEQ_JSON)
    parser.add_argument("--antigen-labels-csv", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--conf-gen-root", type=Path, default=DEFAULT_CONF_GEN_ROOT)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--max-conformers", type=int, default=100)
    parser.add_argument("--materialize-mode", choices=["symlink", "copy", "none"], default="symlink")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    out_root = args.out_root.resolve()
    manifests = out_root / "manifests"
    manifests.mkdir(parents=True, exist_ok=True)

    tierb_rows = read_tsv(args.tierb_pass_manifest)
    seq_meta = load_sequence_meta(args.sequence_json)
    labels = load_antigen_labels(args.antigen_labels_csv)

    target_rows = []
    label_rows = []
    coverage_rows = []
    conformer_rows = []
    missing_sequence = []
    selected_count_distribution = Counter()

    for tier_row in tierb_rows:
        target_id = str(tier_row["target_id"]).lower()
        seq = seq_meta.get(target_id)
        if not seq:
            missing_sequence.append(target_id)
            continue
        label = label_summary(target_id, labels)
        sources, coverage = target_conformer_sources(
            target_id=target_id,
            conf_gen_root=args.conf_gen_root,
            max_conformers=args.max_conformers,
        )
        label_rows.append(label)
        coverage.update(
            {
                "latest_stage_observed": tier_row.get("latest_stage_observed", ""),
                "has_filtered_af3_models": tier_row.get("has_filtered_af3_models", ""),
                "has_af3_inputs": tier_row.get("has_af3_inputs", ""),
                "has_prothunt_output": tier_row.get("has_prothunt_output", ""),
            }
        )
        coverage_rows.append(coverage)
        selected_count_distribution[int(coverage["selected_conformer_count"])] += 1

        target_rows.append(
            {
                **tier_row,
                "target_id": target_id,
                "name": target_id,
                "vh_seq": seq.get("vh_seq", ""),
                "vl_seq": seq.get("vl_seq", ""),
                "antigen_size": seq.get("antigen_size", ""),
                "has_antigen_label": label.get("has_antigen_label", ""),
                "positive_key": label.get("positive_key", ""),
                "selected_conformer_count": coverage["selected_conformer_count"],
                "tierb_dataset_eligible": bool(coverage["selected_conformer_count"]),
            }
        )

        for pick_rank, source in enumerate(sources):
            src = Path(str(source["source_af3_model_path"]))
            relative = Path("conformers") / target_id / f"{pick_rank:03d}_{src.name}"
            dst = out_root / relative
            link_or_copy(src, dst, args.materialize_mode, args.force)
            conformer_rows.append(
                {
                    **source,
                    "pick_rank": pick_rank,
                    "copied_relative_path": str(relative),
                    "materialize_mode": args.materialize_mode,
                }
            )

    target_fields = list(tierb_rows[0].keys()) + [
        "name",
        "vh_seq",
        "vl_seq",
        "antigen_size",
        "has_antigen_label",
        "positive_key",
        "selected_conformer_count",
        "tierb_dataset_eligible",
    ]
    # Preserve order while avoiding duplicate target_id from expansion.
    target_fields = list(dict.fromkeys(["target_id", *target_fields]))
    write_tsv(manifests / "tierb_targets.tsv", target_rows, target_fields)
    write_tsv(
        manifests / "tierb_antigen_labels.tsv",
        label_rows,
        [
            "target_id",
            "has_antigen_label",
            "label_row_count",
            "positive_key",
            "positive_key_count",
            "antigen_names",
            "antigen_species",
            "antigen_types",
            "antigen_chains",
            "short_headers",
            "compounds",
            "organisms",
        ],
    )
    write_tsv(
        manifests / "tierb_conformer_coverage.tsv",
        coverage_rows,
        [
            "target_id",
            "latest_stage_observed",
            "has_filtered_af3_models",
            "has_af3_inputs",
            "has_prothunt_output",
            "af3_result_dir_count",
            "model_cif_count",
            "missing_model_dirs",
            "selected_conformer_count",
            "selected_max_conformers",
        ],
    )
    write_tsv(
        manifests / "tierb_conformers.tsv",
        conformer_rows,
        [
            "target_id",
            "pick_rank",
            "run_id",
            "cycle_id",
            "copied_relative_path",
            "source_af3_model_path",
            "source_af3_result_dir",
            "source_model_kind",
            "source_size_bytes",
            "source_mtime_epoch",
            "materialize_mode",
        ],
    )

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tierb_pass_manifest": str(args.tierb_pass_manifest),
        "sequence_json": str(args.sequence_json),
        "antigen_labels_csv": str(args.antigen_labels_csv),
        "conf_gen_root": str(args.conf_gen_root),
        "out_root": str(out_root),
        "max_conformers": args.max_conformers,
        "materialize_mode": args.materialize_mode,
        "source_writes_performed": False,
        "tierb_targets_requested": len(tierb_rows),
        "targets_with_sequence": len(target_rows),
        "missing_sequence_targets": missing_sequence,
        "targets_with_antigen_label": sum(str(row.get("has_antigen_label")).lower() == "true" for row in label_rows),
        "targets_with_conformers": sum(int(row["selected_conformer_count"]) > 0 for row in coverage_rows),
        "selected_conformer_count": len(conformer_rows),
        "selected_conformer_count_distribution": sorted(selected_count_distribution.items()),
        "outputs": {
            "targets": str(manifests / "tierb_targets.tsv"),
            "conformers": str(manifests / "tierb_conformers.tsv"),
            "antigen_labels": str(manifests / "tierb_antigen_labels.tsv"),
            "coverage": str(manifests / "tierb_conformer_coverage.tsv"),
        },
    }
    (manifests / "source_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_report(out_root, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
