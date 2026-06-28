#!/usr/bin/env python3
"""Build exact selected-conformer MCA cache and export Phase 1.4 embeddings.

This script is intended to run on Ragon under Jerry-owned /project space.
It reads staged copies of the Phase 1.4 medium conformer dataset and Gaeun's
MCA code/checkpoint, then writes only into a run directory supplied by the
caller.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable


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


def copy_code_tree(src: Path, dst: Path) -> None:
    if not dst.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(
            src,
            dst,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".git", "wandb", "logs"),
        )

    for py_path in dst.glob("*.py"):
        text = py_path.read_text()
        text = text.replace("from turtle import position\n", "")
        text = text.replace("sys.path.append('/project/liulab/gkim/code/scripts')\n", "")
        py_path.write_text(text)


def select_smoke_targets(cdr_summary: Path, medium_targets: list[str]) -> list[str]:
    allowed = set(medium_targets)
    rows = []
    for row in read_tsv(cdr_summary):
        if (
            row.get("target_id") in allowed
            and row.get("metric") == "frame_aligned_ca"
            and row.get("region") == "H-CDR3"
        ):
            try:
                rows.append((float(row["pairwise_rmsd_mean"]), row["target_id"]))
            except (KeyError, ValueError):
                continue
    if len(rows) < 3:
        raise RuntimeError(f"Need at least 3 H-CDR3 rows for smoke selection; got {len(rows)}")
    rows.sort()
    picks = [rows[0][1], rows[len(rows) // 2][1], rows[-1][1]]
    out: list[str] = []
    for target in picks:
        if target not in out:
            out.append(target)
    for _, target in rows:
        if len(out) >= 3:
            break
        if target not in out:
            out.append(target)
    return out[:3]


def load_target_ids(dataset_root: Path, cdr_summary: Path, scope: str) -> list[str]:
    target_rows = read_tsv(dataset_root / "manifests" / "medium_targets.tsv")
    medium_targets = [row["target_id"] for row in target_rows]
    if scope == "smoke":
        return select_smoke_targets(cdr_summary, medium_targets)
    if scope == "full":
        return medium_targets
    raise ValueError(f"Unknown scope: {scope}")


def load_target_meta(dataset_root: Path) -> dict[str, dict[str, str]]:
    return {
        row["target_id"]: row
        for row in read_tsv(dataset_root / "manifests" / "medium_targets.tsv")
    }


def grouped_conformer_rows(dataset_root: Path) -> dict[str, list[dict[str, str]]]:
    rows = read_tsv(dataset_root / "manifests" / "medium_conformers.tsv")
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["target_id"], []).append(row)
    for target_rows in grouped.values():
        target_rows.sort(key=lambda r: int(r["pick_rank"]))
    return grouped


THREE_TO_ONE = {
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
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}


def structure_chain_sequences(path: Path) -> dict[str, str]:
    from Bio.PDB import MMCIFParser, PDBParser  # noqa: PLC0415

    parser = MMCIFParser(QUIET=True) if path.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure(path.name, str(path))
    model = list(structure.get_models())[0]
    sequences: dict[str, str] = {}
    for chain in model.get_chains():
        residues = []
        for residue in chain.get_residues():
            if residue.id[0] != " ":
                continue
            if "CA" not in residue:
                continue
            residues.append(THREE_TO_ONE.get(residue.resname, "X"))
        if residues:
            sequences[chain.id] = "".join(residues)
    return sequences


def sequence_match_score(query: str, candidate: str) -> float:
    if not query or not candidate:
        return 0.0
    if query == candidate:
        return 1.0
    if query in candidate or candidate in query:
        return min(len(query), len(candidate)) / max(len(query), len(candidate))

    from Bio import pairwise2  # noqa: PLC0415

    score = pairwise2.align.globalxx(query, candidate, score_only=True)
    return float(score) / float(max(len(query), len(candidate)))


def infer_hl_chain_ids(
    *,
    target: str,
    first_file: Path,
    target_meta: dict[str, str],
) -> tuple[str, str, float, float, dict[str, int]]:
    chain_sequences = structure_chain_sequences(first_file)
    if not chain_sequences:
        raise RuntimeError(f"{target}: no protein chains found in {first_file}")

    vh_seq = target_meta.get("vh_seq", "")
    vl_seq = target_meta.get("vl_seq", "")
    h_scores = sorted(
        ((sequence_match_score(vh_seq, seq), cid) for cid, seq in chain_sequences.items()),
        reverse=True,
    )
    l_scores = sorted(
        ((sequence_match_score(vl_seq, seq), cid) for cid, seq in chain_sequences.items()),
        reverse=True,
    )
    h_score, h_chain = h_scores[0]
    l_score, l_chain = l_scores[0]
    if l_chain == h_chain:
        for candidate_score, candidate_chain in l_scores[1:]:
            if candidate_chain != h_chain:
                l_score, l_chain = candidate_score, candidate_chain
                break
    if h_chain == l_chain:
        raise RuntimeError(f"{target}: could not infer distinct H/L chains from {first_file}")
    if h_score < 0.75 or l_score < 0.75:
        raise RuntimeError(
            f"{target}: weak H/L chain match from {first_file}: "
            f"H={h_chain}:{h_score:.3f}, L={l_chain}:{l_score:.3f}"
        )
    chain_lengths = {cid: len(seq) for cid, seq in chain_sequences.items()}
    return h_chain, l_chain, h_score, l_score, chain_lengths


def build_selected_cache(
    *,
    dataset_root: Path,
    code_dir: Path,
    targets: list[str],
    target_meta: dict[str, dict[str, str]],
    max_conformers: int,
    cache_dir: Path,
    manifest_out: Path,
    chain_mapping_out: Path,
    force: bool,
) -> None:
    sys.path.insert(0, str(code_dir))
    import train_mca_v2 as T  # noqa: PLC0415

    grouped = grouped_conformer_rows(dataset_root)
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / "H").mkdir(exist_ok=True)
    (cache_dir / "L").mkdir(exist_ok=True)
    stale_failed = cache_dir / "failed.json"
    if stale_failed.exists():
        stale_failed.unlink()

    selected_rows: list[dict[str, object]] = []
    chain_mapping_rows: list[dict[str, object]] = []
    completed: list[str] = []
    failed: list[dict[str, object]] = []

    for idx, target in enumerate(targets, start=1):
        target_rows = grouped.get(target, [])[:max_conformers]
        if len(target_rows) != max_conformers:
            raise RuntimeError(
                f"{target}: expected {max_conformers} conformers, found {len(target_rows)}"
            )
        files = [dataset_root / row["copied_relative_path"] for row in target_rows]
        missing = [str(path) for path in files if not path.exists()]
        if missing:
            raise FileNotFoundError(f"{target}: missing copied CIFs: {missing[:5]}")

        h_chain, l_chain, h_score, l_score, chain_lengths = infer_hl_chain_ids(
            target=target,
            first_file=files[0],
            target_meta=target_meta[target],
        )
        chain_mapping_rows.append(
            {
                "target_id": target,
                "source_heavy_chain_id": h_chain,
                "source_light_chain_id": l_chain,
                "heavy_sequence_score": round(h_score, 6),
                "light_sequence_score": round(l_score, 6),
                "source_chain_lengths_json": json.dumps(chain_lengths, sort_keys=True),
            }
        )

        h_shard = cache_dir / "H" / f"{target}.pt"
        l_shard = cache_dir / "L" / f"{target}.pt"
        if h_shard.exists() and l_shard.exists() and not force:
            print(f"[cache] [{idx}/{len(targets)}] {target}: existing shards, skipping")
            completed.append(target)
        else:
            print(
                f"[cache] [{idx}/{len(targets)}] {target}: processing "
                f"{len(files)} copied CIFs",
                flush=True,
            )
            try:
                data_h, original_s_h = T._process_single_chain(
                    [str(path) for path in files],
                    h_chain,
                    max_len=115,
                )
                data_l, original_s_l = T._process_single_chain(
                    [str(path) for path in files],
                    l_chain,
                    max_len=107,
                )
                if int(original_s_h) != len(files) or int(original_s_l) != len(files):
                    raise RuntimeError(
                        f"{target}: processed H={original_s_h}, L={original_s_l}, "
                        f"expected {len(files)}"
                    )
                original_data = {"chains": {"H": data_h, "L": data_l}, "chain_ids": ["H", "L"]}
                T._save_one_shard(str(cache_dir), target, original_data, True, "H")
                completed.append(target)
            except Exception as exc:  # noqa: BLE001
                failed.append({"target_id": target, "error": repr(exc)})
                print(f"[cache/warn] {target}: failed: {exc}", flush=True)
                continue

        for conformer_index, row in enumerate(target_rows):
            selected_rows.append(
                {
                    "target_id": target,
                    "conformer_index": conformer_index,
                    "pick_rank": row["pick_rank"],
                    "copied_relative_path": row["copied_relative_path"],
                    "source_af3_model_path": row.get("source_af3_model_path", ""),
                    "source_heavy_chain_id": h_chain,
                    "source_light_chain_id": l_chain,
                    "sha256": row.get("sha256", ""),
                }
            )

    for chain_id in ["H", "L"]:
        with (cache_dir / chain_id / "index.json").open("w") as handle:
            json.dump({"group_names": completed}, handle, indent=2)

    write_tsv(
        manifest_out,
        selected_rows,
        [
            "target_id",
            "conformer_index",
            "pick_rank",
            "copied_relative_path",
            "source_af3_model_path",
            "source_heavy_chain_id",
            "source_light_chain_id",
            "sha256",
        ],
    )
    write_tsv(
        chain_mapping_out,
        chain_mapping_rows,
        [
            "target_id",
            "source_heavy_chain_id",
            "source_light_chain_id",
            "heavy_sequence_score",
            "light_sequence_score",
            "source_chain_lengths_json",
        ],
    )
    if failed:
        failed_path = cache_dir / "failed.json"
        failed_path.write_text(json.dumps(failed, indent=2))
        raise RuntimeError(f"{len(failed)} targets failed during cache build; see {failed_path}")


def export_embeddings(
    *,
    code_dir: Path,
    ids_txt: Path,
    cache_dir: Path,
    checkpoint: Path,
    emb_dir: Path,
    max_conformers: int,
    device: str,
) -> None:
    emb_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(code_dir / "export_mca_embeddings.py"),
        "--ids-txt",
        str(ids_txt),
        "--cache-path",
        str(cache_dir),
        "--model-pt",
        str(checkpoint),
        "--emb-dir",
        str(emb_dir),
        "--max-conformers",
        str(max_conformers),
        "--device",
        device,
    ]
    print("[export] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(code_dir), check=True)


def verify_embeddings(
    *,
    emb_dir: Path,
    targets: list[str],
    max_conformers: int,
    checkpoint: Path,
    out_path: Path,
) -> None:
    import torch  # noqa: PLC0415

    index = json.loads((emb_dir / "index.json").read_text())
    completed = index.get("group_names", [])
    problems = []
    for target in targets:
        shard_path = emb_dir / f"{target}.pt"
        if target not in completed:
            problems.append({"target_id": target, "error": "missing from index"})
            continue
        if not shard_path.exists():
            problems.append({"target_id": target, "error": "missing shard"})
            continue
        shard = torch.load(shard_path, map_location="cpu")
        if shard.get("source_checkpoint") != str(checkpoint.resolve()):
            problems.append({"target_id": target, "error": "source_checkpoint mismatch"})
        for chain_id in ["H", "L"]:
            reps = shard["chain_representations"][chain_id]["mca_repr"]
            if int(reps.shape[0]) != max_conformers:
                problems.append(
                    {
                        "target_id": target,
                        "chain_id": chain_id,
                        "error": f"n_conformers={int(reps.shape[0])}",
                    }
                )
    payload = {
        "ok": not problems,
        "n_requested_targets": len(targets),
        "n_completed_targets": len(completed),
        "max_conformers": max_conformers,
        "embedding_dir": str(emb_dir),
        "checkpoint": str(checkpoint.resolve()),
        "problems": problems,
    }
    out_path.write_text(json.dumps(payload, indent=2))
    if problems:
        raise RuntimeError(f"Embedding verification failed; see {out_path}")
    print(f"[verify] ok: {len(completed)} targets, {max_conformers} conformers each")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--shared-dir", required=True, type=Path)
    parser.add_argument("--scope", choices=["smoke", "full"], required=True)
    parser.add_argument("--max-conformers", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    start = time.time()
    run_dir = args.run_dir.resolve()
    shared_dir = args.shared_dir.resolve()
    dataset_root = shared_dir / "gaeun_ph_af3_medium_20260626"
    cdr_summary = shared_dir / "aim1_phase1_3_cdr_structural_20260626" / "cdr_diversity_summary.tsv"
    code_src = shared_dir / "mca1000_code_src"
    checkpoint = shared_dir / "checkpoints" / "mca1000_checkpoint.pt"
    code_dir = run_dir / "code" / "mca1000"
    cache_dir = run_dir / "cache" / "selected_mca_cache"
    emb_dir = run_dir / "embeddings" / "mca1000_selected"
    manifests_dir = run_dir / "manifests"
    logs_dir = run_dir / "logs"
    for path in [manifests_dir, logs_dir]:
        path.mkdir(parents=True, exist_ok=True)

    required = [dataset_root, cdr_summary, code_src, checkpoint]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing staged inputs: {missing}")

    copy_code_tree(code_src, code_dir)
    target_meta = load_target_meta(dataset_root)
    targets = load_target_ids(dataset_root, cdr_summary, args.scope)
    ids_txt = manifests_dir / "target_ids.txt"
    ids_txt.write_text("\n".join(targets) + "\n")
    write_tsv(
        manifests_dir / "target_selection.tsv",
        [{"target_id": target, "scope": args.scope} for target in targets],
        ["target_id", "scope"],
    )

    build_selected_cache(
        dataset_root=dataset_root,
        code_dir=code_dir,
        targets=targets,
        target_meta=target_meta,
        max_conformers=args.max_conformers,
        cache_dir=cache_dir,
        manifest_out=manifests_dir / "selected_conformers.tsv",
        chain_mapping_out=manifests_dir / "source_chain_mapping.tsv",
        force=args.force,
    )
    export_embeddings(
        code_dir=code_dir,
        ids_txt=ids_txt,
        cache_dir=cache_dir,
        checkpoint=checkpoint,
        emb_dir=emb_dir,
        max_conformers=args.max_conformers,
        device=args.device,
    )
    verify_embeddings(
        emb_dir=emb_dir,
        targets=targets,
        max_conformers=args.max_conformers,
        checkpoint=checkpoint,
        out_path=run_dir / "verification.json",
    )
    summary = {
        "scope": args.scope,
        "run_dir": str(run_dir),
        "shared_dir": str(shared_dir),
        "n_targets": len(targets),
        "max_conformers": args.max_conformers,
        "elapsed_seconds": round(time.time() - start, 3),
    }
    (run_dir / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print("[done] " + json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
