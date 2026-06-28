#!/bin/bash
#SBATCH --job-name=conffusion-tierb-cdr
#SBATCH --partition=liulab
#SBATCH --nodelist=ragonliu1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/project/liulab/jg1920/conffusion/logs/%x-%j.out
#SBATCH --error=/project/liulab/jg1920/conffusion/logs/%x-%j.err

set -euo pipefail

PY="${PY:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"
SCRIPT="${SCRIPT:-/project/liulab/jg1920/conffusion/scripts/aim1_phase1_cdr_structural.py}"
DATASET_ROOT="${DATASET_ROOT:-/project/liulab/jg1920/conffusion/tierb_vhvl85_dataset_20260627}"
OUT_DIR="${OUT_DIR:-/project/liulab/jg1920/conffusion/tierb_vhvl85_cdr_structural_20260627}"
MAX_TARGETS="${MAX_TARGETS:-}"
MAX_CONFORMERS_PER_TARGET="${MAX_CONFORMERS_PER_TARGET:-100}"
VALIDATE_ONLY="${VALIDATE_ONLY:-0}"

mkdir -p /project/liulab/jg1920/conffusion/logs

echo "host=$(hostname)"
echo "date=$(date -Is)"
echo "job_id=${SLURM_JOB_ID:-none}"
echo "python=$PY"
echo "script=$SCRIPT"
echo "dataset_root=$DATASET_ROOT"
echo "out_dir=$OUT_DIR"
echo "max_targets=${MAX_TARGETS:-all}"
echo "max_conformers_per_target=$MAX_CONFORMERS_PER_TARGET"
echo "validate_only=$VALIDATE_ONLY"

args=(
  "$SCRIPT"
  --dataset-root "$DATASET_ROOT"
  --out-dir "$OUT_DIR"
  --max-conformers-per-target "$MAX_CONFORMERS_PER_TARGET"
)

if [[ -n "$MAX_TARGETS" ]]; then
  args+=(--max-targets "$MAX_TARGETS")
fi

if [[ "$VALIDATE_ONLY" == "1" ]]; then
  args+=(--validate-only)
fi

"$PY" "${args[@]}"

echo "completed=$(date -Is)"
