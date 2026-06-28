#!/bin/bash
#SBATCH --job-name=conffusion-tierb-stage
#SBATCH --partition=liulab
#SBATCH --nodelist=ragonliu1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/project/liulab/jg1920/conffusion/logs/%x-%j.out
#SBATCH --error=/project/liulab/jg1920/conffusion/logs/%x-%j.err

set -euo pipefail

PY="${PY:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"
SCRIPT="${SCRIPT:-/project/liulab/jg1920/conffusion/scripts/manifests/build_tierb_dataset.py}"
OUT_ROOT="${OUT_ROOT:-/project/liulab/jg1920/conffusion/tierb_vhvl85_dataset_20260627}"
MAX_CONFORMERS="${MAX_CONFORMERS:-100}"
MATERIALIZE_MODE="${MATERIALIZE_MODE:-symlink}"

mkdir -p /project/liulab/jg1920/conffusion/logs

echo "host=$(hostname)"
echo "date=$(date -Is)"
echo "job_id=${SLURM_JOB_ID:-none}"
echo "python=$PY"
echo "script=$SCRIPT"
echo "out_root=$OUT_ROOT"
echo "max_conformers=$MAX_CONFORMERS"
echo "materialize_mode=$MATERIALIZE_MODE"
echo "log_hygiene=safe_summary_only"
echo "log_hygiene_note=do_not_print_raw_sequences_fasta_msa_json_or_structure_payloads"

"$PY" "$SCRIPT" \
  --out-root "$OUT_ROOT" \
  --max-conformers "$MAX_CONFORMERS" \
  --materialize-mode "$MATERIALIZE_MODE" \
  --force

echo "completed=$(date -Is)"
