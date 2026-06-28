#!/usr/bin/env bash
#SBATCH --job-name=conffusion-tierb-seq
#SBATCH --partition=liulab
#SBATCH --nodelist=ragonliu1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/project/liulab/jg1920/conffusion/logs/%x-%j.out

set -euo pipefail

PYTHON="${PYTHON:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"
SCRIPT="${SCRIPT:-/project/liulab/jg1920/conffusion/scripts/manifests/build_tierb_sequence_deoverlap.py}"
WORKERS="${WORKERS:-${SLURM_CPUS_PER_TASK:-32}}"

echo "job_id=${SLURM_JOB_ID:-manual}"
echo "host=$(hostname)"
echo "python=$PYTHON"
echo "script=$SCRIPT"
echo "workers=$WORKERS"
echo "started_at=$(date -Is)"
echo "log_hygiene=safe_summary_only"
echo "log_hygiene_note=do_not_print_raw_sequences_fasta_msa_json_or_structure_payloads"

"$PYTHON" "$SCRIPT" --workers "$WORKERS"

echo "finished_at=$(date -Is)"
