#!/bin/bash
#SBATCH --job-name=conff-p15a
#SBATCH --partition=liulab
#SBATCH --nodelist=ragonliu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/project/liulab/jg1920/conffusion/logs/phase15a_%j.out
#SBATCH --error=/project/liulab/jg1920/conffusion/logs/phase15a_%j.err

set -euo pipefail

PY="${PHASE15A_PY:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"
SCRIPT="${PHASE15A_SCRIPT:-/project/liulab/jg1920/conffusion/scripts/phase15a/phase15a_endpoint_audit.py}"
OUT_DIR="${PHASE15A_OUT_DIR:-/project/liulab/jg1920/conffusion/phase15a_20260627_endpoint_audit}"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "host=$(hostname)"
echo "date=$(date -Is)"
echo "job_id=${SLURM_JOB_ID:-none}"
echo "partition=${SLURM_JOB_PARTITION:-none}"
echo "node_list=${SLURM_JOB_NODELIST:-none}"
echo "python=$PY"
echo "script=$SCRIPT"
echo "out_dir=$OUT_DIR"

"$PY" "$SCRIPT" --out-dir "$OUT_DIR"

echo "completed=$(date -Is)"
