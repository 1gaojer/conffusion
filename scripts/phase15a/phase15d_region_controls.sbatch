#!/bin/bash
#SBATCH --job-name=conff-p15d
#SBATCH --partition=liulab
#SBATCH --nodelist=ragonliu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/project/liulab/jg1920/conffusion/logs/phase15d_%j.out
#SBATCH --error=/project/liulab/jg1920/conffusion/logs/phase15d_%j.err

set -euo pipefail

PY="${PHASE15D_PY:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"
SCRIPT="${PHASE15D_SCRIPT:-/project/liulab/jg1920/conffusion/scripts/phase15a/phase15d_region_controls.py}"
OUT_DIR="${PHASE15D_OUT_DIR:-/project/liulab/jg1920/conffusion/phase15d_20260627_region_controls}"
LIMIT_TARGETS="${PHASE15D_LIMIT_TARGETS:-0}"
N_RANDOM_WINDOW_REPS="${PHASE15D_N_RANDOM_WINDOW_REPS:-10}"
N_SHUFFLE_REPS="${PHASE15D_N_SHUFFLE_REPS:-20}"
SEED="${PHASE15D_SEED:-20260627}"

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
echo "limit_targets=$LIMIT_TARGETS"
echo "n_random_window_reps=$N_RANDOM_WINDOW_REPS"
echo "n_shuffle_reps=$N_SHUFFLE_REPS"
echo "seed=$SEED"

"$PY" "$SCRIPT" \
  --out-dir "$OUT_DIR" \
  --limit-targets "$LIMIT_TARGETS" \
  --n-random-window-reps "$N_RANDOM_WINDOW_REPS" \
  --n-shuffle-reps "$N_SHUFFLE_REPS" \
  --seed "$SEED"

echo "completed=$(date -Is)"
