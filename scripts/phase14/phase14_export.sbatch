#!/bin/bash
set -euo pipefail

echo "host=$(hostname)"
echo "date=$(date -Is)"
echo "job_id=${SLURM_JOB_ID:-none}"
echo "node_list=${SLURM_JOB_NODELIST:-none}"
echo "gres=${SLURM_JOB_GRES:-none}"

: "${PHASE14_RUN_DIR:?PHASE14_RUN_DIR is required}"
: "${PHASE14_SHARED_DIR:?PHASE14_SHARED_DIR is required}"
: "${PHASE14_SCOPE:?PHASE14_SCOPE is required}"
: "${PHASE14_MAX_CONFORMERS:?PHASE14_MAX_CONFORMERS is required}"

SCRIPT_DIR="${PHASE14_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY="${PHASE14_PY:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"

export PYTHONUNBUFFERED=1
export WANDB_MODE=disabled
export MPLBACKEND=Agg

echo "python=$PY"
echo "script_dir=$SCRIPT_DIR"
echo "run_dir=$PHASE14_RUN_DIR"
echo "shared_dir=$PHASE14_SHARED_DIR"
echo "scope=$PHASE14_SCOPE"
echo "max_conformers=$PHASE14_MAX_CONFORMERS"

"$PY" "$SCRIPT_DIR/phase14_build_and_export.py" \
  --run-dir "$PHASE14_RUN_DIR" \
  --shared-dir "$PHASE14_SHARED_DIR" \
  --scope "$PHASE14_SCOPE" \
  --max-conformers "$PHASE14_MAX_CONFORMERS" \
  --device cuda

echo "completed=$(date -Is)"
