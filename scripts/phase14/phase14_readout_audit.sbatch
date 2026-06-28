#!/bin/bash
#SBATCH --job-name=phase14-readout-audit
#SBATCH --partition=compute
#SBATCH --nodelist=lambda-scalar
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=1-08:00:00
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

echo "host=$(hostname)"
echo "date=$(date -Is)"
echo "job_id=${SLURM_JOB_ID:-none}"
echo "node_list=${SLURM_JOB_NODELIST:-none}"
echo "cpus=${SLURM_CPUS_PER_TASK:-unset}"

: "${PHASE14_READOUT_RUN_ROOT:?PHASE14_READOUT_RUN_ROOT is required}"
: "${PHASE14_READOUT_DATASET_ROOT:?PHASE14_READOUT_DATASET_ROOT is required}"
: "${PHASE14_READOUT_CDR_SUMMARY:?PHASE14_READOUT_CDR_SUMMARY is required}"
: "${PHASE14_READOUT_CDR_SCRIPT:?PHASE14_READOUT_CDR_SCRIPT is required}"
: "${PHASE14_READOUT_OUT_DIR:?PHASE14_READOUT_OUT_DIR is required}"
: "${PHASE14_READOUT_SCOPE:?PHASE14_READOUT_SCOPE is required}"

if [[ -n "${PHASE14_SCRIPT_DIR:-}" ]]; then
  SCRIPT_DIR="$PHASE14_SCRIPT_DIR"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "$SLURM_SUBMIT_DIR/scripts/phase14/phase14_readout_audit.py" ]]; then
  SCRIPT_DIR="$SLURM_SUBMIT_DIR/scripts/phase14"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
PY="${PHASE14_PY:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"
WORKERS="${PHASE14_READOUT_WORKERS:-12}"

mkdir -p logs "$PHASE14_READOUT_OUT_DIR"

export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PHASE14_TORCH_THREADS="${PHASE14_TORCH_THREADS:-1}"

args=(
  "$SCRIPT_DIR/phase14_readout_audit.py"
  --run-root "$PHASE14_READOUT_RUN_ROOT"
  --dataset-root "$PHASE14_READOUT_DATASET_ROOT"
  --cdr-summary "$PHASE14_READOUT_CDR_SUMMARY"
  --cdr-script "$PHASE14_READOUT_CDR_SCRIPT"
  --out-dir "$PHASE14_READOUT_OUT_DIR"
  --scope-label "$PHASE14_READOUT_SCOPE"
  --workers "$WORKERS"
  --force
)

if [[ -n "${PHASE14_READOUT_TARGET_IDS_FILE:-}" ]]; then
  args+=(--target-ids-file "$PHASE14_READOUT_TARGET_IDS_FILE")
fi

if [[ -n "${PHASE14_READOUT_MAX_TARGETS:-}" ]]; then
  args+=(--max-targets "$PHASE14_READOUT_MAX_TARGETS")
fi

echo "python=$PY"
echo "script_dir=$SCRIPT_DIR"
echo "run_root=$PHASE14_READOUT_RUN_ROOT"
echo "dataset_root=$PHASE14_READOUT_DATASET_ROOT"
echo "cdr_summary=$PHASE14_READOUT_CDR_SUMMARY"
echo "cdr_script=$PHASE14_READOUT_CDR_SCRIPT"
echo "out_dir=$PHASE14_READOUT_OUT_DIR"
echo "scope=$PHASE14_READOUT_SCOPE"
echo "workers=$WORKERS"

"$PY" "${args[@]}"

echo "completed=$(date -Is)"
