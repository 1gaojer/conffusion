#!/bin/bash
set -euo pipefail

echo "host=$(hostname)"
echo "date=$(date -Is)"
echo "job_id=${SLURM_JOB_ID:-none}"
echo "node_list=${SLURM_JOB_NODELIST:-none}"

: "${PHASE14_ANALYSIS_RUN_ROOT:?PHASE14_ANALYSIS_RUN_ROOT is required}"
: "${PHASE14_ANALYSIS_DATASET_ROOT:?PHASE14_ANALYSIS_DATASET_ROOT is required}"
: "${PHASE14_ANALYSIS_CDR_SUMMARY:?PHASE14_ANALYSIS_CDR_SUMMARY is required}"
: "${PHASE14_ANALYSIS_CDR_SCRIPT:?PHASE14_ANALYSIS_CDR_SCRIPT is required}"
: "${PHASE14_ANALYSIS_OUT_DIR:?PHASE14_ANALYSIS_OUT_DIR is required}"
: "${PHASE14_ANALYSIS_SCOPE:?PHASE14_ANALYSIS_SCOPE is required}"

SCRIPT_DIR="${PHASE14_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY="${PHASE14_PY:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"

export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg

args=(
  "$SCRIPT_DIR/phase14_analyze_embeddings.py"
  --run-root "$PHASE14_ANALYSIS_RUN_ROOT"
  --dataset-root "$PHASE14_ANALYSIS_DATASET_ROOT"
  --cdr-summary "$PHASE14_ANALYSIS_CDR_SUMMARY"
  --cdr-script "$PHASE14_ANALYSIS_CDR_SCRIPT"
  --out-dir "$PHASE14_ANALYSIS_OUT_DIR"
  --scope-label "$PHASE14_ANALYSIS_SCOPE"
  --force
)

if [[ -n "${PHASE14_ANALYSIS_TARGET_IDS_FILE:-}" ]]; then
  args+=(--target-ids-file "$PHASE14_ANALYSIS_TARGET_IDS_FILE")
fi

if [[ -n "${PHASE14_ANALYSIS_MAX_TARGETS:-}" ]]; then
  args+=(--max-targets "$PHASE14_ANALYSIS_MAX_TARGETS")
fi

if [[ -n "${PHASE14_ANALYSIS_WORKERS:-}" ]]; then
  args+=(--workers "$PHASE14_ANALYSIS_WORKERS")
fi

echo "python=$PY"
echo "script_dir=$SCRIPT_DIR"
echo "run_root=$PHASE14_ANALYSIS_RUN_ROOT"
echo "dataset_root=$PHASE14_ANALYSIS_DATASET_ROOT"
echo "cdr_summary=$PHASE14_ANALYSIS_CDR_SUMMARY"
echo "cdr_script=$PHASE14_ANALYSIS_CDR_SCRIPT"
echo "out_dir=$PHASE14_ANALYSIS_OUT_DIR"
echo "scope=$PHASE14_ANALYSIS_SCOPE"
echo "workers=${PHASE14_ANALYSIS_WORKERS:-1}"

"$PY" "${args[@]}"

echo "completed=$(date -Is)"
