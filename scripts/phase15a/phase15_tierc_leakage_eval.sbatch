#!/usr/bin/env bash
#SBATCH --job-name=phase15-tierc
#SBATCH --partition=liulab
#SBATCH --nodelist=ragonliu1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

RUN_ROOT="${RUN_ROOT:?set RUN_ROOT}"
SCRIPT_DIR="${SCRIPT_DIR:-$RUN_ROOT/scripts/phase15a}"
PYTHON="${PYTHON:-/project/liulab/jg1920/envs/jerryenv1/bin/python}"

TIERB_ENDPOINT="${TIERB_ENDPOINT:-/external/liulab/jg1920/conffusion/phase15_tierb_step4_full_20260627_233600/endpoint}"
INPUT_BUNDLE="${INPUT_BUNDLE:-/external/liulab/jg1920/conffusion/phase15_tierb_step4_inputs_20260627_233432}"
CHECKPOINT_IDS="${CHECKPOINT_IDS:-/external/liulab/gkim/antigen_prediction/2026.06.19_retrain_mca_1000_confs/split_confs/all_pdb_ids.txt}"
ANTIGEN_LABELS="${ANTIGEN_LABELS:-/external/liulab/gkim/antigen_prediction/datasets/antigen_labels/full_sabdab_relabeled_antigens.csv}"
EMBEDDING_DIR="${EMBEDDING_DIR:-$INPUT_BUNDLE/embeddings/mca1000_selected}"
CDR_ASSIGNMENTS="${CDR_ASSIGNMENTS:-$INPUT_BUNDLE/tierb_vhvl85_cdr_structural_20260627/cdr_numbering_assignments.tsv}"

mkdir -p "$RUN_ROOT"/endpoint "$RUN_ROOT"/phase15c "$RUN_ROOT"/phase15d "$RUN_ROOT"/logs

echo "job_id=${SLURM_JOB_ID:-}"
echo "host=$(hostname)"
echo "run_root=$RUN_ROOT"
echo "python=$PYTHON"
echo "tierb_endpoint=$TIERB_ENDPOINT"
echo "embedding_dir=$EMBEDDING_DIR"
echo "cdr_assignments=$CDR_ASSIGNMENTS"

"$PYTHON" "$SCRIPT_DIR/phase15_tierc_filter_endpoint.py" \
  --endpoint-dir "$TIERB_ENDPOINT" \
  --checkpoint-ids "$CHECKPOINT_IDS" \
  --antigen-labels "$ANTIGEN_LABELS" \
  --out-dir "$RUN_ROOT/endpoint"

"$PYTHON" "$SCRIPT_DIR/phase15c_region_readout.py" \
  --phase15a-run "$RUN_ROOT/endpoint" \
  --embedding-dir "$EMBEDDING_DIR" \
  --cdr-assignments "$CDR_ASSIGNMENTS" \
  --out-dir "$RUN_ROOT/phase15c"

"$PYTHON" "$SCRIPT_DIR/phase15d_region_controls.py" \
  --phase15a-run "$RUN_ROOT/endpoint" \
  --embedding-dir "$EMBEDDING_DIR" \
  --cdr-assignments "$CDR_ASSIGNMENTS" \
  --out-dir "$RUN_ROOT/phase15d"

find "$RUN_ROOT" -maxdepth 2 -type f | sort
