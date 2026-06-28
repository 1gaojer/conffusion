# Aim 1 Phase 1.5 Tier C Leakage-Sensitivity Evaluation

Date: 2026-06-28

## Scope

This run re-scored the Tier B Step 4 endpoint after exact antigen/source-pair
de-overlap against the June 19/21 `1000 conformer` checkpoint universe.

This is a leakage-sensitivity check, not a stable external validation endpoint:
the Tier C filter leaves only 17 self-match-excluded evaluable same-label
queries under the current exact-label retrieval metric.

## Inputs And Execution

Primary input endpoint:

- `/external/liulab/jg1920/conffusion/phase15_tierb_step4_full_20260627_233600/endpoint`

Reused input bundle:

- `/external/liulab/jg1920/conffusion/phase15_tierb_step4_inputs_20260627_233432`
- Embeddings: `embeddings/mca1000_selected`
- CDR assignments:
  `tierb_vhvl85_cdr_structural_20260627/cdr_numbering_assignments.tsv`

Tier C filter inputs:

- Checkpoint IDs:
  `/external/liulab/gkim/antigen_prediction/2026.06.19_retrain_mca_1000_confs/split_confs/all_pdb_ids.txt`
- Antigen/source labels:
  `/external/liulab/gkim/antigen_prediction/datasets/antigen_labels/full_sabdab_relabeled_antigens.csv`

Run:

- Slurm job: `32258`
- Node/partition: `liulab` / `ragonliu1`
- Resources: CPU-only, 8 CPUs, 64 GB
- Runtime: 2026-06-28 07:02:05 to 07:02:25 EDT
- Exit code: `0:0`
- Remote output root:
  `/external/liulab/jg1920/conffusion/phase15_tierc_leakage_20260628_070200`
- Local lightweight copy:
  `figure_inputs/phase15_tierc_leakage_20260628/`

Scripts:

- `scripts/phase15a/phase15_tierc_filter_endpoint.py`
- `scripts/phase15a/phase15_tierc_leakage_eval.sbatch`
- Existing readout/control scripts:
  `scripts/phase15a/phase15c_region_readout.py` and
  `scripts/phase15a/phase15d_region_controls.py`

## Endpoint Summary

Filter rule:

- Keep an endpoint target only if none of its exact normalized
  `antigen_name || antigen_species` pairs occur among the 930 checkpoint IDs.

Observed endpoint:

| endpoint | targets | condition rows |
| --- | ---: | ---: |
| Tier B Step 4 input | 160 | 2,880 |
| Tier C retained | 87 | 1,566 |
| Tier C dropped | 73 | not applicable |

Other checks:

- Checkpoint IDs: 930
- Checkpoint antigen/source pairs: 543
- Missing Tier C labels: 0
- Retained labels with at least two targets: 8
- Self-match-excluded evaluable queries: 17

## Full-Ensemble Readout Result

Metrics below are self-match-excluded and count only the 17 queries with at
least one other retained target sharing the same exact positive label.

| readout | Recall@10 | MRR | median first positive rank |
| --- | ---: | ---: | ---: |
| all CDRs mean+std | 0.471 | 0.283 | 11.0 |
| H-CDR3 | 0.412 | 0.174 | 15.0 |
| all CDRs H/L mean | 0.353 | 0.154 | 14.0 |
| global + H-CDR3 | 0.412 | 0.152 | 16.0 |
| global + all CDRs | 0.412 | 0.136 | 13.0 |
| global H/L mean | 0.353 | 0.130 | 17.0 |

## Controls

Selected Phase 15D full-ensemble controls:

| readout/control | Recall@10 | MRR | median first positive rank |
| --- | ---: | ---: | ---: |
| all CDRs mean+std | 0.471 | 0.283 | 11.0 |
| L-CDR3 | 0.412 | 0.197 | 23.0 |
| H-CDR3 | 0.412 | 0.174 | 15.0 |
| global H/L mean | 0.353 | 0.130 | 17.0 |
| framework H/L mean | 0.294 | 0.131 | 16.0 |

Shuffled-label controls were lower:

| readout | shuffled Recall@10 mean | shuffled MRR mean |
| --- | ---: | ---: |
| global H/L mean | 0.135 | 0.062 |
| H-CDR3 | 0.112 | 0.045 |
| all CDRs H/L mean | 0.144 | 0.061 |
| all CDRs mean+std | 0.135 | 0.064 |
| framework H/L mean | 0.162 | 0.065 |

## Subset Check Under Best Readout

Under `all_cdrs_mean_std`, the sparse Tier C subset does not establish a clear
selector winner:

| condition | Recall@10 | MRR | median first positive rank |
| --- | ---: | ---: | ---: |
| full 100 conformers | 0.471 | 0.283 | 11.0 |
| K64 all-CDR k-center | 0.412 | 0.279 | 15.0 |
| K64 H-CDR3 k-center | 0.412 | 0.262 | 13.0 |
| random K64, best replicate | 0.529 | 0.309 | 10.0 |
| random K64, worst replicate | 0.353 | 0.220 | 14.0 |
| first K64 | 0.529 | 0.287 | 8.0 |
| single first conformer | 0.235 | 0.104 | 42.0 |

## Interpretation

Observed:

- The exact antigen/source-pair de-overlap run completed successfully without
  rebuilding embeddings or CDR assignments.
- The Tier C filter reproduces the prior count: 87 retained targets and 73
  dropped targets.
- CDR/paratope-aware `all_cdrs_mean_std` remains better than global H/L pooling
  on the retained Tier C subset.
- Shuffled-label controls are substantially lower than the true-label
  `all_cdrs_mean_std` result.
- Single-first-conformer performance is weaker, consistent with the earlier
  ensemble-versus-single diagnostic.

Interpretation:

- The leakage-sensitivity result does not collapse to the shuffled-label
  baseline after antigen/source-pair de-overlap.
- The direction of the readout effect is consistent with Tier B: CDR/paratope
  readout is more informative than global H/L pooling.
- The result is too sparse to use as a headline validation metric. With only 17
  evaluable queries, individual labels and targets can move the metric
  substantially.
- This supports showing Tier C as a caveated leakage stress test, not as proof
  of external generalization.

## Recommended Next Step

Use this as supporting evidence for the readout-bottleneck claim, while keeping
Tier B as the main powered endpoint for now. Before claiming cleaner external
validation, change or expand the endpoint so that antigen/source-de-overlapped
targets still have enough same-label positives for stable retrieval metrics.
