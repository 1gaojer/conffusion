# Gaeun Conformer-Ensemble Manifest: Non-1000-Test Targets

Date: 2026-06-27

## Scope

Read-only inventory of Gaeun-owned conformer-generation output roots. The main TSV excludes target IDs that appear in the June 19/21 `1000` checkpoint split test set (`test_pdb_ids.txt`). It still includes flags for overlap with the checkpoint train/all universe.

## Outputs

- Main manifest, excludes 1000-test IDs: `gaeun_conformer_ensembles_non_1000_test_20260627.tsv`
- Stricter companion, excludes all 1000 checkpoint-universe IDs: `gaeun_conformer_ensembles_non_1000_all_20260627.tsv`
- Machine-readable summary: `gaeun_conformer_ensembles_non_1000_test_20260627_summary.json`


## Generated-Output-Only Companions

The original all-stage manifest includes targets that have only AF3 input/template files. For "ensembles already generated," use these generated-output-only TSVs instead:

- Generated outputs, excludes 1000-test IDs: `gaeun_conformer_ensembles_generated_non_1000_test_20260627.tsv` (1087 rows)
- Generated outputs, excludes the full 1000 checkpoint universe: `gaeun_conformer_ensembles_generated_non_1000_all_20260627.tsv` (359 rows)

## Source Roots

- `af3_input_jsons`: `/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results/af3_input_jsons`
- `af3_results`: `/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results/af3_results`
- `filtered_af3_models`: `/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results/filtered_af3_models`
- `prothunt_conf_gen_results`: `/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results/prothunt_output`
- `prothunt_dataset_root`: `/external/liulab/gkim/antigen_prediction/datasets/prothunt_output`

## 1000 Checkpoint Split Counts

- `all`: 930
- `test`: 202
- `train`: 728
- `trained_ids`: 930

## Manifest Counts

- Unique target IDs seen at any scanned stage: 1504
- Main manifest rows, not in 1000 test set: 1302
- Main rows with any generated output beyond input JSONs: 1087
- Main rows with AF3 result directories: 1070
- Main rows with filtered AF3 model target dirs: 616
- Main rows with Protein-Hunter output dirs: 745
- Main rows that are only AF3-input/template stage so far: 215
- Main rows that still overlap the 1000 train set: 728
- Main rows that still overlap the full 1000 checkpoint universe: 728
- Strict non-1000-universe rows: 574

## Stage Counts: Main Manifest

- `af3_inputs_only`: 215
- `af3_results`: 455
- `filtered_af3_models`: 616
- `prothunt_output`: 16

## Interpretation

The main manifest answers the literal request: generated/pipeline target IDs not overlapping the `1000` checkpoint test set. For clean external validation against the June 19 checkpoint, prefer the stricter companion manifest because the main manifest still includes targets overlapping the checkpoint train/all universe.
