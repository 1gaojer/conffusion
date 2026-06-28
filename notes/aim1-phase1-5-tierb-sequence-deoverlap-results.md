# Aim 1 Phase 1.5 Tier B Sequence De-Overlap Results

Date: 2026-06-27

## Purpose

Build a stricter Gaeun-generated conformer endpoint for Conffusion follow-up
experiments.

Tier A used Gaeun's generated conformer outputs that do not overlap the full
June 19/21 `1000` checkpoint PDB universe. Tier B starts from that 359-target
set and removes targets whose antibody variable-region sequence is still too
similar to anything in the checkpoint universe.

## Method

Inputs:

- Tier A candidate manifest:
  `manifests/gaeun_conformer_ensembles_generated_non_1000_all_20260627.tsv`
- Checkpoint-universe IDs:
  `/external/liulab/gkim/antigen_prediction/2026.06.19_retrain_mca_1000_confs/split_confs/all_pdb_ids.txt`
- Sequence source:
  `/external/liulab/gkim/antigen_prediction/datasets/seq_files/full_sabdab_seq_file.json`

Sequence identity rule:

- Compare concatenated `VH+VL` variable-region sequences.
- Use global Needleman-Wunsch through Biopython `PairwiseAligner`.
- Use Gaeun's split-script identity definition:
  `identical_aligned_columns / max(len_a, len_b)`.
- Keep a target only if its maximum `VH+VL` identity to all 930 checkpoint
  targets is below the selected threshold.
- Current selected threshold: `<0.85`.

Run:

- Slurm job: `32217`
- Node/partition: `liulab` / `ragonliu1`
- Script:
  `scripts/manifests/build_tierb_sequence_deoverlap.py`
- Cluster output:
  `/project/liulab/jg1920/conffusion/tierb_sequence_deoverlap_20260627/`

## Results

- Tier A candidates scored: 359
- Checkpoint IDs compared against: 930
- Checkpoint IDs with usable sequences: 930
- Missing candidate sequences: 0
- Tier B pass at `<0.85`: 172
- Tier B fail at `<0.85`: 187
- Runtime: 8.25 seconds with 32 CPU workers

Local outputs:

- Scored manifest:
  `manifests/tierb_sequence_deoverlap_20260627/gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_scored_20260627.tsv`
- Tier B pass manifest:
  `manifests/tierb_sequence_deoverlap_20260627/gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_pass_20260627.tsv`
- Tier B fail manifest:
  `manifests/tierb_sequence_deoverlap_20260627/gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_fail_20260627.tsv`
- Machine-readable summary:
  `manifests/tierb_sequence_deoverlap_20260627/tierb_sequence_deoverlap_summary.json`

Identity distribution:

- Passing targets at `<0.85`: min `0.613445`, median `0.792738`, max `0.848101`
- Failing targets at `<0.85`: min `0.850877`, median `0.991304`, max `1.000000`
- Failing targets with max `VH+VL` identity at least `0.90`: 142
- Failing targets with max `VH+VL` identity at least `0.95`: 108
- Failing targets with max `VH+VL` identity exactly `1.00`: 82

## Interpretation

Tier B is a real tightening step. The 359 exact-PDB-unseen generated targets
still contained many antibodies that were highly sequence-similar to the
checkpoint universe. After applying Gaeun-style `VH+VL` sequence de-overlap,
172 generated-output targets remain at the selected `<0.85` cutoff. The older
`<0.80` files are retained as a stricter archival subset with 100 passing
targets.

This is now the cleanest current Gaeun-generated conformer candidate pool for
checkpoint-unseen testing.

## Caveats

- Tier B is antibody-sequence de-overlapped against the checkpoint universe, not
  antigen-de-overlapped.
- Pass/fail is based on concatenated `VH+VL`. Some passing targets can still
  have high heavy-only or light-only identity annotations, because the official
  rule used here follows the concatenated split criterion.
- This manifest does not by itself guarantee that embeddings or conformers are
  available in the exact format needed for the next retrieval/readout run.

## Next Step

Use the Tier B pass manifest as the current Gaeun-generated test pool. The next
strictening layer should be Tier C antigen/source de-overlap, or a practical
retrieval-readout run on the Tier B pass set if the needed conformer embeddings
already exist or can be exported cleanly.
