# Tier B Sequence De-Overlap (vhvl85)

Date: 2026-06-27

## Definition

Tier B starts from the 359 generated-output targets that do not overlap
the full June 19/21 `1000` checkpoint PDB universe. It keeps only
targets whose concatenated VH+VL variable-region global identity is
below 0.85 to every checkpoint-universe target.

Identity mirrors Gaeun's split script: global Needleman-Wunsch,
`identical_aligned_columns / max(len_a, len_b)`, using concatenated
VH+VL sequences.

## Counts

- Tier A candidates scored: 359
- Tier B pass: 172
- Tier B fail: 187
- Missing identity values: 0

## Outputs

- Scored all candidates: `gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_scored_20260627.tsv`
- Tier B pass manifest: `gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_pass_20260627.tsv`
- Tier B fail manifest: `gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_fail_20260627.tsv`
- Machine-readable summary: `tierb_sequence_deoverlap_vhvl85_summary.json`

## Caveat

Tier B is antibody-sequence de-overlapped against the checkpoint universe.
It does not yet enforce antigen-name, antigen-sequence, or same-antigen
group de-overlap.
