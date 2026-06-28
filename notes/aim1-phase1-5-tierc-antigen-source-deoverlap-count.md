# Aim 1 Phase 1.5 Tier C Antigen/Source De-Overlap Count

Date: 2026-06-28

## Question

How many of the 160 fixed-100 Tier B Step 4 targets remain after
antigen/source de-overlap against the June 19/21 `1000 conformer` checkpoint
universe?

## Data Sources

Step 4 fixed-100 target list:

- `/project/liulab/jg1920/conffusion/phase14_tierb_vhvl85_full100_20260627_215733/run/manifests/target_ids.txt`

Checkpoint universe:

- `/external/liulab/gkim/antigen_prediction/2026.06.19_retrain_mca_1000_confs/split_confs/all_pdb_ids.txt`
- Count used for this audit: 930 unique checkpoint PDB IDs.

Antigen/source labels:

- `/external/liulab/gkim/antigen_prediction/datasets/antigen_labels/full_sabdab_relabeled_antigens.csv`
- Fields used: normalized `antigen_name` and `antigen_species`.

## Rule

Primary Tier C rule used here:

- Drop a Step 4 target if any exact normalized
  `antigen_name || antigen_species` pair appears in any of the 930 checkpoint
  universe IDs.

This is exact label-string de-overlap using Gaeun's relabeled SAbDab antigen
metadata. It is not antigen sequence-identity de-overlap.

## Result

| filter | remain | drop |
| --- | ---: | ---: |
| Tier B Step 4 fixed-100 endpoint | 160 | 0 |
| Tier C exact antigen/source-pair de-overlap | 87 | 73 |

No Step 4 targets were missing antigen/source labels under this label source.

After exact antigen/source-pair de-overlap, the 87 remaining targets span 78
distinct positive labels. Only 8 labels still have at least two targets, covering
17 targets with a same-label neighbor. This means the Tier C target count is
reasonable for a count audit, but the current same-label retrieval endpoint
would become sparse unless the metric or candidate bank is adjusted.

## Sensitivity Checks

| stricter/alternate rule | remain | drop | interpretation |
| --- | ---: | ---: | --- |
| Exact `antigen_name || antigen_species` pair | 87 | 73 | Recommended Tier C count for the current labels. |
| Antigen-name-only overlap | 73 | 87 | Stricter because it ignores source/species differences. |
| Source/species-only overlap | 31 | 129 | Very strict; drops targets just because the organism/source was present. |
| Combined positive-key string overlap | 87 | 73 | Same count as exact pair in this endpoint. |

## Main Dropped Overlap Labels

The exact antigen/source-pair drops are dominated by checkpoint-overlapping
SARS-CoV-2 spike labels:

| overlapping antigen/source label | dropped targets |
| --- | ---: |
| `spike protein s1_1 || severe acute respiratory syndrome coronavirus2` | 20 |
| `spike protein s1_2 || severe acute respiratory syndrome coronavirus2` | 11 |
| `spike protein s1_6 || severe acute respiratory syndrome coronavirus2` | 4 |
| `neuraminidase_5 || influenza a virus` | 4 |
| `hiv fusion peptide residue 512-519 || human immunodeficiency virus 1` | 3 |
| `spike glycoprotein_1 || severe acute respiratory syndrome coronavirus2` | 3 |

## Interpretation

Under exact antigen/source-pair de-overlap, Tier C would keep 87 of the current
160 fixed-100 Step 4 targets. This resolves the first Tier C feasibility count.

However, the same-label retrieval setup becomes sparse after this filter:
most retained labels are singletons. Before rerunning Step 4 metrics on Tier C,
decide whether the metric should remain same-label nearest-neighbor retrieval,
whether the candidate bank should be expanded, or whether the endpoint should
switch to a different antigen/source grouping.
