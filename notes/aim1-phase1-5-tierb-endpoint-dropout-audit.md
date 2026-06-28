# Aim 1 Phase 1.5 Tier B Endpoint Dropout Audit

Date: 2026-06-28

## Question

Why did the Tier B Step 4 endpoint use 160 targets when the `<0.85` Tier B
antibody-sequence-de-overlap pass manifest contains 172 targets?

## Short Answer

The 160-target Step 4 endpoint is exactly the intersection of:

1. targets with 100 selected model CIFs in the Tier B materialized dataset; and
2. targets with CDR structural distance rows from the Tier B CDR structural run.

The 12 dropped targets are explained:

- 8 targets had fewer than 100 selected conformers/model CIFs.
- 4 targets had 100 model CIFs, but all conformers failed the CDR structural
  chain-assignment gate, so no CDR distance or saturation rows were produced.

No additional target loss happened inside the Step 4 readout/control scripts:
`endpoint/target_overlap.tsv`, `endpoint/condition_manifest.tsv`, and
`phase15d/region_extraction_audit.tsv` all contain the same 160 targets.

## Evidence

Local files:

- Tier B pass manifest:
  `liu-lab/conffusion/manifests/tierb_sequence_deoverlap_20260627/gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_pass_20260627.tsv`
- Step 4 endpoint copy:
  `liu-lab/conffusion/figure_inputs/phase15_tierb_step4_20260627/endpoint/target_overlap.tsv`
- Step 4 condition manifest copy:
  `liu-lab/conffusion/figure_inputs/phase15_tierb_step4_20260627/endpoint/condition_manifest.tsv`
- Step 4 CDR extraction audit copy:
  `liu-lab/conffusion/figure_inputs/phase15_tierb_step4_20260627/phase15d/region_extraction_audit.tsv`

Remote read-only checks on Ragon:

- Dataset summary:
  `/project/liulab/jg1920/conffusion/tierb_vhvl85_dataset_20260627/manifests/source_summary.json`
- CDR structural summary:
  `/project/liulab/jg1920/conffusion/tierb_vhvl85_cdr_structural_20260627/run_summary.json`
- Phase 14 embedding/export target list:
  `/project/liulab/jg1920/conffusion/phase14_tierb_vhvl85_full100_20260627_215733/run/manifests/target_ids.txt`

Remote aggregate checks:

- Tier B targets requested: 172.
- Targets with sequence metadata, antigen labels, and at least one conformer:
  172.
- Targets with exactly 100 selected conformers: 164.
- Targets with any CDR distance rows: 167.
- Targets with both exactly 100 selected conformers and CDR distance rows: 160.
- The Phase 14 target list equals that 160-target intersection.

## Dropped Targets

| target | reason | selected conformers | notes |
| --- | --- | ---: | --- |
| `3wlw` | fewer than 100 conformers | 25 | CDR rows exist, but fixed-100 embedding/export excluded it. |
| `5cez` | fewer than 100 conformers | 1 | Not enough conformers for fixed-100 endpoint. |
| `6o3j` | fewer than 100 conformers | 44 | CDR rows exist, but fixed-100 embedding/export excluded it. |
| `7dq7` | fewer than 100 conformers | 89 | CDR rows exist, but fixed-100 embedding/export excluded it. |
| `7eya` | fewer than 100 conformers | 78 | CDR rows exist, but fixed-100 embedding/export excluded it. |
| `7uop` | fewer than 100 conformers | 79 | CDR rows exist, but fixed-100 embedding/export excluded it. |
| `7upk` | fewer than 100 conformers | 60 | CDR rows exist, but fixed-100 embedding/export excluded it. |
| `7upx` | fewer than 100 conformers | 18 | CDR rows exist, but fixed-100 embedding/export excluded it. |
| `2ybr` | CDR chain assignment failed | 100 | All 100 conformers failed: low chain assignment score heavy=0.979, light=0.823. |
| `9ctu` | CDR chain assignment failed | 100 | All 100 conformers failed: low chain assignment score heavy=0.805, light=0.766. |
| `9iw2` | CDR chain assignment failed | 100 | All 100 conformers failed: low chain assignment score heavy=0.752, light=0.750. |
| `9pzv` | CDR chain assignment failed | 100 | All 100 conformers failed: low chain assignment score heavy=0.825, light=0.795. |

## Interpretation

The 160/172 gap is not an unexplained Step 4 readout artifact. It is a
pre-readout eligibility gate caused by the fixed-100 conformer requirement and
the CDR structural assignment requirement.

For the current Phase 1 conclusion, it is fair to describe Step 4 as:

> a 160-target Tier B endpoint consisting of targets from the 172-target Tier B
> pass set that had exactly 100 selected conformers and usable CDR structural
> assignments.

This does not resolve the separate antigen/source de-overlap caveat. Tier B
remains antibody-sequence-de-overlapped against the checkpoint universe, but not
fully antigen/source-de-overlapped.

## Next Decision

If the goal is a fixed-100 readout benchmark, keep the current 160-target set
and document this eligibility gate.

If the goal is maximum Tier B coverage, rerun the endpoint with variable
conformer counts and decide how to handle targets with partial conformer sets
and CDR assignment failures. That would be a different endpoint and should not
be mixed with the current Step 4 benchmark table.
