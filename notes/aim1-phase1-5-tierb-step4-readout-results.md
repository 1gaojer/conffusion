# Aim 1 Phase 1.5 Tier B Step 4 Readout Results

Date: 2026-06-27

## Scope

Step 4 reran the Phase 1.5 CDR/paratope-aware retrieval diagnostic on the
stricter Tier B Gaeun-generated endpoint.

This is cleaner than the original internal 149-target diagnostic because the
endpoint starts from Gaeun-generated conformer outputs that are exact-PDB-unseen
relative to the June 19/21 `1000` checkpoint universe and then antibody
sequence-de-overlapped at `<0.85` concatenated `VH+VL` identity. It is still not
fully external validation because antigen/source de-overlap is not enforced in
the Step 4 metrics.

## Inputs And Execution

Input bundle staged on the `liulab` side:

- `/external/liulab/jg1920/conffusion/phase15_tierb_step4_inputs_20260627_233432`
- 160 embedding shards
- Tier B target IDs
- CDR curves and CDR assignment tables
- Antigen/source labels
- Phase 1.5 scripts copied from Jerry's `conffusion` workspace

Full run:

- Slurm job: `32255`
- Node/partition: `liulab` / `ragonliu1`
- Resources: CPU-only, 8 CPUs, 64 GB
- Runtime: 2026-06-27 23:36:00 to 23:36:17 EDT
- Exit code: `0:0`
- Output root:
  `/external/liulab/jg1920/conffusion/phase15_tierb_step4_full_20260627_233600`

Scripts used or patched for this run:

- Added `scripts/phase15a/phase15_tierb_prepare_endpoint.py`
- Patched `scripts/phase15a/phase15c_region_readout.py` to infer conformer
  count from each shard instead of assuming 128 conformers, and to normalize
  kappa/light-chain CDR assignments to MCA `L` tensors.
- Patched `scripts/phase15a/phase15d_region_controls.py` to use the Tier B
  positive-label field and avoid hardcoded internal-bank wording.

## Output Files

Endpoint preparation:

- `endpoint/condition_manifest.tsv`
- `endpoint/target_overlap.tsv`
- `endpoint/tierb_endpoint_prep_summary.json`

Phase 1.5c readout rerun:

- `phase15c/phase15c_region_readout_report.md`
- `phase15c/region_readout_condition_summary.tsv`
- `phase15c/region_readout_retrieval_results.tsv`
- `phase15c/region_readout_random_rollup.tsv`
- `phase15c/region_readout_delta_vs_global.tsv`
- `phase15c/region_readout_vector_preservation.tsv`

Phase 1.5d controls:

- `phase15d/phase15d_region_control_report.md`
- `phase15d/region_control_summary.tsv`
- `phase15d/region_control_retrieval_results.tsv`
- `phase15d/region_control_random_window_rollup.tsv`
- `phase15d/region_control_shuffled_label_results.tsv`
- `phase15d/region_control_shuffled_label_summary.tsv`
- `phase15d/region_extraction_audit.tsv`

## Endpoint Summary

- Targets with usable endpoint data: 160
- Eligibility audit: the endpoint is the 160-target intersection of the
  172-target Tier B pass set with exactly 100 selected conformers and usable CDR
  structural distance rows. See
  `notes/aim1-phase1-5-tierb-endpoint-dropout-audit.md`.
- Condition rows: 2,880
- Conformers per target in this endpoint: 100
- Random subset replicates: 5
- Labels with at least two targets: 20
- Self-match-excluded evaluable queries in the retrieval summaries: 74

Largest endpoint labels:

| positive label | targets |
| --- | ---: |
| `spike protein s1_1 || severe acute respiratory syndrome coronavirus2` | 20 |
| `spike protein s1_2 || severe acute respiratory syndrome coronavirus2` | 11 |
| `neuraminidase_5 || influenza a virus` | 4 |
| `spike protein s1_6 || severe acute respiratory syndrome coronavirus2` | 4 |
| `hiv fusion peptide residue 512-519 || human immunodeficiency virus 1` | 3 |

## Full-Ensemble Readout Result

Metrics use self-match exclusion and only count queries with at least one other
same-label target in the candidate bank.

| readout | Recall@10 | MRR | median first positive rank |
| --- | ---: | ---: | ---: |
| all CDRs mean+std | 0.622 | 0.356 | 5.0 |
| L-CDR3 | 0.595 | 0.346 | 5.0 |
| L-CDR1 | 0.541 | 0.308 | 7.5 |
| H-CDR3 | 0.541 | 0.248 | 9.0 |
| all CDRs H/L mean | 0.541 | 0.242 | 8.0 |
| global H/L mean | 0.446 | 0.244 | 12.5 |
| framework H/L mean | 0.446 | 0.207 | 13.5 |
| random H-framework windows, mean | 0.372 | 0.184 | 18.8 |

Shuffled-label controls were lower than the true-label CDR-aware result:

| readout | shuffled Recall@10 mean | shuffled MRR mean |
| --- | ---: | ---: |
| global H/L mean | 0.314 | 0.130 |
| H-CDR3 | 0.342 | 0.132 |
| all CDRs mean+std | 0.326 | 0.129 |

## Subset And Compression Check

The strongest readout was `all_cdrs_mean_std`. Under that readout:

| condition | Recall@10 | MRR | median first positive rank | cosine to full |
| --- | ---: | ---: | ---: | ---: |
| full 100 conformers | 0.622 | 0.356 | 5.0 | 1.000 |
| K64 H-CDR3 k-center | 0.622 | 0.335 | 6.5 | 1.000 |
| K64 all-CDR k-center | 0.608 | 0.343 | 6.0 | 1.000 |
| random K64, mean | 0.646 | 0.346 | 5.1 | not applicable |

The global endpoint still showed the expected weak one-conformer baseline:

| condition | Recall@10 | MRR | median first positive rank |
| --- | ---: | ---: | ---: |
| global H/L, full 100 conformers | 0.446 | 0.244 | 12.5 |
| global H/L, first conformer only | 0.311 | 0.106 | 21.0 |

## Interpretation

Observed:

- CDR/paratope-aware readouts again outperform global H/L mean pooling on a
  stricter Tier B endpoint.
- `all_cdrs_mean_std` is the best full-ensemble readout in this run, followed
  closely by L-CDR3.
- Framework-only pooling and random framework-window controls remain weaker.
- Shuffled-label controls drop toward the low-MRR baseline expected from this
  label distribution.
- Structural k-center K64 does not clearly beat random K64 under the strongest
  readout. Random K64 remains competitive and can be slightly better on
  Recall@10 in this endpoint.

Interpretation:

- This strengthens the readout-bottleneck hypothesis: useful antigen-neighbor
  information is more visible when MCA tensors are read out through CDR/paratope
  regions than through global H/L pooling.
- The main current claim should be "CDR/paratope-aware readout is better than
  global pooling," not "structural k-center selection has solved conformer
  compression."
- The selector result points toward adaptive or task-aware selectors, learned
  readouts, or better pooling objectives as the next technical direction.

## Caveats

- Tier B is antibody-sequence-de-overlapped against the checkpoint universe,
  but antigen/source de-overlap is not enforced in these Step 4 metrics.
- Tier C exact antigen/source-pair de-overlap would keep 87 of the 160 fixed-100
  Step 4 targets, but only 17 retained targets have a same-label neighbor under
  the current exact-label retrieval setup. See
  `notes/aim1-phase1-5-tierc-antigen-source-deoverlap-count.md`.
- The endpoint has 160 usable targets from the 172-target `<0.85` Tier B pass
  manifest because 8 pass targets had fewer than 100 selected conformers and 4
  pass targets had 100 conformers but failed the CDR chain-assignment gate.
- The candidate bank is built from the same Tier B endpoint with self-matches
  excluded. It is cleaner than the internal 149-target diagnostic, but still
  not a fully independent external benchmark.
- Label balance is uneven; the largest SARS-CoV-2 spike labels dominate the
  available same-label neighbor structure.

## Gaeun-Facing Takeaway

On the stricter Tier B Gaeun-generated endpoint, CDR/paratope-aware pooling
again retrieves same-antigen neighbors better than global H/L mean pooling. The
effect survives framework and shuffled-label controls, so it looks like a real
readout issue rather than a generic local-window artifact. However, K64
structural k-center is not yet clearly better than random K64, so the next
optimization target may be the readout/selector objective rather than only the
number of generated conformers.

## Next Steps

- Make Tier B Step 4 figures from the new result tables.
- Decide whether to rerun/re-score the sparse Tier C endpoint or change the
  candidate bank/grouping before presenting this as cleaner validation.
- If optimizing compression, test task-aware/adaptive subset selectors rather
  than relying only on geometry-only k-center.
