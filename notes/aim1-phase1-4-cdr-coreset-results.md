# Aim 1 Phase 1.4 Structural CDR Coreset Results

Date: 2026-06-27

This note records the structural-only CDR coreset run for the copied Gaeun
PH/AF3 medium ensemble. The purpose was to test how aggressively the full
teacher ensemble can be compressed while still covering H-CDR3 and all-CDR
geometry.

## Inputs And Execution

Dataset:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626
```

Output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627
```

Job:

```text
32190 phase14_cdr_coresets
```

Execution details:

- CPU-only Slurm job on `ragonliu1`.
- Partition: `liulab`.
- Resources: 8 CPUs, 48 GB RAM, no GPU.
- Runtime: 31:57.
- Exit code: 0:0.
- Stderr was empty.
- Used Jerry-owned inputs and outputs only.

The run reused the Phase 1.3 CDR structural script, so some output filenames
still say `aim1_phase1_3`. The result should be interpreted as the Phase 1.4
structural coreset run.

## Coverage

Run summary:

- Targets analyzed: 149.
- Conformers in manifest: 19,072.
- Assignment-ok conformers: 19,072.
- Mapped conformers used: 19,072.
- Targets with CDR distance matrices: 149.
- K values: 1, 2, 4, 8, 16, 32, 64, 128.
- Random seeds: 0, 1, 2, 3, 4.
- Numbering: AbNumber / ANARCI with IMGT scheme.

Key files:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627/run_summary.json
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627/cdr_saturation_curves.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627/cdr_selection_strategy_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627/cdr_diversity_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627/figures/
```

Current implementation note:

- `cdr_saturation_curves.tsv` contains the selected conformer indices in the
  `selected_indices` column.
- A separate `cdr_coreset_selected_conformers.tsv` manifest has not yet been
  materialized. It can be derived from `cdr_saturation_curves.tsv` when needed
  for embedding or retrieval follow-up.

## Main Result

The structural coreset signal is real, but it is strongest at K=32 to K=64.
Tiny subsets, especially K=8, do not reliably preserve H-CDR3 geometry.

Median across targets of H-CDR3 teacher-to-subset mean-nearest frame-aligned CA
RMSD:

| K | Random | H-CDR3 k-center | Interpretation |
|---:|---:|---:|---|
| 8 | 0.640 A | 0.690 A | k-center is not helpful yet |
| 16 | 0.521 A | 0.496 A | small gain |
| 32 | 0.379 A | 0.321 A | meaningful gain |
| 64 | 0.222 A | 0.143 A | strong gain |

For all-CDR geometry:

| K | Random | all-CDR k-center | Interpretation |
|---:|---:|---:|---|
| 8 | 0.465 A | 0.478 A | no gain |
| 16 | 0.382 A | 0.385 A | no gain |
| 32 | 0.292 A | 0.270 A | modest gain |
| 64 | 0.175 A | 0.139 A | useful gain |

## Threshold View

Target counts below structural coverage thresholds for greedy k-center:

| Region | Threshold | K=8 | K=16 | K=32 | K=64 |
|---|---:|---:|---:|---:|---:|
| H-CDR3 | <= 0.25 A | 25/149 | 39/149 | 58/149 | 113/149 |
| H-CDR3 | <= 0.5 A | 62/149 | 75/149 | 108/149 | 133/149 |
| H-CDR3 | <= 1.0 A | 107/149 | 119/149 | 131/149 | 148/149 |
| all-CDRs | <= 0.25 A | 18/149 | 31/149 | 64/149 | 123/149 |
| all-CDRs | <= 0.5 A | 80/149 | 102/149 | 122/149 | 141/149 |
| all-CDRs | <= 1.0 A | 123/149 | 131/149 | 142/149 | 149/149 |

This suggests that K=32 is already useful for many targets, while K=64 is a
more robust default if the goal is broad CDR structural coverage.

## High-Motion Targets Remain Hard

The most flexible H-CDR3 targets still have meaningful residual uncovered
motion even at K=64:

| Target | H-CDR3 pairwise mean | Random K=64 | k-center K=64 | k-center K=64 max nearest |
|---|---:|---:|---:|---:|
| 8suo | 7.374 A | 1.422 A | 1.226 A | 3.116 A |
| 6pzy | 7.357 A | 1.023 A | 0.930 A | 2.328 A |
| 6ii4 | 6.397 A | 1.006 A | 0.787 A | 1.972 A |
| 8g40 | 6.146 A | 1.072 A | 0.754 A | 1.992 A |
| 8g3v | 6.016 A | 0.962 A | 0.777 A | 2.128 A |

Interpretation:

- Most targets appear compressible.
- High-H-CDR3-motion targets likely need adaptive K, not a single fixed small
  K for every antibody.
- A fixed K=64 is much safer than K=32 for these hard cases, but even K=64
  does not perfectly cover the full H-CDR3 teacher support.

## Scientific Interpretation

Observed:

- The full PH/AF3 medium ensemble is structurally redundant at the CDR level.
- H-CDR3 is the main driver of difficult coverage.
- CDR-aware k-center selection becomes meaningfully better than random at
  K=32 and especially K=64.
- K=8 is too aggressive for structural H-CDR3 coverage.

Interpretation:

- This supports a CDR-aware compression thesis more than a tiny random-subset
  thesis.
- A good first practical target is likely K=32 to K=64 conformers per antibody,
  with adaptive K for high-H-CDR3-diversity targets.
- The compression result is structural only. It should not yet be described as
  preserving MCA/ConFormer behavior or antigen retrieval.

Recommended next step:

- Compare these structural coreset selections against the MCA/ConFormer
  embedding sensitivity run.
- If H-CDR3 k-center subsets preserve embeddings better than random or first-K,
  proceed toward Aim 2 CDR-aware coreset/compression.
- If embeddings do not respond to these structural differences, diagnose model
  readout, pooling, checkpoint, and training signal before optimizing
  conformer generation.

## Status

This result is ready to use as structural evidence for Aim 1 Phase 1.4 with the
caveat that downstream representation and retrieval preservation remain
unproven.
