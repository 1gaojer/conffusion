# Aim 1 Phase 1.2 Antibody Structural Results

Date: 2026-06-26

This note records the antibody-focused follow-up to the first global structural
saturation run. The purpose was to test whether the Phase 1.1 global diversity
remained when measuring only the antibody heavy/light chains.

## Inputs

Dataset:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626
```

Previous global structural output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626
```

Analysis script:

```text
/external/liulab/jg1920/conffusion/scripts/aim1_phase1_antibody_structural.py
```

Tracked source:

```text
scripts/aim1_phase1_antibody_structural.py
```

Execution:

- CPU-only.
- Ran on `ragonliu1`.
- Used Jerry env: `/project/liulab/jg1920/envs/jerryenv1/bin/python`.
- Wrote outputs only under Jerry-owned `/external/liulab/jg1920/...` paths.

## Outputs

Pilot smoke output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_2_pilot_smoke_20260626
```

Medium output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626
```

Key medium files:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/aim1_phase1_2_antibody_structural_report.md
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/run_summary.json
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/chain_role_assignments.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/target_chain_role_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/antibody_diversity_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/antibody_saturation_curves.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/antibody_selection_strategy_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/global_vs_antibody_rmsd.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/figures/saturation_antibody_ca_median_mean_nearest_rmsd.png
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/figures/saturation_heavy_ca_median_mean_nearest_rmsd.png
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626/figures/saturation_light_ca_median_mean_nearest_rmsd.png
```

## Chain Assignment Results

Pilot smoke:

- 20 targets analyzed.
- 2,560 conformers in manifest.
- 2,560 heavy/light assignments succeeded.
- 0 assignment failures.
- 20 targets with antibody distance matrices.

Medium set:

- 149 targets analyzed.
- 19,072 conformers in manifest.
- 19,072 heavy/light assignments succeeded.
- 0 assignment failures.
- 149 targets with antibody distance matrices.

Observed chain pattern:

- Heavy chain was consistently chain `B`.
- Light chain was consistently chain `C`.
- Non-antibody/context chain was generally chain `A`.
- Chain assignment used exact or near-exact VH/VL sequence matching with a
  threshold of 0.95.

## Main Result

The large global structural diversity from Phase 1.1 mostly collapses when
restricted to antibody heavy/light chains.

Medium set global-vs-antibody comparison:

- Median antibody/global mean pairwise RMSD ratio: 0.034.
- Mean antibody/global mean pairwise RMSD ratio: 0.045.
- Ratio range: 0.015 to 0.179.

Interpretation:

- The Phase 1.1 global RMSD signal was real, but it was mostly not antibody-chain
  motion.
- Much of the global diversity likely comes from non-antibody chain/context
  movement or complex-level rearrangement.
- For conformer-generation optimization, antibody-only metrics are the more
  relevant first structural gate.

## Antibody-Only Diversity Snapshot

Medium target-level mean pairwise CA RMSD:

| Metric | Min | Median | Mean | Max |
|---|---:|---:|---:|---:|
| heavy chain CA | 0.159 | 0.563 | 0.680 | 2.873 |
| light chain CA | 0.091 | 0.174 | 0.241 | 2.443 |
| paired heavy+light CA | 0.199 | 0.525 | 0.645 | 3.249 |

Interpretation:

- Antibody-chain structural variation is present but much smaller than the
  global-complex variation.
- Light chain variation is especially small in this first readout.
- The relevant optimization problem may be closer to detecting rare/high-motion
  antibody targets than compressing a broadly diverse antibody ensemble for
  every target.

## Antibody-Only Saturation Snapshot

Metric: median across targets of mean nearest-neighbor RMSD from every conformer
to the selected subset. Lower is better.

Paired heavy+light CA:

| Strategy | K=8 | K=16 | K=32 | K=64 |
|---|---:|---:|---:|---:|
| evenly spaced | 0.285 | 0.235 | 0.185 | 0.114 |
| random | 0.288 | 0.241 | 0.188 | 0.116 |
| greedy k-center | 0.298 | 0.243 | 0.177 | 0.096 |
| first K | 0.286 | 0.239 | 0.189 | 0.117 |

Observed:

- Unlike the global analysis, the simple strategies are very similar under
  antibody-only RMSD.
- Naive `first K` is no longer strongly worse.
- This suggests that the antibody-chain conformers are much more redundant than
  the full global complex appeared.

## Current Scientific Status

Phase 1.2 changes the interpretation of Phase 1.1:

- There is large global conformer diversity.
- The antibody heavy/light chains are comparatively stable.
- The first optimization target should not be "compress a highly diverse
  antibody ensemble" without further qualification.

The next defensible question is:

> Are the small antibody-chain differences still meaningful to Gaeun's
> representation/retrieval model, or are they below the downstream model's
> sensitivity?

## Recommended Next Step

Do not move directly to Phase 2 coreset/compression yet.

Next run should be one of:

1. CDR-H3-focused analysis, if reliable numbering/residue mapping is available.
2. MCA/ConFormer embedding sensitivity, to test whether the small antibody-only
   structural differences change representations.
3. A target-stratified diagnostic that identifies the minority of targets with
   high antibody-only motion.

For thesis framing, the current result is still useful:

- It shows why global structure metrics can be misleading.
- It gives a concrete reason to design antibody-aware conformer-generation
  optimization metrics.
- It suggests that a large fraction of PH/AF3 compute may be redundant for
  antibody-chain geometry, but this needs embedding/downstream confirmation.
