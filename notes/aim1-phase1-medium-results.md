# Aim 1 Phase 1 Medium Structural Results

Date: 2026-06-26

This note records the first CPU-only structural saturation run on Jerry-owned
copies of Gaeun's existing PH/AF3 conformer outputs.

## Inputs

Dataset:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626
```

Key input manifests:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626/manifests/medium_targets.tsv
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626/manifests/medium_conformers.tsv
```

Analysis script:

```text
/external/liulab/jg1920/conffusion/scripts/aim1_phase1_structural.py
```

Tracked source:

```text
scripts/aim1_phase1_structural.py
```

Execution:

- CPU-only.
- Ran on `ragonliu1`.
- Used Jerry env: `/project/liulab/jg1920/envs/jerryenv1/bin/python`.
- No GPU required.
- Wrote outputs only under Jerry-owned `/external/liulab/jg1920/...` paths.

## Outputs

Pilot smoke output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_pilot_smoke_20260626
```

Medium output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626
```

Key medium files:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/aim1_phase1_structural_report.md
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/run_summary.json
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/analysis_catalog.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/target_parse_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/target_diversity_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/saturation_curves.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/selection_strategy_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/parse_failures.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/figures/saturation_median_mean_nearest_rmsd.png
/external/liulab/jg1920/conffusion/aim1_phase1_medium_20260626/figures/target_mean_pairwise_rmsd_hist.png
```

## Parse Results

Pilot smoke:

- 20 targets analyzed.
- 2,560 conformers in manifest.
- 2,560 parse-ok conformers.
- 0 parse failures.
- 20 targets with distance matrices.

Medium set:

- 149 targets analyzed.
- 19,072 conformers in manifest.
- 19,072 parse-ok conformers.
- 0 parse failures.
- 149 targets with distance matrices.

This means the copied CIFs are usable for structural-only Aim 1 Phase 1
experiments.

## Structural Diversity Snapshot

Medium target-level pairwise CA RMSD across common parsed CA atoms:

- Target mean pairwise RMSD range: 8.970 to 19.485 A.
- Median target mean pairwise RMSD: 15.839 A.
- Target median pairwise RMSD range: 8.893 to 18.980 A.
- Median target max pairwise RMSD: 24.011 A.
- Maximum observed target max pairwise RMSD: 44.063 A.

Interpretation:

- The copied ensembles are not trivially identical under this global structural
  readout.
- The first-pass RMSD is global/common-CA over all common chains, so it may
  reflect large complex-level rearrangements or pseudo-binder motion, not only
  antibody CDR movement.
- These numbers should be treated as a screening readout, not as final antibody
  conformational biology.

## Saturation Snapshot

Metric: median across targets of mean nearest-neighbor RMSD from every conformer
to the selected subset. Lower is better.

Medium set:

| Strategy | K=8 | K=16 | K=32 | K=64 |
|---|---:|---:|---:|---:|
| evenly spaced | 11.543 | 9.817 | 7.522 | 4.133 |
| random | 11.555 | 9.993 | 7.763 | 4.636 |
| greedy k-center | 12.092 | 10.321 | 7.892 | 4.113 |
| first K | 12.274 | 11.009 | 9.050 | 5.650 |

Observed:

- The full 128 conformers are substantially more diverse than very small K
  subsets under this global RMSD readout.
- `first K` is consistently worse than evenly spaced/random/k-center at the same
  K, so naive truncation is a poor baseline.
- Even at K=64, the selected subset still leaves roughly 4 A median
  mean-nearest global RMSD to the full ensemble.

## Current Scientific Status

Aim 1 Phase 1 has passed the data-readiness gate:

- A medium ensemble dataset exists in Jerry-owned storage.
- The CIFs parse cleanly.
- The first structural saturation tables and figures exist.
- The run is reproducible from a tracked script.

Aim 1 Phase 1 has not yet answered the final thesis question:

- No CDR-H3-specific analysis yet.
- No VH/VL-only or antibody-only RMSD yet.
- No MCA/ConFormer embedding readout yet.
- No strict-300 disease/retrieval endpoint yet.

## Next Technical Step

Add antibody-focused structural readouts:

- identify chain roles more carefully;
- compute antibody-chain-only RMSD;
- compute VH/VL orientation features;
- add CDR-H3 after numbering or residue mapping is validated.

This should happen before making strong claims about how many conformers are
"enough" for antibody representation.
