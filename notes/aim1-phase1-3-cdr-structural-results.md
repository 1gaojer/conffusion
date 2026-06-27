# Aim 1 Phase 1.3 CDR Structural Results

Date: 2026-06-26

This note records the CDR-focused follow-up to the Phase 1.2 antibody-chain
analysis. The purpose was to test whether whole-VH/VL RMSD was averaging away
meaningful loop diversity, especially in CDR-H3.

## Inputs

Dataset:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626
```

Previous antibody-chain structural output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626
```

Analysis script:

```text
/external/liulab/jg1920/conffusion/scripts/aim1_phase1_cdr_structural.py
```

Tracked source:

```text
scripts/aim1_phase1_cdr_structural.py
```

Execution:

- CPU-only.
- Ran on `ragonliu1`.
- Used Jerry env: `/project/liulab/jg1920/envs/jerryenv1/bin/python`.
- Wrote outputs only under Jerry-owned `/external/liulab/jg1920/...` paths.
- Used AbNumber IMGT numbering with `use_anarcii=True`.
- Used chain/mapping threshold 0.90 because many AF3 CIF antibody chains are
  truncated by a few residues but still unambiguously match the intended VH/VL
  chains.

## Region Extraction Validation

Before the full run, a validation pass was run on the first 10 targets and two
conformers per target:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_validation_20260626
```

The extracted IMGT CDR sequences looked biologically plausible. Example target
`1a2y`:

| Region | Residues | Sequence |
|---|---:|---|
| H-CDR1 | 8 | GFSLTGYG |
| H-CDR2 | 7 | IWGDGNT |
| H-CDR3 | 10 | ARERDYRLDY |
| L-CDR1 | 6 | GNIHNY |
| L-CDR2 | 3 | YTT |
| L-CDR3 | 9 | QHFWSTPRT |

Some CIFs are missing individual CA residues. For example, `1bvk` mapped
H-CDR3 fully but mapped L-CDR3 as 8/9 residues. The full analysis handles this
by using common mapped residues across conformers and records the missing
residues in `cdr_mapping_summary.tsv`.

## Outputs

Full output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626
```

Key files:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/aim1_phase1_3_cdr_structural_report.md
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/run_summary.json
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/cdr_numbering_assignments.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/cdr_mapping_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/cdr_diversity_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/cdr_saturation_curves.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/cdr_selection_strategy_summary.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/cdr_per_residue_variability.tsv
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/whole_vs_cdr_summary.tsv
```

Figures:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626/figures/
```

## Coverage

Full run:

- 149 targets analyzed.
- 19,072 conformers in manifest.
- 19,072 heavy/light assignments succeeded.
- 19,072 mapped conformers used.
- 149 targets with CDR distance matrices.

Mapping completeness caveat:

- H-CDR3 had complete residue mapping for 12,800/19,072 conformer rows.
- L-CDR3 had complete residue mapping for 5,632/19,072 conformer rows.
- This incompleteness is mostly due to missing CA residues in CIFs, not failed
  numbering.
- Diversity metrics use common mapped residues per target/region.

## Main Result

CDR-H3 diversity is substantially larger than whole-antibody VH/VL RMSD
suggested.

Live file comparison across the same 149 targets:

| Metric | Min | Median | Mean | P90 | Max |
|---|---:|---:|---:|---:|---:|
| whole paired VH+VL CA mean pairwise RMSD | 0.199 | 0.525 | 0.645 | 1.171 | 3.249 |
| frame-aligned H-CDR3 CA mean pairwise RMSD | 0.226 | 1.604 | 1.815 | 4.008 | 7.374 |
| frame-aligned all-CDR CA mean pairwise RMSD | 0.245 | 0.948 | 1.186 | 2.371 | 4.328 |

The median H-CDR3 / whole-antibody mean pairwise RMSD ratio was 2.865.

Interpretation:

- Phase 1.2 was correct that the VH/VL framework is comparatively stable.
- That did not mean the ensemble lacks antibody-relevant motion.
- The diversity is concentrated in CDR-H3 and, to a lesser extent, combined
  heavy-chain CDRs.
- Whole VH/VL RMSD is too blunt for evaluating conformer-generation usefulness.

## Region-Level Diversity

Target-level mean pairwise frame-aligned CA RMSD, summarized across targets:

| Region | Min | Median | Mean | P90 | Max |
|---|---:|---:|---:|---:|---:|
| H-CDR1 | 0.125 | 0.196 | 0.310 | 0.681 | 1.564 |
| H-CDR2 | 0.149 | 0.228 | 0.277 | 0.475 | 1.221 |
| H-CDR3 | 0.226 | 1.604 | 1.815 | 4.008 | 7.374 |
| L-CDR1 | 0.105 | 0.172 | 0.318 | 0.653 | 2.061 |
| L-CDR2 | 0.065 | 0.103 | 0.116 | 0.149 | 0.506 |
| L-CDR3 | 0.133 | 0.233 | 0.377 | 0.871 | 2.577 |
| H-CDRs | 0.203 | 1.122 | 1.342 | 2.931 | 5.499 |
| L-CDRs | 0.126 | 0.236 | 0.365 | 0.748 | 1.953 |
| all-CDRs | 0.245 | 0.948 | 1.186 | 2.371 | 4.328 |

Self-aligned loop-shape RMSD showed the same qualitative pattern:

| Region | Median target mean self-aligned RMSD |
|---|---:|
| H-CDR1 | 0.089 |
| H-CDR2 | 0.075 |
| H-CDR3 | 1.057 |
| L-CDR1 | 0.065 |
| L-CDR2 | 0.018 |
| L-CDR3 | 0.114 |
| H-CDRs | 0.997 |
| L-CDRs | 0.162 |
| all-CDRs | 0.894 |

## Saturation Snapshot

Metric: median across targets of mean-nearest frame-aligned CA distance from
all conformers to a selected subset. Lower is better.

H-CDR3:

| Strategy | K=8 | K=16 | K=32 | K=64 |
|---|---:|---:|---:|---:|
| first | 0.695 | 0.533 | 0.384 | 0.227 |
| evenly spaced | 0.654 | 0.506 | 0.368 | 0.220 |
| random | 0.640 | 0.521 | 0.379 | 0.222 |
| greedy k-center | 0.690 | 0.496 | 0.321 | 0.143 |

All CDRs:

| Strategy | K=8 | K=16 | K=32 | K=64 |
|---|---:|---:|---:|---:|
| first | 0.476 | 0.387 | 0.299 | 0.179 |
| evenly spaced | 0.456 | 0.362 | 0.292 | 0.174 |
| random | 0.465 | 0.382 | 0.292 | 0.175 |
| greedy k-center | 0.478 | 0.385 | 0.270 | 0.139 |

Observed:

- H-CDR3 has meaningful residual diversity even when whole VH/VL looks stable.
- K=32 covers much of median H-CDR3 diversity but leaves roughly 0.32-0.38 A
  mean-nearest distance depending on strategy.
- K=64 improves coverage further, especially with greedy k-center.
- Light-chain CDRs saturate much faster and have much smaller diversity.

## High H-CDR3 Motion Targets

Top targets by H-CDR3 mean pairwise frame-aligned CA RMSD:

| Target | Mean | P90 | Max | Common H-CDR3 CA count |
|---|---:|---:|---:|---:|
| 8suo | 7.374 | 10.777 | 17.274 | 20 |
| 6pzy | 7.357 | 11.940 | 15.682 | 20 |
| 6ii4 | 6.397 | 9.765 | 13.949 | 15 |
| 8g40 | 6.146 | 9.233 | 16.656 | 20 |
| 8g3v | 6.016 | 9.175 | 15.711 | 20 |
| 8g30 | 5.962 | 8.776 | 15.891 | 20 |
| 8ez8 | 5.537 | 8.463 | 12.453 | 19 |
| 8g3z | 5.288 | 7.957 | 12.065 | 21 |
| 8g3q | 5.074 | 7.794 | 12.877 | 21 |
| 8g3r | 4.943 | 7.615 | 12.040 | 21 |

## Scientific Status

This result materially updates Aim 1 Phase 1:

- The earlier whole-antibody result should not be interpreted as "no useful
  antibody conformer diversity."
- The correct structural readout is CDR-aware, especially H-CDR3-aware.
- The conformer ensemble may be redundant at the framework level while still
  carrying task-relevant loop diversity.
- The next question is whether MCA/ConFormer embeddings and antigen retrieval
  are sensitive to this H-CDR3 diversity.

## Recommended Next Step

Proceed to a representation-sensitivity test:

1. Select targets spanning low, medium, and high H-CDR3 diversity.
2. Embed conformers with the relevant MCA/ConFormer checkpoint.
3. Test whether conformers that are far apart in H-CDR3 geometry are also far
   apart in embedding space.
4. If yes, run coreset/subsampling on H-CDR3-aware distances and compare
   embedding/retrieval preservation.
5. If no, the downstream model may not currently use the loop diversity, which
   changes the optimization target.

Do not use whole VH/VL RMSD alone as the conformer-generation compression
metric.
