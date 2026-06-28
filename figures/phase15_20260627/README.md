# Phase 1.5 Figure Pack

Generated from Jerry-owned Phase 1.5 outputs.

## Inputs

- Phase 1.5b: `/Users/jerrygao/Research/liu-lab/conffusion/figure_inputs/phase15_20260627/phase15b`
- Phase 1.5c: `/Users/jerrygao/Research/liu-lab/conffusion/figure_inputs/phase15_20260627/phase15c`
- PCA points: `/Users/jerrygao/Research/liu-lab/conffusion/figure_inputs/phase15_20260627/pca_readout_points.tsv`

## Figures

1. `01_readout_condition_heatmap_mrr.png` / `.svg` - MRR summarized across readouts and conformer subset conditions.
2. `02_global_endpoint_condition_dotplot.png` / `.svg` - Full ensemble, subset selectors, first-K controls, and random-control means under global H/L pooling.
3. `03_cdr_readout_delta_vs_global.png` / `.svg` - Full-ensemble readout gains relative to global H/L mean pooling.
4. `04_first_positive_rank_distribution.png` / `.svg` - Distribution of first same-antigen neighbor ranks for full-ensemble readouts.
5. `05_target_condition_failure_heatmap.png` / `.svg` - Target-level reciprocal-rank patterns across global and CDR-aware readouts.
6. `06_selector_vs_random_paired_scatter.png` / `.svg` - Per-target selector reciprocal rank compared with same-budget random-control mean.
7. `07_cosine_preservation_vs_rr_change.png` / `.svg` - Near-identical global vectors can still produce different nearest-neighbor ranks.
8. `08_pca_readout_embeddings.png` / `.svg` - PCA projections of full-ensemble global, H-CDR3, and CDR mean+std vectors.
9. `09_rank_transition_global_to_hcdr3.png` / `.svg` - Target-level rank changes when switching from global H/L pooling to H-CDR3 pooling.
10. `10_pipeline_result_schematic.png` / `.svg` - Minimal schematic tying conformer ensembles, MCA tensors, pooling, and retrieval MRR.

## Caveat

These are internal diagnostic figures from the 149-target Phase 1.5 set. They are not independent external validation figures.
