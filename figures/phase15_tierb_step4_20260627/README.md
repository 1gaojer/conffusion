# Tier B Step 4 Figure Pack

Generated from Jerry-owned Tier B Step 4 outputs.

## Inputs

- Tier B Step 4 tables: `/Users/jerrygao/Research/liu-lab/conffusion/figure_inputs/phase15_tierb_step4_20260627`
- Internal comparison controls: `/Users/jerrygao/Research/liu-lab/conffusion/figure_inputs/phase15d_20260627`
- PCA points: `/Users/jerrygao/Research/liu-lab/conffusion/figure_inputs/phase15_tierb_step4_20260627/pca_readout_points.tsv`

## Figures

1. `01_readout_performance_ladder.png` / `.svg` - Global, framework, CDR-aware, random-window, and shuffled-label readouts ranked by MRR and Recall@10.
2. `02_global_to_cdr_rank_shift.png` / `.svg` - Target-level first-positive ranks when switching from global pooling to all-CDR mean+std pooling.
3. `03_first_correct_neighbor_cdf.png` / `.svg` - Cumulative distribution of first same-antigen neighbor rank for global and CDR-aware readouts.
4. `04_internal_vs_tierb_mrr_heatmap.png` / `.svg` - MRR pattern across readouts in the original internal diagnostic and stricter Tier B endpoint.
5. `05_controls_cdr_vs_framework_shuffle.png` / `.svg` - Focused control comparison showing CDR readouts against framework, random-window, and shuffled-label baselines.
6. `06_compression_honesty_k64.png` / `.svg` - Full 100 conformers versus K64 structural selectors and same-budget random K64 under the strongest readout.
7. `07_endpoint_funnel.png` / `.svg` - How the Gaeun-generated pool narrows into the Tier B Step 4 evaluable retrieval endpoint.
8. `08_per_antigen_readout_heatmap.png` / `.svg` - MRR by antigen label for global pooling, all-CDR mean+std, and L-CDR3.
9. `09_pca_global_vs_cdr_readouts.png` / `.svg` - PCA of actual full-ensemble global and all-CDR mean+std readout vectors, colored by antigen label.
10. `10_readout_bottleneck_schematic.png` / `.svg` - Minimal schematic linking conformer ensembles, MCA tensors, pooling choice, and Tier B retrieval MRR.

## Caveat

Tier B is exact-PDB-unseen and antibody-sequence-de-overlapped against the checkpoint universe, but antigen/source de-overlap is not yet enforced. Treat this as cleaner readout evidence, not final external validation.
