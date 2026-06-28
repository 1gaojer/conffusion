# Aim 1 Phase 1.5 Retrieval Preservation Results

Date: 2026-06-27

## Scope

Phase 1.5 tested whether the copied PH/AF3 conformer ensembles preserve useful
nearest-neighbor antigen retrieval signal under the Gaeun-confirmed June 19
`1000 conformer` MCA checkpoint.

This is an internal diagnostic, not an independent validation claim. The ready
149-target PH/AF3 set substantially overlaps existing reference material, and
the Phase 1.5c candidate bank is built from the same 149 targets being queried
with self-matches excluded.

## Artifacts

Local scripts:

- `scripts/phase15a/phase15a_endpoint_audit.py`
- `scripts/phase15a/phase15b_failure_analysis.py`
- `scripts/phase15a/phase15c_region_readout.py`

Cluster outputs:

- Phase 1.5a endpoint audit:
  `/project/liulab/jg1920/conffusion/phase15a_20260627_175558_endpoint_audit`
- Phase 1.5b failure/random-control analysis:
  `/project/liulab/jg1920/conffusion/phase15b_20260627_failure_analysis`
- Phase 1.5c CDR/paratope-aware readout diagnostic:
  `/project/liulab/jg1920/conffusion/phase15c_20260627_region_readout`
- Phase 1.5c five-target preflight:
  `/project/liulab/jg1920/conffusion/phase15c_20260627_region_readout_preflight`

Key report files:

- `phase15b_failure_analysis_report.md`
- `phase15c_region_readout_report.md`
- `condition_summary.tsv`
- `selector_vs_random_summary.tsv`
- `region_readout_condition_summary.tsv`
- `region_readout_delta_vs_global.tsv`
- `region_readout_random_rollup.tsv`

## Phase 1.5a Endpoint Audit

Observed:

- Confirmed checkpoint path:
  `/external/liulab/gkim/antigen_prediction/2026.06.19_retrain_mca_1000_confs/train_v1/checkpoint.pt`
- Staged checkpoint copy:
  `/project/liulab/jg1920/conffusion/phase14_20260626_2301/shared/checkpoints/mca1000_checkpoint.pt`
- SHA256:
  `a54e938e55de2c2e2253b0f68c8b4fca6fc6090c79a335f7c95f0cccf11f6d74`
- Hash match between confirmed checkpoint and staged copy: yes.
- Phase 1.4 embeddings can be reused for this diagnostic.
- Reference bank:
  `/external/liulab/gkim/antigen_prediction/2026.06.24_hiv_flu_covid_nn_retrieval/outputs_1000/reference_embeddings.npz`
- Reference rows: 931.
- June 24 query count: 278.

Overlap:

| set | count |
| --- | ---: |
| Phase 1.4 target count | 149 |
| Phase 1.4 overlap with reference bank | 144 |
| Phase 1.4 overlap with June 24 queries | 0 |
| Strict-300 count | 300 |
| Strict-300 overlap with June 24 queries | 263 |
| Strict-300 PDB overlap with reference bank | 8 |

Interpretation:

- The checkpoint issue is resolved for this diagnostic.
- The 149-target set is useful for endpoint plumbing and internal readout
  testing.
- The strict-300 set is the cleaner candidate for a future non-overlap-aware
  endpoint once enough conformers/embeddings exist.

## Phase 1.5b Global Endpoint Failure Analysis

This analysis reused the Phase 1.5a global H/L readout and condensed random
controls across five random replicates.

Main metrics used 102 evaluable queries after positive-label and self-match
filtering.

| condition | Recall@10 | MRR | median first positive rank |
| --- | ---: | ---: | ---: |
| full 128 | 0.539 | 0.304 | 9.0 |
| single first conformer | 0.206 | 0.070 | 93.5 |
| H-CDR3 k-center K32 | 0.431 | 0.182 | 20.5 |
| H-CDR3 k-center K64 | 0.461 | 0.231 | 13.5 |
| all-CDR k-center K32 | 0.343 | 0.139 | 21.0 |
| all-CDR k-center K64 | 0.500 | 0.227 | 10.5 |
| random K32 mean | 0.378 | 0.184 | 19.0 |
| random K64 mean | 0.488 | 0.253 | 10.9 |

Selector-vs-random summary:

| selector | K | mean RR minus random | Recall@10 minus random |
| --- | ---: | ---: | ---: |
| H-CDR3 k-center | 32 | -0.003 | +0.053 |
| H-CDR3 k-center | 64 | -0.022 | -0.027 |
| all-CDR k-center | 32 | -0.045 | -0.035 |
| all-CDR k-center | 64 | -0.027 | +0.012 |

Interpretation:

- Full ensemble retrieval is much better than using one conformer.
- K=64 generally beats K=32.
- Under the current global mean H/L endpoint, CDR structural k-center does not
  cleanly beat random K64.
- The global vectors are almost unchanged by subset selection
  (`cosine_to_full` near 0.999 for K64), but retrieval ranks still move. Small
  embedding changes can matter near nearest-neighbor decision boundaries.
- This supports a readout/selector problem more than a simple "CDR k-center is
  enough" story.

## Phase 1.5c CDR/Paratope-Aware Readout Diagnostic

This analysis kept the same per-conformer MCA tensors but changed the pooling
readout. It built an internal 149-target candidate bank using full-ensemble
vectors, then queried it with full/subset/single/random vectors under several
region-aware readouts.

Metrics used 78 evaluable queries after self-match exclusion and positive-label
filtering.

Full-ensemble readout ranking:

| readout | Recall@10 | MRR | median first positive rank |
| --- | ---: | ---: | ---: |
| H-CDR3 only | 0.872 | 0.709 | 1.0 |
| all CDRs mean+std | 0.897 | 0.642 | 1.5 |
| all CDRs H/L mean | 0.808 | 0.614 | 2.0 |
| global + H-CDR3 | 0.859 | 0.596 | 2.0 |
| global + all CDRs | 0.833 | 0.493 | 2.5 |
| global H/L mean | 0.718 | 0.403 | 4.0 |

K64 CDR-aware subset examples:

| readout | subset | Recall@10 | MRR | median first positive rank |
| --- | --- | ---: | ---: | ---: |
| global H/L mean | H-CDR3 k-center K64 | 0.679 | 0.367 | 6.0 |
| global H/L mean | all-CDR k-center K64 | 0.705 | 0.399 | 4.0 |
| H-CDR3 only | H-CDR3 k-center K64 | 0.872 | 0.700 | 1.0 |
| H-CDR3 only | all-CDR k-center K64 | 0.872 | 0.692 | 1.0 |
| all CDRs H/L mean | all-CDR k-center K64 | 0.833 | 0.621 | 1.5 |
| all CDRs mean+std | H-CDR3 k-center K64 | 0.885 | 0.593 | 2.0 |
| all CDRs mean+std | all-CDR k-center K64 | 0.808 | 0.615 | 2.0 |

Random K64 controls under region-aware readouts:

| readout | random K64 Recall@10 mean | random K64 MRR mean |
| --- | ---: | ---: |
| global H/L mean | 0.682 | 0.373 |
| H-CDR3 only | 0.859 | 0.695 |
| all CDRs H/L mean | 0.833 | 0.570 |
| global + H-CDR3 | 0.808 | 0.571 |
| global + all CDRs | 0.782 | 0.504 |
| all CDRs mean+std | 0.841 | 0.672 |

Interpretation:

- CDR-aware readouts recover much stronger antigen-neighbor structure than the
  current global H/L mean readout on this internal diagnostic.
- H-CDR3-only pooling is especially strong in MRR, suggesting that task-relevant
  signal is present in local residue embeddings but diluted by global pooling.
- All-CDR mean+std gives the best Recall@10 in the full-ensemble internal bank,
  suggesting that conformer/region variability features may be useful.
- Random K64 remains competitive with structural CDR k-center. This means the
  readout bottleneck looks clearer than the selector-quality result.
- The result does not prove that the model understands diversity perfectly, nor
  that the generated conformers are physically complete. It shows that the MCA
  tensors contain CDR-local signal that the default endpoint is not fully using.

## Scientific Takeaway

The strongest current thesis-relevant statement is:

> Gaeun's conformer ensembles contain task-relevant CDR-local representation
> signal, but the standard global mean endpoint appears to dilute that signal.
> Before optimizing generation with a new model, we should define or learn a
> CDR/paratope-aware ensemble readout and then test whether smaller conformer
> subsets preserve that readout under leakage-aware evaluation.

Do not claim yet:

- that CDR k-center beats random selection;
- that Phase 1.5 is an independent validation benchmark;
- that the conformer generator itself should be replaced by diffusion;
- that MCA is blind to conformer diversity.

Reasonable claim:

- the bottleneck is plausibly in representation pooling/readout and selector
  alignment, not merely in conformer generation volume.

## Recommended Next Step

Phase 1.5d should turn this into a cleaner decision artifact:

- make one compact table/figure comparing global H/L, H-CDR3, all-CDR, and
  all-CDR mean+std readouts;
- add PCA/UMAP plots for the region-aware vectors if useful for Gaeun's current
  framing;
- repeat the CDR-aware readout on the strict-300 or another de-overlapped set
  once conformers/embeddings are available;
- test task-aware/adaptive subset selectors, because current structural k-center
  does not yet beat random robustly;
- ask Gaeun whether the region-aware retrieval readout matches how she wants to
  frame nearest-neighbor antigen retrieval.

## Figure Pack

Generated 2026-06-27:

```text
figures/phase15_20260627/
```

The folder contains PNG and SVG versions of ten data-focused figures plus
`figure_manifest.tsv`, `figure_inputs.json`, `README.md`, and a
`contact_sheet.png` for quick review.

Figures:

- `01_readout_condition_heatmap_mrr`
- `02_global_endpoint_condition_dotplot`
- `03_cdr_readout_delta_vs_global`
- `04_first_positive_rank_distribution`
- `05_target_condition_failure_heatmap`
- `06_selector_vs_random_paired_scatter`
- `07_cosine_preservation_vs_rr_change`
- `08_pca_readout_embeddings`
- `09_rank_transition_global_to_hcdr3`
- `10_pipeline_result_schematic`

Caveat: these figures visualize the internal Phase 1.5 diagnostic. They should
be shown as readout/endpoint evidence, not as independent external validation.
