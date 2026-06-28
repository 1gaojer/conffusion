# Aim 1 Phase 1.4 Embedding Sensitivity Results

Date: 2026-06-27

This note records the full MCA/ConFormer embedding-sensitivity run on the copied
Gaeun PH/AF3 medium ensemble. The purpose was to test whether the CDR structural
diversity measured in Phase 1.3 and the structural coreset run is visible in the
current frozen MCA representation.

## Inputs And Execution

Embedding analysis output:

```text
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full
```

Embedding source run:

```text
/project/liulab/jg1920/conffusion/phase14_20260626_2301/runs/full
```

Dataset copy:

```text
/project/liulab/jg1920/conffusion/phase14_20260626_2301/shared/gaeun_ph_af3_medium_20260626
```

Embedding directory:

```text
/project/liulab/jg1920/conffusion/phase14_20260626_2301/runs/full/embeddings/mca1000_selected
```

Checkpoint:

```text
/project/liulab/jg1920/conffusion/phase14_20260626_2301/shared/checkpoints/mca1000_checkpoint.pt
```

Job:

```text
32100 phase14_analysis_full
```

Execution details:

- CPU-only analysis job on `lambda-scalar`.
- Runtime: 2:21:45.
- Exit code: 0:0.
- Completed at 2026-06-27 09:43:11 EDT.
- Stderr only contained a Matplotlib deprecation warning.
- Targets analyzed: 149/149.
- Targets with H-CDR3 structural matrix: 149/149.
- Targets with all-CDR structural matrix: 149/149.
- Input validation passed with no missing embedding shards and no manifest
  count mismatches.

Key files:

```text
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/run_summary.json
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/input_validation.json
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/aim1_phase1_4_embedding_sensitivity_report.md
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/cdr_to_embedding_correlations.tsv
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/correlation_summary.tsv
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/subset_strategy_summary.tsv
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/subset_embedding_preservation.tsv
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/order_duplication_controls.tsv
/project/liulab/jg1920/conffusion/phase14_20260627_analysis/full/figures/
```

## Main Result

The current frozen MCA embeddings show, at most, a weak local H-CDR3 sensitivity
signal. The full heavy/light pooled embedding is essentially flat with respect
to H-CDR3 structural distance.

Primary H-CDR3 structural distance versus cosine embedding distance:

| Embedding mode | H-CDR3 stratum | Targets | Median Spearman r | Fraction positive |
|---|---|---:|---:|---:|
| H-CDR3 residues | low | 50 | 0.0229 | 0.74 |
| H-CDR3 residues | medium | 49 | 0.0386 | 0.76 |
| H-CDR3 residues | high | 50 | 0.0509 | 0.80 |
| all CDR residues | low | 50 | 0.0071 | 0.66 |
| all CDR residues | medium | 49 | 0.0166 | 0.69 |
| all CDR residues | high | 50 | 0.0104 | 0.68 |
| full H/L concat | low | 50 | -0.0015 | 0.42 |
| full H/L concat | medium | 49 | -0.0000 | 0.49 |
| full H/L concat | high | 50 | 0.0025 | 0.54 |

Across all 149 targets for H-CDR3 structural distance versus cosine embedding
distance:

| Embedding mode | Median Spearman r | Targets with r > 0.05 | Targets with r > 0.1 | Targets with r < 0 |
|---|---:|---:|---:|---:|
| H-CDR3 residues | 0.0351 | 60/149 | 21/149 | 35/149 |
| all CDR residues | 0.0095 | 13/149 | 2/149 | 48/149 |
| full H/L concat | -0.0001 | 4/149 | 0/149 | 77/149 |
| H global | -0.0009 | 4/149 | 0/149 | 77/149 |
| L global | -0.0044 | 5/149 | 0/149 | 81/149 |

Interpretation:

- The H-CDR3-specific residue embedding contains a small positive signal.
- That signal is weak: even in the high-diversity stratum, the median Spearman
  correlation is only 0.0509.
- The full H/L embedding, which is closer to the representation one would use
  for downstream retrieval, does not meaningfully track H-CDR3 pairwise
  structural distance.

## Subset Preservation

Subset analyses compare smaller selected conformer sets against the full
128-conformer teacher ensemble. Lower cosine distance means the subset's average
embedding is closer to the full ensemble's average embedding.

High-H-CDR3-diversity targets, full H/L mean embedding error:

| K | First | Random | H-CDR3 k-center | all-CDR k-center | Embedding k-center |
|---:|---:|---:|---:|---:|---:|
| 8 | 0.015541 | 0.013282 | 0.014574 | 0.011321 | 0.008930 |
| 16 | 0.006480 | 0.006097 | 0.006166 | 0.005318 | 0.004504 |
| 32 | 0.002645 | 0.002480 | 0.002390 | 0.002415 | 0.002564 |
| 64 | 0.000801 | 0.000774 | 0.000622 | 0.000661 | 0.000967 |

Structural coverage for the same high-diversity targets:

| K | Strategy | H-CDR3 coverage | all-CDR coverage | H/L embedding error |
|---:|---|---:|---:|---:|
| 32 | random | 0.772 A | 0.539 A | 0.002480 |
| 32 | H-CDR3 k-center | 0.695 A | 0.548 A | 0.002390 |
| 32 | all-CDR k-center | 0.697 A | 0.509 A | 0.002415 |
| 64 | random | 0.455 A | 0.329 A | 0.000774 |
| 64 | H-CDR3 k-center | 0.322 A | 0.255 A | 0.000622 |
| 64 | all-CDR k-center | 0.333 A | 0.258 A | 0.000661 |

Interpretation:

- K=32 and K=64 subsets preserve the full H/L mean embedding well.
- H-CDR3 and all-CDR structural k-center selection slightly improve embedding
  preservation at K=32 to K=64 in high-diversity targets.
- The absolute embedding differences are small. This weakens the claim that the
  current global representation strongly rewards CDR-aware structural coverage.
- Embedding k-center is a useful upper-bound/control, but it is circular if the
  endpoint is the same embedding preservation metric.

## Order And Duplication Controls

The mean embedding readout is order-invariant but count/duplicate-sensitive.

Full H/L embedding controls:

| Control | Median cosine distance from full ensemble |
|---|---:|
| Permute conformer order | 0.000000 |
| Repeat first conformer 128 times | 0.105324 |
| Repeat H-CDR3 medoid 128 times | 0.102567 |
| Duplicate first conformer 2x weight | 0.032026 |
| Duplicate first conformer 10x weight | 0.089303 |
| Duplicate first conformer 100x weight | 0.103706 |

H-CDR3 embedding controls:

| Control | Median cosine distance from full ensemble |
|---|---:|
| Permute conformer order | 0.000000 |
| Repeat first conformer 128 times | 0.070444 |
| Repeat H-CDR3 medoid 128 times | 0.061319 |
| Duplicate first conformer 2x weight | 0.019548 |
| Duplicate first conformer 10x weight | 0.058543 |
| Duplicate first conformer 100x weight | 0.069217 |

Interpretation:

- Ordering does not affect the mean embedding readout.
- Repeating one conformer is very different from the full ensemble, so the
  representation is not completely indifferent to ensemble diversity.
- Duplicates or sampling weights can strongly bias the mean embedding. Any
  future subset or generation comparison should control conformer counts and
  duplicate density.

## Scientific Interpretation

Observed:

- Structural CDR diversity is real, especially in H-CDR3.
- The current frozen MCA H-CDR3 residue embeddings weakly reflect that geometry.
- The full H/L pooled embedding barely reflects pairwise H-CDR3 distance.
- CDR-aware structural k-center selection gives strong structural coverage and
  modest embedding-preservation gains at K=32 to K=64.
- One-conformer repeats and duplicate-heavy ensembles move the mean embedding
  much more than K=32 or K=64 subsets.

Interpretation:

- This is a mixed result, not a clean positive result.
- It supports structural CDR-aware compression as a sensible baseline, but it
  does not yet prove that CDR-aware compression preserves task-relevant MCA
  behavior.
- The current global H/L representation may average away or underweight the CDR
  geometric variation that the structural analysis found.
- The local H-CDR3 signal suggests the model is not completely blind to loop
  geometry, but the signal is too weak to claim that generated H-CDR3 diversity
  is strongly used downstream.

Recommended next step:

1. Do not claim downstream sufficiency from this result alone.
2. Use K=32 and K=64 structural coresets as the next controlled subset settings.
3. Test a downstream retrieval or nearest-neighbor endpoint if the antigen bank
   and leakage controls are ready.
4. In parallel, inspect whether a CDR-specific readout, pair representation,
   attention, or cluster-aware pooling is more sensitive than the current full
   H/L mean embedding.
5. Keep duplicate and conformer-count controls in every future comparison.

## Status

This result is ready to discuss as a Phase 1.4 premise check with a conservative
caveat:

> The PH/AF3 ensemble has compressible CDR structural diversity, but the current
> frozen global MCA embedding only weakly sees that diversity. The next decision
> is whether retrieval uses information that this embedding-sensitivity readout
> misses, or whether the model/readout needs to change before conformer-generation
> optimization can be treated as task-preserving.
