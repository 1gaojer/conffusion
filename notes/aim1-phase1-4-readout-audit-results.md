# Aim 1 Phase 1.4 Rich MCA Readout Audit Results

Date: 2026-06-27

This note records the completed richer MCA/ConFormer readout audit after the
initial Phase 1.4 embedding-sensitivity result. The question was whether the
weak signal in the first analysis was a true encoder limitation or a readout
limitation caused by averaging CDR and conformer information too aggressively.

## Inputs And Execution

Output directory:

```text
/project/liulab/jg1920/conffusion/phase14_20260627_readout_audit/full
```

Slurm job:

```text
32194 phase14-readout-audit
```

Execution details:

- CPU-only job on `lambda-scalar`.
- Runtime: 5:23:53.
- Completed: 2026-06-27 16:39:41 EDT.
- Targets requested: 149.
- Targets analyzed: 149.
- Targets with H-CDR3 structural matrix: 149.
- Targets with all-CDR structural matrix: 149.
- Input validation passed.

Inputs:

```text
/project/liulab/jg1920/conffusion/phase14_20260626_2301/runs/full
/project/liulab/jg1920/conffusion/phase14_20260626_2301/shared/gaeun_ph_af3_medium_20260626
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627/cdr_diversity_summary.tsv
/project/liulab/jg1920/conffusion/scripts/aim1_phase1_cdr_structural.py
```

Checkpoint used by the embedding shards:

```text
/project/liulab/jg1920/conffusion/phase14_20260626_2301/shared/checkpoints/mca1000_checkpoint.pt
```

Key output files:

```text
input_validation.json
run_summary.json
region_sanity.tsv
target_readout_summary.tsv
readout_correlations.tsv
readout_correlation_summary.tsv
subset_readout_preservation.tsv
subset_readout_summary.tsv
pair_repr_region_summary.tsv
pair_repr_cross_target_correlations.tsv
aim1_phase1_4_readout_audit_report.md
```

Important limitation:

- The current saved shard format contains final `mca_repr` and `pair_repr`
  only.
- It does not contain layer-wise conformer weights or subset-specific
  `pair_repr`.
- Testing those requires a later GPU re-encode.

## What This Added Beyond The First Embedding Analysis

The first embedding-sensitivity analysis mostly tested mean-pooled readouts:

- H global mean;
- L global mean;
- full H/L global mean;
- H-CDR3 mean;
- all-CDR mean.

This audit added richer summaries:

- flattened H-CDR3 residue embeddings;
- flattened all-CDR residue embeddings;
- mean plus standard-deviation readouts;
- pair-representation region summaries;
- subset preservation for readout-specific and structural k-center choices.

The main point was to test whether CDR signal is present before pooling, even
if it disappears after global averaging.

## Main Observed Result

H-CDR3 structural geometry is more visible when the readout keeps residue-level
H-CDR3 positions, but it is still mostly lost in the full pooled H/L embedding.

H-CDR3 structural distance versus readout distance:

| Readout | Distance | Stratum | Targets | Median Spearman r | Mean Spearman r | Fraction positive |
|---|---|---|---:|---:|---:|---:|
| H-CDR3 flattened residues | cosine | low | 50 | 0.160 | 0.175 | 1.00 |
| H-CDR3 flattened residues | cosine | medium | 49 | 0.227 | 0.250 | 1.00 |
| H-CDR3 flattened residues | cosine | high | 50 | 0.202 | 0.226 | 0.98 |
| H-CDR3 flattened residues | euclidean | low | 50 | 0.197 | 0.214 | 1.00 |
| H-CDR3 flattened residues | euclidean | medium | 49 | 0.301 | 0.321 | 1.00 |
| H-CDR3 flattened residues | euclidean | high | 50 | 0.268 | 0.296 | 1.00 |
| all-CDR flattened residues | cosine | low | 50 | 0.074 | 0.082 | 1.00 |
| all-CDR flattened residues | cosine | medium | 49 | 0.102 | 0.114 | 0.96 |
| all-CDR flattened residues | cosine | high | 50 | 0.090 | 0.102 | 0.96 |
| H-CDR3 mean | cosine | low | 50 | 0.023 | 0.028 | 0.74 |
| H-CDR3 mean | cosine | medium | 49 | 0.039 | 0.048 | 0.76 |
| H-CDR3 mean | cosine | high | 50 | 0.051 | 0.048 | 0.80 |
| H-CDR3 mean+std | cosine | low | 50 | 0.037 | 0.053 | 0.84 |
| H-CDR3 mean+std | cosine | medium | 49 | 0.077 | 0.099 | 0.88 |
| H-CDR3 mean+std | cosine | high | 50 | 0.072 | 0.081 | 0.80 |
| full H/L global mean | cosine | low | 50 | -0.001 | -0.005 | 0.42 |
| full H/L global mean | cosine | medium | 49 | -0.000 | 0.003 | 0.49 |
| full H/L global mean | cosine | high | 50 | 0.003 | -0.000 | 0.54 |

Interpretation:

- The model is not completely blind to H-CDR3 geometry.
- The signal is much stronger when the readout preserves H-CDR3 residue
  positions.
- Simple H-CDR3 averaging weakens the signal.
- Full H/L global averaging nearly erases it.
- This points more toward a pooling/readout bottleneck than a simple "encoder
  has zero CDR information" conclusion.

## Other Structural Metrics

The same pattern held for related structural metrics:

- For H-CDR3, H-CDRs, all CDRs, and paired VH/VL structural distances, the
  strongest cosine readout was usually `h_cdr3_flat`.
- Median Spearman r for `h_cdr3_flat` was roughly 0.19 to 0.23 in medium/high
  strata across those structural metrics.
- Light-chain CDR structural distance had weaker readout correlations; the top
  median cosine correlations were around 0.05 to 0.08.

Interpretation:

- H-CDR3 remains the main region where generated structural diversity is
  reflected in the MCA residue embeddings.
- All-CDR and light-CDR readouts add some signal, but the cleanest signal is
  still H-CDR3-specific.

## Subset Preservation

K=32 and K=64 subsets preserve average embeddings well, but embedding
preservation does not strongly favor structural k-center over random selection.

For H-CDR3 flattened readouts, median equal-mean cosine error by stratum median:

| Strategy | K | Median-of-strata embedding error | Median-of-strata H-CDR3 coverage distance |
|---|---:|---:|---:|
| first K | 8 | 0.0341 | 0.695 |
| random | 8 | 0.0295 | 0.651 |
| H-CDR3 k-center | 8 | 0.0470 | 0.694 |
| all-CDR k-center | 8 | 0.0477 | 0.696 |
| first K | 16 | 0.0156 | 0.533 |
| random | 16 | 0.0140 | 0.526 |
| H-CDR3 k-center | 16 | 0.0219 | 0.496 |
| all-CDR k-center | 16 | 0.0213 | 0.509 |
| first K | 32 | 0.0069 | 0.386 |
| random | 32 | 0.0060 | 0.382 |
| H-CDR3 k-center | 32 | 0.0084 | 0.333 |
| all-CDR k-center | 32 | 0.0081 | 0.333 |
| first K | 64 | 0.0021 | 0.230 |
| random | 64 | 0.0020 | 0.223 |
| H-CDR3 k-center | 64 | 0.0025 | 0.154 |
| all-CDR k-center | 64 | 0.0025 | 0.150 |

Interpretation:

- K=32 and K=64 are still reasonable structural compression baselines.
- Structural k-center improves structural coverage distance, especially by
  K=64.
- But random and first-K are often just as good or slightly better for average
  embedding preservation.
- That is exactly the warning: the mean embedding is easy to preserve even when
  structural coverage is not optimized.
- Therefore, embedding preservation alone is not a sufficient endpoint for
  task-preserving conformer compression.

## Pair Representation Checks

The pair-representation region summary produced a surprising cross-target
pattern:

- H `cdr3_cdrs` `block_scalar_std` versus target H-CDR3 diversity:
  Spearman r = -0.911.
- H `cdr3_cdr3` `block_mean_abs` versus target H-CDR3 diversity:
  Spearman r = -0.874.
- H `cdr3_framework` metrics were also strongly negative.
- H framework-framework metrics were mildly positive.
- Light-chain CDR block metrics were weak by comparison.

Sanity checks:

- H-CDR3 structural diversity and H-CDR3 length are strongly correlated in this
  target set: Spearman r = 0.702.
- The strong negative pair-representation correlations remain large after a
  simple partial Spearman check controlling H-CDR3 length, for example:
  `cdr3_cdrs` `block_scalar_std` stays around -0.872.

Interpretation:

- This is interesting, but it is not yet a task claim.
- These are cross-target pair-representation magnitude summaries, not
  within-target conformer perturbation tests.
- The result may reflect sequence, target-family, checkpoint, normalization, or
  block-summary effects.
- It is worth following up, but should not be presented as proof that the model
  uses conformer diversity.

## Scientific Interpretation

Observed:

- The completed audit validated 149/149 targets.
- H-CDR3 flattened residue embeddings track H-CDR3 geometry better than the
  earlier mean-pooled H-CDR3 readout.
- Full H/L global mean embeddings remain essentially flat.
- Structural k-center improves structural coverage, but random/first subsets
  are competitive for preserving average embeddings.
- Pair-representation summaries show strong cross-target signals, but their
  meaning is unresolved.

Interpretation:

- The first embedding-sensitivity result was not just "MCA sees nothing."
- A more precise conclusion is:

```text
Some H-CDR3 structural signal exists in local MCA residue embeddings, but the
standard global mean-like readout mostly averages it away.
```

- This strengthens the hypothesis that the bottleneck may be pooling/readout,
  not conformer generation alone.
- It also weakens any plan that evaluates conformer compression only by global
  mean embedding preservation.

What this supports:

- CDR/paratope-aware pooling or adapters are worth testing.
- H-CDR3/all-CDR structural k-center remains a good structural baseline.
- Retrieval preservation is still the right downstream endpoint.

What this does not support yet:

- It does not prove K=32 or K=64 preserves antigen retrieval.
- It does not prove Gaeun's conformers are task-useful.
- It does not prove the MCA encoder itself is fully diversity-aware.
- It does not justify diffusion/generative distillation yet.

## Recommended Next Step

Treat this as a readout bottleneck result.

Next technical gate:

1. Use K=32 and K=64 H-CDR3/all-CDR structural coresets as controlled subset
   conditions.
2. Test nearest-neighbor retrieval once the antigen bank, checkpoint, and
   leakage controls are frozen.
3. Run a GPU re-encode only if needed to expose subset-specific `pair_repr`,
   layer-wise conformer weights, or attention-pooling outputs.
4. If retrieval is flat but local CDR readouts remain sensitive, prototype a
   CDR/paratope-aware pooling or small adapter on frozen MCA features.

Current phase conclusion:

> Phase 1.4 now supports a sharper thesis direction: the useful question is not
> only "how many conformers can we remove?" but "which readout preserves
> task-relevant CDR ensemble diversity instead of averaging it away?"
