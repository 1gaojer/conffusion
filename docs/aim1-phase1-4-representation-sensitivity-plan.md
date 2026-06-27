# Aim 1 Phase 1.4 Representation-Sensitivity And CDR-Aware Compression Plan

Date: 2026-06-26

Purpose: define the next step after the CDR-focused structural result. Phase
1.3 showed that Gaeun's copied PH/AF3 ensemble is mostly framework-redundant but
contains substantial CDR-H3 diversity. Phase 1.4 should test whether that
CDR-H3 diversity is visible to MCA/ConFormer representations and, later,
whether it matters for antigen retrieval.

## Current Evidence

Phase 1.2 whole-antibody structural output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_2_antibody_structural_20260626
```

Phase 1.3 CDR structural output:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_3_cdr_structural_20260626
```

Tracked result note:

```text
notes/aim1-phase1-3-cdr-structural-results.md
```

Observed across the same 149 targets and 19,072 copied conformers:

| Metric | Median target mean pairwise RMSD |
|---|---:|
| whole paired VH+VL CA | 0.525 A |
| frame-aligned all-CDR CA | 0.948 A |
| frame-aligned H-CDR3 CA | 1.604 A |

Interpretation:

- Whole VH/VL RMSD is too blunt for conformer-generation optimization.
- The ensemble contains antibody-relevant motion, mainly in H-CDR3.
- Structural diversity alone is not enough; the next question is whether the
  intended downstream model sees and uses that diversity.

## Phase 1.4 Question

Primary question:

> Do MCA/ConFormer embeddings respond to H-CDR3 conformational diversity in the
> PH/AF3 teacher ensemble?

Secondary question:

> Can a smaller CDR-aware subset preserve the full-ensemble representation
> better than random or first-K selection?

This phase is still an Aim 1 premise test. It is not yet diffusion modeling and
not yet a final strict-300 retrieval claim.

## Working Hypotheses

### H1: CDR-H3 Geometry Is Representation-Relevant

If conformers that are far apart in framework-aligned H-CDR3 RMSD are also far
apart in MCA/ConFormer embedding space, then CDR-aware conformer selection is
scientifically justified.

Expected next move if supported:

- proceed to Aim 2 CDR-aware coreset/compression;
- compare K=8/16/32/64 subsets against the full teacher ensemble;
- prioritize H-CDR3/all-CDR coverage over whole VH/VL RMSD.

### H0: CDR-H3 Geometry Is Representation-Invisible

If embeddings barely change despite large H-CDR3 structural differences, the
current model may not exploit the generated conformer diversity.

Expected next move if supported:

- diagnose model architecture, pooling, conformer ordering, and training target;
- avoid claiming conformer-generation optimization helps downstream behavior;
- treat structural diversity as a model-design or training-signal problem, not
  immediately as a generation-cost problem.

## Part A: Build The Phase 1.4 Target Set

Use the copied medium PH/AF3 dataset first:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626
```

Rationale:

- already copied into Jerry-owned space;
- 149 targets;
- 128 conformers per target;
- complete CDR structural summaries;
- no need to wait for strict-300 generation.

Create a target-stratification table with:

- low, medium, and high H-CDR3 diversity targets;
- low and high all-CDR diversity targets;
- targets with high H-CDR3 / whole-antibody RMSD ratio;
- top high-motion examples such as `8suo`, `6pzy`, `6ii4`, `8g40`, and `8g3v`;
- representative low-motion controls.

Recommended pilot:

- 10 low-H-CDR3 targets;
- 10 medium-H-CDR3 targets;
- 10 high-H-CDR3 targets;
- all 128 conformers per target.

Deliverable:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_4_target_strata_YYYYMMDD/
```

Expected files:

- `target_strata.tsv`
- `selected_conformers.tsv`
- `selection_summary.md`

Gate:

- Each selected target must have reliable H-CDR3 and all-CDR mapping.
- Record missing CDR residues; do not silently drop them.

## Part B: Identify The Correct Embedding Source

Before running embedding sensitivity, resolve the exact model and extraction
path.

Needed:

- correct MCA/ConFormer checkpoint;
- code path for frozen embedding extraction;
- expected input format for conformer structures;
- whether the model embeds individual conformers or conformer sets;
- whether the model uses order-sensitive pooling, mean pooling, attention, or
  another aggregation rule;
- whether the model expects full complex, antibody-only chains, or processed
  features.

Known project caveat:

- The older antigen-score comparison used an incorrect checkpoint and is
  provisional.
- Gaeun's 2026-06-24 guidance favors nearest-neighbor retrieval over the older
  closed 300-category classifier framing.

Deliverable:

```text
docs/aim1-phase1-4-embedding-source-audit.md
```

Minimum fields:

- checkpoint path;
- checkpoint provenance;
- embedding script path;
- input schema;
- output tensor fields;
- aggregation behavior;
- GPU/CPU requirements;
- exact command used for a small smoke extraction.

Gate:

- Do not make embedding-sensitivity claims until the checkpoint and extraction
  code are verified live.

## Part C: Run CDR-To-Embedding Sensitivity

For each selected target, compute:

- per-conformer embedding;
- pairwise embedding distances;
- pairwise H-CDR3 RMSD;
- pairwise all-CDR RMSD;
- pairwise whole VH/VL RMSD;
- optional VH/VL orientation distance.

Primary analyses:

- Spearman correlation between H-CDR3 RMSD and embedding distance per target;
- pooled correlation across targets, with target fixed effects or target-level
  summaries;
- comparison against whole VH/VL RMSD and light-CDR RMSD;
- PCA/UMAP plots colored by H-CDR3 cluster or H-CDR3 distance to medoid;
- embedding distance between structural clusters;
- embedding stability of low-H-CDR3 versus high-H-CDR3 targets.

Controls:

- one conformer repeated to the same ensemble size;
- first-K conformers;
- random K conformers;
- greedy H-CDR3 k-center K conformers;
- all-CDR k-center K conformers;
- conformer order permutation if the model consumes sets;
- conformer duplication at 2x, 10x, and 100x if the model consumes sets.

Recommended K values:

```text
K = 1, 2, 4, 8, 16, 32, 64, 128
```

Deliverable:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_4_embedding_sensitivity_YYYYMMDD/
```

Expected files:

- `embedding_manifest.tsv`
- `per_conformer_embeddings.*`
- `embedding_pairwise_summary.tsv`
- `cdr_to_embedding_correlations.tsv`
- `subset_embedding_preservation.tsv`
- `order_duplication_controls.tsv`
- `aim1_phase1_4_embedding_sensitivity_report.md`
- `figures/`

Gate:

- If H-CDR3 distance correlates with embedding distance and CDR-aware subsets
  preserve embeddings better than naive controls, proceed to CDR-aware coreset
  experiments.
- If embedding distance is flat with respect to H-CDR3 distance, diagnose the
  model before optimizing generation.

## Part D: Run Structural CDR Coreset Curves

This can proceed in parallel with embedding-source auditing because it is
CPU-only and uses existing CDR distance outputs.

Goal:

> Determine how many conformers are needed to cover H-CDR3 and all-CDR geometry.

Use selection strategies:

- first K;
- random K, repeated;
- evenly spaced K;
- greedy k-center on H-CDR3 distance;
- greedy k-center on all-CDR distance;
- optional cluster medoids.

Metrics:

- H-CDR3 teacher-to-subset mean-nearest RMSD;
- all-CDR teacher-to-subset mean-nearest RMSD;
- maximum uncovered CDR cluster radius;
- target-level saturation curves;
- target-specific recommended K.

Initial interpretation from Phase 1.3:

- K=32 covers much of median H-CDR3 diversity but leaves roughly 0.32-0.38 A
  mean-nearest frame-aligned H-CDR3 distance.
- K=64 improves coverage further, especially with greedy k-center.
- Light-chain CDRs saturate faster and should not drive the primary budget.

Deliverable:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_YYYYMMDD/
```

Expected files:

- `cdr_coreset_targets.tsv`
- `cdr_coreset_selected_conformers.tsv`
- `cdr_coreset_saturation_summary.tsv`
- `target_specific_k_recommendations.tsv`
- `figures/`

Gate:

- This is structural evidence only. Do not claim downstream sufficiency until
  embedding or retrieval preservation is tested.

## Part E: Add Retrieval Only After The Endpoint Is Frozen

The final Aim 1 test remains downstream retrieval, not structural diversity.

Do this after:

- correct checkpoint is resolved;
- embedding extraction is verified;
- antigen bank and hard negatives are defined;
- strict-300 or an approved pilot endpoint set exists;
- leakage/family split is frozen.

Preferred endpoint:

- nearest-neighbor antigen retrieval over MCA/ConFormer pretrain
  representations.

Metrics:

- Recall@1, Recall@5, Recall@10;
- MRR;
- mean average precision;
- macro averages by antigen group;
- performance under target/family split;
- candidate-bank size and hard-negative construction.

Controls:

- sequence-only;
- one static structure;
- one conformer repeated;
- full PH/AF3 ensemble;
- random K;
- H-CDR3 k-center K;
- all-CDR k-center K.

Deliverable:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_retrieval_YYYYMMDD/
```

Gate:

- Only after retrieval results should the project claim that a reduced ensemble
  preserves task-relevant conformer information.

## Part F: Literature And Baseline Positioning

The current literature makes the next step narrower, not broader.

Use these papers as framing anchors:

- ABodyBuilder4-STEROIDS: mandatory baseline/novelty check for generic
  paired-VH/VL antibody ensemble generation.
- ITsFlexible / ALL-conformations: supports CDR-focused flexibility analysis
  and external CDR-mode validation.
- AlphaFlow and AlphaFlow-Lit: support the idea of amortized sequence-conditioned
  protein ensemble generation.
- BioEmu: raises the evidence bar for physical equilibrium-ensemble claims.
- DiffAb, AbDiffuser, IgDiff, RFantibody: show antibody diffusion/design is an
  active crowded space, so a generic antibody generator is not the first
  defensible contribution.
- DeepAIR: supports sequence/structure immune-receptor benchmarking, but
  downstream benchmark claims must be controlled and leakage-aware.

Working conclusion:

> The thesis niche should be CDR-aware task-preserving compression or
> distillation of Gaeun's PH/AF3 pseudo-bound teacher ensemble, not generic
> antibody conformer generation from scratch.

## Decision Tree

### Outcome 1: H-CDR3 Geometry Tracks Embedding Distance

Meaning:

- The model sees at least some of the structural diversity.

Next:

- run CDR-aware K-subset embedding preservation;
- proceed toward Aim 2 coreset/compression;
- prepare strict-300 pilot once conformers exist.

### Outcome 2: H-CDR3 Geometry Does Not Track Embedding Distance

Meaning:

- The ensemble may be structurally diverse but downstream-invisible.

Next:

- inspect model inputs and pooling;
- test per-residue embeddings or attention rather than only pooled embeddings;
- ask whether training objectives encourage conformer sensitivity;
- pause generation-cost optimization claims.

### Outcome 3: Random Small K Preserves Embeddings And Retrieval

Meaning:

- The ensemble is oversampled enough that simple subsampling works.

Next:

- prioritize reduced pipeline settings and adaptive stopping;
- treat sophisticated generators as unnecessary for the core thesis.

### Outcome 4: CDR-Aware Coresets Beat Random K

Meaning:

- Selection matters.

Next:

- build target-specific selector or k-center/cluster-medoids pipeline;
- test prospective reduced PH/AF3 allocation.

### Outcome 5: No Compact Subset Preserves Signal

Meaning:

- Either rare modes matter or the model is sensitive to many small variations.

Next:

- quantify which targets need large K;
- consider adaptive target-specific budgets;
- delay generative distillation unless it can reproduce rare H-CDR3 modes.

## Immediate Work Order

1. Build `target_strata.tsv` from Phase 1.3 CDR summaries.
2. Audit the MCA/ConFormer checkpoint and embedding extraction path.
3. Run a tiny embedding smoke test on 2-3 targets.
4. Run a 30-target low/medium/high H-CDR3 embedding-sensitivity pilot.
5. In parallel, generate CDR-aware coreset manifests for K=8/16/32/64.
6. Use pilot results to decide whether Aim 2 compression is justified.
7. Repeat on strict-300 only after conformers and endpoint metadata are ready.

## Decisions Or Inputs Needed

- Which MCA/ConFormer checkpoint is the correct one for current evaluation?
- Which script should be treated as the canonical embedding extractor?
- Does the model consume individual conformers, sets of conformers, or
  aggregated conformer features?
- Should the first embedding pilot use Gaeun's copied SAbDab-like medium set,
  strict-300 smoke/pilot outputs when available, or both?
- What non-inferiority margin should count as "preserved" for embedding and
  retrieval behavior?

Recommendation:

- Start with the copied 149-target medium set for embedding sensitivity because
  it is ready now.
- Treat strict-300 as the later retrieval benchmark once conformers and endpoint
  metadata are ready.
- Do not start diffusion or generative distillation until embedding sensitivity
  and CDR-aware coreset baselines are known.
