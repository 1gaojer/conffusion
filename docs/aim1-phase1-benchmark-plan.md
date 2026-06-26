# Aim 1 Phase 1 Benchmark And Data Inventory Plan

Date: 2026-06-26

Purpose: define the first practical benchmark for Aim 1 before doing any
optimization, coreset selection, or generative modeling. Phase 1 should answer:
what data do we already have, what is missing, and what needs to be built so we
can test whether conformer ensembles actually matter for Gaeun's downstream
antigen-specificity/retrieval work.

## Phase 1 Question

Aim 1 is not "train a diffusion model" yet. Aim 1 starts with a stricter
premise test:

> Does a PH/AF3 conformer ensemble provide useful signal beyond sequence-only,
> single-structure, and repeated-single-conformer controls on a leakage-aware
> antigen/disease endpoint?

If this answer is weak or confounded, optimizing conformer generation is not yet
scientifically grounded. If the answer is positive, later phases can ask how
much of the ensemble can be removed or amortized.

## Recommended Benchmark Unit

Use one paired heavy/light antibody target as the benchmark unit.

- One row = one antibody target with VH sequence, VL sequence, metadata, and
  endpoint label or retrieval truth.
- Conformers are nested under that target.
- Splits must happen at target or family level, not conformer level.
- No conformer from the same target should cross train/test boundaries.

This keeps the benchmark aligned with the real question: whether extra
conformers help represent an antibody target better, not whether a model can
memorize near-duplicate conformers from the same target.

## Observed Cluster Evidence

All paths below were inspected read-only.

### 1. Primary Disease-Labeled Target Manifests

Current strict 300-set:

```text
/project/liulab/jg1920/bcr-conformer-runs/model3_iedb_latest_pipeline_inputs_20260624/recommended_strict_300/input/targets.tsv
/project/liulab/jg1920/bcr-conformer-runs/model3_iedb_latest_pipeline_inputs_20260624/recommended_strict_300/input/targets.json
```

Observed:

- 300 target rows.
- Balanced disease labels: 100 Covid19, 100 HIV, 100 InfluenzaA.
- All 300 have paired heavy and light sequences in the TSV.
- All 300 have heavy/light CDR3 and V-gene fields.
- Source databases: CoV-AbDab 100, IEDB 86, PLAbDab 57, LANL CATNAP 55,
  SAbDab/PDB 2.
- 263 rows have existing `selected_structure_path` values that currently exist.
- 37 replacement rows are missing `selected_structure_path`.
- Pipeline JSON contains only `name`, `vh_seq`, and `vl_seq`, which is enough
  for generation but not enough for endpoint analysis by itself.

Interpretation:

- This should be the main Aim 1 Phase 1 disease/retrieval benchmark spine.
- The TSV should be the metadata source of truth.
- The JSON should be treated as the generation input format.
- The 37 missing static structures are not fatal for generation, but they matter
  for static-structure baselines.

Related manifests:

```text
/project/liulab/jg1920/bcr-conformer-runs/model3_iedb_latest_pipeline_inputs_20260624/realtest_278/input/
/project/liulab/jg1920/bcr-conformer-runs/model3_iedb_latest_pipeline_inputs_20260624/supplemented_300/input/
```

Use these for comparison or fallback, but keep `recommended_strict_300` as the
primary benchmark unless Gaeun changes the target set.

### 2. Existing Jerry-Owned Generated Structure Layer

Catalog:

```text
/project/liulab/jg1920/model3-iedb/results/conformer_benchmark_v1/generated_conformer_catalog_balanced100_30992.tsv
```

Observed:

- 2,133 balanced benchmark rows, 711 per class.
- `generated_conformer_source = generated_af3_vh` for all rows.
- 300 rows have status `ok`; 1,833 rows have status `missing_model_cif`.
- 300 selected structure paths exist.
- The old Protein Hunter output root exists:

```text
/project/liulab/jg1920/malid-compare/results/model3_generated_honly_antigen_balanced100_31036/prothunt_output_honly
```

- That root contains 300 target directories.
- Each directory contains one PDB file.
- 263 of those target directory names match the strict 300 query IDs.

Interpretation:

- This is useful for static/single-structure baselines and pipeline plumbing.
- It is not a true PH/AF3 ensemble for Aim 1.
- It can help build the Phase 1 audit machinery before full strict-300
  ensemble generation finishes.

### 3. Gaeun-Owned Existing PH/AF3 Ensemble Outputs

Gaeun source sequence file:

```text
/external/liulab/gkim/antigen_prediction/datasets/seq_files/full_sabdab_seqsim_split_test_seq_file.json
```

Observed:

- 483 records.
- Record fields include `name`, `vh_seq`, `vl_seq`, and `antigen_size`.

Gaeun conformer results root:

```text
/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results
```

Observed top-level evidence:

- `prothunt_output`: 947 target directories.
- `filter_metrics`: 1,099 target directories.
- `filtered_af3_models`: 770 target directories.
- `af3_results`: very large AF3 result tree.
- `af3_input_jsons`: very large AF3 input/template tree.
- `results_boltz`: 2,085 target directories.

Sample filtered targets show large ensembles. For example, targets such as
`1a14`, `1a2y`, `1a3r`, `1acy`, and `1adq` each had roughly 1,600-2,100 kept
CIFs and 2,500 rejected CIFs in the sampled filtered directories.

Interpretation:

- This is the best existing source for studying PH/AF3 ensemble structure,
  counts, filtering, and failure modes.
- It is Gaeun-owned and must remain read-only.
- It is not automatically the primary disease/retrieval benchmark because it is
  SAbDab/PDB-like and does not directly carry the Covid/HIV/Influenza endpoint
  structure from the strict 300-set.
- It can act as a teacher ensemble inventory or structural sanity set while the
  strict 300 generation is pending.

### 4. Current Staged Generation State

Jerry-owned conformer-generation run roots:

```text
/project/liulab/jg1920/bcr-conformer-runs/smoke_3_targets/input/
/project/liulab/jg1920/bcr-conformer-runs/pilot_30_targets/input/
/project/liulab/jg1920/bcr-conformer-runs/smoke_3_targets_phase1_a2/input/
/project/liulab/jg1920/bcr-conformer-runs/smoke_3_targets_phase1_b1/input/
```

Observed/known from current project state:

- Smoke 3 and pilot 30 inputs exist.
- Phase 1 smoke jobs were queued as separate runs.
- Phase 2 and Phase 3 were not submitted.
- The strict 300 full conformer ensemble has not yet been generated with the
  updated pipeline.

Interpretation:

- Phase 1 benchmark design can proceed now.
- Full ensemble performance claims must wait for generated strict-300 ensemble
  outputs or a smaller reviewed pilot ensemble.

## Data We Can Use Now

Use immediately:

- Strict 300 target TSV as the primary target metadata table.
- Strict 300 JSON as the generation-ready sequence input.
- Existing 300 single-structure/Protein Hunter output folders as a
  single-conformer baseline and plumbing test.
- Gaeun's read-only SAbDab-like ensemble outputs for ensemble-count,
  filter-metric, and file-format inventory.
- Existing Model 3 / Mal-ID / DeepAIR comparison artifacts for endpoint design
  context, not as final conformer-ensemble evidence.
- SAbDab2 strict Tier A candidates as a possible later external
  structure-backed set, after leakage/overlap decisions.

Use cautiously:

- The older 278-query antigen-score/reranker results are provisional because
  they need to be rerun with the correct checkpoint and updated framing.
- The 2,133-row generated AF3 VH catalog is useful for coverage diagnostics but
  has only 300 completed selected structures.
- Gaeun's large ensemble outputs are valuable, but they belong to a different
  source distribution than the strict disease-labeled 300-set.

## Data Still Needed

Before Aim 1 Phase 1 can produce a clean premise result, we still need:

1. A frozen endpoint definition.

   Recommended current framing: nearest-neighbor retrieval over MCA pretrain
   representations, with PCA/UMAP inspection, rather than the older 300-category
   classifier framing.

2. Correct checkpoint and embedding source.

   Existing antigen-score results should not be treated as final until the
   checkpoint issue is resolved and the retrieval-style endpoint is rerun.

3. Full or pilot PH/AF3 ensembles for the strict target set.

   The strict 300 manifest exists, but the updated pipeline has not produced the
   strict 300 ensemble yet. A reviewed smoke or pilot ensemble is enough to
   start validating the benchmark machinery.

4. Static baseline coverage for the 37 missing strict-300 rows.

   Options:

   - accept a 263-row static-structure subset for static-baseline comparisons;
   - generate or select static structures for the 37 replacement rows;
   - use `supplemented_300` or `realtest_278` for a baseline-only sensitivity
     check.

5. Family/split metadata.

   We need target-level and family-level split assignments using features such
   as exact H/L pair, heavy sequence, light sequence, V genes, CDR3 clusters,
   and sequence similarity. Existing Model 3 split-family logic can be reused
   conceptually but should be rebuilt or revalidated for the strict 300 target
   table.

6. Antigen group normalization and hard negatives.

   The disease labels are balanced, but retrieval quality depends on the
   antigen truth table: antigen names, antigen groups, source organisms, and
   candidate-negative construction.

7. A conformer catalog with one row per conformer.

   Needed fields include target ID, conformer ID, source phase, PH branch/cycle,
   AF3 seed/sample, file path, keep/reject status, filter reason, quality
   metrics, and endpoint linkage.

## Things To Build

### 1. Aim 1 Phase 1 Inventory Report

Build a read-only inventory script/report that emits:

- target counts by disease and source database;
- sequence completeness;
- antigen metadata completeness;
- static-structure coverage;
- generated-ensemble coverage;
- conformer counts per target;
- kept/rejected counts per target;
- missing or failed target list;
- exact paths to all source manifests and output roots.

This should not copy Gaeun-owned structures. It should record paths and counts.

### 2. Strict 300 Benchmark Manifest

Build a normalized benchmark table from `targets.tsv` plus generation JSON:

- `benchmark_id`
- `query_id`
- `label_class`
- `source_database`
- `source_record_id`
- `antibody_name`
- `heavy_sequence`
- `light_sequence`
- `heavy_cdr3`
- `light_cdr3`
- `heavy_v_gene`
- `light_v_gene`
- `antigen_name`
- `antigen_region`
- `source_organism`
- `selected_structure_path`
- `static_structure_available`
- `generation_input_name`
- `split_family_id`
- `endpoint_group_id`

This becomes the target-level source of truth for Aim 1 Phase 1.

### 3. Conformer Catalog Builder

Build a cataloger for both Jerry-owned and read-only Gaeun-owned outputs.

Required output:

- one row per conformer file;
- link each conformer back to one benchmark target;
- record whether the conformer came from Gaeun's SAbDab-like teacher inventory,
  strict-300 smoke/pilot/full generation, or older single-structure outputs;
- record keep/reject/filter metrics when available;
- record missing targets explicitly instead of silently dropping them.

### 4. Split Builder And Leakage Validator

Build split assignment and checks:

- target-level split;
- exact H/L pair grouping;
- V-gene holdout or stratification;
- CDR3 cluster holdout;
- VH/VL similarity holdout if feasible;
- source-database balance diagnostics;
- disease-label balance diagnostics;
- no conformer-level leakage.

The validator should fail loudly if conformers from one target appear in more
than one split.

### 5. Endpoint Builder

Build the retrieval/classification endpoint wrapper:

- preferred endpoint: nearest-neighbor retrieval over MCA pretrain
  representations;
- secondary endpoint: disease-panel ranking or macro-F1 only as a diagnostic;
- outputs: Hit@K, MRR, macro-F1 where relevant, class/source breakdowns, and
  ranking-bias diagnostics;
- visual checks: PCA/UMAP embeddings, especially HIV antigen-count/diversity
  confounds.

### 6. Baseline Definitions

Predefine baselines before seeing the final result:

- sequence-only;
- one static/predicted structure;
- one conformer;
- repeated one conformer;
- full PH/AF3 ensemble;
- reduced PH/AF3 subsets;
- random conformer subsets;
- simple structure diversity selections;
- optionally an external antibody ensemble baseline such as ABB4-STEROIDS.

The first baseline to beat is not diffusion. It is a cheap, honest subset.

### 7. Phase 1 Readiness Gate

Before moving to Phase 2 coreset/saturation curves, produce:

- target manifest;
- conformer catalog;
- split report;
- endpoint definition;
- baseline list;
- missing-data report;
- smoke or pilot run summary if full strict-300 generation is not ready.

## Recommended Execution Order

1. Freeze the strict 300 target table as the main Aim 1 Phase 1 target universe.
2. Build the inventory report against strict 300, old 300 single-structure
   outputs, and Gaeun's read-only ensemble tree.
3. Build the target-level benchmark manifest.
4. Build the conformer cataloger using path-only inventory first.
5. Define endpoint groups and retrieval truth.
6. Build split assignments and leakage checks.
7. Validate everything on smoke 3 or pilot 30.
8. Only then run full ensemble-vs-baseline comparisons.

## Initial Go/No-Go Logic

Proceed to Aim 1 Phase 2 only if:

- the endpoint is frozen and reproducible;
- target/family leakage checks pass;
- the full or pilot ensemble can be compared to honest controls;
- the ensemble shows some useful signal beyond sequence-only and
  single-conformer controls.

Pause or revise if:

- retrieval results are dominated by source/database shortcuts;
- static structures perform the same as the full ensemble;
- the conformer catalog reveals too many missing or failed targets;
- strict split performance collapses enough that the endpoint is not yet
  interpretable;
- the current checkpoint/embedding source remains unresolved.

## Practical Conclusion

Phase 1 is mostly benchmark engineering and data archaeology. The important
deliverable is not a model. It is a clean, leakage-aware answer to:

> For these antibody targets, do expensive PH/AF3 conformer ensembles add
> downstream value over simpler inputs?

That answer decides whether Conffusion should become a coreset/adaptive
sampling thesis, a generative distillation thesis, or a useful negative result
showing that the current ensemble is not yet the bottleneck.
