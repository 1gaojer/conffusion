# Research Plan

## Core Question

What is the minimal structural distribution needed by a conformer-aware
antibody representation model, and can that distribution be obtained without
thousands of AF3 cofolds per antibody?

## Phase 0: Freeze Scope And Data Contracts

Before any modeling:

- define the primary downstream endpoint;
- freeze a dataset snapshot and split;
- record the exact PH/AF3 pipeline version;
- record all generation parameters;
- define the unit of independence;
- document what information is available at deployment.

Important current version note: project notes and presentation-derived reports
refer to a 25-design PH setup in some contexts, while the local Jerry-owned
wrapper currently exposes defaults of 20 designs, 5 cycles, and 5 AF3 seeds.
Treat design count and total cofold count as versioned parameters until the
exact live Gaeun configuration is verified.

Deliverable:

- `manifest_targets.tsv` with one row per antibody target.
- `manifest_conformers.tsv` with one row per conformer.
- A short provenance report that records data sources, pipeline version, and
  known leakage risks.

## Phase 1: Redundancy And Signal Audit

Goal: determine whether the conformer ensemble contains meaningful structure
for the downstream model.

Analyses:

- count retained conformers per target;
- cluster conformers after framework alignment;
- quantify near-duplicates;
- decompose variation by PH design, cycle, AF3 seed, and AF3 model;
- measure CDR-H3, all-CDR, and VH/VL orientation variability;
- test ConFormer order, duplicate, and multiplicity sensitivity;
- compare full ensemble to sequence-only, static, repeated-single, and reduced
  AF3 controls.

Go/no-go criteria:

- Continue if the full ensemble shows reproducible benefit or has clear
  nonredundant structural modes.
- Pause if the full ensemble does not beat sequence-only or single-structure
  baselines under strict splits.

## Phase 2: Subsampling And Coreset Curves

Goal: estimate the smallest task-sufficient ensemble.

Evaluate subset sizes:

```text
1, 2, 4, 8, 16, 32, 64, 128, 256, 512, full
```

Selection strategies:

- random sampling;
- stratified random sampling by PH branch and AF3 seed;
- one conformer per PH design;
- quality-score selection;
- k-medoids;
- farthest-point or k-center sampling;
- facility-location selection;
- determinantal point process selection;
- embedding-moment matching;
- cluster-weighted medoids.

Deliverables:

- performance-versus-K curves;
- coverage-versus-K curves;
- compute/storage-versus-K curves;
- non-inferiority summary for candidate K values.

Interpretation:

- If random small K works, optimize fixed reduced generation.
- If oracle or coreset small K works but random does not, build a learned
  selector or prototype predictor.
- If no small K works, the full ensemble may be necessary or the downstream
  model may be sensitive to raw multiplicity rather than modes.

## Phase 3: Prospective Compute Reduction

Goal: reduce generation cost before the full archive exists.

Test controlled protocols on untouched targets:

- fewer PH designs;
- fewer PH cycles;
- fewer AF3 seeds;
- fewer AF3 models per seed;
- one AF3 model per PH design followed by selective expansion;
- adaptive stopping based on structural coverage or embedding convergence.

Potential stopping signals:

- no new cluster in recent samples;
- full-ensemble embedding mean/covariance estimate has stabilized;
- downstream prediction has stabilized;
- predicted unseen-mode mass is below threshold;
- marginal cluster coverage per GPU-hour is too low.

Deliverable:

- a prospective reduced-generation protocol with measured GPU-hour, storage,
  and downstream tradeoffs.

## Phase 4: Optional Generative Or Distillation Model

Only proceed if Phases 1-3 show a meaningful target distribution and a useful
small-K endpoint.

Candidate models:

1. Set-of-prototypes transformer:
   - encode paired VH/VL sequence and one anchor structure;
   - use K learned query tokens;
   - predict K CDR-frame medoids and optional weights;
   - train with optimal-transport or Hungarian matching to teacher clusters.

2. CDR-focused residual flow:
   - condition on paired sequence, anchor structure, and CDR masks;
   - generate residue frames or torsions for CDRs;
   - optionally perturb VH/VL orientation;
   - reconstruct side chains with a separate packer or leave side chains out
     unless the downstream model needs them.

3. Representation distillation:
   - predict full-ensemble ConFormer representation from sequence plus a small
     structure set;
   - use when deployment does not require explicit PDB structures.

Preferred result:

```text
32-64 student samples recover most task-relevant structural support and
downstream performance of a 1500+ conformer teacher ensemble at much lower
marginal compute.
```

## 3-Month Milestone

- Frozen manifest and split.
- Redundancy report on 50-100 diverse targets.
- ConFormer sensitivity tests.
- Initial K-curves.
- Baseline comparison against sequence-only, static structure, repeated-single
  conformer, and reduced AF3.
- Written decision on whether generative modeling is justified.

## 6-Month Milestone

- Scaled redundancy and K-curves.
- Coreset method that preserves structural and downstream metrics.
- Prospective reduced-generation pilot.
- End-to-end compute accounting.
- Thesis core is viable without a generator.

## 12-Month Milestone

- Optional learned selector, prototype model, CDR-frame flow, or representation
  distillation model.
- External validation against available multistate, apo/holo, MD, or
  antibody-ensemble baselines.
- Paper- or thesis-quality benchmark and methods writeup.
