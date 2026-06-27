# Research Plan

## Phase 0: Freeze Scope And Evidence

Goal: avoid a moving-target thesis.

Tasks:

- Choose the primary dataset snapshot.
- Record the exact PH/AF3 pipeline version.
- Record the exact generation parameters:
  - PH design count;
  - PH cycle count;
  - AF3 seed count;
  - AF3 model count;
  - contact-conditioning rule;
  - pseudo-binder size rule;
  - filter criteria.
- Build a target-level manifest with one row per conformer.
- Freeze the primary downstream endpoint.
- Freeze train/test splits before model development.

Deliverable:

- `manifest_schema.md` or equivalent data dictionary.
- A short data audit report with target counts, conformer counts, failure modes,
  and known caveats.

Gate:

- Do not start model development until the dataset snapshot, split, and endpoint
  are documented.

## Phase 1: Does The Ensemble Matter?

Goal: test the premise before optimizing.

Detailed working plan: `docs/aim1-phase1-benchmark-plan.md`.

Current follow-up after the CDR structural result:
`docs/aim1-phase1-4-representation-sensitivity-plan.md`.

This phase is Aim 1 in operational form. Keep it organized around three
questions:

1. What is the fair benchmark?
2. Does the full ensemble beat simple controls?
3. Does any structural diversity translate into useful representation or
   retrieval signal?

### 1. Build A Fair Test

Tasks:

- Define the benchmark unit as one antibody/BCR target with sequence,
  conformers, and downstream label or retrieval target.
- Freeze the benchmark set and endpoint before comparing methods.
- Use leakage-safe splits, preferably target-level plus at least one stricter
  family-level split.
- Confirm that conformers from the same antibody target do not appear on both
  sides of train/test.

### 2. Compare Full Ensemble Against Simple Controls

Experiments:

- Full ensemble versus one conformer.
- Full ensemble versus repeated one conformer.
- Full ensemble versus sequence-only.
- Full ensemble versus one static or predicted structure.
- Full ensemble versus reduced AF3-only outputs.
- Full ensemble versus random coordinate perturbations.
- Optional: PH/AF3 versus ABB4-STEROIDS or another antibody ensemble baseline.

### 3. Judge By Useful Outputs

Metrics:

- structural diversity and geometry validity;
- CDR-H3, all-CDR, and VH/VL orientation coverage;
- ConFormer output stability;
- embedding mean and covariance preservation;
- downstream retrieval metrics such as Hit@K and MRR;
- disease-panel ranking or macro-F1 when classification is used;
- ranking-bias diagnostics, especially whether the same antigen groups dominate
  many unrelated queries.

Deliverable:

- premise report: does ensemble variability improve the chosen endpoint, and
  under which split/control conditions?

Go/no-go:

- If the full ensemble does not help, pause generation work and diagnose model,
  endpoint, and leakage.
- If it helps, proceed to compression.

## Phase 2: Saturation And Coreset Curves

Goal: identify the smallest task-sufficient subset.

Evaluate subset sizes:

```text
K = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, full
```

Selection strategies:

- uniform random;
- stratified by PH branch;
- stratified by AF3 seed/model;
- k-medoids;
- k-center/farthest point;
- cluster medoids;
- facility-location objective;
- ConFormer moment matching;
- task-aware selection where appropriate.

Deliverable:

- performance-versus-K curves;
- performance-versus-cost curves;
- recommended default K range.

Gate:

- If random small K works, prioritize reduced pipeline settings and adaptive
  stopping.
- If oracle or learned coresets beat random, prioritize selection.
- If no compact subset works, reconsider whether a small generator is realistic.

## Phase 3: Prospective Cost Reduction

Goal: reduce compute before full generation.

Variance decomposition:

- PH design;
- PH cycle;
- contact conditioning;
- pseudo-binder size;
- AF3 seed;
- AF3 model;
- filter outcome.

Reduced protocols:

- fewer PH designs;
- fewer cycles;
- fewer AF3 seeds;
- fewer AF3 models;
- one AF3 output per PH design before deeper allocation;
- target-specific K;
- branch-specific adaptive allocation.

Stopping criteria:

- no new structural cluster after recent samples;
- embedding mean/covariance has stabilized;
- estimated unseen-mode mass below threshold;
- downstream prediction stable;
- marginal coverage per GPU-hour below threshold.

Deliverable:

- prospective reduced-generation pilot.

Gate:

- Compression is not enough. At least one prospective experiment should show
  that generation itself can be reduced or stopped early.

## Phase 4: Optional Generative Distillation

Goal: test whether a small model can amortize useful teacher modes.

Only proceed after:

- the full ensemble adds value;
- a compact target distribution exists;
- strong non-neural coresets are established;
- privileged conditioning has been separated from blind conditioning.

Candidate models:

- set-of-prototypes transformer predicting K conformer medoids;
- CDR-H3 frame or torsion residual flow;
- all-CDR frame model around an anchor structure;
- ABB4-STEROIDS initialized transport toward PH/AF3 pseudo-bound modes;
- direct full-ensemble representation distillation.

Non-goals:

- full all-atom generation from scratch;
- claiming true equilibrium dynamics;
- treating teacher frequency as physical occupancy.

Deliverable:

- small held-out proof of concept compared against random subsets, coresets,
  reduced AF3, ABB4-STEROIDS, and the full PH/AF3 teacher.

## Tentative Timeline

### First 3 Months

- Freeze data, pipeline version, and endpoint.
- Build manifest.
- Run leakage checks.
- Run premise tests.
- Generate K-curves on a pilot subset.
- Decide whether the thesis is viable as compression/adaptive sampling.

### First 6 Months

- Scale saturation analysis to the full leakage-safe dataset.
- Implement coreset selection.
- Run prospective reduced-generation pilots.
- Quantify end-to-end compute and storage reductions.
- Produce a thesis-core report.

### First 12 Months

- Add generative distillation only if justified.
- Compare against external ensemble baselines.
- Validate on experimental or MD-supported conformational sets where possible.
- Write final thesis and paper-style methods report.
