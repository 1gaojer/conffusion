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

## Phase 1: Define The Retrieval Readout And Test Whether The Ensemble Matters

Goal: define the endpoint and readout before optimizing conformer count.

Detailed working plan: `docs/aim1-phase1-benchmark-plan.md`.

Current follow-up after the CDR structural result:
`docs/aim1-phase1-4-representation-sensitivity-plan.md`.

Current downstream gate after Phase 1.4:
`docs/aim1-phase1-5-retrieval-preservation-plan.md`.

This phase is Aim 1 in operational form. Keep it organized around four
questions:

1. What is the fair benchmark?
2. Which CDR/paratope-aware readout exposes useful ensemble signal?
3. Does the full ensemble beat simple controls under that readout?
4. Does any structural diversity translate into useful representation or
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

### 2. Define And Validate The Readout

Experiments:

- Global H/L pooling versus H-CDR3-only pooling.
- Global H/L pooling versus all-CDR H/L pooling.
- Global H/L pooling versus all-CDR mean plus variance or mean plus standard
  deviation.
- CDR readouts versus framework-only and random framework-window controls.
- True labels versus shuffled-label controls.
- PCA/UMAP inspection for the readout vectors when useful for Gaeun's framing.

### 3. Compare Full Ensemble Against Simple Controls

Experiments:

- Full ensemble versus one conformer.
- Full ensemble versus repeated one conformer.
- Full ensemble versus sequence-only.
- Full ensemble versus one static or predicted structure.
- Full ensemble versus reduced AF3-only outputs.
- Full ensemble versus random coordinate perturbations.
- Optional: PH/AF3 versus ABB4-STEROIDS or another antibody ensemble baseline.

### 4. Judge By Useful Outputs

Metrics:

- structural diversity and geometry validity;
- CDR-H3, all-CDR, and VH/VL orientation coverage;
- CDR/paratope-aware readout stability;
- embedding mean and covariance preservation where they remain informative;
- downstream retrieval metrics such as Hit@K and MRR;
- disease-panel ranking or macro-F1 when classification is used;
- ranking-bias diagnostics, especially whether the same antigen groups dominate
  many unrelated queries.

Deliverable:

- premise report: which readout exposes ensemble value, whether ensemble
  variability improves the chosen endpoint, and under which split/control
  conditions? Current report:
  `docs/aim1-phase1-premise-report.md`.

Go/no-go:

- If the full ensemble does not help, pause generation work and diagnose model,
  endpoint, and leakage.
- If only CDR/paratope readouts help, make readout design part of the core
  thesis before proceeding to compression.
- If it helps under a defensible readout, proceed to task-aware compression.

Current checkpoint status:

- Gaeun confirmed that the correct model checkpoint is the June 19
  `1000 conformer` checkpoint.
- Phase 1.5a verified the exact checkpoint path/hash and used the June 24
  `outputs_1000/reference_embeddings.npz` bank for internal retrieval/readout
  diagnostics.
- The current cleaner Gaeun-generated candidate endpoint is Tier B:
  172 generated targets that are exact-PDB-unseen relative to the full 930-ID
  checkpoint universe and have concatenated `VH+VL` global identity `<0.85` to
  every checkpoint target.
- Tier B pass manifest:
  `manifests/tierb_sequence_deoverlap_20260627/gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_pass_20260627.tsv`
- Remaining leakage control before external validation claims: antigen/source
  de-overlap. Tier C exact antigen/source-pair de-overlap has now been run as a
  leakage-sensitivity check, but only 17 same-label retrieval queries are
  evaluable, so it is not yet a stable external-validation endpoint.

Phase 1 premise decision:

- Continue to Aim 2 under a readout-first framing.
- Treat all-CDR mean+std as the current leading readout, with H-CDR3/L-CDR3,
  all-CDR H/L mean, global H/L, framework, random-window, shuffled-label, and
  single-conformer controls retained.
- Use Tier B Step 4 as the powered development endpoint for now.
- Use Tier C only as a caveated leakage stress test unless the candidate
  bank/grouping is expanded.
- Do not claim solved conformer compression until a selector beats random under
  the fixed CDR/paratope-aware retrieval endpoint.

## Phase 2: Task-Aware Saturation And Coreset Curves

Goal: identify the smallest task-sufficient subset under the fixed readout.

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
- CDR/paratope readout-aware selection;
- retrieval-aware or task-aware selection where appropriate.

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

Goal: reduce compute before full generation while preserving the Aim 1 readout.

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
