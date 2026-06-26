# Conffusion

Conffusion is a documentation-first research workspace for a possible thesis
project on efficient antibody conformer ensemble modeling.

The working idea is not to train a full antibody diffusion model from scratch.
The immediate scientific target is narrower:

> Determine the smallest task-sufficient antibody conformer ensemble, then
> compress, select, or amortize Gaeun's expensive PH/AF3 pseudo-bound ensemble
> while preserving downstream representation and antigen-retrieval value.

## Central Hypothesis

Large Protein Hunter / AlphaFold3 generated antibody conformer ensembles contain
substantial structural and representation-level redundancy. A small,
strategically selected or adaptively generated subset can preserve most
task-relevant conformational coverage and downstream retrieval performance at
much lower compute cost.

This project treats the generated conformers as a computational proposal
distribution, not as a calibrated physical equilibrium ensemble.

## Why This Repo Exists

Gaeun is currently interested in optimizing conformer generation. Jerry is
looking for a thesis idea that connects that practical bottleneck to modern
generative modeling, diffusion, flow matching, and representation learning.

The key question is:

> What is the minimal structural distribution needed by a conformer-aware BCR
> representation model, and can we obtain that distribution without thousands
> of AF3 cofolds per antibody?

## Repository Map

- `docs/proposal.md`: thesis-style proposal with hypothesis and aims.
- `docs/background.md`: scientific context and rationale.
- `docs/research-plan.md`: staged plan, gates, and deliverables.
- `docs/aim1-phase1-benchmark-plan.md`: concrete data inventory and build plan
  for the first Aim 1 benchmark.
- `docs/evaluation.md`: baselines, metrics, controls, and failure modes.
- `docs/implementation-notes.md`: manifest fields, data handling, and future
  implementation notes.
- `docs/questions-for-gaeun-and-sophia.md`: questions to resolve before
  committing.
- `docs/references.md`: papers, methods, and source links.
- `docs/source-notes.md`: notes on the GPT-Pro and Deep Research outputs.
- `docs/model-outputs/`: repo-local copies of the raw GPT-Pro and Deep
  Research markdown outputs.
- `notes/model-report-synthesis.md`: synthesis of the three model reports.
- `notes/dataset-copies.md`: Ragon paths and provenance for Jerry-owned copied
  conformer datasets.
- `notes/`: scratch notes that are not yet evidence.

## Current Position

Green light:

- Ensemble redundancy analysis.
- Task-preserving coreset selection.
- Adaptive stopping and budget allocation for PH/AF3.
- Direct representation distillation as a baseline.

Amber light:

- CDR-focused residual flow or set-of-prototypes generator, only after the
  ensemble benefit and coreset structure are established.

Red light:

- Full all-atom antibody generation from scratch on the current target count.
- Claims that PH/AF3 sample frequencies are physical occupancy probabilities.
- Claims that sequence-only generated ensembles add antigen information unless
  additional deployable conditioning variables are present.

## Immediate Next Step

Before committing to a generative thesis, run a short saturation and leakage
study:

1. Does the full ensemble outperform sequence-only, static-structure, and
   repeated-single-conformer baselines on a leakage-safe endpoint?
2. Can an oracle or simple coreset shrink the ensemble without losing
   structural coverage or downstream performance?
3. Does blind PH/AF3 provide useful modes beyond fast antibody ensemble
   baselines such as ABB4-STEROIDS?

If those answers are positive, a distillation or lightweight generative
extension is justified. If not, the scientifically useful result is that the
current ensemble is redundant, leaky, or not yet exploited by the downstream
model.

## Ownership And Scope

This repo is a Jerry-owned planning and documentation workspace. It should not
modify Gaeun-owned source code, shared generated conformers, cluster jobs, or
collaborator outputs. External assets should be read-only unless copied into a
Jerry-owned workspace with clear provenance.
