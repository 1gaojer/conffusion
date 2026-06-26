# Conffusion

Conffusion is a proposal and planning repository for a possible Jerry-owned
thesis project at the intersection of antibody conformer generation,
conformer-aware representation learning, and generative modeling.

The working idea is not "train Stable Diffusion for antibodies." The sharp
version is:

> Determine the smallest task-sufficient antibody conformer ensemble, then
> compress, select, or amortize Gaeun's expensive Protein Hunter / AF3
> pseudo-bound conformer-generation pipeline while preserving downstream
> retrieval and representation value.

## Current Thesis Hypothesis

Large PH/AF3-generated antibody conformer ensembles contain substantial
structural and representation-level redundancy. A small, strategically selected
or adaptively generated subset can preserve most task-relevant conformational
coverage and downstream antigen-retrieval performance at much lower compute
cost.

This repo treats the generated structures as computational proposal ensembles,
not as calibrated physical equilibrium ensembles.

## Why This Exists

Gaeun's current conformer-generation optimization question is practical and
scientific:

- How many generated conformers are actually needed?
- Which pipeline branches produce new structural modes?
- Does the full ensemble improve downstream retrieval beyond sequence-only,
  static-structure, reduced-AF3, or repeated-single-conformer baselines?
- If the ensemble helps, can its useful distribution be sampled more cheaply?

Jerry's personal interest in diffusion, flow, LoRA-style adaptation, and
generative image/video models motivates the optional modeling direction. The
core thesis, however, should survive even if a learned generator never beats
strong selection and adaptive-sampling baselines.

## Repository Map

- [docs/proposal.md](docs/proposal.md): one-page proposal plus extended aims.
- [docs/background.md](docs/background.md): scientific and modeling rationale.
- [docs/research-plan.md](docs/research-plan.md): phased plan and go/no-go
  gates.
- [docs/evaluation.md](docs/evaluation.md): baselines, metrics, leakage
  controls, and statistical tests.
- [docs/implementation-notes.md](docs/implementation-notes.md): expected data
  manifests, pipeline metadata, and ownership boundaries.
- [docs/questions-for-gaeun-and-sophia.md](docs/questions-for-gaeun-and-sophia.md):
  questions to resolve before committing.
- [docs/references.md](docs/references.md): papers and tools to read.

## Project Boundaries

This is a Jerry-owned planning repository. It should not modify Gaeun-owned
source code, shared outputs, cluster jobs, or collaborator-owned data. Any
future analysis code should read from approved sources and write only to
Jerry-owned paths unless explicit approval says otherwise.

## Current Position

Green light:

- Ensemble redundancy analysis.
- Subsampling and coreset curves.
- ConFormer sensitivity tests.
- Prospective reduced-generation and adaptive-stopping studies.

Amber light:

- CDR-focused flow/diffusion or set-of-prototypes distillation after the
  compression premise survives.

Red light:

- Full all-atom antibody generation from scratch.
- Claims about true antibody dynamics, free energies, or Boltzmann occupancies
  from PH/AF3 sample frequencies alone.

## Immediate Next Step

Freeze one pilot dataset and answer three questions:

1. Does the full ensemble outperform sequence-only, static-structure, and
   repeated-single-conformer controls on a leakage-aware retrieval endpoint?
2. Can an oracle or simple coreset shrink the ensemble while preserving
   structural coverage and downstream behavior?
3. Does blind PH/AF3 outperform a fast antibody-ensemble baseline such as
   ABodyBuilder4-STEROIDS or reduced-seed AF3?

If yes, distillation becomes compelling. If no, the negative result is still
scientifically useful and should redirect the thesis away from generation.
