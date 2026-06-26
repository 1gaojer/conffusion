# References And Source Map

This is a working source map, not a final bibliography. Verify versions,
publication status, and claims before citing in a thesis.

## Project-Local Sources To Verify

These are not tracked in this repo, but they are the first sources to inspect
before making present-tense claims.

- Research dashboard:
  `/Users/jerrygao/Research/research-state.md`
- Gaeun thesis support brief:
  `/Users/jerrygao/Research/projects/gaeun-thesis-support.md`
- BCR conformer pipeline brief:
  `/Users/jerrygao/Research/projects/bcr-conformer-pipeline.md`
- Local conformer pipeline wrapper:
  `/Users/jerrygao/Research/liu-lab/bcr-conformer-pipeline`
- Gaeun source pipeline on cluster, read-only:
  `/external/liulab/gkim/antigen_prediction/datasets/pipeline/conformer_generation`
- Jerry-owned cluster conformer pipeline copy:
  `/project/liulab/jg1920/bcr-conformer-pipeline`

## Canonical Local PDFs

PDFs are stored once in the shared Research paper library:

```text
/Users/jerrygao/Research/papers/
```

Project-specific docs should link to these canonical copies rather than storing
duplicates under project folders.

Conffusion-relevant PDFs:

- `/Users/jerrygao/Research/papers/2026-ABodyBuilder4-STEROIDS-Sampling-Antibody-Conformational-Ensembles.pdf`
- `/Users/jerrygao/Research/papers/2025-ITsFlexible-Antibody-TCR-CDR-Conformational-Flexibility.pdf`
- `/Users/jerrygao/Research/papers/2023-AbDiffuser-Full-Atom-In-Vitro-Functioning-Antibodies.pdf`
- `/Users/jerrygao/Research/papers/2024-AlphaFlow-Lit-Efficient-Protein-Ensemble-Generation.pdf`
- `/Users/jerrygao/Research/papers/2024-AlphaFlow-Flow-Matching-Protein-Ensembles.pdf`
- `/Users/jerrygao/Research/papers/2025-BioEmu-Protein-Equilibrium-Ensembles.pdf`
- `/Users/jerrygao/Research/papers/2022-DiffAb-Antigen-Specific-Antibody-Design-Diffusion.pdf`
- `/Users/jerrygao/Research/papers/2024-FlowAB-Energy-Guided-SE3-Flow-Matching-Antibody-Refinement.pdf`
- `/Users/jerrygao/Research/papers/2023-FrameDiff-SE3-Diffusion-Protein-Backbone-Generation.pdf`
- `/Users/jerrygao/Research/papers/2024-IgDiff-De-Novo-Antibody-Design-SE3-Diffusion.pdf`
- `/Users/jerrygao/Research/papers/2025-CF-Random-AF2-Sequence-Association-Alternative-Conformations.pdf`
- `/Users/jerrygao/Research/papers/2025-Protein-Hunter-Structure-Hallucination-Diffusion-Protein-Design.pdf`
- `/Users/jerrygao/Research/papers/2025-RFantibody-De-Novo-Antibody-Design-RFdiffusion.pdf`
- `/Users/jerrygao/Research/papers/2023-RFdiffusion-De-Novo-Protein-Structure-Function-Design.pdf`
- `/Users/jerrygao/Research/papers/2024-Cfold-Alternative-Protein-Conformation-Prediction.pdf`

## Protein Structure Generative Models

### RFdiffusion

- Watson et al., "De novo design of protein structure and function with
  RFdiffusion."
- Link: https://www.nature.com/articles/s41586-023-06415-8
- Relevance: establishes diffusion for protein backbone generation and
  conditional protein design.
- Caveat: primarily de novo design, not fixed-sequence conformer ensemble
  modeling.

### FrameDiff

- "SE(3) diffusion model with application to protein backbone generation."
- Link: https://arxiv.org/abs/2302.02277
- Relevance: residue-frame diffusion on SE(3), useful representation for CDR
  frame generation.
- Caveat: architecture reference, not an antibody ensemble method.

### AlphaFlow

- "AlphaFold Meets Flow Matching for Generating Protein Ensembles."
- Link: https://arxiv.org/abs/2402.04845
- PMLR: https://proceedings.mlr.press/v235/jing24a.html
- Relevance: direct conceptual precedent for converting static structure
  predictors into sequence-conditioned ensemble samplers.
- Caveat: general protein ensembles, not PH/AF3 pseudo-bound antibody
  distillation.

### AlphaFlow-Lit

- Link: https://openreview.net/forum?id=Z6fPAsu91p
- Related arXiv link reported in model output:
  https://arxiv.org/abs/2407.12053
- Relevance: efficiency-oriented ensemble sampling and amortization precedent.
- Caveat: verify current status and exact speedup before citing.

### BioEmu

- "Scalable emulation of protein equilibrium ensembles with generative deep
  learning."
- Microsoft page:
  https://www.microsoft.com/en-us/research/publication/scalable-emulation-of-protein-equilibrium-ensembles-with-generative-deep-learning/
- Science link:
  https://www.science.org/doi/10.1126/science.adv9817
- PubMed:
  https://pubmed.ncbi.nlm.nih.gov/40638710/
- Relevance: standard for serious equilibrium-ensemble emulation claims.
- Caveat: BioEmu uses MD/static/thermodynamic supervision; PH/AF3 distillation
  is a different and weaker claim.

### ESM3

- "Simulating 500 million years of evolution with a language model."
- Link: https://www.science.org/doi/10.1126/science.ads0018
- Relevance: sequence/structure/function token modeling and possible latent
  structure abstraction.
- Caveat: structure-token resolution may be too coarse for subtle CDR mode
  distinctions.

## Antibody Generative Models

### DiffAb

- "Antigen-Specific Antibody Design and Optimization with Diffusion-Based
  Generative Models."
- Link:
  https://proceedings.neurips.cc/paper_files/paper/2022/hash/3fa7d76a0dc1179f1e98d1bc62403756-Abstract-Conference.html
- Relevance: antibody CDR sequence/structure generation conditioned on antigen.
- Caveat: antibody design, not fixed-sequence conformer ensemble compression.

### AbDiffuser

- "Full-Atom Antibody Design."
- Link:
  https://proceedings.neurips.cc/paper_files/paper/2023/file/801ec05b0aae9fcd2ef35c168bd538e0-Paper-Conference.pdf
- Relevance: full-atom antibody generation with physics-informed structure
  priors.
- Caveat: broader and more complex than needed for fixed-sequence CDR conformer
  selection.

### IgDiff

- Link: https://arxiv.org/abs/2405.07622
- Relevance: antibody-specific adaptation of geometric diffusion.
- Caveat: design/generation prior, not calibrated teacher-ensemble
  distillation.

### RFantibody / Antibody-Specific RFdiffusion

- Link: https://www.nature.com/articles/s41586-025-09721-5
- Relevance: antibody-specific RFdiffusion-style design.
- Caveat: epitope-directed design, not fixed-sequence conformer ensembles.

### ABodyBuilder4-STEROIDS

- Reported title/context: antibody conformational ensembles from paired VH/VL
  sequences using an SE(3) flow model.
- BioRxiv link from GPT-Pro output:
  https://www.biorxiv.org/content/biorxiv/early/2026/04/16/2026.04.14.718378.full.pdf
- Alternate DOI-style link to verify:
  https://www.biorxiv.org/content/10.64898/2026.04.14.718378v1
- Oxford metadata:
  https://www.oqi.ox.ac.uk/publication/2412744/dimensions
- Relevance: mandatory novelty check for generic antibody ensemble generation.
- Caveat: preprint status and benchmark claims need independent verification.

## Teacher-Generation Context

### Protein Hunter

- Reported source:
  https://www.biorxiv.org/content/10.1101/2025.10.10.681530v1
- Alternate link from GPT-Pro with context:
  https://www.biorxiv.org/content/10.1101/2025.10.10.681530v2.full.pdf
- Relevance: part of Gaeun's teacher-generation mechanism.
- Caveat: Protein Hunter is a protein design method, not a calibrated antibody
  dynamics model.

### BindCraft

- Relevance: related structure-prediction-driven binder hallucination context.
- Action: add exact paper/source after verification if needed for thesis text.

## External Validation And Benchmarks

### ALL-conformations

- Link: https://www.nature.com/articles/s42256-025-01131-6
- Relevance: experimentally observed loop conformations and possible external
  support for whether generated CDR modes are observed in PDB-derived data.
- Caveat: PDB-derived and may overlap with SAbDab or training sets; requires
  strict de-overlap.

### AbBiBench

- Relevance: antibody design/evaluation benchmark warning that downstream
  functional evaluation and leakage-aware splits matter.
- Action: add verified primary source link before citing.

## Internal Model Outputs

See `docs/source-notes.md` for how to treat the three attached model reports.

Working rule:

- Deep Research output: conceptually useful, but citation handles need
  reconstruction.
- GPT-Pro with presentation context: best project-specific grounding, but
  presentation-derived claims require live verification.
- GPT-Pro without presentation context: best novelty warning around
  ABB4-STEROIDS, but any presentation-like file-citation claims require
  mapping to actual supplied documents.
