# References And Starting Points

This is a reading list, not a completed literature review.

## Protein Diffusion And Flow Models

- RFdiffusion: de novo protein structure and function design with diffusion.
  https://www.nature.com/articles/s41586-023-06415-8
- FrameDiff: SE(3) diffusion over protein backbone frames.
  https://arxiv.org/abs/2302.02277
- AlphaFlow / ESMFlow: flow matching for protein conformational ensembles.
  https://arxiv.org/abs/2402.04845
- AlphaFlow code.
  https://github.com/bjing2016/alphaflow
- ESM3: multimodal protein sequence, structure, and function generation.
  https://www.science.org/doi/10.1126/science.ads0018

## Protein Ensemble Modeling

- BioEmu: scalable emulation of protein equilibrium ensembles with generative
  deep learning.
  https://www.microsoft.com/en-us/research/publication/scalable-emulation-of-protein-equilibrium-ensembles-with-generative-deep-learning/
- AlphaFlow-Lit: lightweight AlphaFlow-style sampling.
  https://openreview.net/forum?id=Z6fPAsu91p

## Antibody Generative Models

- DiffAb: antigen-specific antibody CDR sequence/structure generation with
  diffusion.
  https://proceedings.neurips.cc/paper_files/paper/2022/hash/3fa7d76a0dc1179f1e98d1bc62403756-Abstract-Conference.html
- IgDiff: SE(3) diffusion for antibody variable-domain design.
  https://arxiv.org/abs/2405.07622
- AbDiffuser: full-atom antibody structure and sequence diffusion.
  https://arxiv.org/abs/2308.05027
- ABodyBuilder4-STEROIDS: antibody conformational ensemble sampling with a
  flow-matching model.
  https://www.biorxiv.org/content/10.64898/2026.04.14.718378v1.full-text

## Related Benchmarking And Guardrails

Topics to expand:

- antibody train/test leakage and family-aware splits;
- PDB/SAbDab overlap and temporal holdouts;
- apo/holo and multistate antibody structure validation;
- CDR-H3 loop clustering and canonical loop definitions;
- retrieval metrics and hard-negative construction;
- non-inferiority testing for compressed ensembles.

## Local Context To Keep In Sync

Relevant Research workspace files:

- `/Users/jerrygao/Research/research-state.md`
- `/Users/jerrygao/Research/projects/gaeun-thesis-support.md`
- `/Users/jerrygao/Research/projects/bcr-conformer-pipeline.md`
- `/Users/jerrygao/Research/projects/thesis-ideas.md`
- `/Users/jerrygao/Research/context/research-environment.md`

Relevant local repo:

- `/Users/jerrygao/Research/liu-lab/bcr-conformer-pipeline`
