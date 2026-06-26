# Model Report Synthesis

This note records the synthesis of three model-generated reports that Jerry
provided while exploring the conffusion idea.

## Reports Reviewed

- `deep-research-report (without gaeun's presentation context).md`
- `gpt-pro (with presentation context).md`
- `gpt-pro (without presentation context).md.md`

The presentation-aware report was the most project-specific because it noticed
details about the PH/AF3 hierarchy, ConFormer pooling, possible contact or
binder-size leakage, and the distinction between post-hoc compression and true
generation-time savings.

The presentation-free GPT-Pro report was most useful for the novelty warning:
ABodyBuilder4-STEROIDS now directly targets antibody conformational ensemble
sampling, so a generic antibody ensemble generator is not the right novelty
claim.

The deep-research report was the cleanest compact summary: distillation and
compression are good thesis framings; claiming true antibody dynamics from
PH/AF3 teacher outputs is weak.

## Synthesis

The reports agree on one central recommendation:

> Make task-preserving ensemble compression and adaptive sampling the thesis
> core. Treat diffusion or flow modeling as a gated extension.

Important preserved ideas:

- PH/AF3 outputs should be described as pseudo-bound proposal ensembles.
- Teacher sample frequencies should not be interpreted as physical occupancies.
- Millions of conformer files are not millions of independent examples.
- Splits must be target-level and family-aware.
- The first baseline to beat is smart subsampling, not diffusion.
- A learned generator is justified only if the full ensemble helps and simple
  selection cannot fully explain the useful signal.

## Proposed One-Sentence Thesis

This project asks how much of a large PH/AF3-generated pseudo-bound antibody
conformer ensemble is actually needed by conformer-aware antigen-retrieval
models, and whether that task-relevant structural distribution can be selected
or sampled much more cheaply.
