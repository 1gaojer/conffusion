# Background And Rationale

## What Conformer Generation Is For

Antibodies are not single rigid objects. Their frameworks are relatively
constrained, but CDR loops, especially CDR-H3, and VH/VL orientation can vary.
For antigen-retrieval or antibody-representation models, the hope is that a
conformer ensemble exposes shape and paratope states that a single predicted
structure misses.

Gaeun's current pipeline uses Protein Hunter and AF3-style cofolding to create
large pseudo-bound conformer ensembles. The intent is practical: produce many
plausible antibody geometries that downstream models can use, then filter and
preprocess them for conformer-aware learning.

## The Distribution Being Modeled

A useful formalization is:

```text
q_teacher(x | s, u)
```

where:

- `x` is an antibody conformation;
- `s` is the paired VH/VL sequence;
- `u` includes PH design branch, cycle, contact-conditioning choices, pseudo
  binder information, AF3 seed/model, templates, MSA state, confidence scores,
  and filtering.

This is different from:

```text
p_physical(x | s, temperature, solvent, binding state, ...)
```

PH/AF3 ensembles are computational proposal ensembles. Their sample counts
reflect a protocol, not thermodynamic occupancies. A mode sampled 100 times by
PH/AF3 is not automatically 100 times more populated in solution than a mode
sampled once.

## Why A Generative Model Is Plausible

The problem is mathematically compatible with diffusion, flow matching, and
transformer-based generation because the output is a conditional distribution
over structures:

```text
p(conformer | antibody sequence, deployable context)
```

or, more practically:

```text
p(conformer cluster | antibody target)
```

The key is the representation. Antibody conformers should not be modeled like
pixels. Reasonable structural variables include:

- residue frames in SE(3);
- backbone torsions;
- CDR loop anchor frames;
- VH/VL orientation;
- side-chain chi torsions if needed;
- learned conformer-cluster or structure-token latents.

The model should encode geometry by construction: global rotation/translation
equivariance, chain continuity, chirality, bond geometry, sterics, loop closure,
and disulfide constraints.

## How This Differs From Image Diffusion And LoRA

The useful analogy from image/video diffusion:

- condition on an input;
- sample from noise;
- generate multiple plausible outputs;
- use guidance or reranking;
- distill a slow sampler into a faster one;
- adapt a pretrained prior.

The non-transferable parts:

- proteins have hard geometric validity constraints;
- conformers are unordered sets, not video frames;
- millions of conformer files are not millions of independent examples;
- sample frequency in a teacher pipeline is not physical occupancy;
- a "style" LoRA analogue may only learn AF3/PH artifacts.

If parameter-efficient adaptation is used, it should be on top of a strong
protein or antibody structural prior. The LoRA-like mechanism is an engineering
strategy, not the scientific thesis.

## Relevant Model Families

General protein generators:

- RFdiffusion shows that diffusion over protein structures can support powerful
  conditional backbone and binder design.
- FrameDiff gives a principled SE(3) diffusion formulation over residue frames.
- Flow-matching variants such as AlphaFlow and ESMFlow show how static
  structure predictors can become sequence-conditioned ensemble samplers.

Protein ensemble models:

- AlphaFlow/ESMFlow are highly relevant because they model conformational
  ensembles rather than single structures.
- BioEmu is a useful reference for what equilibrium-emulation claims require:
  MD data, static structures, experimental observables, and physical validation.

Antibody-specific generators:

- DiffAb, IgDiff, IgFlow, AbDiffuser, and RFantibody-like work show that
  antibody and CDR generation are active, feasible model classes.
- ABodyBuilder4-STEROIDS is especially important because it directly targets
  antibody conformational ensembles using flow matching. That crowds out a
  generic "VH/VL sequence to antibody ensemble" novelty claim.

This project's novelty should therefore be the PH/AF3 pseudo-bound,
task-directed proposal distribution and the downstream question of how much
of that distribution is actually necessary.

## Why Compression Comes Before Generation

Before training a generator, the project must answer:

1. Does the full ensemble beat simple controls?
2. How many nonredundant modes exist?
3. Which pipeline levels create useful diversity?
4. Can simple selection already preserve downstream behavior?
5. Does a fast antibody ensemble model already match PH/AF3?

If simple subsampling works, the practical solution may be adaptive stopping or
branch-aware reduced generation. If oracle coreset selection works but random
selection fails, a learned selector is justified. If PH/AF3 uniquely discovers
task-relevant pseudo-bound modes, then a generative distillation model becomes
scientifically interesting.
