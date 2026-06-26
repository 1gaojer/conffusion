# Background And Rationale

## Problem Setting

Antibodies are not static objects. Their binding-relevant regions, especially
CDR-H3 and sometimes VH/VL orientation, can adopt multiple conformations.
Gaeun's conformer-generation workflow attempts to expose useful antibody
geometries by generating pseudo-bound structures with Protein Hunter and AF3.
Those conformers can then feed conformer-aware models such as MCA/ConFormer for
antigen retrieval or related downstream tasks.

The workflow is expensive. Current versions of the pipeline are parameterized by
PH design count, cycle count, AF3 seed count, AF3 model count, and filters. The
exact arithmetic is version-dependent and must be verified before use. The
important point for this project is that the process can create thousands of
candidate structures per antibody before filtering.

## What Distribution Is Being Modeled?

The generated ensemble should be written as a teacher or proposal distribution:

```text
q_PH/AF3(x | sequence, pipeline settings, pseudo-binder/context, filters)
```

where `x` is an antibody conformation.

This is not the same as a physical equilibrium distribution:

```text
p_physical(x | sequence, temperature, solvent, pH, ligand state, environment)
```

PH/AF3 sample counts do not automatically imply physical occupancies. A mode
sampled 50 times is not necessarily 50 times more populated in solution than a
mode sampled once. Counts may reflect branch settings, pseudo-binder choices,
AF3 seeds, templates, filters, or implementation details.

Therefore the defensible claim is:

> We model, compress, or distill a useful computational proposal ensemble.

The claim to avoid is:

> We recover the true antibody energy landscape.

## Why Compression Comes Before Generation

The first baseline to beat is not a diffusion model. It is smart subsampling of
the existing teacher ensemble.

If random or cluster-aware selection of 16 to 64 conformers preserves downstream
performance, then the key contribution is not a new generator. It is the
discovery that the pipeline is over-sampling redundant modes, plus a procedure
for obtaining a smaller task-sufficient ensemble.

If random subsets fail but oracle coresets work, then a learned selector,
adaptive sampler, or set-of-prototypes model becomes justified.

If even oracle subsets fail, then a small generated ensemble is unlikely to
match the full teacher unless it produces better-than-teacher structures, which
would require independent validation.

## Why This Still Connects To Generative Modeling

The conformer problem is still a generative-modeling problem. Each antibody has
a conditional distribution over possible conformers or conformer clusters.
Diffusion, flow matching, and structure-token transformers are natural tools
for multimodal distributions.

The practical thesis path is:

1. define the useful target distribution;
2. measure how much of it is redundant;
3. learn to select or sample a compact approximation;
4. only then test a lightweight generative model.

The likely useful generative target is local and structured:

- CDR-H3 or all-CDR backbone frames;
- VH/VL rigid-body orientation;
- torsion perturbations around an anchor structure;
- latent conformer cluster prototypes.

Full all-atom antibody generation from scratch is a poor first target because
it spends capacity on framework atoms and side-chain details that may not drive
the downstream task.

## Key Risks

### Pipeline Artifact Learning

A student model may learn PH/AF3 artifacts rather than antibody flexibility. It
could reproduce seed, template, MSA, confidence, pseudo-binder, or filter
patterns.

### Privileged Conditioning

If true antigen contacts, cognate-antigen size, templates, or complex-derived
information are used to generate conformers that are later evaluated on antigen
prediction, the conformer ensemble may carry privileged information unavailable
at deployment.

### Leakage

Conformer-level splitting is invalid. Near-identical conformers from one
antibody target must not appear in both train and test. Splits should group by
paired sequence, CDR-H3 similarity, clonotype, V/J gene where relevant, PDB
provenance, and antigen family for downstream retrieval tasks.

### No Downstream Benefit

The full ensemble may not improve over sequence-only or one-structure baselines.
That would not be a failed experiment. It would show that current conformer
generation is not yet justified for the chosen endpoint.

### Overclaiming Biology

Without independent MD, apo/holo structures, repeated experimental structures,
NMR, HDX, SAXS, or prospective validation, the project should not claim physical
dynamics.

## Working Interpretation

This project is strongest as an engineering-science thesis:

> What information in a huge conformer ensemble is necessary for downstream
> antibody representation learning, and how cheaply can that information be
> preserved?

It is weaker as a pure generative-model thesis:

> Can I train a new antibody diffusion model from a few thousand teacher
> ensembles?

The former directly supports Gaeun's conformer-generation optimization. The
latter is underpowered and crowded by existing protein and antibody generators.
