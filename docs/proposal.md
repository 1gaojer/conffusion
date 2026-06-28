# Proposal

## Working Title

Task-sufficient CDR/paratope ensemble representations for antigen retrieval:
readout-first compression and adaptive sampling of pseudo-bound BCR conformer
ensembles.

## Rationale

Gaeun's conformer-generation workflow uses Protein Hunter and AlphaFold3-style
cofolding to produce large antibody/BCR conformer ensembles. These ensembles are
intended to expose antibody conformational variability that may improve
MCA/ConFormer-style representation learning and downstream antigen retrieval.
The workflow is expensive, version-dependent, and can produce thousands of
candidate structures per antibody before filtering.

The central unresolved question is not whether diffusion models can generate
protein-like structures. That has already been demonstrated by several protein
and antibody generative-model families. The sharper question is whether Gaeun's
large PH/AF3 ensembles contain CDR/paratope-local information that improves
antigen retrieval, and whether that information can be read out, selected, or
generated more cheaply.

This proposal treats the generated conformers as a pseudo-bound computational
proposal ensemble. It does not assume that PH/AF3 sample frequencies are
Boltzmann weights or that the generated structures represent a true physical
equilibrium distribution.

## Central Hypothesis

Large PH/AF3 antibody conformer ensembles contain substantial structural and
representation-level redundancy. The current evidence suggests that useful
antigen-neighbor signal is concentrated in CDR/paratope-local MCA tensor
readouts and can be diluted by global H/L mean pooling. A small, strategically
selected or adaptively generated subset can preserve downstream value only if
the task-relevant readout is defined first.

## Aim 1: Establish The Retrieval Endpoint And CDR/Paratope Readout

The first aim is to establish a leakage-aware antigen-retrieval endpoint and a
CDR/paratope-aware ensemble readout before optimizing conformer count. This aim
is still a benchmark and sanity-check phase, not a generative-modeling phase.
It has three parts.

### Phase 1: Build A Fair Test

Define the benchmark unit as one antibody/BCR target with heavy/light sequence,
generated conformers, CDR/paratope annotations, and a downstream retrieval label
such as antigen group, antigen source, disease panel, or nearest-neighbor truth.
Freeze the benchmark set, endpoint, and split policy before comparing methods.

Splits must be target-level or family-level. No conformers from the same
antibody target, or close antibody family when family information is available,
should appear on both sides of train/test. This phase answers: what exactly is
being tested, and how is leakage being prevented?

### Phase 2: Define And Validate The Readout

For the same conformer tensors, compare:

- global H/L pooling;
- H-CDR3-only readout;
- all-CDR H/L readout;
- all-CDR mean plus variance or mean plus standard deviation;
- framework-only controls;
- random framework-window controls;
- shuffled-label controls;
- sequence-only and static-structure baselines where available.

This phase answers whether useful conformer information is visible only under a
CDR/paratope-aware readout, and whether the effect survives reasonable controls.

### Phase 3: Test Whether The Ensemble Adds Useful Retrieval Signal

For each antibody target, compare:

- full PH/AF3 ensemble;
- sequence-only baseline;
- one predicted or static antibody structure;
- repeated copies of one conformer;
- reduced AF3-only ensemble;
- simple random coordinate perturbations;
- reduced PH/AF3 branches where available.

Evaluation should include three classes of output:

- structural diversity: CDR-H3 diversity, all-CDR diversity, VH/VL orientation
  diversity, and cluster coverage;
- representation behavior: ConFormer embedding preservation, embedding
  stability, and whether ensemble information changes the representation in a
  consistent way;
- downstream performance: nearest-neighbor antigen retrieval, disease-panel
  ranking, Hit@K, MRR, macro-F1 when classification is used, and ranking-bias
  diagnostics.

Structural diversity alone is not sufficient. The full ensemble is useful only
if conformational diversity translates into representation or retrieval signal.

Expected outcome: This aim determines whether conformational variability is
actually being used by the downstream model and which readout exposes it. If the
full ensemble does not beat sequence-only, static-structure, and
repeated-single-conformer controls under a defensible readout, the project
should diagnose architecture, leakage, or unnecessary pipeline complexity before
attempting conformer compression or generative modeling.

## Aim 2: Identify The Smallest Task-Sufficient Ensemble Under The Fixed Readout

The second aim is to build ensemble-size saturation curves after Aim 1 freezes
the readout and endpoint. For each target, evaluate subset sizes:

```text
K = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, full
```

Selection methods should include:

- uniform random sampling, repeated many times;
- stratified sampling across PH design, cycle, AF3 seed, and AF3 model;
- k-medoids or k-center sampling in structural feature space;
- farthest-point sampling;
- cluster medoids;
- representation-aware coreset selection;
- ConFormer moment matching when appropriate;
- retrieval-aware or task-aware selection once the readout is fixed.

Metrics should include CDR-H3 and all-CDR coverage, VH/VL orientation diversity,
cluster recall, teacher-to-subset and subset-to-teacher distances, CDR/paratope
readout preservation, downstream output stability, retrieval metrics, and
compute/storage cost.

Expected outcome: This aim produces a direct answer to Gaeun's practical
optimization question: how far can the ensemble be shrunk before structural or
downstream quality degrades?

## Aim 3: Reduce Generation-Time Cost Prospectively

Post-hoc compression can reduce storage and downstream training cost, but it
does not recover AF3 compute already spent. The third aim is to identify where
new conformational modes actually enter the pipeline and use that information
to reduce generation cost.

Perform hierarchical variance decomposition across:

- PH design;
- PH cycle;
- contact-conditioning choice;
- pseudo-binder size or identity where available;
- AF3 seed;
- AF3 model;
- filter outcome.

Then test reduced or adaptive generation protocols:

- fewer PH designs;
- fewer cycles;
- fewer AF3 seeds;
- fewer AF3 models;
- one AF3 output per PH branch before deeper allocation;
- stop when new cluster discovery saturates;
- allocate additional compute only to branches with high marginal mode yield.

Expected outcome: This aim converts post-hoc insight into prospective compute
savings. A strong result would reduce AF3 calls or retained structures by an
order of magnitude while preserving the Aim 1 readout and downstream performance
within a predefined non-inferiority margin.

## Optional Aim 4: Distill A Small Conformer Prototype Generator

Only after Aims 1-3 show that there is a coherent, useful CDR/paratope readout
and a compact target distribution, train a lightweight generative or
set-prediction model. The recommended first formulation is not full all-atom
antibody generation. It is one of:

- predict 8 to 32 conformer cluster medoids and weights;
- generate CDR backbone frames or torsions around an anchor structure;
- transport a pretrained antibody ensemble prior toward the PH/AF3 pseudo-bound
  teacher distribution;
- distill the full-ensemble CDR/paratope-aware representation directly from
  sequence or a small anchor ensemble.

Expected outcome: This aim tests whether modern generative modeling can
amortize the useful parts of the teacher distribution. It is a second-stage
extension, not the dependency on which the thesis should succeed or fail.

## Success Criteria

A successful core thesis would show that a CDR/paratope-aware ensemble readout
captures useful antigen-retrieval signal, and that a small selected or
adaptively generated ensemble preserves that signal at substantially lower cost.

A defensible positive result:

> A CDR/paratope-aware readout improves retrieval over global pooling, and a
> task-aware adaptive sampler retains that retrieval performance within a
> predefined non-inferiority margin while reducing AF3 calls or retained
> structures by at least an order of magnitude.

A defensible negative result:

> The full ensemble does not outperform strict sequence-only, static-structure,
> or repeated-single-conformer baselines, indicating that current conformer
> generation is redundant, leaky, or not exploited by the downstream model.

Both outcomes are scientifically useful.
