# Proposal

## Working Title

How many antibody conformers are enough? Task-preserving compression and
adaptive sampling of pseudo-bound BCR conformer ensembles for antigen-retrieval
models.

## Rationale

Gaeun's current conformer-generation pipeline produces large antibody/BCR
conformer ensembles using Protein Hunter and AF3-style cofolding. These
ensembles are intended to expose conformational variability that may help
MCA/ConFormer-style representation learning and antigen retrieval. The pipeline
is expensive and can produce thousands of candidate structures per antibody
before filtering. It is not yet clear whether all of this conformational
diversity is necessary, whether it improves downstream performance beyond
sequence or static-structure baselines, or whether most useful signal can be
preserved with a much smaller subset.

This project treats the generated structures as pseudo-bound computational
proposal ensembles. It does not assume that PH/AF3 sample frequencies are
physical occupancies or that the generated set is a true antibody dynamics
ensemble.

## Central Hypothesis

Large PH/AF3-generated antibody conformer ensembles contain substantial
structural and representation-level redundancy, and a small, strategically
selected or adaptively generated subset can preserve most task-relevant
conformational coverage and downstream antigen-retrieval performance at much
lower compute cost.

## Aim 1: Test Whether Full Conformer Ensembles Add Useful Signal

The first aim is to establish whether conformer variability is actually useful
under strict controls. For each antibody target, compare the full PH/AF3
ensemble against:

- sequence-only baselines;
- one static or predicted structure;
- repeated copies of one conformer;
- reduced-seed AF3 without the full PH protocol;
- random coordinate or torsion perturbations;
- shuffled-conformer negative controls.

Evaluation should include structural diversity, frozen ConFormer embedding
preservation, output stability, and downstream nearest-neighbor antigen
retrieval or disease-panel ranking. Splits must be target-level and should be
family-aware to avoid conformer-level leakage.

Expected outcome: this aim determines whether the downstream model uses
meaningful conformational variability. If full ensembles do not outperform
simpler controls, the project should diagnose leakage, architecture, or
unnecessary pipeline complexity before building a generator.

## Aim 2: Identify The Smallest Task-Sufficient Conformer Subset

The second aim is to build ensemble-size saturation curves. Evaluate subset
sizes such as:

```text
K = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, full
```

Compare selection policies:

- uniform random sampling with repeated draws;
- stratified sampling across PH design, cycle, AF3 seed, and AF3 model;
- one conformer per PH design;
- quality-score selection;
- k-medoids or farthest-point sampling in structural space;
- facility-location or determinantal point process selection;
- representation-aware coresets that preserve ConFormer embedding moments.

Metrics should include CDR-H3 and all-CDR RMSD coverage, VH/VL orientation
diversity, cluster recall, teacher-to-subset coverage, subset-to-teacher
precision, embedding mean/covariance preservation, retrieval performance, and
compute/storage cost.

Expected outcome: this aim gives Gaeun a quantitative answer to how far the
ensemble can be shrunk before structural or downstream quality degrades.

## Aim 3: Reduce Generation Cost Prospectively

Post-hoc compression reduces storage and downstream model cost but does not save
the original PH/AF3 compute. The third aim is therefore to analyze where
diversity enters the pipeline and test reduced generation protocols before the
full ensemble is generated.

Decompose variation by:

- PH design;
- PH cycle;
- contact-conditioning choice;
- AF3 seed;
- AF3 model;
- filter status and confidence metrics.

Then test fixed and adaptive protocols:

- fewer PH designs;
- fewer cycles;
- fewer AF3 seeds;
- fewer AF3 models per seed;
- one AF3 prediction per PH design followed by targeted expansion;
- early stopping when cluster coverage, embedding moments, or downstream
  predictions converge.

Expected outcome: this aim determines whether the expensive pipeline can be
shortened before generation, not merely compressed afterward.

## Optional Aim 4: Distill Or Generate A Small Conformer Coreset

Only after Aims 1-3 show that there is a coherent target distribution worth
distilling, train a lightweight student model. The most realistic options are:

- a set-of-prototypes model that predicts 8-32 conformer medoids and weights;
- a CDR-focused residual flow or diffusion model around an anchor structure;
- a representation distillation model that maps sequence plus a small anchor
  set to the full-ensemble ConFormer representation.

The model should generate CDR residue frames, torsions, or cluster prototypes
rather than full all-atom antibody structures from scratch. A pretrained
geometric or antibody-ensemble prior is preferred over training from random
initialization.

## Impact

This project would turn conformer generation from an expensive black-box
preprocessing step into a measurable, optimizable component of antibody
representation learning. It directly supports Gaeun's current optimization
question by defining ensemble-quality readouts, reducing compute burden, and
clarifying whether conformational diversity improves antigen retrieval.

The defensible contribution is not "diffusion for antibodies." It is a rigorous
study of the minimal structural distribution needed by conformer-aware BCR
models, with generative modeling as a gated extension.
