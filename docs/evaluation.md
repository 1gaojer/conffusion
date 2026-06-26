# Evaluation Framework

## Primary Principle

The project should evaluate whether a smaller conformer set preserves useful
information, not whether it produces visually diverse structures.

Every result should be reported with:

- target-level or family-level splits;
- repeated subset trials;
- target-level confidence intervals;
- a predefined non-inferiority margin;
- compute and storage accounting.

## Baselines

### Minimal Baselines

- sequence-only model;
- one static or predicted structure;
- one conformer repeated N times;
- full PH/AF3 teacher ensemble;
- uniform random K-subsets;
- stratified K-subsets.

### Structural Selection Baselines

- k-medoids;
- k-center or farthest-point sampling;
- cluster medoids;
- facility-location objective;
- determinantal point process;
- kernel herding;
- ConFormer moment matching.

### Generation Baselines

- reduced AF3-only;
- reduced PH/AF3 settings;
- one AF3 output per PH design;
- ABB4-STEROIDS or another fast antibody ensemble model if feasible;
- direct full-ensemble representation distillation.

## Structural Metrics

Report coverage and precision separately.

Coverage:

- teacher-to-subset nearest-neighbor distance;
- teacher cluster recall;
- maximum uncovered cluster radius;
- CDR-H3 RMSD coverage;
- all-CDR RMSD coverage;
- VH/VL orientation coverage;
- paratope surface-shape coverage.

Precision:

- subset-to-teacher nearest-neighbor distance;
- generated-to-teacher distance;
- generated-to-experimental distance where available;
- fraction of generated samples near a valid teacher or external mode.

Distributional metrics:

- pairwise RMSD distribution;
- per-residue RMSF;
- torsion distributions;
- pairwise-distance feature MMD or energy distance;
- cluster-weight divergence only when weights are meaningful.

## Geometry And Validity

RMSD alone is insufficient. Measure:

- chain breaks;
- peptide bond geometry;
- bond-length and bond-angle deviations;
- chirality;
- cis/trans peptide state;
- Ramachandran outliers;
- steric clashes;
- CDR loop closure;
- disulfide geometry;
- VH/VL interface clashes;
- side-chain rotamer plausibility;
- confidence scores, with the caveat that they are teacher diagnostics.

Any relaxation, repacking, or revalidation cost must be included in the compute
budget.

## Representation Metrics

Because ConFormer or MCA-style models are intended consumers, evaluate:

- pooled embedding cosine distance;
- mean embedding error;
- covariance or pair-representation error;
- per-residue embedding drift;
- attention-map stability;
- output-logit stability;
- calibration drift;
- prediction stability under conformer order permutation;
- prediction stability under conformer duplication.

If ConFormer embeddings are used for selection, do not use the same embedding
preservation metric as the only success criterion. Include independent
structural and downstream endpoints.

## Downstream Metrics

Preferred endpoint:

- nearest-neighbor antigen retrieval over a defined antigen bank.

Report:

- Recall@1, Recall@5, Recall@10;
- mean reciprocal rank;
- mean average precision;
- candidate-bank size;
- hard-negative construction;
- macro averages across antigen groups;
- performance on unseen sequence families;
- performance on unseen antigen families where possible.

Additional endpoints:

- disease-panel ranking;
- antigen-surface retrieval;
- peptide affinity prediction;
- paratope prediction;
- frozen-model behavior preservation;
- retrained-model behavior at each K.

Closed-set antigen classification should be treated as a debugging or sanity
endpoint unless the split and negative controls are strong.

## Leakage Controls

Required:

- split by antibody target before conformer sampling;
- group paired heavy/light sequence duplicates;
- group near-identical CDR-H3s;
- group clonal or engineered lineages when known;
- group repeated PDB entries of the same antibody;
- control antigen-family similarity for retrieval;
- document overlap with model pretraining sources when possible.

Privileged-conditioning controls:

- compare true-contact, predicted-contact, and random-contact conditioning;
- standardize binder-size rules when testing deployment settings;
- keep privileged teacher outputs as a separately labeled upper bound;
- evaluate final claims using only deployable information.

ConFormer-specific controls:

- permute conformer order;
- replace conformer 0 with different ensemble members;
- duplicate one conformer at 2x, 10x, and 100x;
- hold unique modes fixed while changing multiplicity;
- compare raw conformer-weighted means to cluster-weighted means.

## Compute Metrics

Report:

- PH calls;
- AF3 calls;
- GPU-hours;
- wall-clock time;
- peak memory;
- storage;
- preprocessing and I/O time;
- ConFormer training cost;
- ConFormer inference cost;
- student model training cost;
- student model inference cost;
- break-even number of future targets.

Break-even:

```text
N_break_even = C_train / (C_teacher_per_target - C_student_per_target)
```

## Failure Modes

- full ensemble does not beat sequence/static controls;
- model learns PH/AF3 artifacts;
- conformer-level leakage inflates performance;
- diversity metrics reward broken structures;
- rare useful modes are lost;
- downstream endpoint is too weak or too easy;
- compute savings disappear after relaxation or filtering;
- generated ensembles are overclaimed as physical dynamics;
- thesis scope moves while Gaeun's pipeline changes.
