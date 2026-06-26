# Evaluation Design

## Evaluation Principles

1. Split by target before sampling conformers.
2. Treat conformer files from the same antibody as correlated.
3. Report structural quality, structural coverage, representation preservation,
   downstream performance, and compute cost separately.
4. Define non-inferiority margins before looking at final results.
5. Keep teacher fidelity separate from biological fidelity.

## Baselines

Minimum baselines:

- full PH/AF3 teacher ensemble;
- single static or predicted antibody structure;
- repeated copies of one conformer;
- random subsets of size K;
- stratified subsets of size K by PH design/cycle/AF3 seed;
- quality-score-only selection;
- k-medoids or k-center selection;
- reduced PH design/cycle protocols;
- reduced-seed AF3 without full PH;
- sequence-only model;
- sequence plus one structure;
- shuffled-conformer negative controls;
- fast antibody ensemble model, if deployable, such as ABodyBuilder4-STEROIDS.

Optional baselines:

- BioEmu or AlphaFlow-style general protein ensemble sampling, if compatible
  with antibody chain conventions;
- representation distillation without explicit conformer generation;
- random RMSD-matched perturbations around one anchor.

## Structural Metrics

Geometry and validity:

- chain breaks;
- peptide bond geometry;
- bond and angle outliers;
- Ramachandran outliers;
- chirality;
- cis/trans peptide sanity;
- CDR loop closure;
- disulfide geometry;
- steric clashes;
- side-chain rotamer plausibility;
- VH/VL interface plausibility.

Diversity and coverage:

- framework-aligned CDR-H3 RMSD;
- RMSD for all six CDRs;
- CDR-H3 tip displacement;
- per-residue RMSF;
- backbone torsion distributions;
- pairwise-distance feature diversity;
- VH/VL orientation metrics;
- cluster count and cluster recall;
- teacher-to-subset nearest-neighbor distance;
- subset-to-teacher nearest-neighbor distance;
- maximum mean discrepancy or energy distance in structural feature space.

Important warning: mean pairwise RMSD alone rewards broadness, including
physically bad broadness. Always pair diversity with precision and validity.

## Representation Metrics

Because ConFormer or MCA-style models are expected consumers, evaluate:

- pooled embedding cosine distance;
- mean embedding error;
- covariance or pair-representation error;
- per-residue embedding drift;
- output-logit stability;
- retrieval-score stability;
- attention-map stability where available;
- sensitivity to conformer order and conformer duplication.

Useful controls:

- randomly permute conformer order;
- replace conformer 0 with different ensemble members;
- duplicate one conformer 2x, 10x, and 100x;
- hold unique modes fixed while changing multiplicities;
- compare raw conformer-weighted mean to cluster-weighted mean.

If predictions change when only order or multiplicity changes, some apparent
ensemble-size effects may be architecture artifacts rather than structural
information.

## Downstream Metrics

For antigen retrieval:

- Recall@1, Recall@5, Recall@10;
- mean reciprocal rank;
- mean average precision;
- candidate-bank size;
- macro-averaged performance across antigen groups;
- hard-negative performance;
- calibration and confidence;
- target-family bootstrap confidence intervals.

For disease-panel ranking:

- AUROC and AUPRC;
- macro-averaged metrics by disease;
- patient/cohort split controls where applicable;
- antigen-count and antigen-diversity diagnostics.

Closed-set antigen classification should be treated as a debugging or
sanity-check task unless strict external splits and hard negatives are present.

## Compute Metrics

Report end-to-end cost:

- PH calls;
- AF3 input JSON count;
- AF3 cofold count;
- GPU-hours;
- wall-clock time;
- storage;
- preprocessing time;
- conformer filtering time;
- downstream training and inference time;
- student-model training cost;
- student-model inference cost;
- relaxation or validation cost.

For learned students, report amortization:

```text
break_even_targets =
    training_cost / (teacher_cost_per_target - student_cost_per_target)
```

A fast sampler is not useful if training and validation cost exceed the expected
future savings.

## Leakage Controls

Avoid leakage by grouping:

- all conformers from one antibody target;
- identical paired VH/VL sequences;
- near-identical VH/VL sequences;
- same CDR-H3 sequence and length families;
- clonal relatives where known;
- repeated PDB entries or engineered variants;
- antigen families for retrieval tasks where possible;
- all variants of a shared generation protocol when testing deployment.

Project-specific concern:

- If true antigen contacts, epitope residues, cognate-antigen size, or
  complex-derived templates were used to generate conformers that are later
  evaluated for antigen prediction, label that condition as privileged and keep
  it separate from deployable sequence-only tests.

## Decision Matrix

| Observation | Interpretation | Next step |
|---|---|---|
| Full ensemble does not beat sequence/static controls | Ensemble premise weak or leakage/architecture issue | Do not build generator; diagnose |
| Random small K works | Most conformers redundant | Use reduced fixed generation or random/stratified sampling |
| Coreset works but random does not | Selection matters | Build learned selector or prototype predictor |
| Prospective reduced PH/AF3 works | Compute can be cut directly | Prioritize adaptive generation |
| PH/AF3 beats fast antibody ensemble baseline under blind conditioning | Pseudo-bound proposal may add useful modes | Consider distillation |
| Fast antibody ensemble baseline matches PH/AF3 | Existing model may solve the practical problem | Pivot to benchmark/comparison or hybrid use |
| Generator beats selection baselines | Generative model justified | Expand cautiously |
| Generator loses to simple coreset | Generation not needed yet | Report negative result and use coreset |
