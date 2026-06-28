# Aim 1 Phase 1.4 MCA Diversity-Awareness Interpretation

Date: 2026-06-27

This note records the interpretation after the Phase 1.4 structural CDR
coreset and MCA embedding-sensitivity results. It answers the current question:

> Is diversity awareness relevant to MCA, and should the weak embedding
> sensitivity result change the project direction?

## Bottom Line

Yes, diversity awareness is relevant to MCA, but it should be task-aware rather
than indiscriminate.

The current evidence suggests that Gaeun's PH/AF3 ensemble contains real CDR
structural diversity, especially in H-CDR3, but the current frozen global
MCA/ConFormer readout may underweight or average away much of that diversity.
That is a useful warning sign. It does not prove that MCA is broken, and it
does not prove that Gaeun's conformers are useless. It means the next scientific
question should be:

> Where does CDR conformational signal disappear: the PH/AF3 teacher ensemble,
> the MCA encoder, the pooling/readout, or the downstream retrieval head?

The practical project direction should become:

> Make conformer-aware antibody representations preserve and use
> task-relevant CDR ensemble diversity, instead of merely generating or
> averaging more conformers.

## Observed Evidence

Structural evidence:

- The copied Gaeun PH/AF3 medium set has 149 targets and 19,072 conformers.
- CDR-aware structural coreset analysis showed real compressible CDR diversity.
- H-CDR3 is the main difficult region.
- K=32 and K=64 are the first useful budget range for CDR-aware compression.
- K=8 is too aggressive for reliable H-CDR3 coverage.

Result note:

```text
notes/aim1-phase1-4-cdr-coreset-results.md
```

Embedding evidence:

- H-CDR3 residue embeddings weakly track H-CDR3 structural distance.
- Full H/L pooled embeddings are essentially flat with respect to H-CDR3
  pairwise structural distance.
- K=32 and K=64 subsets preserve average full H/L embeddings well, but
  absolute differences and gains over random are modest.
- One-conformer repeats and duplicate-heavy ensembles move the mean embedding
  much more than K=32/K=64 subset choice, so count/duplicate density matters.
- The richer readout audit found that flattened local H-CDR3 residue embeddings
  recover a stronger but still modest structural signal, while mean-pooled CDR
  and full H/L readouts mostly average it away.

Result note:

```text
notes/aim1-phase1-4-embedding-sensitivity-results.md
notes/aim1-phase1-4-readout-audit-results.md
```

## Why The Result Is Suspicious

It is suspicious because the expensive conformer-generation premise is:

> Multiple generated conformers expose binding-relevant antibody geometries.

If the downstream representation collapses those conformers into a mostly flat
global H/L embedding, then more conformer generation may not improve the
representation or retrieval endpoint. In that case, optimizing generation alone
could miss the bottleneck.

The suspicious object is not the existence of the conformers. The suspicious
object is the current path from ensemble to downstream signal.

## Why The Result Is Not Conclusive

The Phase 1.4 embedding analysis mostly tested simple pooled vectors:

```text
h_global
l_global
hl_concat
h_cdr3
all_cdr_concat
```

The analysis script builds these from mean-pooled residue embeddings:

```text
scripts/phase14/phase14_analyze_embeddings.py
```

That means the result does not prove that all MCA internals are insensitive to
CDR diversity. It proves that the tested global and CDR-mean embedding
distances are weakly sensitive.

MCA may still encode useful ensemble information in:

- per-residue CDR embeddings before global pooling;
- final pair representations;
- layer-wise pair representations;
- conformer attention weights;
- antigen cross-attention outputs;
- downstream retrieval logits or ranks;
- attention-pooling variants rather than arithmetic means.

## Why Diversity Awareness Is Relevant To MCA

MCA is structurally an ensemble model. Its core tensors carry a conformer axis:

```text
mca_repr: [N_conformers, L_residues, D]
pair_repr: [L_residues, L_residues, D_pair]
```

The MCA core also records per-layer conformer weights:

```text
MCA/model.py
```

The outer-product module can compute conformer attention weights, using the
mean conformer representation as a query against all conformers:

```text
MCA/outer_product.py
```

So the architecture is not conceptually blind to ensembles. The issue is that
some downstream heads and extraction paths collapse the conformer axis with
simple mean pooling.

Examples:

- `MCAMeanClassifierHead` averages over residues per conformer, then averages
  over conformers.
- `train_cross_attn_v2.py` defaults to `mca_repr.mean(dim=0)` before antibody
  CDR features interact with antigen surface features.
- Pair and residue heads often aggregate pair or conformer information by
  mean-like operations.

There are also attention-pooling options, such as `ConformerAttnPool` and
`MCAMeanVarAttnClassifierHead`, but Phase 1.4 has not yet shown whether those
recover stronger CDR diversity signal.

## What Diversity Awareness Should Mean

Diversity awareness should not mean "the model reacts to every structural
perturbation." That would be bad. MCA should ignore irrelevant noise.

Useful diversity awareness means preserving ensemble information likely to
matter for antibody-antigen recognition, especially:

- distinct H-CDR3 structural modes;
- all-CDR shape and exposure changes;
- rare conformers that create plausible antigen-compatible geometry;
- VH/VL orientation changes that alter the paratope surface;
- cluster coverage across meaningful CDR modes;
- model uncertainty or conformer weights when they correspond to useful modes.

It should ideally ignore:

- tiny framework jitter;
- AF3/PH duplicate-density artifacts;
- conformer ordering;
- redundant near-identical structures;
- pipeline-specific sampling frequency when frequency is not biologically or
  task calibrated;
- structural variation far from the paratope.

Therefore the right phrase is not generic "diversity awareness." The better
phrase is:

> CDR- and paratope-aware ensemble representation.

## Current Working Hypothesis

Observed:

- The teacher ensemble contains real CDR diversity.
- The tested global frozen MCA embedding weakly sees that diversity.
- MCA has internal components that could, in principle, represent ensemble
  information.
- Several practical readout paths use arithmetic mean pooling.

Hypothesis:

> The current bottleneck may be pooling/readout or training objective, rather
> than conformer generation itself.

Alternative hypotheses to keep open:

1. The MCA encoder does encode CDR diversity, but the tested global embedding
   hides it.
2. The encoder encodes some CDR diversity, but not in a task-useful way.
3. The antigen retrieval head uses signal that the embedding-sensitivity metric
   missed.
4. The PH/AF3 ensemble diversity is structurally real but not useful for the
   downstream endpoint.
5. The current checkpoint or extraction path is not the best representative of
   the intended MCA/ConFormer use case.

## Next Tests

### 1. Retrieval Preservation

Use the same target set and compare:

- full 128-conformer ensemble;
- K=32 H-CDR3 k-center;
- K=64 H-CDR3 k-center;
- K=32 all-CDR k-center;
- K=64 all-CDR k-center;
- random K controls;
- first-K controls;
- one-conformer repeated controls;
- duplicate-heavy controls.

Endpoint:

- nearest-neighbor antigen retrieval, once the antigen bank and leakage controls
  are frozen.

Interpretation:

- If retrieval is also flat, the current downstream model likely does not use
  conformer diversity strongly.
- If retrieval is sensitive while global embeddings are flat, useful signal may
  live in the retrieval head, attention maps, pair features, or local CDR
  residues rather than the global mean embedding.

### 2. MCA Readout Audit

Measure the same structural subset perturbations against richer MCA outputs:

- CDR-only per-residue embeddings;
- `final_pair_repr`;
- CDR-CDR and CDR-framework subblocks of `pair_repr`;
- `conf_weights_list_dict` from each MCA layer;
- attention-pooling outputs;
- mean plus variance/covariance/moment features;
- cluster-weighted means rather than raw conformer-weighted means.

Interpretation:

- If richer readouts track H-CDR3 structure, the encoder may already contain
  signal and the issue is the readout.
- If no readout tracks CDR structure, then the model objective may need to be
  modified.

### 3. Pooling / Adapter Prototype

If the readout audit points to pooling as the bottleneck, test small changes
before retraining the full model:

- CDR-specific pooling instead of full-chain pooling;
- H-CDR3/all-CDR mean plus variance features;
- learned conformer attention over CDR features;
- cluster-balanced pooling to reduce duplicate-density bias;
- pair-representation pooling over CDR/CDR and paratope-relevant blocks;
- a small adapter on top of frozen MCA representations.

### 4. Training-Signal Test

If pooling changes are insufficient, consider a diversity-preserving objective:

- contrastive objective where structurally distinct H-CDR3 conformers are
  separable locally;
- metric-learning objective on CDR structural distance;
- auxiliary loss to predict CDR cluster identity or CDR coverage;
- task-aware loss that rewards preserving retrieval rankings under subset
  compression;
- duplicate-normalized or cluster-weighted ensemble objective.

This should be a later step, not the immediate next action.

## Decision Logic

If attention/CDR/pair readouts recover structural sensitivity:

- focus on readout and retrieval-head design;
- use existing conformers and frozen MCA features first;
- avoid claiming the whole encoder failed.

If retrieval uses diversity even though global embeddings are flat:

- preserve the retrieval-sensitive readout;
- use CDR coresets as a compression baseline;
- evaluate K=32/K=64 against full ensemble.

If all readouts and retrieval are flat:

- current generation may be oversupplying information the model does not use;
- generation-cost optimization may be simple subsampling;
- diversity-aware model training becomes a stronger thesis direction.

If only duplicate/count controls move the embedding:

- raw PH/AF3 sampling frequency is a major confound;
- future comparisons need cluster-balanced or duplicate-normalized pooling.

## How To Say This To Gaeun

Concise technical framing:

> The copied PH/AF3 ensemble has real H-CDR3 structural diversity and is
> structurally compressible around K=32 to K=64. However, the current frozen
> global MCA H/L embedding is nearly flat with respect to H-CDR3 structural
> distance. This suggests the current endpoint may be averaging away some CDR
> ensemble information. I think the next check should be whether retrieval,
> pair representations, conformer attention weights, or CDR-specific readouts
> recover the signal before we claim conformer-generation optimization is
> task-preserving.

What not to overclaim:

- Do not say MCA is blind to conformers.
- Do not say Gaeun's conformers are useless.
- Do not say structural diversity implies antigen specificity.
- Do not say diffusion/generation is now justified.

## Updated Thesis Framing

The strongest current Jerry-owned thesis angle is:

> What information in a large PH/AF3 antibody conformer ensemble is actually
> used by conformer-aware antibody representation models, and can we preserve
> that information with a smaller CDR/paratope-aware ensemble?

This supports Gaeun's conformer-generation optimization more directly than a
generic antibody diffusion model.
