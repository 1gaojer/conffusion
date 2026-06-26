# My judgment

**Yes, this can become a serious thesis—but the strongest thesis is not “train Stable Diffusion for antibodies.”** The scientifically sharp version is:

> **Quantify and remove redundancy from expensive pseudo-bound antibody conformational ensembles, then optionally distill the remaining distribution into a cheaper, CDR-focused conditional sampler.**

I would separate the possible projects this way:

| Project | Judgment |
|---|---|
| Determine the minimum task-sufficient ensemble size | **Strong, feasible, and immediately useful** |
| Develop a task-aware or coverage-aware conformer coreset | **Strong thesis contribution** |
| Learn to allocate Protein Hunter/AF3 sampling budgets prospectively | **Arguably the best compute-saving thesis** |
| Fine-tune a pretrained flow/diffusion model to generate CDR conformers | **Ambitious but defensible as an extension** |
| Train a full-atom antibody ensemble generator from scratch on a few hundred targets | **Too risky for a master’s thesis** |
| Claim the generated ensembles represent true antibody equilibrium dynamics | **Not supported by the current data** |
| Claim ensembles contain new antigen information unavailable from sequence | **Conceptually wrong unless antigen-derived conditioning is used—in which case leakage becomes a concern** |

The April workflow is already a natural distillation target: 25 Protein Hunter designs over 5 cycles produce 125 complexes, followed by 5 AF3 seeds × 5 models for each complex—nominally 3,125 structures per antibody before filtering. fileciteturn0file1 Your reported tree of roughly 1.5 million retained PDBs is therefore not merely “a large dataset”; it is a highly structured, hierarchical Monte Carlo output whose redundancy should be measured before building a new generator.

There is also an unusually favorable detail in the current downstream architecture: the MCA/ConFormer method ultimately mean-pools the conformer dimension for the antibody–antigen cross-attention module, while earlier layers summarize conformer information partly through outer-product statistics. That means preserving the **mean and selected second-order statistics of learned conformer embeddings** may matter more than retaining thousands of raw structures. fileciteturn0file2

---

# The most important conceptual correction

## The ensemble is currently a computational proposal distribution, not a physical distribution

A useful formalization is:

\[
p_{\text{teacher}}(x\mid s)
=
\sum_m p(m\mid s)\,
p_{\text{AF3}}(x\mid s,m)\,
\mathbf{1}\{\text{passes filter}\}/Z,
\]

where:

- \(s\) is the heavy–light antibody sequence;
- \(m\) includes pseudo-binder design, Protein Hunter cycle, contact conditioning, AF3 seed, AF3 model, template, MSA and other protocol choices;
- \(x\) is the resulting antibody conformation.

This is fundamentally different from

\[
p_{\text{physical}}(x\mid s,T,\text{solvent},\text{pH},\text{ligand state},\ldots),
\]

the equilibrium or kinetic conformational distribution one might try to estimate with molecular dynamics or experiments.

The current samples do not have known Boltzmann weights. Twenty-five structures from one PH design do not mean that basin is 25 times more physically occupied than a basin represented by one structure. Cluster populations reflect the sampling protocol as much as antibody thermodynamics.

So the defensible claim is:

> “We model or compress the support and downstream utility of a PH/AF3-generated pseudo-bound conformational proposal ensemble.”

The indefensible claim is:

> “We learn the true conformational distribution of antibodies.”

BioEmu is conceptually closer to the latter because it was explicitly developed to emulate equilibrium ensembles using large-scale MD, static structures and experimental stability information. AlphaFlow and ESMFlow similarly target sequence-conditioned structural ensembles and can be trained on MD-derived distributions. ([science.org](https://www.science.org/doi/pdf/10.1126/science.adv9817?utm_source=chatgpt.com))

Your internal observation that BioEmu explored a narrower CDR-H3 space than PH on several selected antibody groups is interesting, but it does not establish that PH is more physically correct. PH may be recovering useful bound-like modes, or it may simply be broader because pseudo-binders and cofolding force the antibody into many plausible-looking states. The February deck itself shows the important warning: multi-seed AF3 can be very diverse while also lying farther from experimentally observed structure space. fileciteturn0file0

**Diversity is not the same as realism.**

---

## A second conceptual correction: generated ensembles cannot create antigen information from nothing

Suppose antigen identity is \(Y\), antibody sequence is \(S\), and the ensemble is sampled only from the sequence:

\[
E\sim q(E\mid S).
\]

Then the causal structure is approximately

\[
Y \longleftrightarrow S \longrightarrow E.
\]

Conditional on the sequence, the generated ensemble contains no new antigen information:

\[
I(Y;E\mid S)=0.
\]

A geometry-aware representation can still outperform a particular sequence-only model because it provides a useful inductive bias, pretraining signal, regularization or model architecture. But the ensemble has not magically added information that was absent from the sequence.

This matters for interpreting the early antigen-classification results. The February experiment used only 97 PDBs across 20 antigens and reported 1.0 accuracy for MCA attention versus 0.83 for a linear ESM2 baseline. That is encouraging as a pipeline sanity check, but it is not yet evidence that conformational variability carries generalizable antigen-specific signal. fileciteturn0file0 The sequence baseline, architecture, training objective and capacity were not matched, and the task was small enough for family, study, germline or structural provenance shortcuts to dominate.

There is an even more serious issue to resolve. The optimized pipeline states that structurally resolved complexes can use the top three contact residues, while unknown structures use random CDR-H3 contacts; it also says design size may be proportional to a ground-truth binder size. fileciteturn0file1

**You need to determine whether any antigen-derived information was used to generate ensembles that were later used to predict that antigen.**

If true paratope contacts, cognate-antigen size or complex-derived templates enter the generator, then the ensemble is not a sequence-only representation. That may be valid for another task, but it risks direct or indirect label leakage in antigen prediction.

This should be the first question asked before committing to the thesis.

---

# 1. Are antibody conformer ensembles a natural diffusion, flow or transformer target?

**Yes, mathematically they are a natural target.**

The conditional distribution is highly multimodal:

\[
p(x\mid \text{heavy sequence},\text{light sequence},\text{context}),
\]

where different modes may correspond to:

- CDR-H3 loop basins;
- smaller CDR loop rearrangements;
- VH–VL orientation changes;
- framework breathing;
- alternative side-chain rotamers;
- apo-like versus pseudo-bound-like states.

Diffusion and flow matching are suitable because they learn continuous distributions over structure and can generate multiple samples from noise. SE(3)-equivariant models are especially appropriate because global rotation and translation must not affect the answer.

However, I would not initially model the whole antibody as unconstrained Cartesian coordinates. Antibodies have a strong structural prior: the framework is comparatively stable, while much of the useful variation is concentrated in the CDRs and VH–VL orientation.

A more realistic decomposition is:

\[
x =
\left(
x_{\text{anchor framework}},
T_{\text{VH--VL}},
\{\text{CDR residue frames}\},
\{\chi\text{ torsions}\}
\right).
\]

A practical student model could:

1. take the paired sequence and one cheap anchor structure;
2. keep most framework residues fixed or softly restrained;
3. generate CDR backbone frames and VH–VL orientation;
4. reconstruct side chains using sequence-conditioned torsion prediction or a rotamer packer;
5. apply a short geometry repair or relaxation step.

This reduces the learning problem dramatically compared with full-atom generation from unstructured noise.

### Diffusion versus flow matching

Both are plausible. Flow matching may be attractive because inference can use a relatively small number of integration steps. FrameDiff established principled diffusion over residue frames on \(SE(3)\), while FrameFlow adapted that idea to flow matching and reported substantially fewer sampling steps. RFdiffusion demonstrated how a pretrained structure model can be repurposed into a generative denoiser, although its primary purpose is de novo design rather than conformational ensemble prediction. ([nature.com](https://www.nature.com/articles/s41586-023-06415-8?utm_source=chatgpt.com))

### Transformers are not a competing category

“Transformer,” “diffusion” and “flow matching” are not mutually exclusive terms.

- A transformer can be the denoising network inside a diffusion model.
- A transformer can predict the vector field inside a flow model.
- A transformer can instead generate discrete structure tokens, as in ESM3-like masked-token generation.

ESM3 is relevant because it jointly reasons over discrete sequence, structure and function tracks. That makes it interesting for latent conformer representations or coarse structural generation, but its structure-token representation is not automatically calibrated for subtle sub-angstrom CDR differences or ensemble occupancies. ([science.org](https://www.science.org/doi/10.1126/science.ads0018?utm_source=chatgpt.com))

---

# 2. How this differs from image/video diffusion and LoRA

The analogy is useful, but only at the highest level.

## What transfers from image generation

The transferable ideas are:

- conditional generation from noise;
- latent-variable control of multimodal outputs;
- classifier-free or property guidance;
- distillation into faster samplers;
- parameter-efficient adaptation;
- diversity versus fidelity trade-offs;
- latent compression.

## What does not transfer cleanly

### Protein geometry has hard constraints

A slightly malformed hand in an image may still look acceptable. A 1.8 Å peptide bond, broken CDR loop, steric catastrophe or inverted chirality is structurally wrong.

A protein generator must respect:

- covalent connectivity;
- bond lengths and angles;
- stereochemistry;
- chain identity and residue correspondence;
- loop closure;
- disulfide bonds;
- clash avoidance;
- global SE(3) symmetry.

### Your conformers are an unordered set, not a video

A video has temporal ordering and local continuity between frames. The PH/AF3 structures are an ensemble of samples; “cycle 5” is not a later physical time point than “cycle 2.”

A video diffusion model would therefore introduce a false temporal interpretation unless you specifically wanted to model the computational pipeline trajectory. For downstream ensemble representation, the output should be permutation-invariant over conformers.

### The effective dataset size is much smaller than the file count

Image models may see millions of semantically distinct images. Here, 1.5 million structures may represent only several hundred to a thousand paired sequences, with thousands of correlated predictions per sequence.

A good analogy is not “1.5 million independent training images.” It is closer to “1.5 million highly correlated views of 1,000 identities.”

### LoRA is not the generative idea

LoRA is a parameter-efficient way to modify a pretrained network. It works well for image subjects or styles because the base model already knows the broad image manifold.

A protein LoRA would be credible only if the base model already knows:

- protein geometry;
- antibody frameworks;
- CDR loop priors;
- sequence–structure compatibility;
- ideally, conformational ensemble sampling.

Otherwise a small adapter may simply learn the visual equivalent of a “PH/AF3 style”—pipeline artifacts—rather than meaningful antibody flexibility.

I would therefore describe LoRA or adapters as an **engineering strategy for fine-tuning AlphaFlow-, FrameFlow-, BioEmu- or antibody-pretrained models**, not as the scientific thesis itself.

---

# 3. Distilling the pipeline versus learning biophysics

These are three distinct objectives, and the lab should choose explicitly.

## Objective A: teacher-distribution fidelity

“Generate samples statistically resembling the retained PH/AF3 ensemble.”

This is the cleanest distillation problem. The student should reproduce:

- teacher cluster support;
- teacher embedding statistics;
- teacher CDR variability;
- teacher structural quality;
- downstream behavior.

The limitation is that it faithfully reproduces teacher errors and artifacts.

## Objective B: geometric support coverage

“Cover all structurally distinct basins the teacher found, regardless of how frequently it sampled them.”

This may be more useful for retrieval because rare conformers may expose antigen-compatible geometries. It requires cluster-balanced or diversity-weighted training rather than treating raw teacher frequencies as probabilities.

## Objective C: downstream task utility

“Generate or select the structures that preserve antigen retrieval or disease-panel performance.”

This is the most practical objective and may produce a distribution quite different from the teacher distribution.

A task-aware sampler could intentionally overrepresent rare but useful conformers. That is legitimate, but it should be described as **task-optimized proposal generation**, not physical ensemble modeling.

## Improving beyond the teacher

A student trained only to imitate the teacher cannot establish that it is more biophysically correct than the teacher.

It can improve in three narrower senses:

1. **Efficiency:** similar outputs at lower cost.
2. **Regularization:** fewer duplicates and less protocol noise.
3. **Task performance:** optimized with an independent downstream objective.

To claim physical improvement, you need independent supervision such as:

- repeated experimental structures of the same or highly similar sequence;
- bound and unbound antibody pairs;
- ALL-conformations/ITsFlexible groups;
- antibody-focused MD;
- NMR, HDX, SAXS or other ensemble observables;
- prospective experimental validation.

ALL-conformations is useful because it explicitly aggregates observed loop conformations and includes more than 1.2 million loop structures across over 100,000 unique sequences, although crystal structures still have selection and packing biases. ([nature.com](https://www.nature.com/articles/s42256-025-01131-6?utm_source=chatgpt.com))

---

# 4. Relevant model families

## RFdiffusion, FrameDiff and FrameFlow

These establish the core geometric technology:

- **RFdiffusion**: denoises protein structures and supports conditional de novo design.
- **FrameDiff**: diffusion over residue rigid frames with SE(3)-equivariant score modeling.
- **FrameFlow/FoldFlow-style models**: flow matching over structural frames, often with faster sampling.

They are excellent architectural references, but they are primarily trained to generate new backbones, not distributions of conformers for one fixed antibody sequence. ([nature.com](https://www.nature.com/articles/s41586-023-06415-8?utm_source=chatgpt.com))

## AlphaFlow and ESMFlow

These are probably the closest conceptual precedents.

They convert AlphaFold- and ESMFold-like predictors into sequence-conditioned ensemble generators using flow matching. The work explicitly evaluates structural diversity, PDB multistate behavior and MD-derived ensemble observables. ([proceedings.mlr.press](https://proceedings.mlr.press/v235/jing24a.html?utm_source=chatgpt.com))

For your project, an AlphaFlow-like formulation would be:

> fine-tune a pretrained structure predictor so that noisy antibody conformers are mapped along a learned flow toward the PH/AF3 teacher ensemble.

The main practical concern is computational weight. A lightweight or structure-module-only variant would be more realistic than repeatedly executing a full AlphaFold trunk; AlphaFlow-Lit is a relevant efficiency precedent. ([arxiv.org](https://arxiv.org/abs/2407.12053?utm_source=chatgpt.com))

## BioEmu

BioEmu is the strongest reference for the claim that a generative model can amortize expensive protein ensemble generation. It integrates large MD data, static structures and thermodynamic supervision to generate sequence-conditioned equilibrium-like samples efficiently. ([science.org](https://www.science.org/doi/pdf/10.1126/science.adv9817?utm_source=chatgpt.com))

It is a useful baseline and potentially a pretrained initialization. But antibodies, CDR loops and pseudo-bound states may be outside the distributions it models best, and an unbound equilibrium model is not automatically expected to reproduce the PH pseudo-binding protocol.

## ESM3-like structure-token models

These are attractive for:

- compressing a conformer into a sequence of discrete structural tokens;
- predicting latent conformer clusters;
- masked structure completion;
- generating coarse structural alternatives.

Their weakness is resolution and calibration. A token model may distinguish broad folds while merging small but functionally important CDR rearrangements.

## DiffAb

DiffAb jointly generates CDR sequence and structure conditioned on antigen and antibody context. It is directly relevant for its equivariant representation of CDR geometry, but it solves an **antibody design** problem rather than a **fixed-sequence conformational ensemble** problem. ([proceedings.neurips.cc](https://proceedings.neurips.cc/paper_files/paper/2022/hash/3fa7d76a0dc1179f1e98d1bc62403756-Abstract-Conference.html?utm_source=chatgpt.com))

## IgDiff

IgDiff adapts FrameDiff to antibody variable-domain backbone generation and multiple chains. This is highly relevant as an antibody structural prior and as evidence that a general geometric diffusion model can be fine-tuned to the antibody manifold. It still targets de novo backbone generation rather than reproducing conditional conformers of one fixed sequence. ([arxiv.org](https://arxiv.org/abs/2405.07622?utm_source=chatgpt.com))

## AbDiffuser

AbDiffuser jointly generates full-atom antibody structures and sequences using antibody alignment and physics-informed priors. It is relevant for side-chain and full-atom representations, but most of its complexity is unnecessary here because your sequence is already fixed. ([proceedings.neurips.cc](https://proceedings.neurips.cc/paper_files/paper/2023/file/801ec05b0aae9fcd2ef35c168bd538e0-Paper-Conference.pdf?utm_source=chatgpt.com))

## Antibody-specific RFdiffusion work

The recent RFdiffusion antibody-design work shows that an RFdiffusion-style model can be specialized to generate epitope-directed antibody loops while preserving framework context. This supports the feasibility of antibody-specific geometric fine-tuning, but again the task is de novo design, not calibrated conformational sampling. ([nature.com](https://www.nature.com/articles/s41586-025-09721-5?utm_source=chatgpt.com))

## Protein Hunter and BindCraft

Protein Hunter and BindCraft are best viewed as parts of the **teacher-generation mechanism**. BindCraft uses structure prediction weights for binder hallucination, while Protein Hunter iteratively cycles structure hallucination and sequence redesign. Neither is designed to produce thermodynamic antibody ensembles. ([biorxiv.org](https://www.biorxiv.org/content/10.1101/2025.10.10.681530v2.full.pdf?utm_source=chatgpt.com))

### Practical conclusion from the model landscape

The most promising starting points are:

- **AlphaFlow/ESMFlow or BioEmu** for ensemble-generation philosophy;
- **FrameDiff/FrameFlow** for geometric sampling machinery;
- **IgDiff or AbDiffuser** for antibody-specific structural priors;
- **ConFormer** for the task-relevant representation and compression objective.

I would not start from DiffAb unless antigen structures are explicitly part of the conditioning task.

---

# 5. Is the available data enough?

## For ensemble compression: unquestionably yes

Hundreds to thousands of targets and millions of structures are more than sufficient to characterize redundancy, cluster structure, variance sources and downstream saturation.

## For a learned selector or budget allocator: probably yes

A model that predicts which PH design, cycle or AF3 seed is likely to contribute a novel conformer has a much simpler target than full structure generation.

You have a hierarchical dataset:

\[
\text{antibody}
\rightarrow
\text{PH design}
\rightarrow
\text{cycle}
\rightarrow
\text{AF3 seed}
\rightarrow
\text{AF3 model}.
\]

This allows a variance-decomposition analysis:

- How much conformational variation comes from changing PH design?
- How much comes from cycle number?
- How much comes from AF3 seed?
- How much comes from the five models within a seed?
- Which level mostly produces duplicates?

That analysis alone could reveal that, for example, 25 models per design are unnecessary while design diversity is critical—or the reverse. It would directly inform a cheaper generation protocol.

## For CDR-focused generative fine-tuning: plausibly yes

A few hundred targets would be a pilot-scale dataset. Roughly 1,000–3,000 genuinely distinct paired-sequence families could plausibly support antibody-specific fine-tuning of a pretrained flow or diffusion model, especially if the framework is anchored and the output is restricted to CDR frames/torsions.

## For a full-atom model from scratch: probably no

The relevant sample size for sequence-level generalization is the number of independent antibody families, not the number of PDB files.

One thousand sequences with 1,500 correlated conformers each do not provide the same information as 1.5 million distinct sequence–structure examples.

A from-scratch model could easily:

- memorize target identities;
- learn a generic “AF3 perturbation style”;
- fail on unseen CDR-H3 lengths;
- fail on unseen V-gene or clonotype families;
- reproduce protocol metadata rather than antibody mechanics.

## You probably should not train on all 1.5 million structures

First cluster or deduplicate within target. Then use:

- uniform sampling over targets;
- uniform or controlled sampling over structural clusters;
- optional sampling proportional to teacher frequency only when teacher-frequency matching is the explicit goal.

Otherwise targets with higher PH/AF3 success rates dominate training, and dense duplicate modes dominate rare structural basins.

A realistic first generative dataset might contain 32–128 representative conformers per target, producing perhaps \(3\times10^4\) to \(2\times10^5\) meaningful training examples rather than 1.5 million near-duplicates.

---

## Leakage must be controlled at the target level

Never split conformers from one antibody between train and test.

The outer split should group by:

- paired heavy–light sequence;
- clonal lineage where available;
- heavy-chain CDR3 identity and near-identity;
- light-chain identity;
- V/J genes;
- highly similar variable-domain sequences;
- repeated PDB depositions of the same antibody;
- ideally antigen family or deposition time for harder evaluation.

The April deck’s no-V-gene-overlap splits are substantially better than random PDB splits. fileciteturn0file1 But no-V overlap alone does not eliminate:

- similar CDR-H3s;
- shared light chains;
- clonal relatives;
- study-specific artifacts;
- antigen-family shortcuts;
- model-pretraining overlap;
- complex-derived contact leakage.

Recent antibody-design benchmarking is also moving toward unseen-epitope, unseen-antigen-fold and temporal splits rather than random structure splits, reflecting the same concern. ([arxiv.org](https://arxiv.org/abs/2603.13431?utm_source=chatgpt.com))

---

# 6. Which formulation is most realistic?

| Formulation | Feasibility | Scientific value | Main limitation |
|---|---:|---:|---|
| Post-hoc random/cluster subsampling | Very high | Necessary baseline | Does not save original generation cost |
| Task-aware conformer coreset | Very high | High | Can overfit the chosen downstream model |
| Prospective PH/AF3 budget allocation | High | Very high | Requires detailed pipeline metadata and new runs |
| Direct sequence → ensemble embedding distillation | High | High for downstream use | Produces no explicit structures |
| Generate cluster medoids in latent space | Moderate | Moderate–high | Decoder/geometry fidelity remains difficult |
| Generate CDR frames/torsions around an anchor | Moderate | High | Requires geometric model engineering |
| Generate the entire variable domain backbone | Low–moderate | High if successful | Data and compute risk |
| Full-atom sequence-conditioned physical ensemble | Low | Potentially very high | Thesis-scale mismatch and inadequate physical labels |

## My preferred ordering

### First: task-aware coresets

Let \(E_i\) be the complete ensemble for antibody \(i\), and \(S_i^k\subset E_i\) a subset of size \(k\).

Choose \(S_i^k\) to preserve some combination of:

- CDR structural coverage;
- VH–VL orientation coverage;
- ConFormer embedding mean;
- ConFormer embedding covariance;
- downstream predictions.

This is especially well matched to the current architecture because conformer means and outer-product-like statistics are central to ConFormer.

Candidate methods include:

- k-medoids;
- farthest-point sampling;
- facility-location objectives;
- determinantal point processes;
- kernel herding;
- moment matching;
- weighted coresets;
- gradient- or attention-informed selection.

### Second: prospective adaptive sampling

Post-hoc compression reduces storage and downstream cost, but it does not recover the expensive AF3 compute already spent.

The more important system is:

> Given the conformers generated so far, decide whether to stop, run another AF3 seed, sample another PH design or allocate compute to a different target.

This could be framed as:

- sequential experimental design;
- submodular budget allocation;
- uncertainty-guided sampling;
- a multi-armed bandit over PH designs;
- learned prediction of marginal cluster coverage.

A particularly practical protocol would run one AF3 model per PH design first, estimate which designs are generating new modes, and allocate additional seeds/models only to promising or undercovered designs.

### Third: distill the ensemble representation rather than the structures

If the true goal is antigen retrieval, train a student:

\[
g_\theta(\text{heavy sequence},\text{light sequence})
\rightarrow
\text{ConFormer pooled representation}.
\]

This could remove conformer generation entirely at deployment.

It would answer a different scientific question—whether the expensive ensemble acts as a useful training-time teacher—but it may be the highest-return engineering outcome.

### Fourth: CDR-focused generative flow

For explicit structural generation, I would recommend:

- one anchor variable-domain structure;
- generated CDR backbone frames;
- generated VH–VL transform;
- optional side-chain \(\chi\) torsions;
- fixed sequence;
- cluster-balanced teacher targets;
- geometry repair after sampling.

I would not initially generate the full framework or amino-acid identities.

---

# 7. Minimum viable thesis

A strong minimum viable thesis does **not** require training a diffusion model.

## Core research question

> What is the smallest and cheapest antibody conformational ensemble that is statistically non-inferior to the complete PH/AF3 ensemble for structural coverage and downstream antigen-representation tasks?

## Aim 1: quantify where the diversity comes from

Perform hierarchical variance decomposition over:

- PH design;
- cycle;
- contact-conditioning choice;
- AF3 seed;
- AF3 model;
- filter status.

Use invariant structural features such as:

- per-CDR RMSD;
- backbone torsions;
- pairwise distances;
- VH–VL orientation;
- ConFormer embeddings.

This identifies which pipeline dimensions can be cut.

## Aim 2: establish ensemble-size scaling curves

Evaluate:

\[
k\in \{1,2,4,8,16,32,64,128,256,512,\text{full}\}.
\]

For every \(k\), use repeated random subsets and report uncertainty over targets and subset draws.

The key result is not “performance at \(k=64\).” It is a saturation curve showing where additional conformers stop helping.

## Aim 3: show informed selection beats random selection

Compare random, stratified random, cluster medoids, geometric diversity selection and representation-aware selection.

## Aim 4: prospective confirmation

On a smaller set of untouched antibodies, run genuinely reduced generation protocols:

- fewer PH designs;
- fewer cycles;
- fewer AF3 seeds;
- fewer AF3 models;
- adaptive stopping.

This is essential because a post-hoc subset does not prove that the same subset could have been obtained without generating the full ensemble.

## Optional Aim 5: lightweight generative or representation distillation

Only after Aims 1–4 are stable:

- fine-tune a CDR-focused flow model; or
- predict cluster medoids/latent conformer codes; or
- distill the full ensemble embedding directly from sequence.

That is enough for a coherent, publishable master’s thesis even if the generative extension is only a pilot.

---

# 8. Essential baselines and evaluations

## Generation and selection baselines

At minimum:

1. **Uniform random subsampling**, repeated many times.
2. **Stratified random sampling** across PH design, cycle and AF3 seed.
3. **One conformer per PH design**.
4. **Top-confidence AF3 structures** using available confidence measures.
5. **k-medoids or farthest-point sampling**.
6. **Reduced-seed AF3** without Protein Hunter.
7. **Reduced PH designs/cycles plus AF3**.
8. **PH/Boltz structures before AF3 reprediction**, where appropriate.
9. **Single predicted structure**.
10. **A duplicated-single-conformer ensemble** of size \(N\).

The duplicated-conformer baseline is important: if repeating one conformer 1,000 times performs like 1,000 distinct conformers, the model is not using variability.

Additional negative controls should include:

- geometry removed but sequence retained;
- conformers shuffled between sequence-similar antibodies;
- CDR-only versus framework-only geometry;
- mean structure plus random coordinate perturbations;
- a matched-capacity sequence-only model;
- sequence plus one structure;
- sequence plus confidence scores.

## Geometry and physical validity

RMSD alone is insufficient. Measure:

- chain breaks;
- peptide bond geometry;
- bond and angle deviations;
- Ramachandran outliers;
- CDR loop closure;
- steric clashes;
- disulfide geometry;
- VH–VL interface plausibility;
- side-chain rotamer outliers;
- radius of gyration and abnormal collapse;
- AF3 confidence distributions, with the caveat that they are teacher-model diagnostics.

## Structural diversity and coverage

Report separately:

- CDR-H3 RMSD;
- CDR-L3 RMSD;
- other CDR RMSDs;
- per-residue RMSF;
- torsion distributions;
- VH–VL orientation metrics;
- pairwise-distance feature diversity;
- cluster count and cluster recall;
- full-ensemble-to-subset nearest-neighbor distance;
- subset-to-full precision;
- maximum mean discrepancy or energy distance;
- coverage against experimental multistate groups.

For a generated ensemble rather than a subset, use both directions:

\[
\text{teacher}\rightarrow\text{generated}
\]

to measure mode recall, and

\[
\text{generated}\rightarrow\text{teacher}
\]

to measure precision or hallucination.

## Representation preservation

Because ConFormer is the intended consumer, measure:

- pooled embedding cosine distance;
- mean-embedding error;
- covariance or pair-representation error;
- per-residue embedding drift;
- attention-map stability;
- output-logit stability;
- calibration drift.

A very natural coreset objective is:

\[
\left\|
\mu_{\text{full}}-\mu_{\text{subset}}
\right\|^2
+
\lambda
\left\|
\Sigma_{\text{full}}-\Sigma_{\text{subset}}
\right\|_F^2,
\]

where \(\mu\) and \(\Sigma\) are computed from frozen conformer embeddings.

## Downstream evaluation

Use both:

1. **Frozen-model evaluation:** does a subset preserve the behavior of the already trained model?
2. **Retrained-model evaluation:** is a small ensemble sufficient when the model is trained appropriately for that ensemble size?

Relevant tasks include:

- antigen classification;
- antigen retrieval;
- top-\(k\) recall;
- mean reciprocal rank;
- AUROC and AUPRC;
- disease-panel prediction;
- peptide-affinity prediction;
- paratope prediction;
- antigen-surface retrieval;
- calibration and confidence;
- performance on unseen V-gene, CDR-H3 and antigen families.

Use target-level bootstrap confidence intervals and define a non-inferiority margin before examining results.

## Compute evaluation

Report:

- number of PH calls;
- number of AF3 cofolds;
- GPU-hours;
- wall-clock time;
- peak memory;
- storage;
- conformer preprocessing time;
- ConFormer training and inference time;
- amortization cost of training a student generator;
- break-even number of future antibodies.

“Ten times faster at inference” is not compelling if training the student costs more than the AF3 compute it saves for the expected number of targets.

---

# 9. Major failure modes

| Failure mode | Why it is dangerous | How to detect it |
|---|---|---|
| Learning AF3/PH artifacts | Student reproduces seeds, templates or confidence patterns rather than antibody variability | Test on changed protocols, new seed ranges and alternate teachers |
| Ground-truth contact leakage | Antigen-derived contacts enter a supposedly sequence-only task | Regenerate using sequence-only/random-contact protocol |
| Conformer-level leakage | Near-identical structures appear in train and test | Split before conformer generation at paired-sequence/family level |
| V-gene or clonotype shortcuts | Antigen labels correlate with repertoire family | Lineage, CDR-H3 and paired-chain clustering |
| Sequence dominates geometry | ConFormer succeeds without using ensemble variation | Zero-geometry, shuffled-geometry and duplicated-conformer controls |
| Mode collapse | Generator produces valid but narrow structures | Teacher-to-generated coverage and cluster recall |
| Artificial overdiversity | Generator produces broad but unrealistic structures | Generated-to-experimental precision and geometry filters |
| Teacher frequency mistaken for occupancy | Raw cluster counts interpreted thermodynamically | Avoid free-energy language; report support separately from frequency |
| Circular evaluation | ConFormer selects and evaluates the same structures | Use independent structural metrics and separate downstream models |
| No downstream benefit | Full ensemble is not better than sequence or one structure | Run saturation and matched-capacity baselines first |
| Post-hoc compression only | Storage is reduced but AF3 cost is unchanged | Conduct prospective reduced-generation experiments |
| Synthetic structures do not reflect dynamics | Plausible structures are described as physical states | Use “proposal ensemble,” not “equilibrium ensemble” |
| Overclaiming antigen specificity | Geometry is treated as new antigen evidence | Apply the sequence→ensemble information argument and strict controls |
| Unequal target weighting | High-yield targets dominate training | Uniform target sampling and per-target normalization |
| Protocol confounding | Resolved and unresolved antibodies use different conditioning | Stratify or standardize generation protocols |

A particularly important test is **protocol transfer**:

- train the student on conformers generated using one contact-selection rule;
- evaluate against ensembles generated using another;
- test whether its output reflects antibody geometry or the contact-selection procedure.

---

# 10. A compelling but non-overhyped thesis framing

## Recommended title

**Task-aware compression and amortized sampling of antibody conformational ensembles for antigen retrieval**

A more biologically cautious title would be:

**Efficient approximation of pseudo-bound antibody conformational proposal ensembles**

A more algorithmic title:

**Budgeted conformer generation for antibody representation learning**

## Defensible claims

1. PH/AF3 antibody ensembles contain substantial measurable redundancy.
2. A small, intelligently chosen subset can preserve structural and representation-level coverage.
3. Adaptive sampling can reduce AF3 cost while preserving downstream performance.
4. A pretrained geometric generator can potentially be fine-tuned to reproduce selected aspects of the teacher ensemble.
5. Ensemble-derived supervision may improve the inductive bias of antibody representations.

## Claims to avoid

- “We recover the true antibody energy landscape.”
- “Our samples have physical occupancy probabilities.”
- “The ensemble adds antigen information beyond sequence.”
- “More RMSD diversity means better dynamics.”
- “Perfect antigen classification demonstrates general antigen specificity.”
- “A student improves AF3 accuracy” without independent experimental evidence.

The scientifically honest novelty is not merely using diffusion. It is defining **what information in a huge conformer ensemble is actually necessary**, and showing how to preserve it under a fixed computational budget.

---

# 11. Concrete research plan

## First 3 months: establish whether the project is real

### Data and provenance

Build a compact metadata table with one row per conformer:

- antibody target ID;
- paired heavy/light sequence IDs;
- source PDB;
- PH run/design/cycle;
- selected contact residues;
- pseudo-binder length and sequence;
- AF3 seed/model;
- template/MSA regime;
- filter result;
- confidence values;
- file location;
- structural hash or cluster ID.

Convert the relevant coordinates and features from millions of PDB files into a compact format rather than repeatedly parsing PDB text.

### Pilot analysis

Use 50–100 diverse targets spanning:

- CDR-H3 lengths;
- flexible and rigid loops;
- germline families;
- high- and low-yield PH runs;
- different antigen categories where available.

Run:

- duplicate analysis;
- hierarchical variance decomposition;
- random \(k\)-subset curves;
- stratified subset curves;
- k-medoids/farthest-point baselines;
- one-conformer and duplicated-conformer controls;
- frozen ConFormer output preservation.

### Three-month go/no-go criterion

Continue if at least one is true:

- substantial ensemble redundancy exists;
- different pipeline levels contribute identifiable complementary modes;
- informed selection beats random sampling;
- downstream performance saturates at a manageable \(k\);
- a prospective stopping rule appears possible.

Reconsider the project if the full ensemble has no reproducible advantage over one structure or sequence-only baselines.

---

## By 6 months: complete the core thesis

Scale the analysis to all leakage-safe targets.

Develop:

- target-specific adaptive \(k\);
- geometric coreset selection;
- ConFormer moment-preserving selection;
- sequential stopping based on marginal coverage;
- a surrogate predicting which PH designs deserve additional AF3 sampling.

Run prospective generation on a held-out panel with:

- 1 versus 5 AF3 models per seed;
- 1 versus 5 seeds;
- reduced PH design counts;
- reduced cycles;
- adaptive allocation.

Evaluate all structural, downstream and compute outcomes with target-level confidence intervals.

At this point, you should have a complete thesis even without generative modeling.

---

## By 12 months: generative extension

Only proceed if the first six months show that there is a coherent target distribution worth distilling.

### Recommended student model

- initialize from FrameFlow-, AlphaFlow-, BioEmu- or IgDiff-like weights;
- condition on paired sequence and one anchor structure;
- generate CDR residue frames and VH–VL orientation;
- keep framework restrained;
- train on cluster-balanced teacher medoids;
- include frame/FAPE-style loss, torsion loss, clash/closure penalties and representation-distribution matching;
- sample 32–128 structures per target rather than thousands.

Compare:

- student samples;
- random teacher subsets;
- task-aware teacher coresets;
- reduced-seed AF3;
- BioEmu;
- full PH/AF3 teacher ensembles.

A successful generative result would be:

> On held-out antibody families, 32–64 student samples recover most of the task-relevant structural support and downstream performance of a 1,500-conformer teacher ensemble at substantially lower marginal compute.

That would be a strong outcome. It does not require claiming equilibrium dynamics.

---

# 12. Questions to ask Gaeun and Sophia before committing

## Scientific target

1. **Were true antigen contacts, epitope residues, cognate-antigen size or complex-derived templates used when generating any ensembles evaluated on antigen prediction?**
2. Is the desired object a physical ensemble, a broad set of bound-like proposals or simply a useful representation for retrieval?
3. Is preserving rare structural modes more important than reproducing the teacher’s sample frequencies?
4. What result would count as success: storage reduction, AF3-call reduction, faster ConFormer training or improved retrieval?
5. Is there convincing evidence, on a leakage-safe split, that the complete ensemble outperforms sequence-only, one-structure and reduced-AF3 baselines?

## Data and provenance

6. What exactly defines an “ensemble group”: paired sequence, PDB entry, contact-conditioning run, antigen complex or PH run?
7. How many genuinely unique paired heavy–light sequences and clonal families are present?
8. Are PH design, cycle, AF3 seed, AF3 model, templates, confidence metrics and filter reasons stored for every PDB?
9. How many retained conformers are exact or near structural duplicates?
10. Are pseudo-binder structures and sequences retained?
11. Can train/test splits be frozen before any model development, including a temporal or external holdout?
12. Are there experimental multistate, apo/holo or ALL-conformations groups that can serve as independent validation?

## Downstream evaluation

13. Will the downstream ConFormer model be frozen while evaluating compression, retrained for each \(k\), or both?
14. What are the primary endpoints: antigen identity classification, retrieval, epitope-surface matching, peptide affinity or disease panels?
15. Are there matched-capacity sequence-only and single-structure baselines?
16. Can geometry-shuffling and duplicated-conformer negative controls be run?
17. What non-inferiority margin would be considered acceptable?

## Compute and thesis ownership

18. Can you run prospective reduced PH/AF3 experiments, or is only the existing archive available?
19. What is the actual AF3 GPU-hour cost per target and expected number of future targets?
20. Is training or adapting a geometric diffusion/flow model computationally realistic within the lab?
21. Who would maintain the equivariant-model code if the project extends beyond the thesis?
22. Is compression/distillation sufficiently distinct from Gaeun’s PhD contribution for you to own a clear first-author thesis question?
23. Can the data, code and benchmark be published, given AF3 outputs and any project restrictions?
24. Does Sophia prefer a methods thesis with a rigorous negative result, or is a new generative architecture expected?

---

# Final recommendation

**Commit to this direction if the project is defined as task-aware ensemble compression and prospective compute reduction.** That version is scientifically clean, directly useful to Gaeun, and feasible with the data already generated.

Treat a generative model as a **second-stage distillation experiment**, preferably restricted to CDR frames/torsions around an anchor and initialized from a pretrained geometric model.

Do not commit to a thesis whose required claim is that a few hundred PH/AF3 targets are sufficient to learn true antibody dynamics. The current data do not support that interpretation.

The sharpest thesis question is not:

> “Can diffusion generate antibody structures?”

That has largely been answered by existing protein and antibody generation work.

It is:

> **“What is the minimal structural distribution needed by a conformer-aware antibody representation, and can we sample that distribution without thousands of AF3 cofolds?”**

That question is original enough, falsifiable, computationally meaningful and appropriately scoped for a strong master’s thesis.
