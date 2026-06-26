# Advisor’s verdict

**Yes, this can become a serious thesis—but not as “train Stable Diffusion for antibodies.”**

The scientifically sharp version is:

> **Determine the smallest task-sufficient antibody conformer ensemble, then amortize or predict that small ensemble without running the full Protein Hunter/AlphaFold3 pipeline.**

The overly ambitious version is:

> **Train a full-atom antibody generator from scratch on a few hundred or thousand antibodies, claim it learns the true physical conformational distribution, and outperform Protein Hunter, AlphaFold3, molecular dynamics, and existing antibody ensemble models.**

That second version is not a sensible Master’s thesis target.

There is also an important novelty update. As of **June 26, 2026**, a generic “paired VH/VL sequence → antibody conformational ensemble using flow matching” project would overlap directly with the April 2026 **ABodyBuilder4-STEROIDS** preprint. Its authors report training on 4.2 million frames from approximately 136,000 coarse-grained antibody simulations and 83 all-atom simulations. The released implementation takes paired VH/VL sequences and samples ensembles using an SE(3) flow model. Therefore, building another generic antibody ensemble generator from scratch is no longer a compelling novelty claim. ([biorxiv.org](https://www.biorxiv.org/content/biorxiv/early/2026/04/16/2026.04.14.718378.full.pdf))

What remains novel is Gaeun’s **pseudo-bound, task-directed distribution**: Protein Hunter proposes artificial binding partners, then AF3 cofolding perturbs the antibody toward binding-compatible states. That is meaningfully different from learning a general MD-derived antibody ensemble.

My recommendation is therefore:

- **Primary thesis:** ensemble compression, adaptive sampling, and task preservation.
- **Generative extension:** distill the pseudo-bound teacher into a small set of CDR conformer modes, preferably starting from a pretrained antibody ensemble prior such as ABB4-STEROIDS rather than from random initialization.
- **Claim:** approximation of a computational proposal distribution for representation learning—not recovery of true antibody dynamics.

## What is observed, interpreted, and speculative

| Status | Statement |
|---|---|
| **Observed in Gaeun’s materials** | The current optimized pipeline uses 25 designs across 5 cycles, producing 125 Protein Hunter complexes, followed by 5 AF3 seeds × 5 models per design. That implies as many as **3,125 raw cofolds per antibody** before filtering. fileciteturn0file1 |
| **Observed in your tree** | You report hundreds to roughly 1,000 ensemble groups, approximately 1.5 million kept PDBs, and commonly 1,500–2,000 retained conformers per successful ensemble. I am taking those filesystem observations as given rather than independently verified. |
| **Observed in the current representation** | ConFormer processes aligned heavy/light variable-domain conformers with coordinates, atom types, torsions and residue identities; it aggregates ensemble information into pair representations, and the antigen-matching stage ultimately mean-pools the conformer axis. fileciteturn0file0 |
| **Plausible interpretation** | Protein Hunter/AF3 is producing a broad **pseudo-bound proposal ensemble** containing geometries useful to an antigen-retrieval model, but its sample frequencies are determined by design branches, seeds, templates and filters—not thermodynamic probabilities. |
| **Speculation to test** | A lightweight CDR-frame flow model or a set-of-prototypes transformer could reproduce the task-relevant modes with tens rather than thousands of conformers. |
| **Not currently justified** | Calling these samples a Boltzmann ensemble, assigning physical free energies from their frequencies, or claiming they represent real antibody dynamics. |

---

# 1. Are antibody conformer ensembles natural diffusion, flow, or transformer targets?

**Yes.** Mathematically, this is a conditional generative-model problem.

Let:

- \(s\) be the paired antibody sequence;
- \(x\) be an antibody conformation;
- \(u\) contain pipeline settings such as contact conditioning, pseudo-binder size, PH design/cycle, AF3 seed, AF3 model, template and filtering decisions.

The current teacher defines something like

\[
q_{\mathrm{teacher}}(x\mid s,u).
\]

A sequence-only student would instead learn

\[
p_\theta(x\mid s)
\approx
\int q_{\mathrm{teacher}}(x\mid s,u)\,\pi(u\mid s)\,du.
\]

That is perfectly legitimate conditional generative modeling. The problem is that \(\pi(u\mid s)\) is an **algorithmic sampling scheme**, not a physical equilibrium measure.

### Diffusion

Diffusion is natural because it can model multimodal structural distributions and has mature geometric formulations. But conventional diffusion may require many denoising steps, which is awkward when the purpose of the project is computational savings.

### Flow matching

Flow matching is probably the better starting point for this thesis. It learns a vector field transporting a simple prior toward the structural distribution and often needs fewer inference steps than diffusion. AlphaFlow and ESMFlow already demonstrated that a static structure predictor can be converted into a sequence-conditioned ensemble generator through flow matching, including training on MD ensembles. ([arxiv.org](https://arxiv.org/abs/2402.04845))

AlphaFlow-Lit is an especially relevant precedent: it focused on a lightweight structure module and reported approximately 47-fold faster sampling than AlphaFlow while retaining similar performance. That does not prove antibody distillation will work, but it directly supports the feasibility of amortizing an expensive structure pipeline. ([openreview.net](https://openreview.net/forum?id=Z6fPAsu91p))

### Transformer-based generation

“Transformer” is not really an alternative to diffusion or flow. A Transformer can parameterize the diffusion score, flow vector field, latent-mode predictor, or discrete structure-token model.

For your data regime, a transformer may be most useful for **jointly predicting a fixed small set of conformer prototypes**, rather than autoregressively generating every atom:

1. Encode sequence and one reference antibody structure.
2. Introduce \(K\) learned conformer-query tokens.
3. Decode \(K\) CDR backbone/frame prototypes and \(K\) weights.
4. Match the predictions to teacher clusters using optimal transport or Hungarian matching.
5. Train with structural, diversity, closure and task-distillation losses.

That is an **amortized coreset generator**. It is arguably better aligned with “how few conformers are enough?” than an unconstrained stochastic model that must be sampled thousands of times.

### The representation matters more than the generative brand

For a first model, I would use:

\[
\mathrm{SE}(3)^{L_{\mathrm{CDR}}}
\times
\mathbb{T}^{N_{\mathrm{torsion}}},
\]

meaning one rotation/translation frame per flexible residue plus periodic torsion variables.

A practical formulation would:

- hold most of the VH/VL framework fixed;
- model all CDR backbone frames or begin with CDR-H3;
- optionally model a small rigid-body perturbation of VH relative to VL;
- reconstruct backbone atoms analytically;
- add side chains with a deterministic or learned rotamer packer;
- enforce loop closure, chirality, peptide geometry and disulfide constraints.

Pure Cartesian atom noise is possible, but it wastes data learning bond geometry and global SE(3) symmetries that should be encoded by construction.

---

# 2. How this differs from image/video diffusion and LoRA fine-tuning

Your intuition from image generation is useful, but only at the highest level.

### The valid analogy

- Antibody sequence is analogous to a condition or prompt.
- Latent noise selects one of several possible outputs.
- Repeated sampling explores a conditional distribution.
- A pretrained generative prior can be adapted to a narrower domain.
- Guidance or reranking can bias samples toward downstream utility.

### The analogy breaks in several important ways

**Protein coordinates are not pixels.** A slightly wrong pixel may be harmless; a slightly wrong peptide geometry can break chirality, create a chain discontinuity or produce severe steric clashes.

**Global orientation is meaningless.** Rotating an antibody must not change the modeled probability. The model must therefore be invariant or equivariant to SE(3), rather than learning orientation through augmentation alone.

**The ensemble is unordered.** Protein Hunter cycles and AF3 seeds are not a physical time trajectory. Treating them as frames of a video would teach false dynamics. Only genuine MD trajectories have meaningful temporal order, and even then subsampling and correlated frames require care.

**The data are hierarchically correlated.** Two thousand conformers from one antibody are more like two thousand augmentations of the same image subject than two thousand independent concepts. They help estimate that antibody’s local distribution, but they do not teach generalization to two thousand new sequences.

**Sample frequency matters differently.** Image diversity is often judged perceptually. A conformational ensemble ideally needs support, occupancy and possibly thermodynamic calibration. Protein Hunter/AF3 supplies an algorithmic sample count, not an equilibrium probability.

**Guidance changes the distribution.** A task-guided sampler may improve antigen retrieval but no longer approximates the teacher distribution, much less a Boltzmann distribution. That is acceptable, but it must be stated.

### Where LoRA might fit

LoRA is a parameter-efficient adaptation mechanism, not the scientific idea.

A reasonable use would be:

- start from ABB4-STEROIDS, AlphaFlow-like, FrameFlow-like or another accessible geometric ensemble model;
- freeze most of the sequence/structure prior;
- adapt the geometric vector field or final structure module to Gaeun’s PH/AF3 teacher;
- regularize against the original MD/static-structure training objective.

A LoRA applied to an ordinary static predictor does not automatically make it an ensemble model. You still need:

- a noising/interpolation process;
- a conditional generative objective;
- an appropriate geometric representation;
- target-level splits;
- distributional evaluation.

---

# 3. True biophysical distributions versus pipeline distillation

This distinction should be the center of the thesis.

## A. True physical ensemble

The ambitious target is something like

\[
p_{\mathrm{phys}}(x\mid s, T, \text{solvent}, \text{protonation}, \text{glycosylation}, \text{binding state},\ldots).
\]

This distribution has:

- physically meaningful occupancies;
- temperature and environmental dependence;
- free-energy relationships;
- potentially meaningful kinetics if transitions are also modeled.

BioEmu is an example of a model explicitly aimed at equilibrium ensembles. It integrated more than 200 milliseconds of MD, static structures and experimental stability data, and was designed to generate many approximately independent equilibrium samples. ([microsoft.com](https://www.microsoft.com/en-us/research/publication/scalable-emulation-of-protein-equilibrium-ensembles-with-generative-deep-learning/))

ABB4-STEROIDS is now the antibody-specific analogue closest to that ambition, although its current result is a recent preprint and its reported performance should be independently benchmarked. ([biorxiv.org](https://www.biorxiv.org/content/biorxiv/early/2026/04/16/2026.04.14.718378.full.pdf))

## B. Gaeun’s pipeline-induced ensemble

Gaeun’s current process is better viewed as

\[
q_{\mathrm{PH/AF3}}(x\mid s,u),
\]

where \(u\) includes pseudo-antigen design choices and AF3 inference choices.

Protein Hunter itself was developed as a fine-tuning-free **protein design** method using iterative structure/sequence cycling, not as a calibrated antibody dynamics model. Gaeun’s clever contribution is repurposing it as a pose and conformation sampler. ([biorxiv.org](https://www.biorxiv.org/content/10.1101/2025.10.10.681530v1))

The resulting ensemble may be valuable because it explores:

- binding-compatible CDR arrangements;
- artificial induced-fit-like states;
- geometry useful for a representation model;
- broader conformational support than an unbound predictor.

But the number of times a mode appears does not tell you its physical population.

## C. Task-optimal proposal distribution

For antigen retrieval, the most useful target may be neither the true equilibrium distribution nor the raw teacher distribution. It may be

\[
r_{\mathrm{task}}(x\mid s),
\]

a small set of conformers that maximizes retrieval information.

For example, a rare high-energy binding-competent state could be very useful for antigen retrieval even if it has low equilibrium occupancy. Conversely, a dominant low-energy framework fluctuation may add no antigen information.

That makes the following thesis claim defensible:

> “We learn or select a compact set of pseudo-bound antibody conformations that preserves the structural and predictive information of an expensive teacher ensemble.”

It does **not** make this claim defensible:

> “We learn the true antibody conformational landscape.”

## Can a student improve on the teacher?

Only in a qualified sense.

A student trained only on PH/AF3 outputs can improve:

- speed;
- storage;
- geometric smoothness through inductive bias;
- downstream utility through task-aware weighting;
- robustness through regularization.

It cannot establish improved physical truth without external data. To claim physical improvement, you need comparison with MD, repeated experimental structures, NMR/cryo-EM ensembles, apo/holo states or another independent source.

---

# 4. Relevant model families

## Directly relevant ensemble models

### AlphaFlow and ESMFlow

These are the closest general-protein conceptual precedents. They repurpose AlphaFold and ESMFold as flow-matching ensemble generators and distinguish training on static PDB structures from training on MD distributions. They support the basic proposition that an expensive structure predictor can be amortized into a sampler. ([arxiv.org](https://arxiv.org/abs/2402.04845))

### BioEmu

BioEmu is relevant as the standard for what a serious **equilibrium-emulation** claim requires: extensive MD, static structures, experimental observables and distributional validation. It is not merely trained on stochastic predictions from another neural structure pipeline. ([pubmed.ncbi.nlm.nih.gov](https://pubmed.ncbi.nlm.nih.gov/40638710/))

### ABB4-STEROIDS

This is now the mandatory antibody-specific baseline. It directly accepts paired VH/VL sequences and generates conformational ensembles using SE(3) flow matching. Its training scale also demonstrates why a few thousand PH/AF3 targets are not a competitive foundation for a new general-purpose physical antibody ensemble model—but could be enough to adapt a strong prior to a specialized pseudo-bound teacher. ([biorxiv.org](https://www.biorxiv.org/content/biorxiv/early/2026/04/16/2026.04.14.718378.full.pdf))

## Relevant geometric generator architectures

### FrameDiff and FrameFlow

FrameDiff formulated diffusion over residue rigid frames in SE(3), while FrameFlow translated this idea to flow matching. These are strong architectural templates for CDR backbone generation because rotations, translations and global symmetries are handled explicitly. They were developed for de novo backbone generation, not fixed-sequence conformer distributions. ([arxiv.org](https://arxiv.org/abs/2302.02277))

### RFdiffusion

RFdiffusion demonstrated powerful conditional protein-backbone generation from noise, including target-, motif- and symmetry-conditioned design. Its importance here is architectural and conceptual. It does not by itself learn equilibrium conformers of a fixed antibody sequence. ([nature.com](https://www.nature.com/articles/s41586-023-06415-8))

## Antibody design models

### DiffAb

DiffAb jointly diffuses CDR sequence and structure in the context of an antigen using equivariant networks. It is relevant for antigen-conditioned CDR generation, but its goal is antibody design, not sampling the conformational distribution of one fixed antibody sequence. ([proceedings.neurips.cc](https://proceedings.neurips.cc/paper_files/paper/2022/hash/3fa7d76a0dc1179f1e98d1bc62403756-Abstract-Conference.html))

### IgDiff and IgFlow

IgDiff adapts SE(3) backbone diffusion to paired antibody variable domains and conditional CDR design. IgFlow applies flow matching to related antibody-generation problems. These offer useful antibody-specific architectural priors, but they optimize designability and novelty rather than teacher-ensemble fidelity. ([arxiv.org](https://arxiv.org/html/2405.07622v1))

### AbDiffuser

AbDiffuser jointly generates full-atom antibody structures and sequences using equivariant, physics-informed diffusion. It shows that atomistic antibody generation is possible, but its data and objective concern design—not producing calibrated multiple conformations for a fixed sequence. ([proceedings.neurips.cc](https://proceedings.neurips.cc/paper_files/paper/2023/file/801ec05b0aae9fcd2ef35c168bd538e0-Paper-Conference.pdf))

### RFantibody

RFantibody specializes RFdiffusion for framework- and epitope-conditioned CDR design, followed by sequence design and structural filtering. It is highly relevant to antibody design but not a direct conformer-ensemble model. ([nature.com](https://www.nature.com/articles/s41586-025-09721-5))

## ESM3-like structure-token models

ESM3 is a multimodal model over protein sequence, structure and function. Structure-token models are attractive for:

- compressing conformers;
- predicting conformer clusters;
- indexing large structural sets;
- generating latent structural modes with a Transformer.

However, the reconstruction error of the tokenizer must be substantially smaller than the CDR movements you care about. If a tokenizer introduces roughly the same geometric error as the distinction between two CDR modes, it destroys the signal. ESM3-like tokens are therefore more attractive for latent mode modeling and retrieval than as the first coordinate generator. ([science.org](https://www.science.org/doi/10.1126/science.ads0018))

---

# 5. Is the available training data enough?

## The central rule

**Millions of conformers do not equal millions of independent training examples.**

For generalization to unseen antibodies, the effective sample size is closer to:

- the number of distinct paired VH/VL sequences;
- the number of clonotypes or sequence families;
- the number of distinct CDR-H3 sequence/length families;
- the number of independently generated structural modes.

The 1.5 million PDBs are valuable for learning within-target variation, but they do not compensate for only hundreds or thousands of independent targets.

## My feasibility judgment

| Objective | Few hundred targets | ~1,000 targets | A few thousand targets |
|---|---:|---:|---:|
| Post hoc conformer selection | More than enough | Yes | Yes |
| Adaptive stopping/budget allocation | Plausible | Strong | Strong |
| Predict full-ensemble embedding | Plausible | Strong | Strong |
| Predict latent conformer clusters | Risky but testable | Plausible | Good with pretraining |
| CDR-H3 residual flow from a pretrained prior | Proof of concept | Plausible | Credible |
| Full paired-Fv backbone flow from scratch | Weak | High risk | Still high risk |
| Full-atom antibody generator from scratch | No | No | Probably underpowered |
| True equilibrium model competitive with ABB4/BioEmu | No | No | No |

These are qualitative judgments, not universal thresholds. The correct way to answer the data question is a **target-level learning curve**:

\[
N_{\mathrm{targets}}\in
\{50,100,250,500,1000,\ldots\}
\]

with a fixed, family-held-out test set.

## Required sampling strategy

Do not sample PDBs uniformly from the 1.5 million files. That lets targets with 2,000 retained conformers dominate targets with 100.

Use hierarchical batches:

1. sample an antibody target uniformly;
2. sample a structural cluster or PH design branch;
3. sample one conformer from that mode or branch.

For a generator, normalize losses per target, not per file.

## Deduplication

Before model training, estimate how many independent modes exist after:

- framework alignment;
- CDR-specific RMSD clustering;
- torsion clustering;
- VH/VL orientation clustering;
- near-duplicate removal within each AF3 seed/design branch.

It would not surprise me if 2,000 retained PDBs collapsed to tens of meaningful modes. That is exactly what the first thesis experiment should discover.

## Leakage controls

At minimum, group all of the following into the same split:

- identical paired VH/VL sequences;
- sequence-near-identical antibodies;
- clonally related antibodies;
- the same CDR-H3 sequence, especially with similar light chains;
- multiple PDBs of the same antibody;
- engineered variants and affinity-maturation lineages;
- all conformers generated from one target;
- repeated complexes or antigen constructs.

A V-gene-disjoint split is not sufficient by itself. Gaeun’s proposed no-V-gene-overlap splits are better than random splitting, but sequence similarity and lineage leakage can remain. fileciteturn0file1

For downstream antigen retrieval, also control:

- antigen sequence and structure similarity;
- antigen family;
- epitope overlap;
- candidate-bank composition;
- antigen size and fold;
- patient, study and cohort for disease datasets.

A temporal split based on deposition date is particularly valuable, provided all model-training cutoffs and database snapshots are documented.

## ALL-conformations and ITsFlexible

ALL-conformations contains 1.2 million loop structures representing more than 100,000 unique sequences and is a valuable source of experimentally observed loop-state support. But it is PDB-derived and contains antibody structures related to SAbDab, so it is not automatically independent of a SAbDab-derived training set. Exact and near-sequence overlap must be removed. It is also better suited to validating whether modes are experimentally observed than to assigning physical occupancy probabilities. ([nature.com](https://www.nature.com/articles/s42256-025-01131-6))

---

# A major project-specific leakage concern

The April pipeline description says that, for structurally resolved complexes, Protein Hunter is conditioned using the **top three known contact residues**, and that designed binder size is proportional to the **ground-truth binder size**. Unknown structures instead use random CDR-H3 contacts. fileciteturn0file1

This creates a possible privileged-information pathway:

- true paratope/contact positions come from the cognate complex;
- antigen or binder size may correlate strongly with antigen class;
- resolved training examples and sequence-only deployment examples receive different conditioning;
- the generated conformations may encode information unavailable for a new BCR.

That does not necessarily invalidate the current work, but it must be isolated.

Essential controls are:

1. Use random contact conditioning and a fixed binder-size distribution for **all** antibodies.
2. Compare true-contact, predicted-contact and random-contact conditions.
3. Train a classifier to predict contact-conditioning scheme or binder size from the resulting antibody ensemble. High accuracy would demonstrate a pipeline artifact.
4. Evaluate retrieval only using information available at deployment.
5. Keep the privileged teacher as a separate upper bound rather than mixing it with the blind pipeline.

The preliminary 20-antigen classification result—MCA accuracy 1.0 versus 0.83 for an ESM2 linear head—is promising, but it should be treated as a debugging result, not evidence of general antigen inference, until these controls and strict splits are complete. fileciteturn0file2

---

# 6. Which model formulation is most realistic?

## My ranking

### 1. Post hoc conformer selection and coreset construction — **best first project**

Given a full ensemble \(E\), choose a subset \(S_K\) that preserves:

- structural mode coverage;
- ConFormer pair/ensemble statistics;
- downstream predictions;
- quality.

Useful methods include:

- random and stratified sampling;
- k-medoids;
- greedy k-center or farthest-point sampling;
- facility-location objectives;
- determinantal point processes;
- kernel herding to preserve embedding means;
- cluster-weighted medoids.

This answers how many conformers are needed and establishes an oracle upper bound for any future generator.

Its limitation is that post hoc compression saves storage and training cost but **not initial AF3 compute**.

### 2. Hierarchical budget allocation and adaptive stopping — **highest likely practical value**

The pipeline has a clear hierarchy:

\[
\text{PH design}
\rightarrow
\text{cycle}
\rightarrow
\text{AF3 seed}
\rightarrow
\text{AF3 model}.
\]

Decompose structural variance across those levels. Determine where new modes actually originate.

For example:

- If most diversity comes from different PH designs, one AF3 prediction per design may be enough.
- If the five AF3 models within a seed are nearly redundant, remove them.
- If later PH cycles contribute little new coverage, stop early.
- If some antibodies saturate at 16 samples and others need 128, predict a target-specific budget.

A sequential policy can generate batches and stop when:

- no new cluster has appeared recently;
- embedding mean/covariance has converged;
- estimated unseen-mode mass is below a threshold;
- downstream prediction is stable;
- marginal coverage per GPU-hour falls below a threshold.

This directly reduces generation cost and may be more valuable than a generative model.

### 3. Direct representation distillation — **essential baseline**

Train a cheap model to predict the full-ensemble ConFormer representation from:

- paired sequence;
- one static structure;
- or a very small conformer set.

If the only use is antigen retrieval, this may solve the practical problem without generating PDBs.

Scientifically, it is less general because it is tied to the current encoder. But it must be included: otherwise a structure generator might be solving a harder problem than the downstream task requires.

### 4. Latent conformer clusters or amortized coreset prediction — **best generative thesis formulation**

For each teacher ensemble:

1. framework-align conformers;
2. cluster CDR geometries;
3. choose cluster medoids and weights;
4. encode medoids into CDR frames, torsions or a continuous structural latent;
5. train a conditional model to predict the small set of medoids and weights.

This changes the target from “generate thousands of samples” to:

> “Predict the 8–32 structural modes needed to represent the teacher.”

A set-prediction Transformer with learned queries and optimal-transport matching is a good formulation. It is likely more sample-efficient than full coordinate diffusion and directly controls output budget.

### 5. Residual CDR frame/torsion flow — **best stochastic extension**

Condition on:

- paired VH/VL sequence;
- a reference structure;
- CDR masks and loop anchor frames;
- optional deployable contact predictions;
- a latent mode variable.

Generate:

- CDR residue frames;
- backbone torsions;
- optionally VH/VL orientation;
- later, side-chain torsions.

I would not start from completely isotropic whole-antibody noise. Start from a reference structure plus calibrated noise or from an ABB4-generated proposal. This makes the problem a learned structural transport or refinement problem.

An especially strong model would be:

\[
x_0\sim p_{\mathrm{ABB4}}(x\mid s),\qquad
x_1=T_\theta(x_0,s,z),
\]

where \(T_\theta\) transports an MD-informed antibody prior toward the PH/AF3 pseudo-bound distribution. That is much more data-efficient and novel than rebuilding ABB4.

### 6. Full variable-domain generation — **future work**

Jointly generating all VH/VL backbone frames may be justified if ablations show that:

- framework motion;
- VH/VL orientation;
- non-CDR residues

are needed for retrieval. Otherwise it wastes capacity modeling mostly stable regions.

### 7. Full-atom paired-Fv generation from noise — **do not make this the primary thesis**

It adds:

- side-chain rotamers;
- hydrogen-bond geometry;
- disulfide constraints;
- clashes;
- atom naming and completeness;
- a much larger loss space.

The current downstream model must first demonstrate that such atomistic detail provides information beyond backbone/CDR geometry.

---

# A critical ConFormer audit before generative modeling

The current methods document raises two immediate representation-level questions.

First, ConFormer’s outer-product update computes conformer weights using **conformer 0 as the query**. Unless conformer 0 is a deliberately defined canonical reference, this can introduce ensemble-order sensitivity.

Second, the downstream antigen-matching stage arithmetic-mean-pools conformer embeddings. Conformer duplicates and the teacher’s arbitrary sampling density can therefore affect the representation, even when they do not add structural information. fileciteturn0file0

Run these tests immediately:

1. Randomly permute conformer order many times.
2. Replace conformer 0 with different ensemble members.
3. Duplicate one conformer 2×, 10× and 100×.
4. Hold unique modes fixed while changing their multiplicities.
5. Compare a cluster-weighted mean with a raw conformer-weighted mean.
6. Compare repeated copies of one central conformer with a genuinely diverse ensemble.
7. Vary \(N\) without changing structural support.

If predictions materially change, some apparent “ensemble-size effects” may actually be ordering or multiplicity artifacts. Fixing this may yield larger gains than training a generator.

---

# 7. Minimum viable thesis

## Recommended thesis aims

### Aim 1: Quantify redundancy and information saturation

Create a clean target-level dataset manifest containing:

- paired sequences and sequence-family identifiers;
- source PDB and antigen;
- PH design and cycle;
- contact-conditioning scheme;
- pseudo-binder length and identity;
- AF3 seed and model;
- template and MSA settings;
- filter outcome and quality scores;
- final conformer cluster.

Measure:

- number of conformers;
- number of nonredundant modes;
- variance contribution of each pipeline level;
- mode saturation as samples are added.

### Aim 2: Build a task-preserving ensemble coreset

Evaluate \(K\) in:

\[
K\in\{1,2,4,8,16,32,64,128,256,512\}.
\]

Compare random, stratified and structural-selection methods. Use the full ensemble as the teacher upper bound and repeat random sampling many times.

Evaluate both:

- frozen ConFormer at inference, isolating the effect of shrinking \(N\);
- retraining ConFormer with small ensembles, testing whether the whole system can adapt.

### Aim 3: Reduce generation-time cost

A post hoc subset alone does not save AF3 compute. Therefore run a factorial or staged budget study varying:

- number of PH designs;
- number of PH cycles;
- AF3 seeds per design;
- models per seed.

Train an adaptive stopping or acquisition policy if simple fixed budgets are suboptimal.

### Optional Aim 4: Predict the coreset or augment a small anchor ensemble

Only after Aims 1–3 succeed:

- predict teacher cluster medoids in latent/frame space;
- or train a CDR residual flow;
- or augment 2–8 PH/AF3 anchor structures with cheap generated conformers.

## Three decisive experiments

These should determine whether the generative extension is justified.

### Experiment 1: Does variability matter?

Compare:

- full ensemble;
- single best conformer;
- repeated copies of the same conformer;
- mean structure;
- sequence-only model;
- single static structure;
- RMSD-matched random perturbations.

If the full ensemble does not reliably beat these on strict splits, there is no reason to distill it.

### Experiment 2: Is there a compact oracle subset?

Select the best \(K\)-member subset after seeing the full ensemble.

- If random \(K=16\) works, use random or stratified sampling.
- If oracle selection works but random does not, learned selection is justified.
- If even the oracle subset fails at the desired \(K\), a \(K\)-sample generator cannot reproduce the full teacher unless it creates better-than-teacher structures.

### Experiment 3: Does PH/AF3 beat ABB4-STEROIDS?

Compare blind, sequence-only ensembles from:

- ABB4-STEROIDS;
- multi-seed AF3;
- reduced Protein Hunter/AF3;
- full Protein Hunter/AF3.

If ABB4 matches structural coverage and downstream retrieval at much lower cost, the correct conclusion may be to use ABB4—not to build another generator.

If PH/AF3 clearly outperforms ABB4 under blind conditioning, that supports the hypothesis that pseudo-bound sampling contributes task-relevant states and provides a strong rationale for distillation.

---

# 8. Essential baselines and evaluation

## Generation and selection baselines

At minimum:

1. Single static structure.
2. Repeated copies of one conformer.
3. Random subsampling, with many repetitions.
4. Sampling stratified by PH design and AF3 seed.
5. One AF3 model per PH design.
6. Reduced PH designs and cycles.
7. Reduced-seed unbound AF3.
8. Quality-score-only selection.
9. k-medoids or k-center selection.
10. Facility-location or DPP selection.
11. Full teacher ensemble.
12. ABB4-STEROIDS.
13. Direct full-ensemble embedding distillation.

BioEmu or AlphaFlow can be additional general-protein comparators if they support the antibody lengths and chain conventions reliably, but ABB4 is the more important direct comparator.

## Structural diversity and fidelity

Do not report only mean pairwise RMSD. That metric rewards physically absurd outliers.

Use both **coverage** and **precision**:

- teacher-to-subset nearest-neighbor distance;
- subset-to-teacher nearest-neighbor distance;
- fraction of teacher clusters represented;
- fraction of generated samples lying near a teacher or external reference mode;
- cluster-weight Jensen–Shannon divergence, only when weights are considered meaningful;
- energy distance, MMD or Wasserstein distance in an independent structural representation.

Antibody-specific descriptors should include:

- framework-aligned CDR-H3 RMSD;
- RMSD and torsion distributions for all six CDRs;
- CDR-H3 tip displacement;
- per-residue RMSF;
- backbone \(\phi,\psi,\omega\);
- loop-contact maps;
- canonical-loop cluster occupancy for non-H3 loops;
- VH/VL orientation;
- paratope surface shape and solvent exposure;
- known apo-to-holo state coverage.

Embedding diversity is useful, but avoid circularity. If selection uses ConFormer embeddings, do not declare success solely because the selected subset preserves those same embeddings.

## Physical and geometric validity

Measure:

- bond-length and bond-angle deviations;
- peptide-chain breaks;
- chirality;
- cis/trans peptide geometry;
- Ramachandran violations;
- steric clashes;
- disulfide geometry;
- CDR loop closure;
- VH/VL interface clashes;
- side-chain rotamer plausibility;
- energy before and after a short standardized minimization.

Any relaxation or repacking cost must be included in the compute budget. A “fast” generator that requires expensive AF3 revalidation or long molecular minimization for every sample may not provide real savings.

## External conformational validation

Separate two questions:

1. **Teacher fidelity:** does the student reproduce PH/AF3?
2. **Biological fidelity:** does it cover experimental or MD-supported conformations?

Possible external sets:

- strict-overlap-filtered ALL-conformations;
- repeated structures of identical or nearly identical antibody sequences;
- apo/holo antibody pairs;
- a small held-out all-atom MD set;
- NMR or cryo-EM multi-state examples;
- ABB4’s available MD benchmark data.

Teacher fidelity alone is insufficient for claims about dynamics.

## Downstream antigen retrieval

For a real retrieval problem, report:

- Recall@1, @5, @10;
- mean reciprocal rank;
- mean average precision;
- candidate-bank size;
- performance against fold-, size- and surface-matched hard negatives;
- macro-averaged metrics across antigens;
- calibration and confidence;
- target-family bootstrap confidence intervals.

Closed-set 20-class classification is not equivalent to general antigen retrieval.

Essential downstream controls include:

- ESM or another sequence-only baseline;
- static-structure baseline;
- ensemble with shuffled or repeated conformers;
- ensemble-summary statistics only;
- full versus compressed ensemble;
- blind versus privileged PH conditioning.

For disease-panel prediction, split by patient and cohort, not merely by antibody sequence.

## Compute savings

Report separately:

- PH calls;
- AF3 calls;
- GPU-hours;
- wall-clock time;
- storage;
- preprocessing and I/O;
- generator training cost;
- inference cost;
- relaxation/filter cost.

For a distilled model, report the amortization break-even point:

\[
N_{\mathrm{break-even}}
=
\frac{C_{\mathrm{train}}}
{C_{\mathrm{teacher/target}}-C_{\mathrm{student/target}}}.
\]

## Statistical testing

Use:

- target-level or family-level bootstrap intervals;
- repeated random-subset trials;
- hierarchical analysis respecting targets and branches;
- predefined non-inferiority margins;
- performance-versus-cost Pareto curves.

Do not interpret “not statistically significantly worse” as equivalence. Define how much retrieval or coverage loss is acceptable before looking at the result.

---

# 9. Major failure modes

## 1. Learning Protein Hunter or AF3 artifacts

The student may reproduce template, seed, MSA, pseudo-binder or filter signatures rather than antibody flexibility.

**Control:** predict pipeline metadata from ensemble representations; use pipeline-version holdouts; randomize nonbiological settings.

## 2. Distilling privileged antigen information

Known contacts and ground-truth binder size may leak deployability-adjacent information.

**Control:** blind conditioning for all final tests; privileged teacher as a separately labeled upper bound.

## 3. Conformer-level train/test leakage

Randomly splitting 1.5 million PDBs would produce nearly identical structures on both sides.

**Control:** split by antibody family before any conformer sampling.

## 4. SAbDab family leakage

V-gene separation alone may leave highly similar CDRs, clonotypes, engineered variants or repeated PDB entries.

**Control:** paired-sequence clustering, CDR-H3 clustering, lineage grouping and temporal tests.

## 5. Mode collapse

A generator can produce beautiful central conformers while missing rare binding-compatible modes.

**Control:** teacher-to-student coverage, cluster recall and multiple-sample best-of-\(K\) metrics.

## 6. Diversity gaming

A model can improve pairwise RMSD by generating clashes or broken loops.

**Control:** report precision and physical validity alongside diversity.

## 7. No downstream benefit

The full ensemble may not outperform sequence or one structure under strict splits.

**Interpretation:** this would challenge the central ensemble hypothesis, not merely indicate a failed generator.

## 8. Current pooling does not exploit multimodality

Mean pooling or conformer-0 dependence may make the network sensitive to counts but insensitive to meaningful modes.

**Control:** permutation, duplication, repeated-conformer and cluster-weighted pooling tests.

## 9. Synthetic conformers are not real dynamics

PH/AF3 samples are not time-resolved and have no thermodynamic weighting.

**Control:** use “pseudo-bound proposal ensemble” language; external MD/experimental evaluation for any physical claim.

## 10. The student cannot reproduce hidden conditioning

If the teacher depends heavily on true contacts, binder size or branch metadata and the student sees only sequence, it may average incompatible modes.

**Control:** preserve deployable conditioning variables, introduce an explicit latent mode, or generate a task-oriented coreset rather than matching raw frequencies.

## 11. Compute savings disappear downstream

Full-atom reconstruction, minimization, AF3 validation or quality filtering may cost more than generation.

**Control:** end-to-end accounting.

## 12. Antigen-specificity overclaiming

An antibody may have incomplete specificity labels, unknown cross-reactivities or sequence-family correlations with antigen.

**Control:** call the task retrieval or ranking within a specified antigen bank; include hard negatives and uncertainty; do not imply universal cognate-antigen identification.

## 13. Moving-target thesis

If Gaeun’s conformer pipeline, ConFormer architecture and downstream datasets change continuously, your benchmarks will become uninterpretable.

**Control:** freeze a data release, pipeline version, split and primary metric before model development.

---

# 10. A compelling but defensible thesis framing

## Recommended title

**How many antibody conformers are enough? Task-preserving compression and amortized distillation of pseudo-bound conformational ensembles for antigen retrieval**

Alternative:

**Budgeted sampling of antibody conformational ensembles for conformer-aware antigen representation learning**

For a more generative version:

**Amortized conformer coreset prediction for pseudo-bound antibody ensembles**

## Defensible thesis claims

1. Quantify redundancy and hierarchical sources of structural diversity in a large PH/AF3 antibody ensemble dataset.
2. Establish the smallest ensemble size that preserves structural mode coverage and downstream retrieval.
3. Develop an adaptive or learned selection method that reduces generation or storage cost.
4. Optionally predict a small conformer coreset or distill pseudo-bound CDR modes from a pretrained antibody ensemble prior.
5. Compare pseudo-bound teacher ensembles with MD-trained antibody ensemble generation.

## Claims to avoid

- “We recover the true antibody free-energy landscape.”
- “AF3 seeds are samples from a Boltzmann distribution.”
- “Our model simulates antibody dynamics.”
- “Conformational diversity proves antigen specificity.”
- “Millions of conformers provide millions of independent examples.”
- “Our model improves biological realism” without independent physical or experimental evidence.

## The most interesting scientific question

The project can be deeper than simple compression:

> **Do pseudo-bound PH/AF3 ensembles contain antigen-retrieval information that is absent from fast sequence-conditioned MD-like antibody ensembles?**

That question separates Gaeun’s approach from ABB4-STEROIDS and generic ensemble generation. Possible outcomes are all informative:

- PH/AF3 wins: pseudo-bound sampling adds useful task-specific states.
- ABB4 wins: physical/intrinsic flexibility is enough; PH/AF3 may be unnecessary.
- Both are equivalent: downstream representation is robust to how diversity is generated.
- Neither beats a static or sequence baseline: current ensemble benefit was likely leakage, sample-size or architecture related.

---

# 11. Concrete research plan

## First 3 months: establish whether the project premise survives

### Month 1: audit and freeze

- Build a complete data manifest.
- Count unique paired sequences, clonotypes and structural families.
- Identify every use of true contacts, antigen size or other privileged information.
- Freeze pipeline version, data snapshot and downstream split.
- Cluster and deduplicate conformers.
- Run ConFormer permutation, conformer-0 and duplication tests.
- Install and benchmark ABB4-STEROIDS on a small clean target set.

### Month 2: saturation and variance decomposition

- Compute structural variance attributable to PH design, cycle, AF3 seed and model.
- Run \(K\)-curves for random and stratified subsampling.
- Implement k-medoids, k-center and facility-location coresets.
- Compare cluster-uniform, branch-uniform and conformer-uniform weighting.
- Evaluate structural coverage and geometry.

### Month 3: downstream gates

- Freeze ConFormer and test full versus compressed ensembles.
- Compare sequence-only, static, repeated-single, random-noise and ABB4 baselines.
- Test blind versus privileged PH conditioning.
- Define non-inferiority and compute-savings criteria.
- Produce the first performance–cost Pareto curve.

### Three-month go/no-go decisions

- **Full ensemble does not beat static/repeated:** do not build a generator; investigate representation and leakage.
- **Random small \(K\) works:** focus on early stopping and reduced pipeline branching.
- **Oracle coreset works but random does not:** build a learned selector or prototype predictor.
- **PH/AF3 beats ABB4 under blind tests:** proceed to pseudo-bound distillation.
- **ABB4 matches PH/AF3:** pivot to benchmarking, representation distillation or hybrid use rather than training a new ensemble model.

## First 6 months: complete a defensible Master’s thesis core

- Expand evaluation to all strict train/test families.
- Train an adaptive target-specific budget predictor.
- Develop a two-stage selector: choose PH designs before AF3, then select final conformers.
- Retrain ConFormer using compressed ensembles.
- Include direct full-ensemble representation distillation.
- Validate on strict-overlap-filtered experimental and MD sets.
- Quantify end-to-end compute and storage reductions.
- Write a paper-quality methods and benchmark draft.

A strong six-month result would be:

> “A cluster-aware adaptive sampler retains downstream performance within a predefined non-inferiority margin while reducing AF3 calls and retained structures by an order of magnitude.”

That is already a serious thesis.

## First 12 months: add the generative component

- Build teacher clusters and target cluster medoids.
- Train a set-of-prototypes model or CDR-frame residual flow.
- Start with CDR-H3, then test all CDRs and VH/VL orientation.
- Compare direct PH/AF3 distillation with ABB4-initialized transport.
- Evaluate generated structural support, validity and downstream utility.
- Test hybrid strategies:
  - 2–8 AF3 anchors plus generated augmentation;
  - ABB4 proposal plus pseudo-bound reranking;
  - predicted cluster medoids plus deterministic atom reconstruction.
- Run target-level learning curves.
- Complete external validation and thesis writing.

The generator should be a gated extension. The thesis must remain successful even if the generator underperforms the coreset baseline.

---

# 12. What to ask Gaeun and Sophia

## Questions for Gaeun

1. **What exactly is an independent “target” in the output tree?** Is it a unique paired sequence, a PDB ID, an antibody–antigen complex or a pipeline run?

2. **How many unique paired VH/VL sequences and clonotypes remain after deduplication?**

3. **Are PH design, cycle, contact residues, binder size, AF3 seed/model, template, confidence and filter outcomes retained for every PDB?**

4. **Were true antigen-contact residues or antigen-derived binder sizes used for train and test examples in the reported antigen-prediction experiments?**

5. **What information will be available for a new BCR at deployment?** Sequence only, predicted structure, predicted paratope, or known antigen candidates?

6. **What is the actual optimization objective?**
   - preserve structural coverage;
   - preserve ConFormer embeddings;
   - preserve antigen retrieval;
   - preserve physical plausibility;
   - or some weighted combination?

7. **Does the complete ensemble materially outperform a repeated single conformer and a strict sequence-only baseline on the current hardest split?**

8. **Is conformer 0 a canonical reference, or is it arbitrary?**

9. **How sensitive are ConFormer outputs to ensemble ordering, duplicate conformers and variable \(N\)?**

10. **What does the current filtering remove, and are rejected structures available?** Rejected examples could train a useful quality model.

11. **Which component is actually expensive: PH generation, AF3 inference, storage, ConFormer training or all of them?**

12. **Can one dataset snapshot, pipeline version and downstream benchmark be frozen for the thesis?**

13. **Which downstream task is mature enough to be the primary endpoint—closed-set antigen classification, epitope retrieval, antigen-surface retrieval or disease-panel prediction?**

14. **Would a negative result—showing that 16 random conformers or one static structure is enough—be considered scientifically acceptable?**

## Questions for Sophia

1. **Is the lab’s priority computational cost reduction, a new generative model, or evidence that conformational ensembles improve antigen prediction?** These are different projects.

2. **What is the minimum contribution expected for the Master’s thesis if the generative model fails?**

3. **How should ownership be divided so this does not overlap ambiguously with Gaeun’s central PhD contribution?**

4. **Can the thesis be structured with compression/adaptive sampling as the guaranteed core and generation as an optional extension?**

5. **What GPU and AF3 budget is actually available for controlled ablations?**

6. **Is there support for maintaining a frozen dataset release, reproducible manifest and strict family-level split?**

7. **What claim language is acceptable regarding “conformations,” “pseudo-bound states,” “dynamics” and “antigen specificity”?**

8. **Would the lab support comparison against ABB4-STEROIDS even if it shows the existing PH/AF3 pipeline is unnecessary for some tasks?**

9. **What publication venue or audience is intended: structural bioinformatics, ML methodology, antibody engineering or computational immunology?**

10. **Can everyone agree now on a three-month go/no-go gate before committing the thesis to full generative modeling?**

11. **Will the downstream model remain stable long enough to evaluate the thesis, or is ConFormer itself still undergoing major redesign?**

12. **Are there institutional or licensing constraints on distributing AF3-derived structures, model weights or the resulting distilled model?**

---

# Final recommendation

### Green light

**Task-preserving compression, branch-aware sampling and adaptive stopping of PH/AF3 antibody ensembles.**

This directly addresses Gaeun’s question, is falsifiable, has strong non-neural baselines, can produce meaningful negative results and does not depend on a large generative model succeeding.

### Amber light

**A CDR-only residual flow or amortized set-of-conformer-prototypes model, preferably initialized from ABB4-STEROIDS or another pretrained geometric prior.**

This is scientifically interesting once you demonstrate that:

- the full ensemble adds genuine downstream value;
- a small oracle coreset exists;
- PH/AF3 provides useful modes beyond ABB4;
- the evaluation is blind to true antigen contacts and sizes.

### Red light

**Full all-atom antibody generation from noise, trained from scratch on the current targets, presented as a model of true biophysical dynamics.**

That is both too broad and now directly crowded by ABB4-STEROIDS.

My strongest advice is to commit only after a short saturation-and-leakage study answers three questions: **Does the ensemble actually help, can an oracle shrink it, and does pseudo-bound PH/AF3 beat ABB4 under blind conditioning?** If those answers are yes, the distillation thesis is compelling. If not, the rigorous negative result and adaptive compression work are still more valuable than forcing a diffusion model onto the problem.
