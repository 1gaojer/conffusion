# Conffusion

Conffusion is a documentation-first research workspace for a possible thesis
project on efficient antibody conformer ensemble modeling.

The working idea is not to train a full antibody diffusion model from scratch.
The immediate scientific target is narrower and now readout-first:

> Define a task-sufficient CDR/paratope-aware ensemble representation for
> antigen retrieval, then compress, select, or amortize Gaeun's expensive
> PH/AF3 pseudo-bound ensemble while preserving that representation and
> downstream retrieval value.

## Central Hypothesis

Large Protein Hunter / AlphaFold3 generated antibody conformer ensembles contain
substantial structural and representation-level redundancy, but the current
evidence says the main bottleneck is not only conformer count. CDR/paratope-local
signal is visible in MCA tensors and in retrieval, while global H/L mean pooling
can dilute that signal. A small, strategically selected or adaptively generated
subset is useful only if it preserves the right CDR/paratope-aware readout.

This project treats the generated conformers as a computational proposal
distribution, not as a calibrated physical equilibrium ensemble.

## Why This Repo Exists

Gaeun is currently interested in optimizing conformer generation. Jerry is
looking for a thesis idea that connects that practical bottleneck to modern
generative modeling, diffusion, flow matching, and representation learning.

The key question is:

> What CDR/paratope-aware ensemble representation is task-sufficient for
> antigen retrieval, and how cheaply can we generate or select it?

## Repository Map

- `docs/proposal.md`: thesis-style proposal with hypothesis and aims.
- `docs/background.md`: scientific context and rationale.
- `docs/research-plan.md`: staged plan, gates, and deliverables.
- `docs/aim1-phase1-premise-report.md`: Phase 1 synthesis report and go/no-go
  interpretation for moving into readout-fixed subset selection.
- `docs/aim1-phase1-benchmark-plan.md`: concrete data inventory and build plan
  for the first Aim 1 benchmark.
- `docs/aim1-phase1-4-representation-sensitivity-plan.md`: next-step plan after
  the CDR result, focused on whether MCA/ConFormer embeddings use H-CDR3
  diversity.
- `docs/aim1-phase1-5-retrieval-preservation-plan.md`: downstream retrieval and
  readout-preservation plan after the Phase 1.4 representation diagnostics.
- `docs/evaluation.md`: baselines, metrics, controls, and failure modes.
- `docs/implementation-notes.md`: manifest fields, data handling, and future
  implementation notes.
- `docs/questions-for-gaeun-and-sophia.md`: questions to resolve before
  committing.
- `docs/references.md`: papers, methods, and source links.
- `docs/source-notes.md`: notes on the GPT-Pro and Deep Research outputs.
- `docs/model-outputs/`: repo-local copies of the raw GPT-Pro and Deep
  Research markdown outputs.
- `manifests/`: Jerry-owned candidate manifests, including the current Tier B
  Gaeun-generated checkpoint-unseen pass set.
- `notes/model-report-synthesis.md`: synthesis of the three model reports.
- `notes/dataset-copies.md`: Ragon paths and provenance for Jerry-owned copied
  conformer datasets.
- `notes/aim1-phase1-medium-results.md`: first CPU-only structural saturation
  run on the medium copied PH/AF3 dataset.
- `notes/aim1-phase1-2-antibody-structural-results.md`: antibody-chain follow-up
  showing how much global diversity remains on heavy/light chains.
- `notes/aim1-phase1-3-cdr-structural-results.md`: CDR-focused follow-up showing
  that H-CDR3 carries substantially more diversity than whole VH/VL RMSD.
- `notes/aim1-phase1-5-tierb-sequence-deoverlap-results.md`: current Tier B
  sequence-de-overlapped Gaeun-generated test-pool result.
- `notes/aim1-phase1-5-tierb-step4-readout-results.md`: Tier B Step 4
  CDR-aware readout rerun and controls on the stricter Gaeun-generated endpoint.
- `figures/phase15_tierb_step4_20260627/`: ten PNG/SVG Tier B Step 4 figures,
  a contact sheet, figure manifest, and provenance files.
- `scripts/aim1_phase1_structural.py`: reproducible CPU-only parser,
  pairwise-RMSD, and saturation analysis script.
- `scripts/aim1_phase1_antibody_structural.py`: CPU-only heavy/light chain
  assignment and antibody-only RMSD/saturation analysis.
- `scripts/aim1_phase1_cdr_structural.py`: CPU-only IMGT CDR numbering,
  framework-aligned CDR RMSD, and CDR saturation analysis.
- `scripts/manifests/`: manifest-building helpers, including Tier B sequence
  de-overlap against the June 19/21 `1000` checkpoint universe.
- `notes/`: scratch notes that are not yet evidence.

## Current Position

Green light:

- CDR/paratope-aware ensemble readout analysis.
- Leakage-aware retrieval endpoint construction.
- Task-aware coreset and subset selection after the readout is fixed.
- Adaptive stopping and budget allocation for PH/AF3.
- Direct representation distillation as a baseline.

Amber light:

- CDR-focused residual flow or set-of-prototypes generator, only after the
  readout, endpoint, and task-aware subset behavior are established.

Red light:

- Full all-atom antibody generation from scratch on the current target count.
- Claims that PH/AF3 sample frequencies are physical occupancy probabilities.
- Claims that sequence-only generated ensembles add antigen information unless
  additional deployable conditioning variables are present.

## Immediate Next Step

Before committing to a generative thesis, finish the readout, saturation, and
leakage study around a clean endpoint:

1. Which CDR/paratope-aware readout should be treated as the retrieval endpoint?
2. Does the full ensemble outperform sequence-only, static-structure, and
   repeated-single-conformer baselines under that readout?
3. Can task-aware or adaptive subset selection shrink the ensemble without
   losing retrieval value?
4. Does blind PH/AF3 provide useful modes beyond fast antibody ensemble
   baselines such as ABB4-STEROIDS?

If those answers are positive, a distillation or lightweight generative
extension is justified. If not, the scientifically useful result is that the
current ensemble is redundant, leaky, or not yet exploited by the downstream
model.

Current cleaner Gaeun-generated candidate endpoint:

- Tier B pass manifest:
  `manifests/tierb_sequence_deoverlap_20260627/gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_pass_20260627.tsv`
- Definition: exact-PDB unseen relative to the full 930-ID checkpoint universe,
  then concatenated `VH+VL` global sequence identity `<0.85` to every checkpoint
  target.
- Remaining caveat: Tier B does not enforce antigen/source de-overlap. The
  Tier C leakage-sensitivity rerun enforces exact antigen/source-pair
  de-overlap but is sparse under the current same-label retrieval metric.

Current Tier B readout result:

- Step 4 usable endpoint: 160 targets, 100 conformers per target, 74 evaluable
  self-match-excluded retrieval queries.
- all-CDR mean+std outperformed global H/L pooling
  (Recall@10 0.622 / MRR 0.356 versus Recall@10 0.446 / MRR 0.244).
- L-CDR3 was also strong; framework/random-window and shuffled-label controls
  were weaker.
- Structural K64 k-center did not clearly beat random K64, so the current
  strongest claim is CDR/paratope-aware readout improvement rather than solved
  conformer compression.
- Result note:
  `notes/aim1-phase1-5-tierb-step4-readout-results.md`
- Figure pack:
  `figures/phase15_tierb_step4_20260627/`

Current Tier C leakage-sensitivity result:

- Exact antigen/source-pair de-overlap kept 87 of the 160 Tier B Step 4 targets
  and left 17 self-match-excluded evaluable same-label retrieval queries.
- CPU-only job `32258` completed on `liulab` / `ragonliu1` in 20 seconds.
- all-CDR mean+std remained above global H/L pooling on the retained subset
  (Recall@10 0.471 / MRR 0.283 versus Recall@10 0.353 / MRR 0.130), and
  shuffled-label controls were lower.
- Interpretation: useful as a leakage stress test, but too sparse to treat as a
  headline external-validation endpoint.
- Result note:
  `notes/aim1-phase1-5-tierc-leakage-eval-results.md`

Current Phase 1 premise decision:

- Phase 1 supports continuing under a readout-first framing.
- Strongest claim: CDR/paratope-aware readouts expose useful retrieval signal
  that global H/L pooling dilutes.
- Main caveat: structural k-center has not yet clearly beaten random K64 under
  the fixed retrieval readout, and Tier C is too sparse for a headline external
  validation claim.
- Premise report:
  `docs/aim1-phase1-premise-report.md`

Current redesigned aims:

1. Establish a leakage-aware antigen-retrieval endpoint and CDR/paratope-aware
   ensemble readout.
2. Identify the smallest task-sufficient ensemble or subset under that fixed
   readout.
3. Reduce PH/AF3 generation cost prospectively using the readout and subset
   evidence.

## Ownership And Scope

This repo is a Jerry-owned planning and documentation workspace. It should not
modify Gaeun-owned source code, shared generated conformers, cluster jobs, or
collaborator outputs. External assets should be read-only unless copied into a
Jerry-owned workspace with clear provenance.
