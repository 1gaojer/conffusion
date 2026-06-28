# Aim 1 Phase 1 Premise Report

Date: 2026-06-28

## Purpose

Phase 1 asked whether the Conffusion project has a real premise before moving
into task-aware compression or prospective PH/AF3 cost reduction.

The premise was:

> Gaeun's PH/AF3 conformer ensembles may contain CDR/paratope-local signal that
> matters for antigen-neighbor retrieval, but the signal must be read out and
> evaluated carefully before optimizing conformer count.

This report synthesizes the completed Phase 1 evidence. It is not a model
development report and not an external-validation claim.

## Bottom Line

Phase 1 supports continuing, but it changes the center of the project.

The strongest supported claim is:

> The current bottleneck is the ensemble readout and evaluation endpoint, not
> simply the number of generated conformers.

Observed:

- Global structure metrics were misleading because most large global RMSD came
  from non-antibody or complex-level motion.
- Antibody framework motion was much smaller, but H-CDR3 and combined CDRs
  retained meaningful structural diversity.
- Local H-CDR3 MCA residue embeddings contain some structural signal, but
  standard global H/L mean pooling mostly erases it.
- CDR/paratope-aware retrieval readouts repeatedly outperformed global H/L
  pooling.
- Framework-only, random framework-window, and shuffled-label controls were
  weaker.
- On the cleaner Tier B endpoint, all-CDR mean+std again beat global H/L
  pooling.
- On Tier C antigen/source-de-overlap, the same direction remained, but only 17
  same-label queries were evaluable, so the result is a leakage stress test
  rather than a stable validation endpoint.

Phase 1 does not yet support the stronger claim that structural k-center
selection has solved conformer compression. Random K64 remains competitive in
the retrieval endpoint.

## Evidence Ladder

| Step | Question | Main observation | Implication |
| --- | --- | --- | --- |
| Phase 1.1 global structural screen | Are copied PH/AF3 conformer files usable, and are they structurally diverse? | 149 targets and 19,072 conformers parsed cleanly; global mean pairwise RMSD was large. | Data were usable, but global RMSD was only a screening readout. |
| Phase 1.2 antibody-only structural screen | Does global diversity remain on antibody H/L chains? | Antibody/global RMSD ratio was small; paired H/L median mean pairwise RMSD was 0.525 A. | Global RMSD mostly measured non-antibody or complex-level motion. |
| Phase 1.3 CDR structural screen | Was whole-VH/VL RMSD averaging away loop diversity? | H-CDR3 median mean pairwise frame-aligned RMSD was 1.604 A; H-CDR3/whole-antibody ratio was 2.865. | The relevant structural signal is CDR-aware, especially H-CDR3. |
| Phase 1.4 structural coreset | Can CDR geometry be covered by fewer conformers? | K=32 and K=64 H-CDR3/all-CDR k-center improved structural coverage; K=8 was too aggressive. | CDR-aware compression is plausible structurally, with adaptive K likely needed. |
| Phase 1.4 embedding sensitivity | Does frozen MCA global embedding track CDR geometry? | H-CDR3 residue embeddings had weak positive signal; full H/L global pooling was essentially flat. | MCA is not completely blind, but global pooling dilutes local CDR signal. |
| Phase 1.4 richer readout audit | Is weak signal an encoder problem or a readout problem? | Flattened H-CDR3 residue readouts tracked geometry more clearly; global H/L remained flat. | The likely bottleneck is readout/pooling. |
| Phase 1.5 internal retrieval | Does a CDR-aware readout improve same-label retrieval? | Internal 149-target readout: H-CDR3 and all-CDR readouts beat global H/L pooling. | CDR-local signal is visible in retrieval-like behavior. |
| Phase 1.5d controls | Is the CDR result just any local-window effect? | Framework-only, random framework-window, and shuffled-label controls were weaker. | The effect is region-specific enough to pursue. |
| Tier B Step 4 | Does the result persist on a stricter Gaeun-generated endpoint? | 160 targets, 74 evaluable queries; all-CDR mean+std Recall@10 0.622 / MRR 0.356 versus global H/L Recall@10 0.446 / MRR 0.244. | The readout effect persisted on a cleaner powered endpoint. |
| Tier C leakage check | Does the result collapse after antigen/source de-overlap? | 87 retained targets, 17 evaluable queries; all-CDR mean+std Recall@10 0.471 / MRR 0.283 versus global H/L Recall@10 0.353 / MRR 0.130. | The direction survived as a leakage stress test, but the endpoint is too sparse for headline validation. |

## What Phase 1 Supports

### Supported Claim 1: CDR/paratope-aware readout matters

The most consistent result is that CDR/paratope-aware readouts expose more
same-label retrieval signal than global H/L pooling.

Evidence:

- Internal Phase 1.5c: H-CDR3 and all-CDR readouts beat global H/L pooling.
- Phase 1.5d: framework-only and random framework-window controls were weaker.
- Tier B Step 4: all-CDR mean+std remained stronger than global H/L pooling on
  a stricter generated endpoint.
- Tier C: the direction remained after exact antigen/source-pair de-overlap,
  although the evaluable set was small.

Interpretation:

- The MCA tensors likely contain useful local signal.
- The default global H/L pooling is not the right endpoint for deciding whether
  conformer ensembles are useful.

### Supported Claim 2: CDR structural diversity exists, but it is not global

The early structural analysis changed the premise from "large global ensemble
diversity" to "localized CDR/paratope diversity."

Evidence:

- Global conformer RMSD was large.
- Antibody H/L RMSD was much smaller.
- H-CDR3 diversity was much larger than whole-antibody RMSD suggested.

Interpretation:

- Whole-complex or whole-antibody structural metrics are too blunt.
- CDR-aware structure metrics should remain the structural baseline.

### Supported Claim 3: K=32 to K=64 is a reasonable structural baseline

K=32 and K=64 are useful structural compression baselines, especially for
H-CDR3/all-CDR geometry.

Evidence:

- Structural k-center improved CDR coverage at K=32 and especially K=64.
- K=8 did not reliably preserve H-CDR3 geometry.
- High-H-CDR3-motion targets remained hard even at K=64.

Interpretation:

- Fixed tiny K is not defensible.
- Adaptive K is probably necessary if compression becomes the main task.

## What Phase 1 Does Not Yet Support

### Not Supported: "Structural k-center solves compression"

Structural k-center improved structural coverage, but in retrieval it did not
clearly beat random K64.

Current interpretation:

- Geometry-only selection is not enough.
- Aim 2 should test task-aware, readout-aware, or adaptive selectors.

### Not Supported: "Tier C proves external generalization"

Tier C is cleaner for leakage, but sparse for the current exact-label retrieval
metric.

Observed:

- 87 targets retained after exact antigen/source-pair de-overlap.
- Only 17 retained targets had a same-label neighbor.

Current interpretation:

- Tier C is useful as a leakage stress test.
- A better-powered external endpoint needs an expanded bank, different grouping,
  or a different retrieval target definition.

### Not Supported: "A diffusion/generative model is justified now"

Phase 1 supports readout and endpoint work first. It does not yet justify a
new generator.

Current interpretation:

- A generative or distillation model remains gated on Aim 2 and Aim 3 evidence.
- The next technical work should stabilize the endpoint and selector objective.

## Decision For The Project

Proceed to Aim 2 only under a readout-first framing:

1. Treat all-CDR mean+std and H-CDR3/L-CDR3-style CDR readouts as the leading
   endpoints to test, with global H/L pooling kept as a control.
2. Keep Tier B as the main powered endpoint for immediate development.
3. Use Tier C as a caveated leakage-sensitivity result, not the main score.
4. Do not claim solved conformer compression until a selector beats random
   under the fixed CDR/paratope-aware retrieval endpoint.
5. Prioritize task-aware/adaptive subset selection over geometry-only k-center.

## Recommended Aim 2 Gate

Aim 2 should start with a fixed endpoint:

- primary readout: all-CDR mean+std;
- secondary readouts: H-CDR3, L-CDR3, all-CDR H/L mean;
- controls: global H/L mean, framework H/L mean, random framework windows,
  shuffled labels, single-first conformer, random K;
- current powered endpoint: Tier B Step 4;
- leakage caveat endpoint: Tier C exact antigen/source-pair de-overlap.

Aim 2 should answer:

> Can any subset selector preserve the fixed CDR/paratope-aware retrieval
> behavior better than random at the same K?

If yes, continue toward task-aware compression and prospective cost reduction.
If no, the project should focus on readout design, pooling, and endpoint
construction rather than conformer-generation optimization.

## Source Notes

Detailed evidence lives in:

- `notes/aim1-phase1-medium-results.md`
- `notes/aim1-phase1-2-antibody-structural-results.md`
- `notes/aim1-phase1-3-cdr-structural-results.md`
- `notes/aim1-phase1-4-cdr-coreset-results.md`
- `notes/aim1-phase1-4-embedding-sensitivity-results.md`
- `notes/aim1-phase1-4-readout-audit-results.md`
- `notes/aim1-phase1-5-retrieval-preservation-results.md`
- `notes/aim1-phase1-5d-region-control-results.md`
- `notes/aim1-phase1-5-tierb-step4-readout-results.md`
- `notes/aim1-phase1-5-tierc-leakage-eval-results.md`
