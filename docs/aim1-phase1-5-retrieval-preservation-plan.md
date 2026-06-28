# Aim 1 Phase 1.5 Retrieval Preservation Plan

Date: 2026-06-27

Purpose: define the next Aim 1 gate after Phase 1.4. Phase 1.4 showed that
Gaeun's copied PH/AF3 ensemble contains real CDR structural diversity and that
MCA1000 sees some H-CDR3 signal in local residue embeddings, while the usual
full H/L global mean readout mostly averages that signal away.

Phase 1.5 should test whether that conformer diversity matters for the
downstream endpoint Gaeun currently cares about: nearest-neighbor antigen
retrieval over MCA/ConFormer representations.

## Current Status Entering Phase 1.5

Observed Phase 1.4 evidence:

- Structural CDR coreset analysis supports K=32 to K=64 as the first useful
  H-CDR3/all-CDR compression range.
- K=8 is too aggressive for reliable H-CDR3 coverage.
- Local H-CDR3 residue embeddings show a modest structural signal.
- Full H/L global mean embeddings are essentially flat with respect to H-CDR3
  structural variation.
- The richer readout audit supports a pooling/readout bottleneck hypothesis,
  not a downstream retrieval claim.

Confirmed guidance entering this phase:

- Gaeun confirmed that the correct checkpoint is the June 19
  `1000 conformer` checkpoint.
- At Phase 1.5 entry, the exact checkpoint path/hash and matching retrieval or
  antigen bank still needed live verification before any retrieval result could
  be called clean.

Previously used local checkpoint copy:

```text
/project/liulab/jg1920/conffusion/phase14_20260626_2301/shared/checkpoints/mca1000_checkpoint.pt
```

Phase 1.5a must verify whether this file is the confirmed June 19 checkpoint or
whether embeddings need to be re-exported.

Execution update, 2026-06-27:

- Phase 1.5a verified that the staged Phase 1.4 checkpoint copy matches the
  Gaeun-confirmed June 19 `1000 conformer` checkpoint.
- Phase 1.5b ran the global endpoint failure/random-control analysis.
- Phase 1.5c ran a CDR/paratope-aware readout diagnostic using saved MCA tensors.
- Results are recorded in
  `notes/aim1-phase1-5-retrieval-preservation-results.md`.
- Main result: the current global H/L mean retrieval endpoint shows full
  ensemble benefit but does not make CDR k-center clearly better than random;
  CDR-aware readouts recover much stronger internal antigen-neighbor structure.
- Follow-up Step 4 reran the CDR-aware readout and controls on a stricter Tier B
  Gaeun-generated endpoint. Results are recorded in
  `notes/aim1-phase1-5-tierb-step4-readout-results.md`.

## Phase 1.5 Question

Primary question:

> Does the full PH/AF3 conformer ensemble, or a smaller CDR-aware subset,
> preserve or improve nearest-neighbor antigen retrieval compared with cheap
> controls?

Secondary question:

> If retrieval is insensitive to CDR-aware structural subsets, is that because
> conformer diversity is unimportant, or because the current global MCA readout
> washes out local CDR signal?

## Hypotheses

### H1: CDR-Aware Subsets Preserve Retrieval

If K=32 or K=64 H-CDR3/all-CDR subsets recover most of the full-ensemble
retrieval behavior and beat naive controls, then Phase 2 coreset/compression is
justified.

Expected next move:

- materialize reusable coreset manifests;
- scale retrieval preservation to the strict-300 endpoint when conformers are
  ready;
- test prospective reduced generation or adaptive stopping.

### H0: Retrieval Is Flat Across Conditions

If full, random K, structural K, and repeated-one controls perform similarly,
then the current retrieval endpoint may not use conformer diversity strongly.

Expected next move:

- diagnose whether retrieval uses global H/L mean embeddings that average away
  local CDR signal;
- test CDR/paratope-aware pooling or adapter readouts before optimizing
  generation;
- avoid claiming task-preserving conformer compression.

### H2: Full Ensemble Helps But K=32/K=64 Loses Signal

If the full ensemble helps but compact subsets fail, then either rare conformer
modes matter or the current selection metrics miss task-relevant modes.

Expected next move:

- inspect target-specific failures;
- test adaptive K;
- compare structural k-center against readout-aware or retrieval-aware
  selectors.

## Conditions To Compare

Use the same target-level query set wherever possible.

Required conditions:

- full ensemble;
- K=64 H-CDR3 k-center;
- K=64 all-CDR k-center;
- K=32 H-CDR3 k-center;
- K=32 all-CDR k-center;
- random K=32 and K=64 controls, repeated across multiple seeds;
- first K=32 and K=64 controls;
- one conformer repeated to the full ensemble size;
- single/static structure baseline where available;
- sequence-only baseline if the endpoint supports it.

Optional later conditions:

- CDR-flattened local readout;
- CDR/paratope-aware mean plus variance readout;
- pair-representation region readout;
- attention-pooling readout, if exposed by a later GPU re-encode.

## Metrics

Primary retrieval metrics:

- Recall@1;
- Recall@5;
- Recall@10;
- MRR;
- mean average precision, if multiple positives exist per query.

Diagnostics:

- per-antigen-group performance;
- per-disease performance for strict-300-style disease labels;
- query-to-bank self-match handling;
- antigen group dominance or shortcut behavior;
- split/family leakage checks;
- performance difference between full, structural K, random K, first K, and
  repeated-one controls.

## Phase 1.5a: Endpoint Feasibility And Freeze

Phase 1.5a is the immediate next step. It should not be a large science run.
It should make the endpoint reproducible and decide whether a Phase 1.5 pilot
is meaningful.

### 1. Verify The Confirmed Checkpoint

Goal:

- identify the exact June 19 `1000 conformer` checkpoint path;
- compute or record its hash, modification time, and size;
- compare it against the local Phase 1.4 checkpoint copy;
- record whether current Phase 1.4 embeddings can be reused.

Deliverable:

```text
docs/aim1-phase1-5a-endpoint-audit.md
```

Minimum fields:

- confirmed checkpoint path;
- hash/size/mtime;
- whether it matches the local copied checkpoint;
- embedding export path to reuse, if matched;
- re-export requirement, if not matched.

### 2. Audit Retrieval Bank And Target Overlap

Goal:

- find the matching antigen/retrieval bank for the confirmed checkpoint;
- inspect required fields and representation format;
- count overlap with the 149-target copied PH/AF3 medium set;
- count overlap with strict-300 targets, if metadata are available;
- identify whether each query has a usable antigen truth label or group;
- flag train/test/retrieval-bank leakage risks.

Key question:

> Is the 149-target medium set usable for retrieval preservation, or is it only
> usable for endpoint plumbing?

If overlap is too small, use Phase 1.5a as a plumbing audit and defer the main
retrieval claim to strict-300 or another endpoint set.

### 3. Freeze The Retrieval Definition

Record:

- query representation;
- candidate bank representation;
- distance metric;
- positive label definition;
- negative candidate policy;
- self-match exclusion policy;
- split or family-level leakage policy;
- metric table schema.

The default recommendation is cosine nearest-neighbor retrieval over the
confirmed MCA1000 representation, with explicit self-match exclusion and
target-level leakage reporting.

### 4. Materialize Subset Manifests

Build target-level condition manifests for:

- full 128 conformers;
- K=64 H-CDR3 k-center;
- K=64 all-CDR k-center;
- K=32 H-CDR3 k-center;
- K=32 all-CDR k-center;
- K=32/K=64 random controls with seeds;
- K=32/K=64 first-K controls;
- repeated-one controls.

Use existing Phase 1.4 CDR coreset outputs where possible:

```text
/external/liulab/jg1920/conffusion/aim1_phase1_4_cdr_coresets_20260627
```

### 5. Run A Tiny Retrieval Smoke

Use a small number of targets first, for example 5 to 10, only after the
checkpoint and retrieval bank are verified.

Smoke goals:

- confirm query embeddings load;
- confirm candidate bank embeddings load;
- confirm ranking runs end to end;
- confirm metrics are nonempty;
- confirm self-match/leakage policy is applied;
- confirm output tables are interpretable.

Smoke deliverable:

```text
/project/liulab/jg1920/conffusion/phase15a_YYYYMMDD_retrieval_smoke/
```

Expected files:

- `endpoint_audit.json`;
- `target_overlap.tsv`;
- `condition_manifest.tsv`;
- `retrieval_smoke_results.tsv`;
- `retrieval_smoke_summary.md`.

## Compute Requirements

Phase 1.5a is CPU-first.

CPU-only tasks:

- checkpoint path/hash audit;
- retrieval-bank inventory;
- target overlap counting;
- split/leakage checks;
- subset manifest materialization;
- nearest-neighbor ranking from already-exported embeddings;
- metric aggregation and report writing.

GPU is needed only if:

- the confirmed June 19 checkpoint does not match the existing Phase 1.4
  embedding export and embeddings must be re-exported;
- strict-300 conformers need fresh MCA/ConFormer embedding inference;
- subset-specific model outputs must be recomputed instead of averaged from
  saved per-conformer embeddings;
- layer-wise conformer weights, attention-pooling outputs, or subset-specific
  `pair_repr` are required.

Practical default:

- Start Phase 1.5a as CPU-only.
- Escalate to a GPU job only after the audit shows that current embeddings are
  missing, stale, or insufficient for the endpoint.

## Phase 1.5 Pilot Gate

Proceed from 1.5a to a Phase 1.5 pilot only if:

- the confirmed checkpoint is resolved to an exact path/hash;
- query and bank embeddings can be loaded reproducibly;
- there are enough target overlaps for a meaningful pilot;
- positive labels and self-match exclusion are defined;
- a tiny retrieval smoke produces nonempty metrics.

Recommended first pilot:

- use the copied 149-target PH/AF3 medium set if overlap is adequate;
- otherwise use it only for pipeline plumbing and wait for strict-300 or a
  better matching endpoint set for scientific claims.

## Interpretation Rules

Do not overclaim:

- Retrieval preservation is required before claiming task-preserving conformer
  compression.
- Embedding preservation alone is not enough.
- A small 149-target diagnostic may be useful but is not automatically a clean
  external validation set.
- If the retrieval bank overlaps training data, report it as an internal
  diagnostic rather than independent validation.

Strong result:

- full ensemble beats one/repeated/static controls;
- K=32 or K=64 CDR-aware subsets preserve most full-ensemble retrieval;
- CDR-aware subsets beat random or first-K controls under the same budget.

Weak or negative result:

- all conditions perform similarly;
- random K matches structural K;
- retrieval is dominated by antigen/source shortcuts;
- results collapse under leakage-aware filtering.

## Decisions Needed Before Execution

Default recommendation:

- Start Phase 1.5a with the copied 149-target PH/AF3 medium set because it is
  ready now and already has structural coresets and MCA1000 embeddings.

Decision to confirm:

- Should Phase 1.5a use the 149-target copied medium set as the first endpoint
  feasibility audit, with strict caveats about overlap and leakage?

Decision status:

- Executed with the 149-target copied medium set as an internal diagnostic.
- Keep strict caveats: this is endpoint plumbing and readout evidence, not an
  independent validation result.
- Next decision is whether to package the CDR-aware readout result for Gaeun and
  repeat it on strict-300 or another cleaner non-overlapped target set once
  conformers/embeddings exist.

## Strict Gaeun-Generated Endpoint Update

Date: 2026-06-27

Tier A:

- Read-only inventory of Gaeun-owned generated conformer outputs produced a
  359-target candidate set that excludes the full 930-ID June 19/21 `1000`
  checkpoint PDB universe.
- Manifest:
  `manifests/gaeun_conformer_ensembles_generated_non_1000_all_20260627.tsv`

Tier B:

- Tier B filters the 359 Tier A targets by Gaeun-style concatenated `VH+VL`
  global sequence identity against all 930 checkpoint-universe IDs.
- Rule: keep only targets with max `VH+VL` identity `<0.85`.
- Result: 172 pass, 187 fail, zero missing identity values.
- Pass manifest:
  `manifests/tierb_sequence_deoverlap_20260627/gaeun_conformer_ensembles_generated_non_1000_all_tierB_vhvl85_pass_20260627.tsv`
- Result note:
  `notes/aim1-phase1-5-tierb-sequence-deoverlap-results.md`

Use Tier B as the current Gaeun-generated checkpoint-unseen, antibody-sequence
de-overlapped candidate pool. It is not yet antigen/source de-overlapped, so it
should still be described as cleaner than the internal 149-target diagnostic,
not as final external validation.

## Tier B Step 4 Result Update

Date: 2026-06-27

Step 4 reran the Phase 1.5c/1.5d readout and controls on the Tier B endpoint.

- Usable endpoint: 160 targets from the 172-target Tier B pass manifest.
- Conformers per target: 100.
- Evaluable self-match-excluded retrieval queries: 74.
- Full-ensemble all-CDR mean+std readout: Recall@10 0.622, MRR 0.356.
- Full-ensemble global H/L mean readout: Recall@10 0.446, MRR 0.244.
- L-CDR3 was also strong: Recall@10 0.595, MRR 0.346.
- Framework/random-window and shuffled-label controls were weaker.
- Structural K64 k-center remained close to the full result but did not clearly
  beat random K64.

Interpretation:

- The CDR/paratope-aware readout finding survives on the stricter Tier B
  endpoint.
- The result supports a readout/pooling bottleneck more strongly than it
  supports geometry-only conformer compression.
- Tier C antigen/source de-overlap is still needed before presenting this as
  clean external validation.

Immediate follow-ups:

- Make Tier B Step 4 figures.
- Check why Step 4 used 160 of the 172 Tier B pass targets.
- Decide whether Tier C antigen/source de-overlap is required before sharing as
  a stronger validation endpoint.
