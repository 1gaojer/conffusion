# Conffusion Aim Redesign

Date: 2026-06-28

## Reason For Update

The original project framing emphasized "how many antibody conformers are
enough?" and treated conformer compression as the central next move after the
premise test.

The completed Phase 1.4 and Phase 1.5 diagnostics changed the emphasis:

- structural CDR diversity exists, especially in H-CDR3;
- local CDR residue embeddings show signal that global H/L mean pooling can
  wash out;
- Tier B Step 4 showed CDR/paratope-aware readouts outperform global H/L pooling
  on a stricter Gaeun-generated endpoint;
- structural K64 k-center did not clearly beat random K64 under the strongest
  current readout.

Therefore the project should be readout-first before it is compression-first.

## Updated North Star

Define a task-sufficient CDR/paratope-aware ensemble representation for antigen
retrieval, then compress, select, or amortize Gaeun's expensive PH/AF3
pseudo-bound ensemble while preserving that representation and downstream
retrieval value.

## Updated Aims

1. Establish a leakage-aware antigen-retrieval endpoint and CDR/paratope-aware
   ensemble readout.
2. Identify the smallest task-sufficient ensemble or subset under that fixed
   readout.
3. Reduce PH/AF3 generation cost prospectively using the readout and subset
   evidence.

Optional later extension:

- Distill a small conformer prototype generator only after the readout, endpoint,
  and subset behavior are established.

## Current Interpretation

The strongest current claim is CDR/paratope-aware readout improvement, not solved
geometry-only conformer compression. The next decision is whether Tier C
antigen/source de-overlap is needed before presenting Tier B as cleaner
validation, then whether task-aware/adaptive selectors improve over random under
the fixed readout.
