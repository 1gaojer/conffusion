# Aim 1 Phase 1.5d Region-Control Results

Date: 2026-06-27

## Scope

Phase 1.5d tested whether the Phase 1.5c CDR-aware retrieval improvement is
specific to CDR/paratope regions, or whether any smaller/local region would
look similarly good.

This is still an internal 149-target diagnostic using saved MCA tensors. It is
not an independent validation endpoint.

## Execution

Cluster scripts:

- `scripts/phase15a/phase15d_region_controls.py`
- `scripts/phase15a/phase15d_region_controls.sbatch`

Cluster outputs:

- `/project/liulab/jg1920/conffusion/phase15d_20260627_region_controls`

Local copied outputs:

- `figure_inputs/phase15d_20260627/`

Jobs:

- Initial preflight: `32214`
- Corrected preflight: `32215`
- Full diagnostic: `32216`

The initial preflight caught a light-chain mapping issue: kappa light chains
were recorded as `chain_type=K` in the CDR assignment table, while MCA tensor
keys use `L`. Phase 1.5d normalizes light-chain assignments through
`chain_role=light -> L` before pooling.

## Region Extraction Check

The corrected preflight showed plausible CDR ranges on the first five targets.
Example for `1a2y`:

| region | positions |
| --- | --- |
| H-CDR1 | 26-33 |
| H-CDR2 | 51-57 |
| H-CDR3 | 96-105 |
| L-CDR1 | 27-32 |
| L-CDR2 | 50-52 |
| L-CDR3 | 89-97 |

In the full audit, all intended CDR/framework/random-window rows were non-empty
except two `L-CDR2` rows (`5d70`, `5d71`). This should be treated as a minor
region-assignment caveat, not a full-run failure.

## Main Result

Metrics used 78 evaluable queries after self-match exclusion and positive-label
filtering.

| readout | Recall@10 | MRR | median first positive rank |
| --- | ---: | ---: | ---: |
| H-CDR3 | 0.872 | 0.709 | 1.0 |
| all CDRs mean+std | 0.897 | 0.690 | 1.0 |
| L-CDR1 | 0.795 | 0.644 | 1.0 |
| all CDRs H/L mean | 0.859 | 0.613 | 2.0 |
| H all CDRs | 0.846 | 0.594 | 2.0 |
| global H/L mean | 0.718 | 0.403 | 4.0 |
| framework H/L mean | 0.718 | 0.391 | 5.0 |
| random H-framework windows, mean over 10 | 0.577 | 0.295 | 8.0 |

Shuffled-label controls over 20 replicates were much lower:

| readout | shuffled Recall@10 mean | shuffled MRR mean |
| --- | ---: | ---: |
| global H/L mean | 0.348 | 0.141 |
| H-CDR3 | 0.344 | 0.133 |
| all CDRs H/L mean | 0.337 | 0.133 |
| all CDRs mean+std | 0.350 | 0.135 |
| framework H/L mean | 0.352 | 0.139 |

## Interpretation

Observed:

- H-CDR3 and all-CDR mean+std remain much stronger than global H/L pooling.
- Framework-only pooling does not improve over global pooling.
- Random framework windows matched to H-CDR3 length are weaker than true CDR
  regions.
- Shuffled labels collapse toward chance-like behavior for this label
  distribution.
- Light-chain signal is not empty after fixing the chain normalization; L-CDR1
  is surprisingly strong in this internal set.

Interpretation:

- This strengthens the CDR/paratope-readout hypothesis. The result is less
  likely to be just "small local window beats global pooling."
- The useful antigen-neighbor signal appears concentrated in biologically
  plausible antigen-facing regions, especially H-CDR3 and combined CDR
  variability features.
- The readout bottleneck remains the strongest current explanation: MCA tensors
  contain useful local signal, while global H/L pooling dilutes it.

Important caveat:

- This is still an internal diagnostic on the 149-target bank. It should be
  repeated on strict-300 or a de-overlapped benchmark before being treated as an
  external validation result.

## Gaeun-Facing Takeaway

Using the same MCA tensors and the same internal retrieval setup, CDR-aware
pooling substantially improves same-antigen nearest-neighbor retrieval compared
with global H/L mean pooling. Framework-only and random framework-window
controls do not reproduce the effect, and shuffled labels drop strongly. This
suggests the improvement is tied to CDR/paratope-local representation signal,
not merely to using a smaller region.

## Next Step

Repeat the CDR-aware readout diagnostic on a stricter or de-overlapped target
set once conformers/embeddings are available. In parallel, treat H-CDR3 and
all-CDR mean+std as the most promising readouts for testing conformer subset
preservation.
