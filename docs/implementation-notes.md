# Implementation Notes

This file is a scaffold for future implementation. The repo currently contains
planning documents only.

## Expected Inputs

Target-level manifest fields:

- `target_id`;
- `vh_seq`;
- `vl_seq`;
- source database or PDB ID;
- antigen identity where available;
- antigen family where available;
- split assignment;
- sequence-family or clonotype group;
- CDR-H3 sequence and length;
- pipeline version.

Conformer-level manifest fields:

- `target_id`;
- `conformer_id`;
- PDB/CIF path;
- PH design index;
- PH cycle;
- contact-conditioning scheme;
- selected contact residues, if any;
- pseudo-binder length and identifier, if available;
- AF3 input JSON path;
- AF3 seed;
- AF3 model index;
- confidence scores;
- filter status;
- filter rejection reason, if rejected;
- structural hash;
- cluster ID.

## Known Local Context

Current Research workspace notes identify the conformer-generation optimization
question as: define a downstream ensemble-quality readout, then test how much
the Prot Hunt / AF3 ensemble can be subsampled while preserving conformational
variability and downstream prediction quality.

Known Jerry-owned conformer pipeline planning path:

```text
/Users/jerrygao/Research/liu-lab/bcr-conformer-pipeline
```

Known Jerry-owned cluster-side planning path:

```text
/project/liulab/jg1920/bcr-conformer-pipeline
```

Gaeun-owned source paths and generated archives should remain read-only unless
Jerry gives explicit approval and ownership is clear.

## Data Handling

Large structural files should not be committed here. This repo should track:

- manifests;
- scripts;
- small summaries;
- plots;
- reports;
- exact source paths and checksums where appropriate.

It should not track:

- bulk PDB/CIF conformers;
- AF3 output archives;
- model weights;
- private collaborator datasets;
- secrets or credentials.

## Analysis Modules To Build Later

Possible future scripts:

- `build_manifest.py`: construct target and conformer manifests from approved
  output roots.
- `cluster_conformers.py`: framework-align and cluster conformers by CDR and
  VH/VL features.
- `subsample_curves.py`: compute K-curves for random, stratified, and coreset
  policies.
- `representation_sensitivity.py`: test ConFormer order, duplication, and
  multiplicity sensitivity.
- `compute_accounting.py`: summarize PH/AF3 calls, runtime, GPU-hours, and
  storage.
- `adaptive_policy_pilot.py`: simulate or evaluate early-stopping policies.

Do not implement these until input manifests and allowed data roots are clear.

## Suggested Directory Layout For Future Work

```text
conffusion/
  docs/
  scripts/
  configs/
  manifests/
  reports/
  figures/
```

Bulk outputs should go under ignored `outputs/` or an external Jerry-owned
storage path, with summaries copied back into `reports/`.

## Reproducibility Bar

A result is "working" only when:

- the command completed;
- expected outputs exist;
- input manifest and code version are recorded.

A result is "reproduced" only when:

- a rerun or independent check matches closely enough;
- input, split, code version, and parameters are recorded.

A result is "ready to show" only when:

- caveats are written down;
- figure/table/report paths exist;
- no known blocker would make the result misleading.
