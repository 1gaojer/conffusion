# Conffusion Agent Instructions

This repo is a Jerry-owned planning workspace for a possible thesis project on
efficient antibody conformer ensemble modeling.

## Operating Rules

- Keep the project grounded in Gaeun's conformer-generation optimization need.
- Treat Gaeun-owned code, data, outputs, jobs, and shared lab assets as
  read-only.
- Do not edit, delete, move, relaunch, cancel, or otherwise control
  collaborator-owned files, outputs, jobs, or runs.
- Do not store secrets, tokens, private keys, or broad credentials in this
  repo.
- Separate observed evidence from interpretation, hypothesis, and speculation.
- Do not claim that PH/AF3 conformers are physical equilibrium samples unless
  external physical or experimental validation supports that claim.
- Prefer small, controlled evaluations before model-building.
- For present-tense pipeline claims, verify the live repo or run configuration.
- For Python work in Liu Lab repos, use the shared configured lab environment
  rather than creating repo-local virtual environments.

## Scientific Guardrails

- The core thesis should be ensemble compression, adaptive sampling, and
  task-preserving distillation.
- Generative modeling is a gated extension, not the first deliverable.
- The first baseline to beat is smart subsampling, not diffusion.
- Train/test splits must be target-level or family-level, never conformer-level.
- Downstream antigen retrieval claims require leakage-aware baselines and hard
  negatives.

## Biological Payload And Log Hygiene

- Treat sequence manifests, conformer manifests, antigen/source tables,
  embedding-export metadata, PDB/mmCIF-derived tables, and generated target
  records as biological payloads when they contain or can expose raw sequences.
- Do not print raw `vh_seq`, `vl_seq`, antigen sequences, FASTA/MSA records,
  full JSON records, or residue dumps in stdout/stderr, notes, or chat unless
  Jerry explicitly asks and the direct payload is necessary.
- Default diagnostics should use safe summaries: target IDs, row counts,
  source labels, split labels, sequence lengths, conformer counts, missing-field
  names, checksums, metrics, and paths.
- When adding analysis scripts or Slurm wrappers, make logs path/count/metric
  oriented. Write sequence-bearing outputs to files and refer to those files by
  path plus checksum or manifest row ID.
- Avoid `cat`, broad `head`, broad `tail`, broad `grep`, `sed`, or unfiltered
  `jq` on sequence-bearing files. Prefer scripts that explicitly drop sequence
  columns before printing.

## Documentation Style

- Keep docs concise but specific.
- Record exact paths, run IDs, manifests, and versions when they become known.
- Label unverified claims from model reports or presentations as "reported" or
  "to verify".
- Keep papers and source links in `docs/references.md`.
