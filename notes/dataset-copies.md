# Dataset Copies

This note tracks Jerry-owned copies of conformer data made from Gaeun-owned
read-only source trees. The copied datasets are not stored in this Git repo;
they live on Ragon under Jerry-owned `/external/liulab/jg1920/...` paths.

## Gaeun PH/AF3 Pilot Copy

Created: 2026-06-26

Destination:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_pilot_20260626
```

Source paths, treated read-only:

```text
/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results/filter_metrics
/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results/af3_results
/external/liulab/gkim/antigen_prediction/datasets/seq_files/full_sabdab_seqsim_split_test_seq_file.json
```

Selection policy:

- top 20 sequence-backed targets by kept conformer count;
- 128 kept conformers per target;
- conformers selected evenly across sorted kept rows by `run_id`, `cycle_id`,
  and source path.

Verified contents:

- 20 selected targets;
- 2,560 copied CIF files;
- 20 selected filter metric CSVs;
- no `.partial` directory left after validation;
- manifest records `source_writes_performed: false`.

Key manifests:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_pilot_20260626/manifests/source_summary.json
/external/liulab/jg1920/conffusion/gaeun_ph_af3_pilot_20260626/manifests/pilot_targets.tsv
/external/liulab/jg1920/conffusion/gaeun_ph_af3_pilot_20260626/manifests/pilot_conformers.tsv
```

## Gaeun PH/AF3 Medium Copy

Created: 2026-06-26

Destination:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626
```

Source paths, treated read-only:

```text
/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results/filter_metrics
/external/liulab/gkim/antigen_prediction/datasets/conf_gen_results/af3_results
/external/liulab/gkim/antigen_prediction/datasets/seq_files/full_sabdab_seqsim_split_test_seq_file.json
```

Selection policy:

- all sequence-backed targets with at least 128 kept conformers;
- 128 kept conformers copied per target;
- conformers selected evenly across sorted kept rows by `run_id`, `cycle_id`,
  and source path.

Verified contents:

- created on `ragonliu1`;
- 149 selected targets;
- 128 conformers per target;
- 19,072 copied CIF files;
- 149 selected filter metric CSVs;
- copied CIF bytes recorded in manifest: 5,240,943,122;
- apparent destination size by `du -sh`: 2.7G;
- selected targets had 154 to 2,294 kept conformers in the source metrics,
  median 1,968;
- no `.partial` directory left after validation;
- manifest records `source_writes_performed: false`.

Key manifests:

```text
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626/manifests/source_summary.json
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626/manifests/medium_targets.tsv
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626/manifests/medium_conformers.tsv
/external/liulab/jg1920/conffusion/gaeun_ph_af3_medium_20260626/manifests/all_filter_metrics_index.tsv
```

## Caveats

- These are SAbDab/PDB-like PH/AF3 ensemble copies from Gaeun's existing output
  tree, not the strict-300 Covid/HIV/Influenza benchmark.
- The copied conformers should be used for Aim 1 Phase 1 manifest/catalog
  development, structural diversity checks, and coreset/subsampling experiments.
- Do not use them to claim disease-specific antigen retrieval performance
  without a separate endpoint mapping and leakage-aware benchmark.
