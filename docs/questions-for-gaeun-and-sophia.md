# Questions For Gaeun And Sophia

## Questions For Gaeun

1. What exactly defines one independent target in the conformer output tree:
   paired sequence, PDB entry, antibody-antigen complex, or pipeline run?
2. How many unique paired VH/VL sequences remain after deduplication?
3. How many clonal families, CDR-H3 families, and antigen families are present?
4. Which PH/AF3 parameters are current for the pipeline version we should use:
   number of designs, cycles, AF3 seeds, and models per seed?
5. Are PH design, cycle, contact residues, pseudo-binder size, AF3 seed/model,
   template state, confidence metrics, and filter outcomes retained per PDB?
6. Were true antigen contacts, epitope residues, cognate-antigen size, or
   complex-derived templates used for any ensembles later evaluated on antigen
   prediction or retrieval?
7. What information will be available for a new BCR at deployment: sequence
   only, predicted structure, predicted paratope, known antigen candidates, or
   known cognate antigen?
8. Does the complete ensemble materially outperform sequence-only, static
   structure, and repeated-single-conformer baselines on the current hardest
   split?
9. Is conformer 0 a canonical reference or arbitrary?
10. How sensitive are ConFormer outputs to conformer ordering, duplicated
    conformers, and variable ensemble size?
11. What does the current filtering remove, and are rejected structures
    available for quality-model training?
12. Which component is actually expensive: PH generation, AF3 inference,
    filtering, storage, ConFormer training, or all of them?
13. Can one dataset snapshot, pipeline version, split, and downstream benchmark
    be frozen for this thesis?
14. Which downstream endpoint is mature enough to be primary: nearest-neighbor
    retrieval, antigen-surface retrieval, closed-set antigen classification,
    disease-panel ranking, or something else?
15. What non-inferiority margin would be acceptable for a compressed ensemble?
16. Are there experimental multistate, apo/holo, MD, or ALL-conformations
    groups that can serve as independent validation?
17. Would a rigorous negative result, such as "16 random conformers are enough"
    or "one static structure is enough," be scientifically acceptable?

## Questions For Sophia

1. Is the lab priority compute reduction, a new generative model, or evidence
   that conformer ensembles improve antigen retrieval?
2. What is the minimum contribution expected if the generative model fails?
3. Can the thesis be structured with compression/adaptive sampling as the core
   and generation as an optional extension?
4. How should ownership be divided so this does not overlap ambiguously with
   Gaeun's PhD contribution?
5. What GPU and AF3 budget is available for prospective reduced-generation
   ablations?
6. Can we freeze one dataset snapshot, pipeline version, split, and primary
   metric for thesis evaluation?
7. What claim language is acceptable around "conformations", "pseudo-bound
   states", "dynamics", and "antigen specificity"?
8. Should the project explicitly benchmark against ABodyBuilder4-STEROIDS or
   other fast antibody ensemble models, even if they make part of the PH/AF3
   pipeline unnecessary?
9. What venue or audience matters most: structural bioinformatics, ML methods,
   antibody engineering, or computational immunology?
10. Is a three-month go/no-go gate acceptable before committing thesis effort
    to a diffusion or flow model?
11. Are there licensing or collaboration constraints on publishing AF3-derived
    conformer summaries, distilled model weights, or benchmark splits?
12. If the project shows that an external antibody ensemble model is sufficient,
    would that be acceptable as a thesis direction or should the thesis remain
    centered on Gaeun's PH/AF3 pipeline?

## Questions For Jerry

1. Is the personal thesis goal primarily a methods thesis, a co-scientist
   support contribution, or a bridge between both?
2. How much of the project should emphasize generative modeling as intellectual
   identity versus compute reduction as lab utility?
3. What result would feel satisfying even if the learned generator never wins?
