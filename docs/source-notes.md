# Source Notes

This file records how the three model outputs should be treated as sources.

## Copied Model Outputs

The raw outputs have been copied into this repository so they are not dependent
on transient Codex attachment paths.

### Deep Research Output Without Presentation Context

Path:

```text
docs/model-outputs/deep-research-without-presentation-context.md
```

Assessment:

- Contains many citation markers.
- Most are ChatGPT internal handles such as `turn12view0` rather than usable
  bibliographic links.
- Useful as a compact conceptual synthesis.
- Not sufficient as a source list without reconstructing the cited papers.

Best-use claims:

- distillation/compression is stronger than a true-dynamics claim;
- effective sample size is target count, not conformer file count;
- leakage and downstream validation are central;
- smart subsampling is the first baseline to beat.

### GPT-Pro Output With Presentation Context

Path:

```text
docs/model-outputs/gpt-pro-with-presentation-context.md
```

Assessment:

- Contains direct literature links and file-citation markers.
- The file-citation markers refer to the presentation/context material supplied
  during that model run.
- Most useful for project-specific issues:
  - PH/AF3 pipeline arithmetic as reported in the presentation context;
  - ConFormer mean-pooling / order-sensitivity concerns;
  - possible true-contact or binder-size leakage;
  - distinction between post-hoc compression and prospective compute savings.

Important caveat:

- Presentation-derived pipeline numbers must be rechecked against the live
  pipeline before being treated as current. The local Jerry wrapper currently
  defaults to `NUM_DESIGNS=20`, `NUM_CYCLES=5`, and `AF3_NUM_SEEDS=5`, while
  presentation/context notes may describe a 25-design variant.

### GPT-Pro Output Without Presentation Context

Path:

```text
docs/model-outputs/gpt-pro-without-presentation-context.md
```

Assessment:

- Contains direct literature links and some file-citation markers.
- The direct links are more usable than the Deep Research citation handles.
- Best contribution is the novelty warning around ABB4-STEROIDS and the
  recommendation to avoid a generic antibody ensemble generator thesis.
- Any presentation-like claims in this file should be treated cautiously unless
  mapped back to an actual supplied file.

Best-use claims:

- ABB4-STEROIDS makes generic paired VH/VL-to-ensemble generation less novel;
- PH/AF3 pseudo-bound distillation remains a distinct niche;
- first three decisive tests should be:
  1. does the full ensemble help;
  2. can an oracle shrink it;
  3. does blind PH/AF3 beat an antibody ensemble baseline.

## Source Handling Rule

For planning:

- Use model outputs as idea generators and checklists.

For thesis claims:

- Use primary papers, live pipeline files, generated manifests, and verified
  experiment outputs.

For Gaeun/Sophia discussions:

- Label model-output-derived points as hypotheses or questions unless verified.
