## Overview

`amalgkit finalize` exports final per-species tables from filtered metadata. It is Python-only and is the current place for optional batch correction.

## Input

- `--metadata`: usually `wsfilter/metadata.tsv` or `csfilter/metadata.tsv`
- `--input_dir`: `merge` or `cstmm` output directory

If `--metadata inferred` is used, `finalize` automatically picks the newest of:

- `out_dir/wsfilter/metadata.tsv`
- `out_dir/csfilter/metadata.tsv`
- otherwise `out_dir/metadata/metadata.tsv`

## Batch-correction backends

`--batch_effect_alg` supports:

- `no`
- `sva`
- `ruvseq`
- `combatseq`
- `latent_glm`

All backends are implemented in Python in current releases.

## Example

No batch correction:

```bash
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg no
```

Nonnegative latent-factor correction:

```bash
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg latent_glm \
    --latent_family nb \
    --latent_k auto
```

## Main outputs

Top level:

- `finalize/metadata.tsv`
- `finalize/finalize_exclusion.pdf`

Per species:

- `<Species>_metadata.tsv`
- `<Species>_expression.tsv`
- `<Species>_expression_uncorrected.tsv`
- `<Species>_sample_group_mean.tsv`
- `<Species>_sample_group_mean_uncorrected.tsv`
- `<Species>_tau.tsv`
- `<Species>_correlation_statistics.tsv`
- `<Species>_batch_effect_summary.tsv`
- `<Species>_curation_round_summary.tsv`
- `<Species>_curation_final_summary.tsv`
- `<Species>_batch_compare_<alg>.pdf`
- `<Species>_tau_hist_<alg>.pdf`

## Useful options

General:

- `--norm`
- `--clip_negative`
- `--maintain_zero`
- `--seed`

SVA:

- `--sva_nsv`
- `--sva_B`
- `--sva_B_auto_max`

RUVSeq:

- `--ruvseq_control_genes`
- `--ruvseq_k`
- `--ruvseq_k_max`

latent_glm:

- `--latent_family poisson|nb`
- `--latent_k INT|auto`
- `--latent_k_max INT`
- `--latent_max_iter INT`
- `--latent_tol FLOAT`
