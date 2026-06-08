## Overview

`amalgkit finalize` exports final per-species expression tables from metadata and merged abundance data.

It is also the current location for optional batch-effect correction. All current backends are Python implementations.

## Inputs

`finalize` expects:

- metadata, usually from `wsfilter/metadata.tsv` or `csfilter/metadata.tsv`
- abundance tables from `merge` or `cstmm`

If `--input_dir inferred` is used, AMALGKIT reads:

```text
out_dir/cstmm if it exists, otherwise out_dir/merge
```

If `--metadata inferred` is used, AMALGKIT can reuse the newest prior filter metadata when available:

```text
out_dir/wsfilter/metadata.tsv
out_dir/csfilter/metadata.tsv
out_dir/metadata/metadata.tsv
```

For reproducible runs, pass the metadata path explicitly:

```bash
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg no
```

## Batch-Correction Backends

`--batch_effect_alg` supports:

- `no`
- `sva`
- `ruvseq`
- `combatseq`
- `latent_glm`

Examples:

```bash
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg sva
```

```bash
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg latent_glm \
    --latent_family nb \
    --latent_k auto
```

## Main Outputs

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

## General Options

| Option | Default | Use |
| --- | --- | --- |
| `--norm` | `log2p1-fpkm` | expression transformation before optional batch correction |
| `--clip_negative` | `yes` | clip negative corrected values to zero |
| `--maintain_zero` | `yes` | preserve input zero values after correction |
| `--seed` | `auto` | random seed for stochastic steps |

## Backend-Specific Options

SVA:

- `--sva_nsv`
- `--sva_B`
- `--sva_B_auto_max`

RUVSeq:

- `--ruvseq_control_genes`
- `--ruvseq_k`
- `--ruvseq_k_max`
- `--ruvseq_control_top_n`
- `--ruvseq_min_controls`

latent_glm:

- `--latent_family poisson|nb`
- `--latent_k INT|auto`
- `--latent_k_max INT`
- `--latent_max_iter INT`
- `--latent_tol FLOAT`

Backend selectors:

- `--sva_backend python`
- `--combatseq_backend python`
- `--ruvseq_backend python`

## Next Steps

Check outputs after a full run:

```bash
amalgkit sanity --out_dir ./ --check finalize
```
