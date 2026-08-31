## Overview

`amalgkit wsfilter` performs within-species outlier filtering. It writes filtered metadata, excluded-row summaries, and per-species QC plots.

`wsfilter` is Python-only.

## Inputs

`wsfilter` expects:

- metadata from the matching `merge`/`cstmm` output or a previous filter step
- abundance tables from `merge` or `cstmm`

If `--input_dir inferred` is used, AMALGKIT reads:

```text
out_dir/cstmm if it exists, otherwise out_dir/merge
```

With `--metadata inferred`, AMALGKIT uses the last successful filter recorded
in `filter_metadata_state.json`, then the selected input directory's metadata.
Legacy workspaces prefer `csfilter` over `wsfilter` with a warning, not by file
modification time. See [metadata and normalization](https://github.com/kfuku52/amalgkit/wiki/Metadata-and-normalization)
for the complete precedence rules and long-read settings.

For the first filter pass, an explicit metadata path is often clearest:

```bash
amalgkit wsfilter --out_dir ./ --input_dir ./cstmm --metadata ./cstmm/metadata.tsv
```

If CSTMM was skipped, use `--input_dir ./merge --metadata ./merge/metadata.tsv`.
Do not use the original pre-CSTMM metadata with CSTMM counts: it lacks the
original library sizes needed to preserve TMM during FPKM calculation.
For Oarfish runs, add `--norm log2p1-none` to every downstream command.

## Basic Use

```bash
amalgkit wsfilter --out_dir ./
```

Restrict to selected sample groups:

```bash
amalgkit wsfilter \
    --out_dir ./ \
    --sample_group leaf,root,flower
```

## Main Outputs

- `wsfilter/metadata.tsv`
- `wsfilter/excluded.tsv`
- `wsfilter/wsfilter_exclusion.pdf`
- `wsfilter/<Species>/<Species>_within_group_correlation_no.pdf`
- `wsfilter/<Species>/<Species>_tau_histogram_no.pdf`

## Useful Options

| Option | Default | Use |
| --- | --- | --- |
| `--mapping_rate` | `0.2` | mapping-rate cutoff |
| `--dist_method` | `pearson` | distance/correlation method |
| `--norm` | `log2p1-fpkm` | expression transformation before filtering |
| `--margin_threshold` | `0.0` | robust-margin threshold |
| `--robust_z_threshold` | `-2.5` | robust z-score threshold |
| `--one_outlier_per_iter` | `no` | remove at most one outlier per group/project per iteration |
| `--plot_intermediate` | `no` | write intermediate filtering plots |

## Chaining

Use only within-species filtering:

```bash
amalgkit wsfilter --out_dir ./
amalgkit finalize \
    --out_dir ./ \
    --metadata ./wsfilter/metadata.tsv \
    --batch_effect_alg no
```

Continue to cross-species filtering:

```bash
amalgkit wsfilter --out_dir ./
amalgkit csfilter \
    --out_dir ./ \
    --metadata ./wsfilter/metadata.tsv \
    --dir_busco ./busco
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg latent_glm
```
