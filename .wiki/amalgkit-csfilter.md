## Overview

`amalgkit csfilter` performs cross-species outlier filtering using ortholog information. It writes filtered metadata, excluded-row summaries, and cross-species QC plots.

`csfilter` is Python-only.

## Inputs

`csfilter` expects:

- filtered or unfiltered metadata
- abundance tables from `merge` or `cstmm`
- one ortholog source

Provide one of:

- `--dir_busco`
- `--orthogroup_table`

If `--input_dir inferred` is used, AMALGKIT reads:

```text
out_dir/cstmm if it exists, otherwise out_dir/merge
```

## Basic Use

After `wsfilter`:

```bash
amalgkit csfilter \
    --out_dir ./ \
    --metadata ./wsfilter/metadata.tsv \
    --dir_busco ./busco
```

Without `wsfilter`:

```bash
amalgkit csfilter \
    --out_dir ./ \
    --metadata ./metadata/metadata.tsv \
    --dir_busco ./busco
```

Using an orthogroup table:

```bash
amalgkit csfilter \
    --out_dir ./ \
    --metadata ./wsfilter/metadata.tsv \
    --orthogroup_table ./Orthogroups.tsv
```

## Main Outputs

- `csfilter/metadata.tsv`
- `csfilter/excluded.tsv`
- `csfilter/csfilter_exclusion.pdf`
- `csfilter/csfilter_overview.pdf`
- `csfilter/csfilter_unaveraged_pca_PC12.pdf`
- `csfilter/csfilter_heatmap.pdf`
- `csfilter/csfilter_within_group_cor.pdf`
- `csfilter/csfilter_outlier_scatter.pdf`

## Useful Options

| Option | Default | Use |
| --- | --- | --- |
| `--missing_strategy` | `em_pca` | missing-value handling before dimensionality reduction |
| `--norm` | `log2p1-fpkm` | expression transformation used for temporary tables |
| `--margin_threshold` | `0.0` | robust-margin threshold |
| `--robust_z_threshold` | `-2.5` | robust z-score threshold |
| `--sample_group` | all groups | comma-separated sample groups to include |
| `--sample_group_color` | automatic | colors for selected sample groups |

## Chaining

```bash
amalgkit wsfilter --out_dir ./
amalgkit csfilter \
    --out_dir ./ \
    --metadata ./wsfilter/metadata.tsv \
    --dir_busco ./busco
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg sva
```
