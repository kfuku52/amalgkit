## Overview

`amalgkit csfilter` performs cross-species outlier filtering using ortholog information and writes filtered metadata plus cross-species QC plots. It is Python-only.

## Required input

- filtered or unfiltered `metadata.tsv`
- `merge` or `cstmm` output directory via `--input_dir`
- one of:
  - `--dir_busco`
  - `--orthogroup_table`

## Example

```bash
amalgkit csfilter --out_dir ./ --dir_busco ./busco
```

## Main outputs

- `csfilter/metadata.tsv`
- `csfilter/excluded.tsv`
- `csfilter/csfilter_exclusion.pdf`
- `csfilter/csfilter_overview.pdf`
- `csfilter/csfilter_unaveraged_pca_PC12.pdf`
- `csfilter/csfilter_heatmap.pdf`
- `csfilter/csfilter_within_group_cor.pdf`
- `csfilter/csfilter_outlier_scatter.pdf`

## Useful options

- `--missing_strategy em_pca|nipals|row_mean`
- `--margin_threshold`
- `--robust_z_threshold`

## Typical chaining

```bash
amalgkit wsfilter --out_dir ./
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg sva
```
