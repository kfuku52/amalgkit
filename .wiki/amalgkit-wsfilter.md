## Overview

`amalgkit wsfilter` performs within-species outlier filtering and writes filtered metadata plus per-species QC plots. It is Python-only.

## Input

- `--metadata`: metadata table from `amalgkit metadata` or `amalgkit integrate`
- `--input_dir`: `merge` or `cstmm` output directory

If `--metadata inferred` is used, `wsfilter` automatically picks the newest of:

- `out_dir/wsfilter/metadata.tsv`
- `out_dir/csfilter/metadata.tsv`
- otherwise `out_dir/metadata/metadata.tsv`

## Example

```bash
amalgkit wsfilter --out_dir ./
```

## Main outputs

- `wsfilter/metadata.tsv`
- `wsfilter/excluded.tsv`
- `wsfilter/wsfilter_exclusion.pdf`
- `wsfilter/<Species>/<Species>_within_group_correlation_no.pdf`
- `wsfilter/<Species>/<Species>_tau_histogram_no.pdf`

## Useful options

- `--mapping_rate`
- `--dist_method`
- `--margin_threshold`
- `--robust_z_threshold`
- `--plot_intermediate yes|no`

## Typical chaining

```bash
amalgkit wsfilter --out_dir ./
amalgkit finalize --out_dir ./ --metadata ./wsfilter/metadata.tsv --batch_effect_alg no
```

With cross-species filtering:

```bash
amalgkit wsfilter --out_dir ./
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg latent_glm
```
