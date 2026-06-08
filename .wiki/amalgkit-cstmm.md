## Overview

`amalgkit cstmm` applies cross-species TMM normalization using single-copy genes. It is Python-only in current releases.

This step is optional. A common workflow is:

```bash
amalgkit merge
amalgkit cstmm        # optional
amalgkit wsfilter
amalgkit csfilter     # optional
amalgkit finalize
```

## Input choices

Provide one of:

- `--dir_busco`
- `--orthogroup_table`

For single-species TMM, pass:

```bash
amalgkit cstmm --out_dir ./ --orthogroup_table ""
```

## Example

```bash
amalgkit cstmm --out_dir ./ --dir_busco ./busco
```

## Main outputs

- `cstmm/cstmm_multispecies_busco_table.tsv`
- `cstmm/cstmm_orthogroup_genecount.tsv`
- `cstmm/cstmm_exclusion.pdf`
- `cstmm/cstmm_normalization_factor_scatter.pdf`
- `cstmm/cstmm_normalization_factor_histogram.sample_group.pdf`
- `cstmm/cstmm_normalization_factor_histogram.scientific_name.pdf`
- `cstmm/cstmm_mean_expression_boxplot.pdf`
- `cstmm/<Species>/<Species>_cstmm_counts.tsv`
- `cstmm/<Species>/<Species>_eff_length.tsv`

## Notes

- `--tmm_backend` defaults to `python`
- `cstmm` itself does not require R
