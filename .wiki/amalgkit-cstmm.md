## Overview

`amalgkit cstmm` applies cross-species TMM normalization using single-copy genes. It is optional and Python-only.

Use it after `merge` when you want ortholog-aware normalization before `wsfilter`, `csfilter`, or `finalize`.

```text
merge -> cstmm -> wsfilter -> csfilter -> finalize
```

## Inputs

`cstmm` expects:

- selected metadata
- per-species count tables from `merge`
- one ortholog source

Provide one of:

- `--dir_busco`
- `--orthogroup_table`

`--dir_count inferred` reads:

```text
out_dir/merge
```

## Examples

Using BUSCO tables:

```bash
amalgkit cstmm --out_dir ./ --dir_busco ./busco
```

Using an orthogroup table:

```bash
amalgkit cstmm \
    --out_dir ./ \
    --orthogroup_table ./Orthogroups.tsv
```

For single-species TMM:

```bash
amalgkit cstmm --out_dir ./ --orthogroup_table ""
```

## Main Outputs

- `cstmm/cstmm_multispecies_busco_table.tsv`
- `cstmm/cstmm_orthogroup_genecount.tsv`
- `cstmm/cstmm_exclusion.pdf`
- `cstmm/cstmm_normalization_factor_scatter.pdf`
- `cstmm/cstmm_normalization_factor_histogram.sample_group.pdf`
- `cstmm/cstmm_normalization_factor_histogram.scientific_name.pdf`
- `cstmm/cstmm_mean_expression_boxplot.pdf`
- `cstmm/<Species>/<Species>_cstmm_counts.tsv`
- `cstmm/<Species>/<Species>_eff_length.tsv`

## Useful Options

| Option | Use |
| --- | --- |
| `--metadata` | metadata table; inferred from `out_dir/metadata/metadata.tsv` |
| `--dir_count` | merge output directory; inferred from `out_dir/merge` |
| `--dir_busco` | directory of per-species BUSCO tables |
| `--orthogroup_table` | OrthoFinder-style table, such as `Orthogroups.tsv` or `N0.tsv` |
| `--redo yes/no` | rerun even if output exists |
| `--tmm_backend python` | current TMM backend |

## Next Steps

```bash
amalgkit wsfilter --out_dir ./
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv
```
