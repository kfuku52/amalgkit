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
- `cstmm/metadata.tsv`
- `cstmm/<Species>/<Species>_cstmm_counts.tsv`
- `cstmm/<Species>/<Species>_eff_length.tsv`

## Normalization Metadata

`cstmm/metadata.tsv` records three related columns for each retained sample:

| Column | Meaning |
| --- | --- |
| `tmm_library_size` | Sum of the uncorrected estimated counts for all targets in the corresponding `merge` output. This is the library size supplied when estimating the TMM factor. |
| `tmm_normalization_factor` | Cross-species TMM factor estimated from the selected single-copy orthologs. |
| `tmm_effective_library_size` | Product of `tmm_library_size` and `tmm_normalization_factor`; this is the effective library size for a raw-count model. |

For sample `i`, the relationship is:

```text
tmm_effective_library_size[i]
    = tmm_library_size[i] * tmm_normalization_factor[i]
```

The per-species `cstmm_counts.tsv` files already contain counts divided by the
TMM normalization factor. They are used together with the original
`tmm_library_size` by the normal amalgkit `wsfilter -> csfilter -> finalize`
workflow. Do not apply `tmm_normalization_factor` to those corrected counts a
second time.

## Using CSTMM Normalization in edgeR

Although CSTMM output is used directly by downstream amalgkit commands, the
normalization can also be passed to external programs. edgeR fits
negative-binomial models to uncorrected counts using library-size offsets, so
CSTMM normalization is supplied through the effective library size rather than
by altering the downstream count matrix.

To connect CSTMM to edgeR, start from uncorrected estimated counts in the
`merge` output (or HOG counts obtained by summing those uncorrected counts), and
supply `tmm_effective_library_size` as the library size. Align metadata and
count columns explicitly; the example below assumes the count columns are run
IDs.

```r
library(edgeR)

cstmm_metadata <- read.delim(
    "cstmm/metadata.tsv",
    check.names = FALSE
)

sample_index <- match(colnames(hog_counts), cstmm_metadata$run)
if (anyNA(sample_index)) {
    stop("Some count columns do not match cstmm metadata run IDs")
}

effective_library_size <-
    cstmm_metadata$tmm_effective_library_size[sample_index]
if (any(!is.finite(effective_library_size) | effective_library_size <= 0)) {
    stop("CSTMM effective library sizes must be finite and positive")
}

dge <- DGEList(
    counts = hog_counts,
    group = group,
    lib.size = effective_library_size,
    norm.factors = rep(1, length(effective_library_size))
)

keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = TRUE]

design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
result <- glmQLFTest(fit, coef = 2)
```

Do not call `calcNormFactors()` in this workflow: that would estimate a new
normalization from the downstream count matrix instead of using CSTMM. Also do
not use `keep.lib.sizes = FALSE`, because it would replace the effective
library sizes with post-filter column sums. If count columns use custom sample
names rather than run IDs, construct and validate an explicit mapping before
creating the `DGEList`.

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
