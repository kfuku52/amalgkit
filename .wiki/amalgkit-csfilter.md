## Overview

`amalgkit csfilter` performs cross-species outlier filtering using ortholog information. It writes filtered metadata, excluded-row summaries, and cross-species QC plots.

`csfilter` is Python-only.

## Inputs

Per-species tau summaries accept `--tau_unit` and `--tau_balance_projects` and use linear arithmetic means. Cross-species averaged correlation plots retain their transformed-scale run means. See [Tau and replicate aggregation](Tau-and-replicate-aggregation.md).

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

Without `wsfilter`, starting directly from CSTMM:

```bash
amalgkit csfilter \
    --out_dir ./ \
    --input_dir ./cstmm \
    --metadata ./cstmm/metadata.tsv \
    --dir_busco ./busco
```

If CSTMM was skipped, use `--input_dir ./merge --metadata ./merge/metadata.tsv`.
Metadata inference follows the last successful filter state, then the selected
input directory; it does not use file modification times. For Oarfish runs,
add `--norm log2p1-cpm`. See [metadata and normalization](https://github.com/kfuku52/amalgkit/wiki/Metadata-and-normalization).

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
- `csfilter/csfilter_run_pca_pc12_pre_correction.pdf`
- `csfilter/csfilter_run_pca_pc12_post_correction.pdf`
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

## Reference and support sensitivity analysis

Filtering remains automatic. The existing defaults are unchanged; these
alternatives require validation for the intended dataset.

| Option | Default | Use |
| --- | --- | --- |
| `--small_group_policy` | `margin_fallback` | `retain` keeps groups with fewer than three finite margins. |
| `--min_common_genes` | `0` | Minimum finite ortholog pairs for every group correlation; `0` disables the extra cutoff, other values must be at least `2`. |
| `--reference_exclusion` | `run` | `species` excludes the target species from both same-group and other-group references. |
| `--robust_z_scope` | `sample_group` | `species_group` computes robust z within each species and sample group. |

Run exclusion already removes the evaluated sample from its within-group
reference. Species exclusion additionally removes all other runs of that
species; it does not make the remaining species phylogenetically independent.
Species-group z can preserve a consistent species-specific shift, but can also
miss species-wide technical problems. With species-group z, the small-group
policy applies to each species/group, not to the pooled group.

Correlations use observed values before PCA imputation. Each comparison retains
its own gene set. If an enabled support cutoff fails for the within-group
reference or any competing group, the margin is unavailable and correlation
filtering retains the run. A missing species-excluded reference is not replaced
with a reference containing the target species. Missing margins are not proof
of good quality. `--missing_strategy` controls visualization, not these scores.

`metadata.tsv` and `excluded.tsv` include `cs_within_common_genes`,
`cs_min_nongroup_common_genes`, `cs_min_common_genes`, and
`cs_reference_exclusion`; metadata also records `cs_robust_z_scope` and
`cs_small_group_policy`, `cs_margin_threshold`, and `cs_robust_z_threshold`.
Previously excluded rows retain their scoring evidence and settings; active
rows receive newly computed values, including missing scores when unscoreable.
Counts are finite pairs with the aggregated reference,
not independent species or effective sample sizes. csfilter makes one scoring
pass. A robust z threshold is not a calibrated p-value.

Compare settings using the same explicit pre-filter metadata and distinct output
directories, then recompute downstream finalization and summaries if exclusions
change. See [filter validation](https://github.com/kfuku52/amalgkit/wiki/Filter-validation).
