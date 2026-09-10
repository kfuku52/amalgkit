## Overview

`amalgkit wsfilter` performs within-species outlier filtering. It writes filtered metadata, excluded-row summaries, and per-species QC plots.

`wsfilter` is Python-only.

## Inputs

Tau panels use linear arithmetic means and accept `--tau_unit` and `--tau_balance_projects`. These settings affect tau summaries, not the correlation-based exclusion decision. See [Tau and replicate aggregation](Tau-and-replicate-aggregation.md).

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
For Oarfish runs, add `--norm log2p1-cpm` to every downstream command.

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

## Reference and support sensitivity analysis

Filtering remains automatic. The existing defaults are unchanged; alternative
policies are opt-in because retention and detection can trade off.

| Option | Default | Use |
| --- | --- | --- |
| `--small_group_policy` | `margin_fallback` | Use `retain` to keep groups with fewer than three finite margins; this is insufficient evidence, not a quality pass. |
| `--min_common_genes` | `0` | Require this many finite gene pairs for every group comparison; `0` disables the additional cutoff. Values other than `0` must be at least `2`. |
| `--reference_exclusion` | `run` | `bioproject` excludes the target project from both same-group and other-group references. Complete project labels are required (`not_provided` is missing). |
| `--max_filter_iterations` | `None` | Set `1` for a single pass, or another positive integer for a bounded number of rounds. Omission repeats until stable. |

Each comparison still uses its own observed gene pairs. If an enabled support
cutoff fails for the within-group reference or any other-group reference, the
margin is unavailable and that sample is retained by correlation filtering.
Other filters, including mapping-rate filtering, still apply. Excluding a
project can leave no reference and does not silently fall back to run exclusion.

`metadata.tsv` and `excluded.tsv` include `ws_within_common_genes`,
`ws_min_nongroup_common_genes`, `ws_min_common_genes`, and
`ws_reference_exclusion`. Counts describe finite pairs with the aggregated
reference, not the number of independent observations. Metadata also records
`ws_small_group_policy`, `ws_filter_iterations`, `ws_max_filter_iterations`,
and `ws_filter_stop_reason`. Removed samples retain their removal-round scores and scoring settings, including
`ws_margin_threshold` and `ws_robust_z_threshold`. Previously excluded samples
are not relabelled by the mapping-rate filter. On rescoring an active run, an
unavailable score clears any earlier finite score in the final metadata.
When `--one_outlier_per_iter yes` is used, candidates are considered by increasing
margin (run ID breaks ties), enforcing at most one per group and per known
project in the same round.

For a comparison, start each run from the same explicit **pre-filter** metadata
in separate output directories. Reusing inferred filtered metadata prevents
previously removed samples from being reconsidered. Recompute downstream
finalization and summaries when the retained sample set changes.

See [filter validation](https://github.com/kfuku52/amalgkit/wiki/Filter-validation)
for the evaluation design, limitations, and adoption decision.
