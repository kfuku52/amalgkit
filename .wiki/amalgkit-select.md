## Overview

`amalgkit select` applies `select_rules.tsv` to metadata and marks rows for downstream analysis.

It updates fields such as:

- `sample_group`
- `exclusion`
- `is_qualified`
- `is_sampled`

Current releases use one TSV rule file. The former config-directory workflow has been removed.

## Basic Use

Create a starter rule file, edit it if needed, then run `select`:

```bash
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
amalgkit select --out_dir ./
```

With an explicit rule path:

```bash
amalgkit select \
    --out_dir ./ \
    --metadata ./metadata/metadata.tsv \
    --select_rules_tsv ./select_rules.tsv
```

If `--select_rules_tsv inferred` is used, AMALGKIT reads:

```text
out_dir/select_rules.tsv
```

## Rule Stages

`select_rules.tsv` uses stage rows to describe selection behavior:

| Stage | Role |
| --- | --- |
| `parameter` | set values such as `min_nspots`, `max_sample`, `sample_group`, and `sampling_strategy` |
| `aggregate` | append sparse metadata fields into a target field |
| `normalize` | assign normalized values such as `sample_group` from matching metadata text |
| `exclude` | set `exclusion` for unwanted library types, keywords, or other patterns |
| `control` | keep control-like samples within a project/group scope |
| `filter` | apply numeric or missing-value filters |
| `dedup` | remove redundant records such as repeated BioSample entries |
| `validate` | collect rule-validation hints |

Available bundled rule sets are `base`, `test`, `plantae`, and `vertebrate`.

```bash
amalgkit dataset --list
amalgkit dataset --rule_set plantae --out_dir ./ --overwrite yes
```

## Species-Wise Batch Mode

For large species lists, `select` can operate on per-species metadata written by `metadata --species_tsv`.

```bash
amalgkit metadata \
    --out_dir ./work \
    --species_tsv ./work/species.tsv \
    --entrez_email example@email.com
amalgkit select \
    --out_dir ./work/batch \
    --species_tsv ./work/species.tsv \
    --metadata_specieswise_dir ./work/metadata_specieswise
```

Useful batch options:

| Option | Default |
| --- | --- |
| `--metadata_specieswise_dir` | `dirname(out_dir)/metadata_specieswise` |
| `--summary_tsv` | `out_dir/select_summary.tsv` |
| `--queue_tsv` | `out_dir/select_queue.tsv` |
| `--manifest_tsv` | `out_dir/external_manifest.tsv` |
| `--batch_label` | `basename(out_dir)` |

Batch mode writes:

- `select_summary.tsv`
- `select_queue.tsv`
- `external_manifest.tsv`
- manifest sidecars for `all_tissues_ge30`, `all_tissues_ge3`, `all_tissues_ge1`, and `any_tissues_ge1`

## Main Outputs

For regular mode:

- updated metadata under `metadata/`

For batch mode:

- summary, queue, and manifest TSVs under `--out_dir`
- species-specific selected metadata in batch workspaces

Rows with `exclusion != no` or `is_sampled != yes` are skipped by downstream commands.

## Next Steps

```bash
amalgkit getfastq --out_dir ./
amalgkit quant --out_dir ./ --build_index yes
```
