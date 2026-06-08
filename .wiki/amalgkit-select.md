## Overview

`amalgkit select` marks metadata rows for downstream analysis using a single `select_rules.tsv` file. It updates metadata columns such as `sample_group`, `exclusion`, `is_qualified`, and `is_sampled`.

The old config-directory workflow has been replaced. Current releases read one TSV rule file through `--select_rules_tsv` or, by default, `out_dir/select_rules.tsv`.

## Example command

```bash
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
amalgkit select --out_dir ./
```

With an explicit rule file:

```bash
amalgkit select --out_dir ./ --select_rules_tsv ./select_rules.tsv
```

## Rule file location

If `--select_rules_tsv inferred` is used, AMALGKIT reads:

```text
out_dir/select_rules.tsv
```

Create one with:

```bash
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
```

Available bundled rule sets include `base`, `test`, `plantae`, and `vertebrate`.

## Rule stages

`select_rules.tsv` uses stage rows to describe selection behavior:

- `parameter`: set values such as `min_nspots`, `max_sample`, `sample_group`, and `sampling_strategy`
- `aggregate`: append sparse metadata fields into a target field
- `normalize`: assign normalized `sample_group` values from matching metadata text
- `exclude`: set `exclusion` for unwanted library types or keywords
- `control`: keep control-like samples within a project/group scope
- `filter`: apply numeric or missing-value filters
- `dedup`: remove redundant records such as repeated BioSample entries
- `validate`: collect hints used by rule validation

## Batch selection

For large species lists, `select` can operate in native batch mode:

```bash
amalgkit metadata --out_dir ./work --species_tsv ./species.tsv
amalgkit select --out_dir ./work/batch --species_tsv ./species.tsv
```

Batch mode reads per-species metadata from `metadata_specieswise/` and writes:

- `select_summary.tsv`
- `select_queue.tsv`
- `external_manifest.tsv`
- manifest sidecars such as `external_manifest.all_tissues_ge30.tsv`

## Main outputs

- updated metadata under `metadata/`
- selection summaries and manifests in batch mode
- rows with `exclusion != no` or `is_sampled != yes` are skipped by downstream commands
