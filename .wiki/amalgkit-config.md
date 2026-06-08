## Legacy Command

`amalgkit config` has been removed from the current CLI.

Selection configuration now lives in one `select_rules.tsv` file.

## Current Workflow

Export a bundled rule set:

```bash
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
```

Edit `select_rules.tsv`, then run:

```bash
amalgkit select --out_dir ./
```

See [amalgkit select](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select) for the active workflow.

## Migration Guide

| Former file | Current `select_rules.tsv` stage |
| --- | --- |
| `group_attribute.config` | `aggregate` rows |
| `control_term.config` | `control` rows |
| `exclude_keyword.config` | `exclude` rows |

The current rule file is a TSV with explicit rule IDs, stages, priorities, columns, patterns, actions, outcomes, and optional parameter rows.

## Bundled Rule Sets

```bash
amalgkit dataset --list
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
amalgkit dataset --rule_set plantae --out_dir ./plant_run --overwrite yes
```

Available rule sets are:

- `base`
- `test`
- `plantae`
- `vertebrate`
