## Legacy command

`amalgkit config` has been removed from the current CLI.

Use `amalgkit dataset --rule_set ...` to export a bundled `select_rules.tsv`, then run `amalgkit select`.

```bash
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
amalgkit select --out_dir ./
```

## Migration guide

| Former file | Current `select_rules.tsv` stage |
| --- | --- |
| `group_attribute.config` | `aggregate` rows |
| `control_term.config` | `control` rows |
| `exclude_keyword.config` | `exclude` rows |

The current rule file is a TSV with explicit rule IDs, stages, priorities, columns, patterns, actions, outcomes, and optional parameter rows. See [amalgkit select](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select) for the active workflow.

## Bundled rule sets

```bash
amalgkit dataset --list
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
amalgkit dataset --rule_set plantae --out_dir ./plant_run --overwrite yes
```

Available rule sets are currently `base`, `test`, `plantae`, and `vertebrate`.
