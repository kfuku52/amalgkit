## Overview

`amalgkit dataset` extracts bundled datasets, initializes empty workspaces, and exports bundled `select_rules.tsv` rule sets.

Use it before `metadata` and `select` when you want a ready-made workspace or a starting rule file.

## Quick Commands

```bash
amalgkit dataset --list
amalgkit dataset --name init --out_dir ./work
amalgkit dataset --name yeast --out_dir ./demo
amalgkit dataset --rule_set base --out_dir ./work --overwrite yes
```

## Workspace Scaffold

`--name init` creates an empty workspace:

```bash
amalgkit dataset --name init --out_dir ./work
```

Main files and directories:

- `species.tsv`
- `organ_terms.tsv`
- `select_rules.tsv`
- `WORKSPACE_README.md`
- `fasta/`
- `private_fastq/`
- `downloads/`
- `metadata/`
- `metadata_specieswise/`

Use this mode for species-wise metadata collection or private FASTQ projects.

## Bundled Dataset

`--name yeast` extracts a compact tutorial dataset:

```bash
amalgkit dataset --name yeast --out_dir ./test_run
```

It includes:

- small FASTA files for *Saccharomyces cerevisiae* and *Schizosaccharomyces pombe*
- precomputed BUSCO full tables
- a yeast-oriented `select_rules.tsv`
- standard workspace directories

## Rule Sets

`--rule_set` writes a bundled `select_rules.tsv` without extracting a dataset:

```bash
amalgkit dataset --rule_set plantae --out_dir ./plant_run --overwrite yes
```

Available rule sets are:

- `base`
- `test`
- `plantae`
- `vertebrate`

Run `amalgkit dataset --list` to show available datasets and rule sets in the installed version.

## Output Behavior

By default, `dataset` does not overwrite existing files. Add `--overwrite yes` when intentionally replacing `select_rules.tsv` or an existing scaffold file.

## Next Steps

After `--name init`:

```bash
amalgkit metadata \
    --out_dir ./work \
    --species_tsv ./work/species.tsv \
    --entrez_email example@email.com
amalgkit select --out_dir ./work/batch --species_tsv ./work/species.tsv
```

After `--rule_set base`:

```bash
amalgkit metadata \
    --out_dir ./work \
    --entrez_email example@email.com \
    --search_string '("type rnaseq"[Filter])'
amalgkit select --out_dir ./work
```

After `--name yeast`, continue with [Tutorial 1](https://github.com/kfuku52/amalgkit/wiki/Tutorial-1).
