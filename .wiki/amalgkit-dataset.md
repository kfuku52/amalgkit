# amalgkit dataset

## Overview

`amalgkit dataset` extracts bundled test datasets for quick pipeline checks and tutorials.
It also initializes empty workspaces and exports bundled `select_rules.tsv` rule sets.

## Usage

```bash
amalgkit dataset --list
amalgkit dataset --name init --out_dir ./work
amalgkit dataset --name yeast --out_dir ./test_run
amalgkit dataset --rule_set base --out_dir ./work --overwrite yes
```

## Workspace scaffold: `init`

`amalgkit dataset --name init --out_dir ./work` creates a starter workspace:

- `species.tsv`
- `organ_terms.tsv`
- `select_rules.tsv`
- `WORKSPACE_README.md`
- `fasta/`
- `private_fastq/`
- `downloads/`
- `metadata/`
- `metadata_specieswise/`

## Available dataset: `yeast`

The bundled yeast dataset contains:

- FASTA files for *Saccharomyces cerevisiae* and *Schizosaccharomyces pombe*
- precomputed BUSCO full tables
- a yeast-oriented `select_rules.tsv` for `amalgkit select`

## Extracted structure

```text
out_dir/
|-- fasta/
|-- busco/
|-- downloads/
|-- metadata/
|-- metadata_specieswise/
|-- private_fastq/
`-- select_rules.tsv
```

## Rule sets

Bundled rule sets can be exported without extracting a dataset:

```bash
amalgkit dataset --rule_set base --out_dir ./work --overwrite yes
amalgkit dataset --rule_set plantae --out_dir ./plant_run --overwrite yes
```

Available rule sets are listed by `amalgkit dataset --list`.

## Example pipeline after extraction

```bash
amalgkit dataset --name yeast --out_dir ./test_run
amalgkit metadata --out_dir ./test_run --entrez_email your@email.com --search_string '("ERP109456"[Accession] AND "Saccharomyces cerevisiae"[Organism]) OR ("SRP565465"[Accession] AND "Schizosaccharomyces pombe"[Organism])'
amalgkit select --out_dir ./test_run
amalgkit getfastq --out_dir ./test_run
amalgkit quant --out_dir ./test_run
amalgkit merge --out_dir ./test_run
amalgkit cstmm --out_dir ./test_run --dir_busco ./test_run/busco
amalgkit wsfilter --out_dir ./test_run
amalgkit csfilter --out_dir ./test_run --metadata ./test_run/wsfilter/metadata.tsv --dir_busco ./test_run/busco
amalgkit finalize --out_dir ./test_run --metadata ./test_run/csfilter/metadata.tsv --batch_effect_alg no
```
