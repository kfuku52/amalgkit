# amalgkit dataset

## Overview

`amalgkit dataset` extracts bundled test datasets for quick pipeline checks and tutorials.

## Usage

```bash
amalgkit dataset --list
amalgkit dataset --name yeast --out_dir ./test_run
```

## Available dataset: `yeast`

The bundled yeast dataset contains:

- FASTA files for *Saccharomyces cerevisiae* and *Schizosaccharomyces pombe*
- precomputed BUSCO full tables
- yeast-optimized config files for `amalgkit select`

## Extracted structure

```text
out_dir/
├── fasta/
├── busco/
└── config/
```

## Example pipeline after extraction

```bash
amalgkit dataset --name yeast --out_dir ./test_run
amalgkit config --out_dir ./test_run --config base --overwrite yes
cp ./test_run/config/*.config ./test_run/config_base/
amalgkit metadata --out_dir ./test_run --entrez_email your@email.com --search_string '("ERP109456"[Accession] AND "Saccharomyces cerevisiae"[Organism]) OR ("SRP565465"[Accession] AND "Schizosaccharomyces pombe"[Organism])'
amalgkit select --out_dir ./test_run --config_dir ./test_run/config_base
amalgkit getfastq --out_dir ./test_run
amalgkit quant --out_dir ./test_run
amalgkit merge --out_dir ./test_run
amalgkit cstmm --out_dir ./test_run --dir_busco ./test_run/busco
amalgkit wsfilter --out_dir ./test_run
amalgkit csfilter --out_dir ./test_run --metadata ./test_run/wsfilter/metadata.tsv --dir_busco ./test_run/busco
amalgkit finalize --out_dir ./test_run --metadata ./test_run/csfilter/metadata.tsv --batch_effect_alg no
```
