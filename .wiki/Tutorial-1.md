# Tutorial 1: Yeast Pipeline

## Overview

This tutorial runs the current Python-only AMALGKIT workflow with the bundled yeast assets:

- *Saccharomyces cerevisiae*
- *Schizosaccharomyces pombe*

The flow is:

```text
dataset -> metadata -> select -> getfastq -> quant -> merge -> cstmm -> wsfilter -> csfilter -> finalize
```

## Prerequisites

- AMALGKIT installed
- `fasterq-dump` from `sra-tools >= 3`
- `fastp` unless you pass `--fastp no`
- `kallisto` for the short-read quantification step

You do not need R.

## 0. Create a working directory

```bash
mkdir ./amalgkit_tutorial
cd ./amalgkit_tutorial
```

## 1. Extract the bundled dataset

```bash
amalgkit dataset --name yeast --out_dir ./
```

This creates:

```text
fasta/
busco/
downloads/
metadata/
metadata_specieswise/
private_fastq/
select_rules.tsv
```

The yeast dataset includes small FASTA files, precomputed BUSCO tables, and a yeast-oriented `select_rules.tsv`.

## 2. Review selection rules

```bash
less select_rules.tsv
```

`amalgkit select` now reads `select_rules.tsv` directly. The old config-directory workflow is no longer part of the current CLI.

## 3. Retrieve metadata

```bash
amalgkit metadata \
    --out_dir ./ \
    --entrez_email your@email.com \
    --search_string '("ERP109456"[Accession] AND "Saccharomyces cerevisiae"[Organism]) OR ("SRP565465"[Accession] AND "Schizosaccharomyces pombe"[Organism])'
```

## 4. Set `sample_group`

Yeast metadata often needs a manual fix so that `sample_group` reflects genotype or condition.

```bash
python - <<'PY'
import pandas as pd
df = pd.read_csv('metadata/metadata.tsv', sep='\t')
if 'genotype' in df.columns:
    x = df['genotype'].fillna('').astype(str).str.lower().str.strip()
    x = x.str.replace(r'[\s.(),/]+', '_', regex=True)
    x = x.str.replace(r'_+', '_', regex=True).str.strip('_')
    df['sample_group'] = x
df.to_csv('metadata/metadata.tsv', sep='\t', index=False)
print(sorted(df['sample_group'].dropna().astype(str).unique().tolist()))
PY
```

## 5. Select samples

```bash
amalgkit select --out_dir ./
```

This writes the selected metadata back under `metadata/` and records selection/exclusion labels according to `select_rules.tsv`.

## 6. Download FASTQ files

```bash
amalgkit getfastq \
    --out_dir ./ \
    --max_bp 70000000 \
    --threads 2
```

## 7. Quantify expression

```bash
amalgkit quant \
    --out_dir ./ \
    --build_index yes \
    --clean_fastq no \
    --threads 2
```

With `--quant_backend auto`, AMALGKIT uses kallisto for short-read runs and oarfish for long-read runs. This tutorial uses kallisto.

## 8. Merge abundance tables

```bash
amalgkit merge --out_dir ./
```

Main result:

- `merge/<Species>/<Species>_est_counts.tsv`
- `merge/<Species>/<Species>_eff_length.tsv`
- `merge/<Species>/<Species>_tpm.tsv`

## 9. Run cross-species TMM normalization

```bash
amalgkit cstmm \
    --out_dir ./ \
    --dir_busco ./busco
```

Main result:

- `cstmm/<Species>/<Species>_cstmm_counts.tsv`

## 10. Run within-species filtering

```bash
amalgkit wsfilter --out_dir ./
```

Main result:

- `wsfilter/metadata.tsv`
- `wsfilter/excluded.tsv`

## 11. Run cross-species filtering

```bash
amalgkit csfilter \
    --out_dir ./ \
    --metadata ./wsfilter/metadata.tsv \
    --dir_busco ./busco
```

Main result:

- `csfilter/metadata.tsv`
- `csfilter/excluded.tsv`
- `csfilter/*.pdf`

## 12. Export final tables

For a minimal tutorial run:

```bash
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg no
```

If you want batch correction, current releases also support:

- `--batch_effect_alg sva`
- `--batch_effect_alg ruvseq`
- `--batch_effect_alg combatseq`
- `--batch_effect_alg latent_glm`

## Key finalize outputs

Per species:

- `<Species>_expression.tsv`
- `<Species>_expression_uncorrected.tsv`
- `<Species>_sample_group_mean.tsv`
- `<Species>_tau.tsv`
- `<Species>_batch_effect_summary.tsv`

Top level:

- `finalize/metadata.tsv`
- `finalize/finalize_exclusion.pdf`

## Optional sanity check

```bash
amalgkit sanity --out_dir ./ --all
```

If the sanity report records failed targets, inspect the dry-run plan before executing reruns:

```bash
amalgkit rerun --out_dir ./ --dry_run yes
```
