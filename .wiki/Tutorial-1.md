# Tutorial 1: Yeast Pipeline

## Overview

This tutorial runs a compact AMALGKIT workflow with the bundled yeast assets for:

- *Saccharomyces cerevisiae*
- *Schizosaccharomyces pombe*

The tutorial follows the current Python-only CLI:

```text
dataset -> metadata -> select -> getfastq -> quant -> merge -> cstmm -> wsfilter -> csfilter -> finalize
```

The BUSCO and filtering steps are included so the tutorial demonstrates the cross-species path. For a minimal single-species workflow, `merge -> finalize` is enough after quantification.

## Prerequisites

- AMALGKIT installed
- `fasterq-dump` from `sra-tools >= 3`
- `fastp`, unless you pass `--fastp no`
- `kallisto` for short-read quantification

R is not required.

## 0. Create a Workspace

```bash
mkdir ./amalgkit_tutorial
cd ./amalgkit_tutorial
```

## 1. Extract the Yeast Assets

```bash
amalgkit dataset --name yeast --out_dir ./
```

This creates a small working layout:

```text
fasta/
busco/
downloads/
metadata/
metadata_specieswise/
private_fastq/
select_rules.tsv
```

The bundled FASTA files are intentionally small tutorial assets, not full production transcriptomes.

## 2. Review Selection Rules

```bash
less select_rules.tsv
```

`amalgkit select` reads `select_rules.tsv` directly. The former config-directory workflow is no longer used by the current CLI.

## 3. Retrieve Metadata

```bash
amalgkit metadata \
    --out_dir ./ \
    --entrez_email your@email.com \
    --search_string '("ERP109456"[Accession] AND "Saccharomyces cerevisiae"[Organism]) OR ("SRP565465"[Accession] AND "Schizosaccharomyces pombe"[Organism])'
```

Main output:

- `metadata/metadata.tsv`

## 4. Review `sample_group`

`sample_group` is the biological grouping column used by `select`, `wsfilter`, `csfilter`, and `finalize`. Yeast metadata often benefits from normalizing this field from genotype or condition columns.

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

## 5. Select Samples

```bash
amalgkit select --out_dir ./
```

`select` applies `select_rules.tsv` and updates metadata columns such as `exclusion`, `is_qualified`, and `is_sampled`.

## 6. Download and Prepare FASTQ Files

```bash
amalgkit getfastq \
    --out_dir ./ \
    --max_bp 70000000 \
    --threads 2
```

Main outputs are written under `getfastq/`, with run-level statistics summarized in `getfastq_stats.tsv`.

## 7. Quantify Expression

```bash
amalgkit quant \
    --out_dir ./ \
    --build_index yes \
    --clean_fastq no \
    --threads 2
```

With `--quant_backend auto`, AMALGKIT uses kallisto for short-read runs and oarfish for long-read runs. This tutorial uses kallisto.

## 8. Merge Per-Run Abundance Tables

```bash
amalgkit merge --out_dir ./
```

Main per-species outputs:

- `merge/<Species>/<Species>_est_counts.tsv`
- `merge/<Species>/<Species>_eff_length.tsv`
- `merge/<Species>/<Species>_tpm.tsv`

## 9. Run Cross-Species TMM Normalization

```bash
amalgkit cstmm \
    --out_dir ./ \
    --dir_busco ./busco
```

Main per-species output:

- `cstmm/<Species>/<Species>_cstmm_counts.tsv`

If you skip this step, downstream filters and `finalize` infer `merge/` as the input directory.

## 10. Run Within-Species Filtering

```bash
amalgkit wsfilter --out_dir ./
```

Main outputs:

- `wsfilter/metadata.tsv`
- `wsfilter/excluded.tsv`
- `wsfilter/wsfilter_exclusion.pdf`

## 11. Run Cross-Species Filtering

```bash
amalgkit csfilter \
    --out_dir ./ \
    --metadata ./wsfilter/metadata.tsv \
    --dir_busco ./busco
```

Main outputs:

- `csfilter/metadata.tsv`
- `csfilter/excluded.tsv`
- `csfilter/csfilter_exclusion.pdf`
- `csfilter/*.pdf`

## 12. Export Final Tables

For a minimal tutorial run:

```bash
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg no
```

Batch-correction choices are:

- `--batch_effect_alg no`
- `--batch_effect_alg sva`
- `--batch_effect_alg ruvseq`
- `--batch_effect_alg combatseq`, after installing optional `inmoose`
- `--batch_effect_alg latent_glm`

## Final Outputs

Top level:

- `finalize/metadata.tsv`
- `finalize/finalize_exclusion.pdf`

Per species:

- `<Species>_metadata.tsv`
- `<Species>_expression.tsv`
- `<Species>_expression_uncorrected.tsv`
- `<Species>_sample_group_mean.tsv`
- `<Species>_tau.tsv`
- `<Species>_batch_effect_summary.tsv`

## Optional Sanity Check

```bash
amalgkit sanity --out_dir ./ --all
```

Preview reruns before executing them:

```bash
amalgkit rerun --out_dir ./ --dry_run yes
```
