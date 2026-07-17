## Start Here

AMALGKIT builds RNA-seq expression datasets from public SRA records and private FASTQ files. The current CLI is Python-only: the main pipeline no longer requires R, `Rscript`, or R packages.

Good first pages:

- [Installation and dependencies](https://github.com/kfuku52/amalgkit/wiki/Installation-and-dependencies)
- [Tutorial 1](https://github.com/kfuku52/amalgkit/wiki/Tutorial-1)
- [Parallel processing](https://github.com/kfuku52/amalgkit/wiki/Parallel-processing)
- [GitHub releases](https://github.com/kfuku52/amalgkit/releases)

## Workflow Map

The usual public-SRA workflow is:

```text
dataset -> metadata -> select -> getfastq -> quant -> merge -> finalize
```

Cross-species projects often add ortholog-aware normalization and filtering after `merge`:

```text
merge + BUSCO/orthogroups -> cstmm -> wsfilter -> csfilter -> finalize
```

Private FASTQ projects start by integrating local files into metadata:

```text
integrate -> getfastq -> quant -> merge -> finalize
```

Validation and recovery are separate maintenance steps:

```text
sanity -> rerun
```

## Typical Public-SRA Run

```bash
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
amalgkit metadata \
    --out_dir ./ \
    --entrez_email example@email.com \
    --search_string '("platform illumina"[Properties]) AND ("type rnaseq"[Filter])'
amalgkit select --out_dir ./
amalgkit getfastq --out_dir ./
amalgkit quant --out_dir ./ --build_index yes
amalgkit merge --out_dir ./
amalgkit finalize --out_dir ./ --batch_effect_alg no
```

Add BUSCO-based cross-species steps when you have compatible transcriptome FASTA files or BUSCO tables:

```bash
amalgkit busco --out_dir ./ --lineage eukaryota_odb12
amalgkit cstmm --out_dir ./ --dir_busco ./busco
amalgkit wsfilter --out_dir ./
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg latent_glm
```

## Active Commands

| Command | Role |
| --- | --- |
| [`amalgkit dataset`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-dataset) | extract bundled datasets, initialize workspaces, and export `select_rules.tsv` |
| [`amalgkit metadata`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata) | retrieve SRA metadata by Entrez query or species-wise batch files |
| [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate) | add local FASTQ files to AMALGKIT metadata |
| [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select) | apply `select_rules.tsv` and mark rows for downstream processing |
| [`amalgkit getfastq`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-getfastq) | download public SRA data or stage private FASTQ files |
| [`amalgkit quant`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant) | quantify transcript abundance with kallisto or oarfish |
| [`amalgkit merge`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-merge) | merge per-run quantification into per-species tables |
| [`amalgkit busco`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-busco) | create BUSCO tables for ortholog-aware downstream steps |
| [`amalgkit cstmm`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-cstmm) | run cross-species TMM normalization using single-copy genes |
| [`amalgkit wsfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-wsfilter) | perform within-species outlier filtering |
| [`amalgkit csfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csfilter) | perform cross-species outlier filtering |
| [`amalgkit finalize`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-finalize) | export final expression tables and optional batch-corrected outputs |
| [`amalgkit sanity`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-sanity) | check expected outputs in a workspace |
| [`amalgkit rerun`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-rerun) | rerun failed targets recorded by `sanity` |

## Legacy Command Migration

These commands existed in older AMALGKIT releases but are not part of the current CLI:

| Old command | Current workflow |
| --- | --- |
| [`amalgkit config`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-config) | use `amalgkit dataset --rule_set ...` to write `select_rules.tsv`, then run `amalgkit select` |
| [`amalgkit curate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-curate) | run `amalgkit wsfilter`, optionally `amalgkit csfilter`, then `amalgkit finalize` |
| [`amalgkit csca`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csca) | use `amalgkit csfilter` for cross-species filtering and `amalgkit finalize` for final tables |
