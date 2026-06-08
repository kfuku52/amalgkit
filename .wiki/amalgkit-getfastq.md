## Overview

`amalgkit getfastq` downloads public SRA objects, extracts FASTQ files with `fasterq-dump`, optionally runs `fastp`, and can run MMseqs2-based rRNA or contaminant filters.

Inputs can come from:

- metadata produced by `amalgkit metadata`, `amalgkit integrate`, or `amalgkit select`
- one explicit BioProject, BioSample, or run accession through `--id`
- a file of accessions through `--id_list`

## Examples

Use inferred `out_dir/metadata/metadata.tsv`:

```bash
amalgkit getfastq --out_dir ./
```

Use an explicit metadata table:

```bash
amalgkit getfastq --out_dir ./ --metadata /PATH/TO/metadata.tsv
```

Generate FASTQ files for one accession:

```bash
amalgkit getfastq --out_dir ./ --id DRR461654
```

## Two-step FASTQ extraction

When `--max_bp` is set, AMALGKIT can perform a compensatory extraction workflow for SRA-derived runs:

1. First-round extraction targets the requested size.
2. A second round can compensate for reads lost during `fasterq-dump`, `fastp`, rRNA filtering, or contaminant filtering.

This feature is for public SRA-derived data, not private FASTQ files.

## Interpreting `getfastq_stats.tsv`

`getfastq_stats.tsv` contains both count and base metrics. For paired-end libraries, count columns can use different units depending on the stage:

- `num_dumped`, `num_written`, `num_rrna_in/out`, and `num_contam_in/out`: spot counts
- `num_fastp_in/out`: read counts reported by `fastp`
- `bp_*`: total bases

For stage-by-stage removal fractions, `bp_*` columns are the safest values to compare across all stages.

## Download providers

By default, `getfastq` tries enabled public providers and moves on when one provider is unavailable or throttled:

- NCBI cloud objects: `--ncbi yes`
- AWS: `--aws yes`
- GCP: `--gcp yes`
- ENA: `--ena yes`
- DDBJ for DRA accessions: `--ddbj yes`

Each provider has an optional shared concurrency limit such as `--ncbi_download_max_concurrency`, `--aws_download_max_concurrency`, `--ena_download_max_concurrency`, and `--ddbj_download_max_concurrency`. These limits use `--download_lock_dir`, which defaults to `out_dir/downloads/locks`.

## Optional filters

```bash
amalgkit getfastq --out_dir ./ --rrna_filter yes
amalgkit getfastq --out_dir ./ --contam_filter yes --contam_filter_rank superkingdom
amalgkit getfastq --out_dir ./ --filter_order rrna,contam,fastp
```

`--rrna_filter yes` uses MMseqs2 with a SILVA database. `--contam_filter yes` uses MMseqs2 taxonomy against a downloadable database such as UniRef90.

## Local FASTQ files

If you use local FASTQ files that are not available in public SRA storage, first run [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate). `getfastq` can then process those local entries through the same filtering and staging logic as public runs.

## Parallel processing

`--batch` processes one selected metadata row by one-based index:

```bash
amalgkit getfastq --out_dir ./ --batch 3
```

SLURM example:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --array=1-3

amalgkit getfastq \
    --out_dir ./ \
    --threads "$SLURM_CPUS_PER_TASK" \
    --batch "$SLURM_ARRAY_TASK_ID"
```

For multi-process downloads on shared storage, keep `--download_dir` and `--download_lock_dir` shared across jobs.
