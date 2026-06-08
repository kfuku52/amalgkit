## Overview

`amalgkit getfastq` turns selected metadata rows into processed FASTQ files. For public data, it downloads SRA objects and extracts FASTQ with `fasterq-dump`. For private data prepared by `integrate`, it stages local FASTQ files into the same workflow.

Optional processing includes:

- `fastp`
- MMseqs2 rRNA filtering
- MMseqs2 contaminant filtering

## Inputs

`getfastq` can start from:

- metadata from `metadata`, `integrate`, or `select`
- one BioProject, BioSample, or run accession through `--id`
- a file of accessions through `--id_list`

Use inferred metadata:

```bash
amalgkit getfastq --out_dir ./
```

Use an explicit metadata table:

```bash
amalgkit getfastq --out_dir ./ --metadata ./metadata/metadata.tsv
```

Generate FASTQ files for one accession:

```bash
amalgkit getfastq --out_dir ./ --id DRR461654
```

## Download Providers

By default, public downloads try enabled providers and continue to the next provider when one is unavailable or throttled.

| Provider | Option |
| --- | --- |
| NCBI cloud objects | `--ncbi yes` |
| AWS | `--aws yes` |
| GCP | `--gcp yes` |
| ENA | `--ena yes` |
| DDBJ for DRA accessions | `--ddbj yes` |

Provider concurrency caps use `--download_lock_dir`, which defaults to `out_dir/downloads/locks`.

| Provider cap | Meaning |
| --- | --- |
| `--ncbi_download_max_concurrency` | cap concurrent NCBI downloads across processes |
| `--aws_download_max_concurrency` | cap concurrent AWS downloads across processes |
| `--gcp_download_max_concurrency` | cap concurrent GCP downloads across processes |
| `--ena_download_max_concurrency` | cap concurrent ENA downloads across processes |
| `--ddbj_download_max_concurrency` | cap concurrent DDBJ downloads across processes |

Set a cap to `0` or `auto` to disable throttling.

## FASTQ Processing

Common options:

| Option | Default | Use |
| --- | --- | --- |
| `--layout single/paired/auto` | `auto` | choose library layout |
| `--max_bp` | very large | target number of bases to extract |
| `--min_read_length` | `25` | minimum read length forwarded through processing |
| `--fastp yes/no` | `yes` | run `fastp` |
| `--remove_sra yes/no` | `yes` | remove downloaded SRA files after extraction |
| `--remove_tmp yes/no` | `yes` | remove temporary files |

When `--max_bp` is set, AMALGKIT can run a compensatory two-step extraction for public SRA-derived data:

1. first-round extraction targets the requested size
2. a second round compensates for reads lost during extraction or filtering

This compensation is for public SRA-derived runs, not private FASTQ files.

## Optional MMseqs2 Filters

```bash
amalgkit getfastq --out_dir ./ --rrna_filter yes
amalgkit getfastq --out_dir ./ --contam_filter yes --contam_filter_rank superkingdom
amalgkit getfastq --out_dir ./ --filter_order rrna,contam,fastp
```

`--rrna_filter yes` uses MMseqs2 with a SILVA database.

`--contam_filter yes` uses MMseqs2 taxonomy against a database such as UniRef90.

`--filter_order` accepts comma or `>` separators, such as `fastp,rrna,contam` or `rrna>contam>fastp`.

## Output Statistics

`getfastq_stats.tsv` contains count and base metrics. For paired-end libraries, count columns can use different units depending on stage:

- `num_dumped`, `num_written`, `num_rrna_in/out`, and `num_contam_in/out`: spot counts
- `num_fastp_in/out`: read counts reported by `fastp`
- `bp_*`: total bases

For stage-by-stage removal fractions, compare `bp_*` columns.

## Array Jobs

`--batch` processes one selected metadata row by one-based index:

```bash
amalgkit getfastq --out_dir ./ --batch 3
```

SLURM example:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --array=1-100

amalgkit getfastq \
    --out_dir ./ \
    --threads "$SLURM_CPUS_PER_TASK" \
    --batch "$SLURM_ARRAY_TASK_ID"
```

For multi-process downloads on shared storage, keep `--download_dir` and `--download_lock_dir` shared across jobs.

## Next Steps

```bash
amalgkit quant --out_dir ./ --build_index yes
amalgkit sanity --out_dir ./ --check getfastq
```
