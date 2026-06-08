## Overview

`amalgkit integrate` detects local FASTQ files and appends them to AMALGKIT metadata so they can be processed by the same downstream pipeline as SRA-derived runs.

## Example

```bash
amalgkit integrate --fastq_dir /PATH/TO/FASTQ --out_dir ./ --metadata ./metadata/metadata.tsv
```

Without `--metadata`, standalone mode writes a private FASTQ metadata table.

## FASTQ discovery

Files in `--fastq_dir` are scanned recursively. The first subdirectory below `--fastq_dir` is parsed as `scientific_name`, with underscores converted to spaces.

Example:

```text
/data/private_fastq/Homo_sapiens/brain/Sample1_1.fq.gz
/data/private_fastq/Homo_sapiens/brain/Sample1_2.fq.gz
```

This assigns `scientific_name` to `Homo sapiens`. If the same FASTQ basename appears under multiple species directories, AMALGKIT prefixes the generated run ID with the species directory token.

## FASTQ naming

- paired-end: `Sample1_1.fq.gz` and `Sample1_2.fq.gz`
- single-end: `Sample1.fq.gz`

Supported extensions:

- `.fastq`
- `.fastq.gz`
- `.fq`
- `.fq.gz`

## Metadata fields to review

`integrate` cannot infer every biological label. After creating metadata, review at least:

- `scientific_name`
- `sample_group`
- `exclusion`

`sample_group` is the grouping column used by `wsfilter`, `csfilter`, and `finalize`.

## Example simplified metadata columns

| scientific_name | sample_group | run | read1_path | read2_path | exclusion | lib_layout | private_file |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Mus musculus | brain | MMB1 | /path/MMB1_1.fq.gz | /path/MMB1_2.fq.gz | no | paired | yes |
| Mus musculus | liver | MML1 | /path/MML1_1.fq.gz | /path/MML1_2.fq.gz | no | paired | yes |

## Downstream workflow

```bash
amalgkit integrate --fastq_dir ./private_fastq --out_dir ./
amalgkit getfastq --out_dir ./
amalgkit quant --out_dir ./
amalgkit merge --out_dir ./
amalgkit finalize --out_dir ./
```
