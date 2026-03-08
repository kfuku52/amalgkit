## Overview

`amalgkit integrate` detects local FASTQ files and appends them to AMALGKIT metadata so they can be processed by the same downstream pipeline as SRA-derived runs.

## Example

```bash
amalgkit integrate --fastq_dir /PATH/TO/FASTQ --out_dir ./ --metadata ./metadata/metadata.tsv
```

## FASTQ naming

Files in `--fastq_dir` must be stored directly in that directory.

- paired-end: `Sample1_1.fq.gz` and `Sample1_2.fq.gz`
- single-end: `Sample1.fq.gz`

Supported extensions:

- `.fastq`
- `.fastq.gz`
- `.fq`
- `.fq.gz`

## Metadata fields to fill manually

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
amalgkit integrate ...
amalgkit getfastq ...
amalgkit quant ...
amalgkit merge ...
amalgkit finalize ...
```
