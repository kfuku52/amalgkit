## Overview

`amalgkit metadata` queries the NCBI SRA and writes the metadata table used by downstream AMALGKIT commands.

Run it in one of two modes:

- one Entrez query with `--search_string`
- species-wise batch queries with `--species_tsv`

## Single Query

```bash
amalgkit metadata \
    --out_dir ./ \
    --entrez_email example@email.com \
    --search_string '("platform illumina"[Properties]) AND ("type rnaseq"[Filter])'
```

Main output:

- `metadata/metadata.tsv`

## Species-Wise Batch Query

Start from the scaffold files created by `dataset --name init`:

```bash
amalgkit dataset --name init --out_dir ./work
amalgkit metadata \
    --out_dir ./work \
    --species_tsv ./work/species.tsv \
    --organ_terms_tsv ./work/organ_terms.tsv \
    --mode title_union \
    --entrez_email example@email.com
```

`--species_tsv` must contain a `scientific_name` column. `--organ_terms_tsv` can provide:

- `sample_group`
- semicolon-separated `title_terms`

`--mode title_union` and `--mode title_split` use those title terms to build organ- or tissue-aware queries.

## Useful Options

| Option | Use |
| --- | --- |
| `--search_string` | one Entrez query |
| `--species_tsv` | species-wise batch queries |
| `--mode base/title_union/title_split` | query construction mode for species-wise runs |
| `--organ_terms_tsv` | sample group and title-term file |
| `--species_limit` | process only the first N species from `--species_tsv` |
| `--merge yes/no` | merge per-query outputs into a species-level table in batch mode |
| `--resolve_names yes/no` | resolve scientific names through NCBI taxonomy IDs |
| `--ncbi_metadata_max_concurrency` | throttle NCBI metadata requests across shared processes |

## Columns to Review

Review at least these columns before running `select`:

- `scientific_name`
- `sample_group`
- `exclusion`
- `bioproject`
- `run`

`scientific_name` groups runs into species-level processing for `merge`, `cstmm`, `wsfilter`, `csfilter`, and `finalize`.

`sample_group` is the main biological grouping column used by `select`, `wsfilter`, `csfilter`, and `finalize`.

Rows with `exclusion != no` are skipped by downstream commands.

## Example Rows

| scientific_name | sample_group | bioproject | run | exclusion |
| --- | --- | --- | --- | --- |
| Nepenthes minima | flower_bud | PRJDB15224 | DRR461729 | no |
| Nepenthes gracilis | root | PRJDB15224 | DRR461727 | no |

## Next Steps

Export or review `select_rules.tsv`, then run:

```bash
amalgkit select --out_dir ./
```

For local FASTQ files, use [amalgkit integrate](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate) before `getfastq`.
