## Overview

`amalgkit metadata` queries the NCBI SRA and writes `metadata.tsv`, the main sample sheet used by downstream AMALGKIT commands.

## Example

Single query:

```bash
amalgkit metadata \
    --out_dir ./ \
    --entrez_email example@email.com \
    --search_string '("platform illumina"[Properties]) AND ("type rnaseq"[Filter])'
```

Species-wise batch query:

```bash
amalgkit dataset --name init --out_dir ./work
amalgkit metadata \
    --out_dir ./work \
    --species_tsv ./work/species.tsv \
    --organ_terms_tsv ./work/organ_terms.tsv \
    --mode title_union \
    --entrez_email example@email.com
```

`--species_tsv` requires a `scientific_name` column. `--organ_terms_tsv` can provide `sample_group` and semicolon-separated `title_terms` for `title_union` and `title_split` query modes.

## Important columns to review manually

- `scientific_name`
- `sample_group`
- `exclusion`
- `bioproject`

## Why `sample_group` matters

`sample_group` is the main biological grouping column used by:

- `amalgkit select`
- `amalgkit wsfilter`
- `amalgkit csfilter`
- `amalgkit finalize`

If multiple samples should be treated as the same biological group, make sure they share the same `sample_group` value.

## Why `scientific_name` matters

Samples with the same `scientific_name` are merged and processed together in species-level steps such as `merge`, `cstmm`, `wsfilter`, and `finalize`.

## Why `exclusion` matters

Rows with `exclusion != no` are skipped by downstream processing.

## Example simplified columns

| scientific_name | sample_group | bioproject | run | exclusion |
| --- | --- | --- | --- | --- |
| Nepenthes minima | flower_bud | PRJDB15224 | DRR461729 | no |
| Nepenthes gracilis | root | PRJDB15224 | DRR461727 | no |

## Related commands

- [amalgkit dataset](https://github.com/kfuku52/amalgkit/wiki/amalgkit-dataset)
- [amalgkit select](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select)
- [amalgkit integrate](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate)
- [amalgkit finalize](https://github.com/kfuku52/amalgkit/wiki/amalgkit-finalize)
