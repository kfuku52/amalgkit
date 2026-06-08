## Legacy Command

`amalgkit curate` has been removed from the current CLI.

The former curation workflow is now split into explicit filtering and export steps.

## Current Workflow

```bash
amalgkit wsfilter --out_dir ./
amalgkit csfilter \
    --out_dir ./ \
    --metadata ./wsfilter/metadata.tsv \
    --dir_busco ./busco
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg no
```

## Migration Guide

| Former role of `curate` | Current command |
| --- | --- |
| within-species outlier filtering | `amalgkit wsfilter` |
| cross-species outlier filtering | `amalgkit csfilter` |
| final per-species table export | `amalgkit finalize` |
| optional batch-effect removal | `amalgkit finalize --batch_effect_alg ...` |

`wsfilter` and `csfilter` each write filtered `metadata.tsv` and `excluded.tsv`. Run one or both filters, then run `finalize` once from the filtered metadata you want to use.

## Minimal Replacement

If you only need final expression tables and do not want filtering:

```bash
amalgkit finalize --out_dir ./ --metadata ./metadata/metadata.tsv --batch_effect_alg no
```

If you want within-species filtering only:

```bash
amalgkit wsfilter --out_dir ./
amalgkit finalize --out_dir ./ --metadata ./wsfilter/metadata.tsv --batch_effect_alg no
```
