## Legacy Command

`amalgkit csca` has been removed from the current CLI.

Use `amalgkit csfilter` for cross-species filtering and QC plots, then use `amalgkit finalize` for final expression tables.

## Current Workflow

```bash
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

| Former role of `csca` | Current workflow |
| --- | --- |
| cross-species correlation/outlier review | `amalgkit csfilter` |
| BUSCO-based ortholog input | `amalgkit csfilter --dir_busco ./busco` |
| OrthoFinder-style orthogroup input | `amalgkit csfilter --orthogroup_table PATH` |
| final exported tables | `amalgkit finalize` |

`csfilter` accepts either BUSCO full tables or an orthogroup table and writes:

- `csfilter/metadata.tsv`
- `csfilter/excluded.tsv`
- `csfilter/csfilter_exclusion.pdf`
- cross-species QC PDFs

See [amalgkit csfilter](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csfilter) and [amalgkit finalize](https://github.com/kfuku52/amalgkit/wiki/amalgkit-finalize) for the active command pages.
