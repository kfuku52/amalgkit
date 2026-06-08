## Legacy command

`amalgkit csca` has been removed from the current CLI.

Use `amalgkit csfilter` for cross-species filtering and QC plots, then use `amalgkit finalize` for final expression tables.

```bash
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg no
```

## Migration guide

| Former role of `csca` | Current workflow |
| --- | --- |
| cross-species correlation/outlier review | `amalgkit csfilter` |
| BUSCO-based ortholog input | `amalgkit csfilter --dir_busco ./busco` |
| OrthoFinder-style orthogroup input | `amalgkit csfilter --orthogroup_table PATH` |
| final exported tables | `amalgkit finalize` |

`csfilter` accepts either BUSCO full tables or an orthogroup table and writes `csfilter/metadata.tsv`, `csfilter/excluded.tsv`, and cross-species QC PDFs.
