## Getting started

AMALGKIT is a Python-only toolkit for transcriptome amalgamation across studies and species. Current releases do not require R or any R packages.

- [Installation and dependencies](https://github.com/kfuku52/amalgkit/wiki/Installation-and-dependencies)
- [Parallel processing](https://github.com/kfuku52/amalgkit/wiki/Parallel-processing)
- [Tutorial 1](https://github.com/kfuku52/amalgkit/wiki/Tutorial-1)

## Active commands

- [`amalgkit dataset`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-dataset): extract bundled test datasets
- [`amalgkit metadata`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata): retrieve and organize SRA metadata
- [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate): add local FASTQ files to metadata
- [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select): mark samples for downstream analysis using `select_rules.tsv`
- [`amalgkit getfastq`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-getfastq): download or prepare FASTQ files
- [`amalgkit quant`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant): quantify transcript abundance with kallisto or oarfish
- [`amalgkit merge`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-merge): merge per-run abundance tables by species
- [`amalgkit busco`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-busco): prepare BUSCO tables for downstream ortholog-based steps
- [`amalgkit cstmm`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-cstmm): cross-species TMM normalization using single-copy genes
- [`amalgkit wsfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-wsfilter): within-species outlier filtering
- [`amalgkit csfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csfilter): cross-species outlier filtering
- [`amalgkit finalize`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-finalize): export final expression tables, summaries, and batch-corrected outputs
- [`amalgkit sanity`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-sanity): check expected pipeline outputs
- [`amalgkit rerun`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-rerun): rerun failed targets recorded by `sanity`

## Typical pipeline

```bash
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
amalgkit metadata --out_dir ./ --search_string 'YOUR_QUERY'
amalgkit select --out_dir ./
amalgkit getfastq --out_dir ./
amalgkit quant --out_dir ./
amalgkit merge --out_dir ./
amalgkit cstmm --out_dir ./ --dir_busco ./busco        # optional
amalgkit wsfilter --out_dir ./                         # optional
amalgkit csfilter --out_dir ./ --dir_busco ./busco      # optional
amalgkit finalize --out_dir ./
```

`finalize` supports Python backends for `no`, `sva`, `ruvseq`, `combatseq`, and `latent_glm`.

## Legacy command migration

| Old command | Current workflow |
| --- | --- |
| `amalgkit config` | `amalgkit dataset --rule_set ...` writes `select_rules.tsv`; edit that file before `amalgkit select` |
| `amalgkit curate` | run `amalgkit wsfilter`, optionally `amalgkit csfilter`, then `amalgkit finalize` |
| `amalgkit csca` | use `amalgkit csfilter` for cross-species filtering and `amalgkit finalize` for final tables and plots |
