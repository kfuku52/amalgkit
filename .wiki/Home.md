## Getting started

AMALGKIT is a Python-only toolkit for transcriptome amalgamation across studies and species. Current releases do not require R or any R packages.

- [Installation and dependencies](https://github.com/kfuku52/amalgkit/wiki/Installation-and-dependencies)
- [Parallel processing](https://github.com/kfuku52/amalgkit/wiki/Parallel-processing)
- [Tutorial 1](https://github.com/kfuku52/amalgkit/wiki/Tutorial-1)

## Active commands

- [`amalgkit dataset`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-dataset): extract bundled test datasets
- [`amalgkit metadata`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata): retrieve and curate SRA metadata
- [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate): add local FASTQ files to metadata
- [`amalgkit config`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-config): create config templates for sample selection
- [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select): mark samples for downstream analysis
- [`amalgkit getfastq`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-getfastq): download or prepare FASTQ files
- [`amalgkit quant`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant): quantify transcript abundance with kallisto
- [`amalgkit merge`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-merge): merge per-run abundance tables by species
- [`amalgkit busco`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-busco): prepare BUSCO tables for downstream ortholog-based steps
- [`amalgkit cstmm`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-cstmm): cross-species TMM normalization using single-copy genes
- [`amalgkit wsfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-wsfilter): within-species outlier filtering
- [`amalgkit csfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csfilter): cross-species outlier filtering
- [`amalgkit finalize`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-finalize): export final expression tables, summaries, and batch-corrected outputs
- [`amalgkit sanity`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-sanity): check expected pipeline outputs

## Typical pipeline

```bash
amalgkit metadata
amalgkit select
amalgkit getfastq
amalgkit quant
amalgkit merge
amalgkit cstmm          # optional
amalgkit wsfilter       # optional
amalgkit csfilter       # optional
amalgkit finalize
```

`finalize` supports Python backends for `no`, `sva`, `ruvseq`, `combatseq`, and `latent_glm`.

## Legacy pages

The pages below are kept only as historical reference for older workflows. They are not exposed by the current CLI.

- [`amalgkit curate` (legacy)](https://github.com/kfuku52/amalgkit/wiki/amalgkit-curate)
- [`amalgkit csca` (legacy)](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csca)
