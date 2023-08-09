![](logo/logo_amalgkit_large.png)

## Overview
**AMALGKIT** ([/əm`ælgkit/](http://ipa-reader.xyz/?text=%C9%99m%60%C3%A6lgkit&voice=Joanna)) is a toolkit to integrate RNA-seq data from [the NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra) and from private fastq files to generate unbiased cross-species transcript abundance dataset for a large-scale evolutionary gene expression analysis.

## Installation
```
# Installation with pip
pip install git+https://github.com/kfuku52/amalgkit

# This should show complete options
amalgkit -h
```

## Functions
See [Wiki](https://github.com/kfuku52/amalgkit/wiki) for details.

- [`amalgkit metadata`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata): NCBI SRA metadata retrieval

- [`amalgkit config`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-config): Creating a series of config files for the metadata selection

- [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select): Selecting SRA entries for analysis

- [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate): Appending local fastq info to a metadata table

- [`amalgkit getfastq`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-getfastq): Generating fastq files

- [`amalgkit quant`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant): Transcript abundance estimation

- [`amalgkit merge`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-merge): Generating transcript abundance tables

- [`amalgkit cstmm`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-cstmm): Cross-species TMM normalization using single-copy genes

- [`amalgkit curate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-curate): Automatic removal of outlier samples and unwanted biases

- [`amalgkit csca`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csca): Generating plots with cross-species correlation analysis

- [`amalgkit sanity`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-sanity): Checking the integrity of AMALGKIT input and output files

## Citation
Although **AMALGKIT** supports novel unpublished functions, some functionalities including metadata curation, expression level quantification, and further curation steps have been described in this paper, in which we reported the transcriptome amalgamation of 21 vertebrate species.

Fukushima K*, Pollock DD*. 2020. Amalgamated cross-species transcriptomes reveal organ-specific propensity in gene expression evolution. Nature Communications 11: 4459 (DOI: 10.1038/s41467-020-18090-8) [open access](https://www.nature.com/articles/s41467-020-18090-8)

## Licensing
**amalgkit** is BSD-licensed (3 clause). See [LICENSE](LICENSE) for details.
