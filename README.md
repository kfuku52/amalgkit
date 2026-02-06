![](logo/logo_amalgkit_large.png)

[![Run Tests](https://github.com/kfuku52/amalgkit/actions/workflows/tests.yml/badge.svg)](https://github.com/kfuku52/amalgkit/actions/workflows/tests.yml)
[![GitHub release](https://img.shields.io/github/v/tag/kfuku52/amalgkit?label=release)](https://github.com/kfuku52/amalgkit/releases)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/amalgkit.svg)](https://anaconda.org/bioconda/amalgkit)
[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)](https://github.com/kfuku52/amalgkit)
[![Platforms](https://img.shields.io/conda/pn/bioconda/amalgkit.svg)](https://anaconda.org/bioconda/amalgkit)
[![Downloads](https://img.shields.io/conda/dn/bioconda/amalgkit.svg)](https://anaconda.org/bioconda/amalgkit)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Overview
**AMALGKIT** ([/əm`ælgkit/](http://ipa-reader.xyz/?text=%C9%99m%60%C3%A6lgkit&voice=Joanna)) is a toolkit to integrate RNA-seq data from [the NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra) and from private fastq files to generate unbiased cross-species transcript abundance dataset for a large-scale evolutionary gene expression analysis.

![](logo/flowchart_00.png)

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

- [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate): Appending local fastq info to a metadata table

- [`amalgkit config`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-config): Creating a series of config files for the metadata selection

- [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select): Selecting SRA entries for analysis

- [`amalgkit getfastq`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-getfastq): Generating fastq files

- [`amalgkit quant`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant): Transcript abundance estimation

- [`amalgkit merge`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-merge): Generating transcript abundance tables

- [`amalgkit cstmm`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-cstmm): Cross-species TMM normalization using single-copy genes

- [`amalgkit curate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-curate): Automatic removal of outlier samples and unwanted biases

- [`amalgkit csca`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csca): Generating plots with cross-species correlation analysis

- [`amalgkit sanity`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-sanity): Checking the integrity of AMALGKIT input and output files

- [`amalgkit dataset`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-dataset): Extracting bundled test datasets

## Citation
Although **AMALGKIT** supports novel unpublished functions, some functionalities including metadata curation, expression level quantification, and further curation steps have been described in this paper, in which we reported the transcriptome amalgamation of 21 vertebrate species.

Fukushima K*, Pollock DD*. 2020. Amalgamated cross-species transcriptomes reveal organ-specific propensity in gene expression evolution. Nature Communications 11: 4459 (DOI: 10.1038/s41467-020-18090-8) [open access](https://www.nature.com/articles/s41467-020-18090-8)

## Licensing
**amalgkit** is BSD-licensed (3 clause). See [LICENSE](LICENSE) for details.
