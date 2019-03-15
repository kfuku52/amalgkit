## Overview
**amalgkit** contains tools to amalgamate RNA-seq data from diverse research projects to enable a large-scale evolutionary gene expression analysis.

## Dependency
* [python 3](https://www.python.org/)
* [biopython](https://biopython.org/)
* [numpy](https://github.com/numpy/numpy)
* [pandas](https://github.com/pandas-dev/pandas)

## Installation
```
# Installation with pip
pip install git+https://github.com/kfuku52/amalgkit

# This should show complete options
amalgkit -h 
```

## SRA metadata curation
**amalgkit metadata** is a subcommand that fetches and curates metadata from the [NCBI's SRA database](https://www.ncbi.nlm.nih.gov/sra). This program needs many config files to enable a tailored metadata curation. See `/amalgkit/config/test/`. Currently, the config files are available only for RNA-seq data from vertebrate organs. To get a fairly good metadata for other taxa/tissues, you would have to extensively edit the config files. 

#### Test run

```
mkdir -p amalgkit_out; cd $_

svn export https://github.com/kfuku52/amalgkit/trunk/config

config_dir="./config/test"

amalgkit metadata \
--config_dir ${config_dir} \
--work_dir . \
--entrez_email 'aaa@bbb.com' # Use your own email address. Don't worry, you won't get spam messages.
```

If you get a network connection error, simply rerun the same analysis. The program will resume the analysis using intermediate files in `--work_dir`.

#### Output
* **metadata_01_raw_YYYY_MM_DD-YYYY_MM_DD.tsv**: This table is a reformatted version of SRA metadata in the xml format.
* **metadata_02_grouped_YYYY_MM_DD-YYYY_MM_DD.tsv**: Similar attributes (columns) are grouped into a few categories according to .config settings.
* **metadata_03_curated_YYYY_MM_DD-YYYY_MM_DD.tsv**: A variety of curation steps are applied according to .config settings. Data unsuitable for evolutionary gene expression analysis such as those from miRNA-seq are marked `No` in the `is_qualified` column. There are particular samples which have been intensively sequenced (e.g., livers of *Bos taurus*). Those samples can be subsampled by the `--max_sample` option and excluded data are marked `No` in the `is_sampled` column.
* **pivot_\*.tsv**: "species x tissue" pivot tables.

## What comes next?
After metadata curation, expression level quantification and further curations had been done in the [paper](https://www.biorxiv.org/content/10.1101/409888v1) where we described the transcriptome amalgamation. The downstream analyses will be added to **amalgkit** as sub-commands in future. Meanwhile, unpackaged scripts we used in the paper are available in `/amalgkit/util/`. For example, **kallisto_20180207.sh** is the quantification step we performed immediately after the metadata curation.

## Licensing
**amalgkit** is BSD-licensed (3 clause). See [LICENSE](LICENSE) for details.