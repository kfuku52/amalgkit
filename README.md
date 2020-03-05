## Overview
**amalgkit** contains tools to amalgamate RNA-seq data from diverse research projects to enable a large-scale evolutionary gene expression analysis.

## Dependency
* [python 3](https://www.python.org/)
* [biopython](https://biopython.org/)
* [numpy](https://github.com/numpy/numpy)
* [pandas](https://github.com/pandas-dev/pandas)
* [lxml](https://lxml.de/)

## Installation
```
# Installation with pip
pip install git+https://github.com/kfuku52/amalgkit

# This should show complete options
amalgkit -h
```

## `amalgkit metadata` – SRA metadata curation
**amalgkit metadata** is a subcommand that fetches and curates metadata from the [NCBI's SRA database](https://www.ncbi.nlm.nih.gov/sra). This program needs many config files to enable a tailored metadata curation. See `/amalgkit/config/test/`. Currently, the config files are available only for RNA-seq data from vertebrate organs. To get a fairly good metadata for other taxa/tissues, you would have to extensively edit the config files. 

#### Subcommand dependency
- Nothing

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

## `amalgkit getfastq` – Generate assembly-ready fastq
**amalgkit getfastq** takes a BioProject/BioSample/SRA ID as input and generates RNA-seq fastq files for transcriptome assembly. In the assembly process, the more RNA-seq libraries you include, the more transcripts you get. However, it's often computationally challenging to get an assembly from overwhelming amount of data. **amalgkit getfastq** can automatically subsample RNA-seq reads from different libraries. The amount of data you need (specified by `--max_bp`) depends on many factors including the assembly program you use. See [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146062) for example.

#### Subcommand dependency
- [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump) for `--pfd yes` (default)
- [fastp](https://github.com/OpenGene/fastp) for `--fastp yes` (default)

#### Test run
```
mkdir fastq_files

amalgkit getfastq \
--entrez_email 'aaa@bbb.com' \
--id 'PRJDB4514' \
--threads 2 \
--work_dir ./fastq_files \
--max_bp '75,000'
```
## `amalgkit quant` - quantification of RNAseq data
**amalkit quant** quantifies abundances of transcripts from RNAseq data using Kallisto. All required input and intermediary files are assumed to be in the working directory (default `./`).

#### Subcommand dependencies
- [kallisto](https://pachterlab.github.io/kallisto/)

#### Other specifics
- Needs fastq files (single end or paired end) for quantification, ideally processed by `amalgkit getfastq`, but should be able to handle custom data as well.
- Needs a reference file (usually a fasta file of cdna sequences) for index building, if `--build_index yes` (default), OR an index file if `--build_index no`
- `--index` is either the name given to the index file (default: `id_name.idx`) for index building (optional in this case), or index file if `build_oindex no`
- results are stored in `results_quant`

#### Usage example
#### Contents of working directory:
- `SRR8819967_1.amalgkit.fastq.gz`
- `SRR8819967_2.amalgkit.fastq.gz`
- `arabidopsis_thaliana.fasta` <- this is a reference genome

```
amalgkit quant \
--id SRR8819967 \
--index arabidopsis_thaliana.idx \
--ref arabidopsis_thaliana.fasta \
--work_dir ./fastq_files
```

#### Output
* **SRR8819967_abundance.h5**: bootstrap results in `h5dump` format
* **SRR8819967_run_info.json**: contains run info
* **SRR8819967_abundance.tsv**: contains target_id, lentgh, eff_length, est_counts and tpm in human readable .tsv


## `amalgkit curate` - transcriptome curation

### Subcommand dependencies
- [R](https://www.r-project.org), with various libraries:
    - Biobase
    - pcaMethods
    - colorspace
    - RColorBrewer
    - sva
    - MASS
    - NMF
    - dendextend
    - amap
    - pvclust
    - Rtsne
    - vioplot

### Other specifics
- needs transcriptome file, containing expression data for a single species across multiple runs/conditions/tissues
- needs metadata file, containing SRA runs to be curated

### Usage example

```
amalgkit curate\
--infile transcriptome.tsv
--metafile metadata.tsv
--dist_method 'pearson'
--tissues brain liver heart embryo
--work_dir './'
```
#### Output



## What comes next?
After metadata curation, expression level quantification and further curations had been done in the [paper](https://www.biorxiv.org/content/10.1101/409888v1) where we described the transcriptome amalgamation. The downstream analyses will be added to **amalgkit** as sub-commands in future. Meanwhile, unpackaged scripts we used in the paper are available in `/amalgkit/util/`. For example, **kallisto_20180207.sh** is the quantification step we performed immediately after the metadata curation.

## Licensing
**amalgkit** is BSD-licensed (3 clause). See [LICENSE](LICENSE) for details.