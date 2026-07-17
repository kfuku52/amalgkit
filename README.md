![AMALGKIT logo](https://raw.githubusercontent.com/kfuku52/amalgkit/master/logo/logo_amalgkit_large.png)

[![Run Tests](https://github.com/kfuku52/amalgkit/actions/workflows/tests.yml/badge.svg)](https://github.com/kfuku52/amalgkit/actions/workflows/tests.yml)
[![GitHub release](https://img.shields.io/github/v/tag/kfuku52/amalgkit?label=release)](https://github.com/kfuku52/amalgkit/releases)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/amalgkit.svg)](https://anaconda.org/bioconda/amalgkit)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue)](https://github.com/kfuku52/amalgkit)
[![Platforms](https://img.shields.io/conda/pn/bioconda/amalgkit.svg)](https://anaconda.org/bioconda/amalgkit)
[![Downloads](https://img.shields.io/conda/dn/bioconda/amalgkit.svg)](https://anaconda.org/bioconda/amalgkit)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
**AMALGKIT** (/əm`ælgkit/) is a toolkit to integrate RNA-seq data from [the NCBI SRA database](https://www.ncbi.nlm.nih.gov/sra) and from private fastq files to generate unbiased cross-species transcript abundance dataset for a large-scale evolutionary gene expression analysis.

The README intentionally keeps the workflow summary text-based. The historical flowchart has been removed from this page because it drifts out of date faster than the CLI and wiki documentation.

## Installation
```bash
# Install the latest GitHub version with pip
pip install git+https://github.com/kfuku52/amalgkit

# Or install the packaged Bioconda version
mamba install -c bioconda amalgkit

# Show top-level commands
amalgkit -h

# Show command-specific help
amalgkit help metadata
```

AMALGKIT supports Linux and macOS with Python 3.11 or later. The Bioconda
package can lag behind the latest GitHub release; run `amalgkit --version` when
reproducing an analysis.

`amalgkit getfastq` requires `fasterq-dump` from `sra-tools >= 3` on `PATH`.
If you manage external tools separately, install it explicitly, for example:

```bash
mamba install -c conda-forge -c bioconda "sra-tools>=3"
```

Commands such as `getfastq`, `quant`, and `busco` use additional external
bioinformatics tools. See [Installation and dependencies](https://github.com/kfuku52/amalgkit/wiki/Installation-and-dependencies)
for the command-by-command dependency table.

## Commands
See [Wiki](https://github.com/kfuku52/amalgkit/wiki) for detailed examples and option descriptions.

- [`amalgkit metadata`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata): NCBI SRA metadata retrieval

- [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate): Appending local fastq info to a metadata table

- [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select): Selecting SRA entries for analysis

- [`amalgkit getfastq`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-getfastq): Generating fastq files

- [`amalgkit quant`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant): Transcript abundance estimation with auto-selected kallisto/oarfish backend

- [`amalgkit merge`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-merge): Generating transcript abundance tables

- [`amalgkit busco`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-busco): Generating BUSCO tables for downstream cstmm/csfilter

- [`amalgkit cstmm`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-cstmm): Cross-species TMM normalization using single-copy genes

- [`amalgkit wsfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-wsfilter): Within-species outlier filtering (`metadata.tsv` + `excluded.tsv` + `wsfilter_exclusion.pdf` + `wsfilter/<Species>/<Species>_*.pdf`)

- [`amalgkit csfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csfilter): Cross-species outlier filtering (`metadata.tsv` + `excluded.tsv` + `csfilter_exclusion.pdf` + `csfilter/*.pdf`)

- [`amalgkit finalize`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-finalize): Export final tables from filtered metadata (with optional batch-effect removal)

- [`amalgkit sanity`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-sanity): Checking the integrity of AMALGKIT input and output files

- [`amalgkit rerun`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-rerun): Rerunning failed sanity targets from `sanity_report.json` and writing `rerun_manifest.json`

- [`amalgkit dataset`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-dataset): Extracting bundled test datasets

Legacy commands from earlier AMALGKIT releases have been replaced:

- `amalgkit config` -> `amalgkit dataset --rule_set ...` plus `select_rules.tsv`
- `amalgkit curate` -> `amalgkit wsfilter`, `amalgkit csfilter`, and `amalgkit finalize`
- `amalgkit csca` -> `amalgkit csfilter` and downstream `amalgkit finalize` outputs

## Typical Workflows
### Initialize an empty workspace
```bash
amalgkit dataset --name init --out_dir ./work
```

This writes `WORKSPACE_README.md` plus starter `species.tsv`, `organ_terms.tsv`, and `select_rules.tsv`.

### Metadata to merged quantification tables
```bash
# 1. Retrieve metadata from SRA
amalgkit metadata --search_string 'vertebrata[Organism] AND liver'

# 2. Export/edit select rules, then select runs
amalgkit dataset --out_dir ./ --rule_set base --overwrite yes
amalgkit select --out_dir ./

# 3. Optionally append private FASTQ files to metadata
amalgkit integrate --out_dir ./ --fastq_dir ./private_fastq

# 4. Download/process FASTQ, quantify, and merge per-species abundance tables
amalgkit getfastq --out_dir ./
amalgkit quant --out_dir ./
amalgkit merge --out_dir ./
```

### Cross-species normalization and filtering
```bash
# Prepare single-copy ortholog tables
amalgkit busco --out_dir ./ --lineage eukaryota_odb12

# Cross-species TMM normalization
amalgkit cstmm --out_dir ./ --dir_busco ./busco

# Metadata filtering and final export
amalgkit wsfilter --out_dir ./
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg no
```

## Split Filtering Workflow
`wsfilter` and `csfilter` are decoupled filters that output `metadata.tsv`, `excluded.tsv`, exclusion summary PDF, and species PDFs (without a `plots/` directory).  
Run one or both in any order, then export tables once with `finalize`.
When `--metadata inferred` is used in these commands, the latest filter metadata (`wsfilter/metadata.tsv` or `csfilter/metadata.tsv`) is auto-detected.

```bash
# Example: wsfilter -> csfilter -> finalize
amalgkit wsfilter --out_dir ./
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg no
```

## Bundled Demo Data
AMALGKIT ships with an empty workspace scaffold (`init`) and a small bundled dataset for smoke testing and examples. The `yeast` dataset uses small BUSCO-focused test FASTAs rather than full gene sets, so its BUSCO completeness is intentionally modest.

```bash
amalgkit dataset --list
amalgkit dataset --name init --out_dir ./work
amalgkit dataset --name yeast --out_dir ./demo
```

## Citation
Although **AMALGKIT** supports novel unpublished functions, some functionalities including metadata curation, expression level quantification, and further curation steps have been described in this paper, in which we reported the transcriptome amalgamation of 21 vertebrate species.

Fukushima K*, Pollock DD*. 2020. Amalgamated cross-species transcriptomes reveal organ-specific propensity in gene expression evolution. Nature Communications 11: 4459 (DOI: 10.1038/s41467-020-18090-8) [open access](https://www.nature.com/articles/s41467-020-18090-8)

## Licensing
**amalgkit** is MIT-licensed. See [LICENSE](LICENSE) for details.
