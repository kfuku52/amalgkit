## Installation

AMALGKIT runs on Python 3.9 or later.

```bash
mamba install -c bioconda amalgkit
amalgkit -h
```

For the latest GitHub revision:

```bash
pip install git+https://github.com/kfuku52/amalgkit
amalgkit help metadata
```

## Runtime Model

Current AMALGKIT releases are Python-only. The main pipeline does not require `R`, `Rscript`, or R packages, including the former merge, curation, and batch-correction stages.

Python package dependencies are installed with AMALGKIT. Important runtime libraries include `numpy`, `pandas`, `scipy`, `matplotlib`, `statsmodels`, `inmoose`, `biopython`, and `ete4`.

## External Tools

Some commands call external bioinformatics tools. Install only the tools needed for the workflow you run.

| Tool | Used by | Required when |
| --- | --- | --- |
| [sra-tools / fasterq-dump](https://github.com/ncbi/sra-tools) | `getfastq` | public SRA extraction |
| [SeqKit](https://github.com/shenwei356/seqkit) | `integrate`, `getfastq` | FASTQ statistics and compression workflows |
| [fastp](https://github.com/OpenGene/fastp) | `getfastq` | `--fastp yes`, which is the default |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | `getfastq` | `--rrna_filter yes` or `--contam_filter yes` |
| [kallisto](https://github.com/pachterlab/kallisto) | `quant` | short-read quantification |
| [oarfish](https://github.com/COMBINE-lab/oarfish) | `quant` | long-read quantification |
| [BUSCO](https://busco.ezlab.org/) | `busco` | `--tool busco` |
| [compleasm](https://github.com/huangnengCSU/compleasm) | `busco` | `--tool compleasm`, or `--tool auto` when compleasm is selected |

Example environment for a short-read public-SRA workflow:

```bash
mamba create -n amalgkit -c conda-forge -c bioconda \
    amalgkit "sra-tools>=3" seqkit fastp kallisto
mamba activate amalgkit
```

Add optional tools as needed:

```bash
mamba install -c conda-forge -c bioconda mmseqs2 busco compleasm
```

## Commands Without Extra Executables

Beyond AMALGKIT's Python dependencies, these commands do not require separate command-line programs:

- `amalgkit dataset`
- `amalgkit metadata`
- `amalgkit select`
- `amalgkit merge`
- `amalgkit cstmm` when BUSCO or orthogroup inputs already exist
- `amalgkit wsfilter`
- `amalgkit csfilter`
- `amalgkit finalize`
- `amalgkit sanity`
- `amalgkit rerun`, except for rerun targets that invoke `getfastq`, `quant`, or `busco`

## Selection Rules

The former `amalgkit config` command has been replaced by `select_rules.tsv`.

```bash
amalgkit dataset --rule_set base --out_dir ./ --overwrite yes
amalgkit select --out_dir ./
```

Available bundled rule sets are `base`, `test`, `plantae`, and `vertebrate`. Use `amalgkit dataset --list` to inspect available datasets and rule sets.

## Batch Correction

`amalgkit finalize` provides Python implementations for all current batch-correction choices:

- `--batch_effect_alg no`
- `--batch_effect_alg sva`
- `--batch_effect_alg ruvseq`
- `--batch_effect_alg combatseq`
- `--batch_effect_alg latent_glm`

Example:

```bash
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg latent_glm
```
