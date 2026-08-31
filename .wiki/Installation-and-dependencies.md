## Installation

AMALGKIT requires Python 3.11 or later; CI covers Python 3.11–3.14 on Linux,
with additional macOS smoke tests. These pages describe the current default
branch, which can be newer than GitHub Releases or Bioconda packages.

```bash
mamba install -c conda-forge -c bioconda --strict-channel-priority amalgkit
amalgkit -h
```

Channel order and strict priority follow the
[Bioconda installation guide](https://bioconda.github.io/index.html#with-conda).

For the latest GitHub revision:

```bash
pip install --upgrade git+https://github.com/kfuku52/amalgkit
amalgkit help metadata
```

## Runtime Model

Current AMALGKIT releases are Python-only. The main pipeline does not require `R`, `Rscript`, or R packages, including the former merge, curation, and batch-correction stages.

Python package dependencies are installed with AMALGKIT. Important runtime libraries include `numpy`, `pandas`, `scipy`, `matplotlib`, `statsmodels`, `scikit-learn`, and `biopython`.

AMALGKIT handles NCBI taxonomy lookups with its built-in SQLite backend. The NCBI taxonomy dump is downloaded automatically when a local taxonomy database is not available.

Existing taxonomy caches created by ETE4 are supported without installing ETE4. AMALGKIT reads the ETE4 `taxa.sqlite` schema directly and reuses a colocated `taxdump.tar.gz` without modifying either file. For commands with `--download_dir`, the cache location remains `<download_dir>/ete_taxonomy`; for the default cache, AMALGKIT also detects ETE4's `~/.local/share/ete` directory. The optional ETE4 `taxa.sqlite.traverse.pkl` file is not needed by AMALGKIT.

## Core External Tools

Some commands call external bioinformatics tools. Install only the tools needed for the workflow you run.

| Tool | Used by | Required when |
| --- | --- | --- |
| [sra-tools / fasterq-dump](https://github.com/ncbi/sra-tools) | `getfastq` | required at startup for every invocation, including private-only runs; use `sra-tools >= 3` |
| [SeqKit](https://github.com/shenwei356/seqkit) | `integrate`, `getfastq` | required by `getfastq`; `integrate` can fall back to Python FASTQ statistics when SeqKit is unavailable |
| [fastp](https://github.com/OpenGene/fastp) | `getfastq` | `--fastp yes`, which is the default |
| [kallisto](https://github.com/pachterlab/kallisto) | `quant` | short-read quantification |

Example environment for a short-read public-SRA workflow:

```bash
mamba create -n amalgkit -c conda-forge -c bioconda \
    --strict-channel-priority amalgkit "sra-tools>=3" seqkit fastp kallisto
mamba activate amalgkit
```

## Optional Dependencies

Install these only when using the corresponding feature:

| Dependency | Required when |
| --- | --- |
| [inmoose](https://inmoose.readthedocs.io/) | `amalgkit finalize --batch_effect_alg combatseq` |
| [oarfish](https://github.com/COMBINE-lab/oarfish) | `quant --quant_backend oarfish`, or `--quant_backend auto` when long-read metadata selects Oarfish |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | `amalgkit getfastq --rrna_filter yes` or `--contam_filter yes` |
| [BUSCO](https://busco.ezlab.org/) | `amalgkit busco --tool busco` |
| [compleasm](https://github.com/huangnengCSU/compleasm) | `amalgkit busco --tool compleasm`, or `--tool auto` when compleasm is selected |

Example:

```bash
mamba install -c conda-forge -c bioconda --strict-channel-priority inmoose oarfish mmseqs2 busco compleasm
```

For Oarfish outputs, also choose a compatible downstream normalization; the
default FPKM is not defined for this backend. See
[metadata and normalization](https://github.com/kfuku52/amalgkit/wiki/Metadata-and-normalization).

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
- `--batch_effect_alg combatseq`, which requires optional `inmoose`
- `--batch_effect_alg latent_glm`

Example:

```bash
amalgkit finalize \
    --out_dir ./ \
    --metadata ./csfilter/metadata.tsv \
    --batch_effect_alg latent_glm
```
