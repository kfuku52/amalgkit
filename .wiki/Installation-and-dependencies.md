## Installation

AMALGKIT runs on Python 3.9 or later.

```bash
pip install git+https://github.com/kfuku52/amalgkit
amalgkit -h
```

## Runtime model

Current releases are Python-only. `R`, `Rscript`, and R packages are not required for `merge`, `cstmm`, `wsfilter`, `csfilter`, or `finalize`.

## Python packages installed with AMALGKIT

| Package | Used for | Required when |
| --- | --- | --- |
| [numpy](https://github.com/numpy/numpy) | core array operations | always |
| [pandas](https://github.com/pandas-dev/pandas) | metadata and table I/O | always |
| [scipy](https://scipy.org/) | statistics and linear algebra | always |
| [matplotlib](https://matplotlib.org/) | PDF and figure generation | always |
| [statsmodels](https://www.statsmodels.org/) | GLM fitting for `ruvseq` and `latent_glm` backends | always |
| [inmoose](https://inmoose.readthedocs.io/) | `combatseq` backend | always |
| [biopython](https://biopython.org/) | sequence utilities | always |
| [ete4](https://github.com/etetoolkit/ete) | taxonomy utilities used by `getfastq` | always |

## External tools

| Tool | Used in step(s) | Required when |
| --- | --- | --- |
| [SeqKit](https://github.com/shenwei356/seqkit) | `amalgkit integrate`, `amalgkit getfastq` | recommended for `integrate`; used by `getfastq` workflows |
| [fasterq-dump](https://github.com/ncbi/sra-tools) | `amalgkit getfastq` | required for public SRA download |
| [fastp](https://github.com/OpenGene/fastp) | `amalgkit getfastq` | required when `--fastp yes` (default) |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | `amalgkit getfastq` | required when `--rrna_filter yes` or `--contam_filter yes` |
| [kallisto](https://github.com/pachterlab/kallisto) | `amalgkit quant` | required |
| [BUSCO](https://busco.ezlab.org/) | `amalgkit busco` | required when `--tool busco` |
| [compleasm](https://github.com/huangnengCSU/compleasm) | `amalgkit busco` | required when `--tool compleasm` |

## Commands with no extra external tool requirement

Beyond the Python packages above, these commands do not require additional executables:

- `amalgkit dataset`
- `amalgkit metadata`
- `amalgkit select`
- `amalgkit config`
- `amalgkit merge`
- `amalgkit cstmm` if BUSCO or orthogroup inputs are already available
- `amalgkit wsfilter`
- `amalgkit csfilter`
- `amalgkit finalize`
- `amalgkit sanity`

## Notes on batch correction

`amalgkit finalize` provides Python implementations for:

- `--batch_effect_alg sva`
- `--batch_effect_alg ruvseq`
- `--batch_effect_alg combatseq`
- `--batch_effect_alg latent_glm`

Example:

```bash
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg latent_glm
```
