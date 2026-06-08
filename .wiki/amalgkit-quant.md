## Overview

`amalgkit quant` estimates transcript abundance from `getfastq` outputs. The current CLI supports:

- `--quant_backend kallisto` for short-read RNA-seq
- `--quant_backend oarfish` for long-read RNA-seq
- `--quant_backend auto` to choose from metadata

## Examples

Use auto backend selection:

```bash
amalgkit quant --out_dir ./ --threads 8
```

Use existing indices:

```bash
amalgkit quant --out_dir ./ --index_dir ./index
```

Build missing indices from FASTA files:

```bash
amalgkit quant --out_dir ./ --fasta_dir ./fasta --build_index yes
```

Force long-read quantification:

```bash
amalgkit quant \
    --out_dir ./ \
    --quant_backend oarfish \
    --oarfish_seq_tech ont-cdna
```

## Reference FASTA and indices

AMALGKIT expects one reference transcriptome FASTA per species when building indices. If `metadata.tsv` contains *Mus musculus*, AMALGKIT searches for a FASTA file prefixed with `Mus_musculus` under `--fasta_dir`.

Generated index suffixes depend on the backend:

- kallisto: `.idx`
- oarfish: `.mmi`

Shared index build locks prevent concurrent batch jobs from building the same species index at the same time. Tune lock waiting with `--index_lock_poll` and `--index_lock_timeout`.

## Backend options

- `--kallisto_options`: extra shell-style options passed to `kallisto quant`
- `--oarfish_options`: extra shell-style options passed to `oarfish`
- `--oarfish_seq_tech`: long-read sequencing technology preset

## Parallel processing

`--batch` processes one selected metadata row by one-based index:

```bash
amalgkit quant --out_dir ./ --batch 3
```

SLURM example:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --array=1-3

amalgkit quant \
    --out_dir ./ \
    --threads "$SLURM_CPUS_PER_TASK" \
    --batch "$SLURM_ARRAY_TASK_ID"
```

## Output files

Typical per-run outputs include:

- `<RUN>_abundance.tsv`: target ID, length, effective length, estimated counts, and TPM
- `<RUN>_run_info.json`: quantification run information
- backend-specific auxiliary files such as kallisto HDF5 output
