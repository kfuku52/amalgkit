## Overview

`amalgkit quant` estimates transcript abundance from `getfastq` outputs.

Supported backends:

| Backend | Use |
| --- | --- |
| `--quant_backend auto` | choose from metadata |
| `--quant_backend kallisto` | short-read RNA-seq |
| `--quant_backend oarfish` | long-read RNA-seq |

## Basic Use

Auto-select the backend using already prepared indices:

```bash
amalgkit quant --out_dir ./ --threads 8
```

Build missing indices from FASTA files:

```bash
amalgkit quant \
    --out_dir ./ \
    --fasta_dir ./fasta \
    --build_index yes
```

Use an existing index directory:

```bash
amalgkit quant --out_dir ./ --index_dir ./index
```

Force long-read quantification:

```bash
amalgkit quant \
    --out_dir ./ \
    --quant_backend oarfish \
    --oarfish_seq_tech ont-cdna
```

## Reference FASTA and Indices

When `--build_index yes` is set, AMALGKIT expects one reference transcriptome FASTA per species under `--fasta_dir`.

`--build_index` defaults to `no`; `quant` does not download reference FASTAs.
Use `--build_index yes` for a new workspace after supplying the references, or
provide existing indices compatible with the chosen backend. Automatic backend
selection also requires Oarfish when metadata identifies a long-read run.

If metadata contains `Mus musculus`, AMALGKIT searches for a FASTA file prefixed with `Mus_musculus`.

Accepted FASTA suffixes include:

- `.fa`
- `.fasta`
- `.fa.gz`
- `.fasta.gz`

Generated index suffixes depend on the selected backend:

- kallisto: `.idx`
- oarfish: `.mmi`

Shared index-build locks prevent concurrent batch jobs from building the same species/backend index. Tune waiting with:

- `--index_lock_poll`
- `--index_lock_timeout`

## Backend Options

| Option | Use |
| --- | --- |
| `--kallisto_options` | extra shell-style options passed to `kallisto quant` |
| `--oarfish_options` | extra shell-style options passed to `oarfish` |
| `--oarfish_seq_tech` | long-read sequencing technology preset |
| `--clean_fastq yes/no` | remove processed FASTQ files after successful quantification |

`--oarfish_seq_tech auto` infers ONT/PacBio subtype from metadata when possible.

Extra backend options cannot override AMALGKIT-managed inputs, output paths,
indices, thread counts, layouts, or sequencing-technology flags.

## Array Jobs

`--batch` processes one run by one-based index after `is_sampled=yes` filtering,
preserving metadata row order (unlike the species batches used by filters):

```bash
amalgkit quant --out_dir ./ --batch 3
```

SLURM example:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --array=1-100

amalgkit quant \
    --out_dir ./ \
    --threads "$SLURM_CPUS_PER_TASK" \
    --batch "$SLURM_ARRAY_TASK_ID"
```

## Main Outputs

Typical per-run outputs include:

- `<RUN>_abundance.tsv`
- `<RUN>_run_info.json`
- backend-specific auxiliary files

`<RUN>_abundance.tsv` contains target ID, length, effective length, estimated counts, and TPM.

For kallisto these are effective-length-normalized TPM values. For Oarfish,
`length` is the annotated transcript length, `eff_length=1` is a placeholder,
and `tpm` represents abundance per million without length correction (a supplied
backend `tpm` column is preserved; otherwise it is computed from counts).
Run-info records `length_model=none`, which `merge` propagates to a per-species
quant-model sidecar. Do not apply FPKM to these runs. For CSTMM plus Oarfish,
use `--norm log2p1-none` in `wsfilter`, `csfilter`, and `finalize`; TPM is also
incompatible with CSTMM. See [metadata and normalization](https://github.com/kfuku52/amalgkit/wiki/Metadata-and-normalization).

When selection columns are populated, only rows with `exclusion == no` and
`is_sampled == yes` are quantified. A run is complete only after its abundance
table and run-info JSON pass schema, row, finite-number, nonnegative-value, and
pseudoalignment-range checks. `--redo yes` builds in a staging directory and
replaces an existing valid result only after the new result validates.

With the default `--clean_fastq yes`, owned FASTQ entries are replaced by
`.safely_removed` markers after successful quantification. For private gzip
inputs linked by `getfastq`, only those managed links are retired; source files
are never removed. Cleanup failures restore the links/files and prior markers.
If restoration itself fails, the error names the preserved recovery directory.

## Next Steps

```bash
amalgkit merge --out_dir ./
amalgkit sanity --out_dir ./ --check quant
```
