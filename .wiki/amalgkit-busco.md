## Overview

`amalgkit busco` creates per-species BUSCO tables used by ortholog-aware downstream steps.

Typical downstream uses:

- `amalgkit cstmm`
- `amalgkit csfilter`

## Basic Use

`--lineage` is required:

```bash
amalgkit busco --out_dir ./ --lineage eukaryota_odb12
```

Choose the underlying tool explicitly if needed:

```bash
amalgkit busco --out_dir ./ --tool busco --lineage eukaryota_odb12
amalgkit busco --out_dir ./ --tool compleasm --lineage eukaryota_odb12
```

With `--tool auto`, AMALGKIT prefers compleasm when available and falls back to BUSCO.

## Input FASTA

By default, AMALGKIT reads transcriptome FASTA files from:

```text
out_dir/fasta
```

File names should start with the species name using underscores:

```text
Arabidopsis_thaliana.v1.fa.gz
```

For one explicit FASTA:

```bash
amalgkit busco \
    --out_dir ./ \
    --fasta ./fasta/Arabidopsis_thaliana.fa.gz \
    --species "Arabidopsis thaliana" \
    --lineage eukaryota_odb12
```

## Useful Options

| Option | Use |
| --- | --- |
| `--tool auto/busco/compleasm` | choose BUSCO implementation |
| `--lineage` | lineage dataset, such as `eukaryota_odb12` |
| `--fasta_dir` | directory of species FASTA files |
| `--fasta` and `--species` | process one explicit FASTA |
| `--tool_args` | pass additional arguments to the selected tool |
| `--download_dir` and `--download_lock_dir` | share lineage/download locks across jobs |

## Main Outputs

- `busco/<Species>_busco.tsv`
- `busco/busco_completeness.pdf`
- summary files used by `cstmm` and `csfilter`

## Existing BUSCO Tables

If you already have compatible per-species BUSCO full tables, you can skip this command and pass the table directory directly:

```bash
amalgkit cstmm --out_dir ./ --dir_busco ./busco
amalgkit csfilter --out_dir ./ --dir_busco ./busco
```
