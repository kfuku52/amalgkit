## Overview

`amalgkit busco` prepares per-species BUSCO tables used by ortholog-based downstream steps.

Typical downstream uses:

- `amalgkit cstmm`
- `amalgkit csfilter`

## Examples

Using automatic tool selection:

```bash
amalgkit busco --out_dir ./ --tool auto --lineage eukaryota_odb12
```

Using BUSCO:

```bash
amalgkit busco --out_dir ./ --tool busco --lineage eukaryota_odb12
```

Using compleasm:

```bash
amalgkit busco --out_dir ./ --tool compleasm --lineage eukaryota_odb12
```

With `--tool auto`, AMALGKIT prefers compleasm when available and falls back to BUSCO.

## Input FASTA

By default, AMALGKIT reads FASTA files from `out_dir/fasta`. File names should start with the species name using underscores, such as:

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

## Main outputs

- `busco/<Species>_busco.tsv`
- `busco/busco_completeness.pdf`
- summary files needed by `cstmm` and `csfilter`

## Notes

- `busco` itself is a wrapper step; the required external executable depends on `--tool`.
- If you already have compatible BUSCO full tables, skip this command and pass `--dir_busco` directly to downstream commands.
- `--download_dir` and `--download_lock_dir` can be shared by parallel jobs to avoid duplicate lineage downloads.
