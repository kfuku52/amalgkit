## Overview

`amalgkit merge` consolidates per-run `quant` outputs into per-species abundance tables.

The merged tables are the standard input for:

- `amalgkit cstmm`
- `amalgkit wsfilter`
- `amalgkit csfilter`
- `amalgkit finalize`

`merge` is Python-only in current releases.

## Basic Use

Use inferred metadata from `out_dir/metadata/metadata.tsv`:

```bash
amalgkit merge --out_dir ./
```

Use an explicit metadata table:

```bash
amalgkit merge \
    --out_dir ./ \
    --metadata ./metadata/metadata.tsv
```

## Inputs

`merge` expects:

- selected metadata
- completed `quant` outputs for selected runs

Rows with `exclusion != no` are ignored.

Each abundance table must contain data rows, unique nonempty `target_id` values,
and finite nonnegative effective lengths, counts and TPM. Invalid inputs stop
the merge with the run, path and offending column in the error; previously
published merge outputs remain intact. Validation shares the quant/sanity
contract and does not reread abundance tables.

Target identifiers remain text throughout quant, merge, orthology joins and
normalization: `0001`, `1` and `NA` are distinct IDs. Empty cells are missing
identifiers; numeric missing values are still rejected in abundance columns.

## Main Outputs

For each species:

- `merge/<Species>/<Species>_est_counts.tsv`
- `merge/<Species>/<Species>_eff_length.tsv`
- `merge/<Species>/<Species>_tpm.tsv`
- `merge/<Species>/<Species>_quant_model.tsv` (run, backend, and length model)

`merge/metadata.tsv` is the metadata handoff to CSTMM or downstream filters.
Oarfish uses `length_model=none`: its effective-length entries are unit
placeholders, and its `tpm` values have no length correction. Preserve the
quant-model sidecar and use `--norm log2p1-none` for Oarfish in downstream
filters/finalization. See [metadata and normalization](https://github.com/kfuku52/amalgkit/wiki/Metadata-and-normalization)
for backend-specific scales and CSTMM compatibility. The examples below use
short-read, effective-length inputs.

Summary PDFs are also produced at the merge level, including mapping, library-layout, duplication, insert-size, and expression summary plots.

## Parallel Options

`merge` accepts the shared CPU-budget options:

- `--threads`
- `--internal_jobs`
- `--internal_cpu_budget`

Most users should set only `--threads`.

## Next Steps

For a simple workflow:

```bash
amalgkit finalize --out_dir ./ --batch_effect_alg no
```

For cross-species normalization and filtering:

```bash
amalgkit busco --out_dir ./ --lineage eukaryota_odb12
amalgkit cstmm --out_dir ./ --dir_busco ./busco
amalgkit wsfilter --out_dir ./
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv
```
