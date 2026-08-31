# Metadata, abundance models, and batch selection

This page describes the current default-branch CLI. Use matching counts and
metadata throughout a workflow; a metadata path is part of the analysis input.

## Which metadata is read?

An explicit `--metadata PATH` takes precedence over inference.

| Command | With `--metadata inferred` |
| --- | --- |
| `getfastq`, `quant`, `merge` | `out_dir/metadata/metadata.tsv` |
| `cstmm` | `dir_count/metadata.tsv`, normally `out_dir/merge/metadata.tsv` |
| `wsfilter`, `csfilter`, `finalize` | Last successful filter recorded in `out_dir/filter_metadata_state.json`; otherwise the selected input directory's `metadata.tsv` |

For filters and finalization, `--input_dir inferred` chooses `out_dir/cstmm`
if present, otherwise `out_dir/merge`. When the state file is absent or invalid,
existing `csfilter/metadata.tsv` takes precedence over `wsfilter/metadata.tsv`,
with a compatibility warning. Modification times do not determine the order.
If no prior filter metadata exists, AMALGKIT uses the selected input directory's
metadata; an invalid state still produces a warning. Review that warning before
continuing or specify the intended metadata explicitly.

To start from a particular stage, specify **both** paths; changing only
`--input_dir` does not disable prior filter metadata inference:

```bash
amalgkit wsfilter --out_dir ./ --input_dir ./cstmm --metadata ./cstmm/metadata.tsv
```

For raw merged counts instead, use `--input_dir ./merge --metadata ./merge/metadata.tsv`.
The resulting filter metadata retains the normalization information needed by
the next filter and `finalize`.

`integrate` does not replace the default metadata file. Explicitly pass its
output to `getfastq`, `quant`, and `merge`; see the
[private FASTQ workflow](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate).
After `merge`, continue with `merge/metadata.tsv`, then `cstmm/metadata.tsv` or
its filtered descendants.

## Counts, lengths, and normalization

| Quantification model | `abundance.tsv` / merged table meaning | Downstream normalization |
| --- | --- | --- |
| kallisto, `length_model=effective` | Estimated counts, effective lengths, length-normalized TPM | FPKM or TPM with raw `merge` counts; FPKM with CSTMM counts |
| oarfish, `length_model=none` | Estimated counts; `length` is the reference annotation length, `eff_length` is a unit placeholder; `tpm` is count abundance per million, without length correction | Use `--norm log2p1-none` with CSTMM counts; raw `merge` counts also support CPM-like `*-tpm` |

Oarfish's unit effective length does **not** make FPKM biologically defined.
AMALGKIT rejects FPKM for this model. The adapter preserves a supplied `tpm`
column after validation; otherwise it computes `count / sum(count) * 1e6`.
Do not interpret that column as kallisto's length-normalized TPM.

`merge/<Species>/<Species>_quant_model.tsv` records each run's backend and length
model, and `cstmm` copies it forward. Keep these sidecars with the tables.
Legacy inputs without a sidecar are treated as effective-length data; that
compatibility behavior is not appropriate for Oarfish data.

The default `--norm` for all three downstream commands is `log2p1-fpkm`.
For a long-read workflow after `merge`, explicitly use:

```bash
amalgkit cstmm --out_dir ./ --dir_busco ./busco
amalgkit wsfilter --out_dir ./ --metadata ./cstmm/metadata.tsv --norm log2p1-none
amalgkit csfilter --out_dir ./ --metadata ./wsfilter/metadata.tsv --dir_busco ./busco --norm log2p1-none
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv --batch_effect_alg no --norm log2p1-none
```

`*-none` uses estimated counts (already divided by TMM factors after CSTMM),
with only the requested log transformation. This defines the input scale; it
does not guarantee that samples from different library protocols are comparable.

CSTMM and TPM cannot be combined: dividing by each corrected library's column
sum would cancel the TMM factor. For the same reason, FPKM from CSTMM counts
requires the **original** `tmm_library_size` in `cstmm/metadata.tsv` or a filtered
descendant. AMALGKIT rejects missing or invalid library sizes instead of
recomputing them from corrected counts. See
[CSTMM normalization metadata](https://github.com/kfuku52/amalgkit/wiki/amalgkit-cstmm#normalization-metadata).

## What does `--batch` select?

| Commands | One-based `--batch N` unit |
| --- | --- |
| `getfastq`, `quant` | Nth run after `is_sampled=yes` filtering, in metadata table order |
| `wsfilter`, `csfilter`, `finalize` | Nth nonempty scientific name in sorted unique order; all its metadata rows are selected before downstream exclusions |

Species batch jobs are not equivalent to one joint cross-species analysis.
In particular, splitting `csfilter` by species removes its cross-species context.
Do not run concurrent filter/finalize batch jobs into the same output directory:
their published metadata and output directories are shared. Prefer internal
parallelism with `--threads`, or separate output directories and explicit input
paths for independent analyses. See
[Parallel processing](https://github.com/kfuku52/amalgkit/wiki/Parallel-processing).
