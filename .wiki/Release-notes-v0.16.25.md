## Release Notes for v0.16.25

These notes summarize the CLI line merged into `master` through PR #166.

## Major Workflow Changes

- AMALGKIT is Python-only for the main pipeline. R scripts and R smoke tests have been removed.
- `amalgkit config` has been removed. Selection behavior now lives in one `select_rules.tsv` file.
- `amalgkit curate` has been replaced by the split workflow `wsfilter -> csfilter -> finalize`.
- `amalgkit csca` has been replaced by `csfilter` plus downstream `finalize` outputs.
- `amalgkit rerun` has been added to rerun failed targets recorded by `sanity_report.json`.

## Workspace and Selection

- `amalgkit dataset --name init` creates a scaffold with `species.tsv`, `organ_terms.tsv`, `select_rules.tsv`, and standard workspace directories.
- `amalgkit dataset --rule_set base|test|plantae|vertebrate` exports bundled selection rules to `out_dir/select_rules.tsv`.
- `amalgkit select` reads `select_rules.tsv` through `--select_rules_tsv` and supports species-wise batch mode with summary, queue, and manifest outputs.

## Metadata, FASTQ, and Quantification

- `metadata` supports one-off Entrez search strings and species-wise batch queries from `--species_tsv`.
- `integrate` scans local FASTQ directories and appends private runs to AMALGKIT metadata.
- `getfastq` supports NCBI, AWS, GCP, ENA, and DDBJ provider fallbacks.
- Shared download locks and provider concurrency caps help large multi-process runs.
- Optional MMseqs2 rRNA and contaminant filters can be ordered with `--filter_order`.
- `quant` auto-selects kallisto for short-read runs and oarfish for long-read runs when `--quant_backend auto`.
- Quant index builds use shared locks to avoid duplicate species/backend index builds under array jobs.

## Filtering and Final Output

- `merge` consolidates run-level abundance files into per-species count, effective-length, and TPM tables.
- `cstmm` performs cross-species TMM normalization with BUSCO or orthogroup inputs.
- `wsfilter` performs within-species outlier filtering and writes filtered metadata plus QC PDFs.
- `csfilter` performs cross-species outlier filtering with BUSCO or orthogroup inputs.
- `finalize` exports final per-species expression tables and supports Python backends for `no`, `sva`, `ruvseq`, `combatseq`, and `latent_glm`.

## Validation and Recovery

- `sanity` checks `getfastq`, `index`, `quant`, `merge`, `busco`, and `finalize` outputs.
- `rerun` resolves missing or failed targets from sanity reports and writes `rerun_manifest.json` unless disabled with `--manifest none`.

## Migration Summary

| Old command | Current replacement |
| --- | --- |
| `amalgkit config` | `amalgkit dataset --rule_set ...` and `select_rules.tsv` |
| `amalgkit curate` | `amalgkit wsfilter`, `amalgkit csfilter`, `amalgkit finalize` |
| `amalgkit csca` | `amalgkit csfilter` and `amalgkit finalize` |
