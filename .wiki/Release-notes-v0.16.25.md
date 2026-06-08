## Draft release notes for v0.16.25

These notes summarize the `kfdevel` line prepared for merge into `master`.

## Major workflow changes

- AMALGKIT is now Python-only for the main pipeline. R scripts and R smoke tests have been removed.
- `amalgkit config` has been removed. Selection rules now live in one `select_rules.tsv` file.
- `amalgkit curate` has been replaced by the split workflow `wsfilter -> csfilter -> finalize`.
- `amalgkit csca` has been replaced by current cross-species filtering and finalization outputs.
- `amalgkit rerun` was added to rerun failed targets from `sanity_report.json`.

## Selection workflow

- `amalgkit dataset --rule_set base|test|plantae|vertebrate` exports bundled `select_rules.tsv`.
- `amalgkit dataset --name init` creates an empty workspace scaffold.
- `amalgkit select` supports species-wise batch mode and writes summary, queue, and manifest tables.

## Download and quantification

- `getfastq` supports NCBI, AWS, GCP, ENA, and DDBJ provider fallbacks.
- Shared download locks and provider concurrency caps help large parallel runs.
- Optional MMseqs2 rRNA and contaminant filters can be ordered with `--filter_order`.
- `quant` can auto-select kallisto for short-read data and oarfish for long-read data.
- Quant index builds use shared locks to avoid duplicate builds under array jobs.

## Filtering and final output

- `wsfilter` performs within-species filtering and writes filtered metadata plus PDFs.
- `csfilter` performs cross-species filtering with BUSCO or orthogroup inputs.
- `finalize` exports final per-species expression tables and supports Python backends for `no`, `sva`, `ruvseq`, `combatseq`, and `latent_glm`.

## Validation and recovery

- `sanity` checks getfastq, index, quant, merge, busco, and finalize outputs.
- `rerun` resolves missing or failed targets from sanity reports and writes `rerun_manifest.json`.

## Migration summary

| Old command | Current replacement |
| --- | --- |
| `amalgkit config` | `amalgkit dataset --rule_set ...` and `select_rules.tsv` |
| `amalgkit curate` | `amalgkit wsfilter`, `amalgkit csfilter`, `amalgkit finalize` |
| `amalgkit csca` | `amalgkit csfilter` and `amalgkit finalize` |
