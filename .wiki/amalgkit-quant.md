## Overview

`amalgkit quant` estimates transcript abundance from `getfastq` outputs.

Supported backends:

| Backend | Use |
| --- | --- |
| `--quant_backend auto` | choose from metadata |
| `--quant_backend kallisto` | short-read RNA-seq |
| `--quant_backend oarfish` | long-read RNA-seq |

For GSA inputs, read-count checks can be deferred until getfastq measures the
original FASTQs. See [GSA native inputs](./GSA-native-inputs) for the measured
metadata snapshot and downstream handoff.

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

If metadata contains `Mus musculus`, use a FASTA stem matching `Mus_musculus`,
for example `Mus_musculus.fa.gz`. The legacy stem
`Mus_musculus_for_kallisto_index` is also accepted. Extra assembly/version suffixes
are not matched: rename downloaded references to the species stem and retain
release/accession information separately. Keep only one matching reference per
species.

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

## Single-end fragment lengths

For single-end kallisto, supply the mean and standard deviation (SD) of the
sequenced **insert/fragment**, in bp. Read length, `spot_length`, summed mate
lengths and a fastp insert-size peak are not estimates of this mean/SD pair.
Instrument measurements should describe the insert after accounting for adapters.
The [kallisto manual](https://pachterlab.github.io/kallisto/manual) requires both
parameters in single-end mode; its 200/20 example is not a universal measurement.

| Option | Default | Meaning |
| --- | --- | --- |
| `--fragment_length_file` | `None` | Run-specific TSV; highest priority |
| `--fragment_length_mean` | `None` | Common mean in bp; requires common SD |
| `--fragment_length_sd` | `None` | Common SD in bp; requires common mean |
| `--fragment_length_policy` | `assume` | `assume` warns and fills missing components; `error` stops |

Sources are selected as a pair in this order:

1. A row in `--fragment_length_file`.
2. Run-specific metadata `fragment_length_mean` / `fragment_length_sd`.
3. The common CLI mean/SD pair.
4. Metadata `nominal_length` / `nominal_sdev`, or a separate `mean_insert_size`
   candidate when nominal length is absent.

After selecting a source, `assume` fills only missing mean with **200 bp** and
missing SD with **20 bp**. For example, mean 150 with missing SD becomes 150/20;
mean 350 with missing SD becomes 350/20. Known values below 200 and independent
SD values are preserved. Every supplemented run emits a warning naming the
assumption. The SD is never silently calculated as 10% of the mean. Use `error`
when assumptions are unacceptable. These defaults are assumptions, not measured
values or a claim of accurate quantification; assess sensitivity when using them.

Non-numeric values, infinities, zero/negative values, ranges and conflicting
values are errors, not missing-data fallbacks. Blank metadata cells are missing.
SRA nominal fields describe a submitted **expected paired-library insert size**
and its SD, not a verified run-level measurement. `mean_insert_size` remains
separate: disagreement with `nominal_length` requires an explicit override;
its mean cannot borrow an unassociated `nominal_sdev`.
During XML import, SAMPLE attributes named `nominal_length`, `nominal_sdev` or
`fragment_length_*` are kept as `sample_attribute_*` fields for review. They do
not overwrite EXPERIMENT library statistics or become run-specific overrides.
`getfastq` preserves raw fragment fields so invalid text cannot turn into a
missing-value assumption during read-statistic initialization.

Common values for runs without explicit fragment metadata:

```bash
amalgkit quant --out_dir ./ --fragment_length_mean 150 --fragment_length_sd 7
```

For mixed libraries, create `fragments.tsv` with literal tab separators:

```tsv
run	fragment_length_mean	fragment_length_sd	source	source_detail
SRR000001	150	7	measured	Library L1; insert assay after adapter subtraction
SRR000002	250	35	user	Library L2; specified from protocol P2
```

```bash
amalgkit quant --out_dir ./ --fragment_length_file fragments.tsv --fragment_length_policy error
```

File rows require both numeric values and nonempty `source` / `source_detail`.
Run IDs must be unique and present in the input metadata; the file may cover a
subset. Validation precedes batch selection so one full file can serve all array
jobs. File values override metadata explicitly; common CLI values do not override
run-specific fragment metadata. For metadata, use `fragment_length_source` and
`fragment_length_source_detail`; an omitted source means `user`. Non-user sources
such as `measured` require a detail explaining the measurement and library.
AMALGKIT records these declarations but does not independently certify them.
Components are not silently assembled from different source pairs.

`--kallisto_options` still rejects `-l`/`-s` and their long forms; use the dedicated
options above. Paired-end kallisto estimates its distribution from paired reads.
Common options are ignored with a warning for paired-end/Oarfish inputs; a file
row specifically targeting either is an error. The actual FASTQ layout governs
this decision, including runs whose metadata says paired but whose input is single.

## Fragment provenance and re-quantification

Each new kallisto `*_run_info.json` contains `amalgkit_fragment_length` schema 1:
actual and metadata layout, selected mean/SD with supplied values and sources,
assumed components, policy, extra kallisto options, executed command and kallisto version. Paired-end
entries identify kallisto estimation without inventing a measured mean/SD.
The abundance TSV schema is unchanged. Keep the run-info files with merged and
normalized results to retain this provenance.

Completed outputs with this provenance are reused only when the backend, current
input layout, resolved distribution, provenance and extra kallisto options match.
For example, changing `--single-overhang` requires re-quantification. When input
files are gone after cleanup, an unchanged metadata layout retains the recorded
layout correction. Changes require `--redo yes`, including when
FASTQs have already been cleaned up. Legacy outputs remain readable; default
reuse warns that their fragment provenance is unknown. Explicit new parameters,
run-specific fragment metadata or strict policy cannot certify a legacy single-end output
and require re-quantification. Merely updating AMALGKIT does not recalculate it.

To apply corrected parameters, retain or restore the same processed FASTQs, run
`quant` with the chosen parameters and `--redo yes`, then rebuild affected merge
and downstream normalization/filter/final outputs. Restore safely removed inputs
with `getfastq` first if necessary; keep extraction/filter settings unchanged.
Changed normalization inputs can affect other runs sharing normalization factors.
An unchanged compatible reference index can be reused. A generic `rerun` does
not recover historical common CLI or file overrides. Under the run lock, it
checks surviving run-info against current settings and stops rather than silently
replacing those settings with defaults, even when abundance output is damaged.
Pass the chosen settings to `quant` explicitly when this check fails. If no
run-info survived, historical parameters cannot be recovered automatically.

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
