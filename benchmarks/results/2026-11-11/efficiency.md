# Pipeline efficiency comparison

The repository's revision benchmark compared the pre-session revision
`4b0d5ec` with `9ee929c` on Apple M2 Max (macOS arm64, Python 3.14.7,
NumPy 2.5.3, pandas 3.0.6, SciPy 1.18.1). Each workload used one
native thread, one warmup pair and three measured pairs in alternating
revision order. The result arrays were identical in all six cases.

| Workload | Wall time (s) | CPU time (s) | Peak RSS (MiB) |
| --- | ---: | ---: | ---: |
| wsfilter | 0.291 → 0.287 | 0.290 → 0.281 | 253.0 → 256.1 |
| csfilter | 0.302 → 0.157 | 0.300 → 0.157 | 113.7 → 110.1 |
| Imputation | 0.045 → 0.036 | 0.044 → 0.035 | 180.8 → 167.3 |
| Tissue specificity | 0.136 → 0.071 | 0.134 → 0.071 | 154.1 → 154.6 |
| TMM rounds | 0.320 → 0.275 | 0.318 → 0.272 | 127.3 → 120.9 |
| TPM transformation | 0.031 → 0.015 | 0.031 → 0.015 | 259.5 → 256.7 |

These are medians of kernel measurements. RSS includes process startup and
input construction. The wsfilter timings and the small RSS increases in
wsfilter and tissue specificity do not establish a meaningful change for
those metrics in this workload. The sizes and generated cases do not represent
every species, missingness pattern, or end-to-end workflow.

To reproduce the comparison from a checkout at `9ee929c`:

```bash
baseline_dir=$(mktemp -d)
results_dir=$(mktemp -d)
git archive 4b0d5ec amalgkit | tar -x -C "$baseline_dir"
python benchmarks/benchmark_recent_changes.py --baseline "$baseline_dir" \
    --cases wsfilter csfilter imputation tau tmm_rounds tpm_transform \
    --genes 15000 --samples 120 --orthologs 3500 --repeats 3 \
    --output "$results_dir"
```

See the [quantification benchmarks](README.md) for the separate measurements
on real sequencing inputs and known spike-in concentrations.

## MMseqs FASTQ staging

On the same machine, three measured runs after warmup compared
`write_mmseqs_rrna_query_chunk` before and after `27cea1a`. The measurements
include opening, reading, and writing the query and replay chunks, but exclude
MMseqs search and final FASTQ rewriting. Decompressed records and their order
matched. MMseqs 18.8cc5c also produced identical hits and retained reads in a
350-spot smoke test with the SEQC ERCC reference, in both layouts.

| Input prefix | Wall time (s) | Temporary files (MB) |
| --- | ---: | ---: |
| 60,000 SG-NEx HEYA8 direct RNA reads, replicate 1 run 1 | 26.85 → 2.12 | 71.25 → 37.27 |
| 100,000 [SEQC BGI ERCC pairs](https://www.ebi.ac.uk/ena/browser/view/SRR896989) | 16.04 → 6.76 | 25.11 → 26.41 |
| 100,000 synthetic pairs | 14.04 → 8.56 | 44.39 → 44.71 |

The single-end query replaces a separate copy of the same reads and uses
gzip level 3. Paired chunks still need separate replay files; gzip level 6
avoids most of the temporary disk growth observed with level 3. Final
FASTQ compression is unchanged.

## CSTMM imputation sensitivity

We also compared TMM factors before and after masking counts simulated from
the [Fukushima et al. ortholog expression profiles](https://doi.org/10.17632/3vcstwdbrn.3).
We used 1,377 orthogroups from 21 species, with either liver alone or brain,
heart, liver, and testis. Within each sample, the relative weights were
proportional to `max(2**expression - 1, 0)`. For six seeds (709–714), we drew
library sizes from a lognormal distribution (median 1 million, log SD 0.6).
We drew the fraction allocated to orthogroups uniformly from 2.5% to 12%,
then drew Poisson counts using that fraction and the relative weights.
We hid 20% of cells independently or hid each
orthogroup within a species with probability 20%. With one sample per
species, both masks were identical.

The table reports median RMSE on log2 round-two TMM factors against the
unmasked factors. EM-PCA used library-size scaling and 50 iterations.

| Input and missingness | Rank 4 | Rank 5 | Rank 8 |
| --- | ---: | ---: | ---: |
| Liver, either mask | 0.074 | 0.067 | 0.065 |
| Four organs, independent cells | 0.079 | 0.085 | 0.081 |
| Four organs, species blocks | 0.083 | 0.082 | 0.080 |

The higher ranks did not consistently improve factor recovery. None of
these deliberately masked simulations met the default convergence tolerance,
even for the tested rank 4 and 5 runs at 150 iterations. The comparison
explicitly inspected intermediate unconverged estimates; the standard
CSTMM failure policy would reject them. Neither the rank nor the
convergence policy was changed.
