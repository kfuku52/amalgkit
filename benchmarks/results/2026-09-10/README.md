# September 2026 performance comparison

The September 8–10 changes increase some analysis costs substantially. The
largest measured runtime increase is CSTMM's default all-pairs reference audit;
tau aggregation also takes more time and memory. These are component timings,
not a claim about the elapsed time of an entire AMALGKIT pipeline.

## Revisions and environment

- Baseline: `7ac5250`, immediately before the September 8–10 change sequence.
- Current: `93d5f52`, including the log-scale-zero and stale batch-diagnostic fixes.
- Same local Python environment for both source trees: Python 3.11.15,
  NumPy 2.4.6, pandas 2.3.3, SciPy 1.17.1; macOS 26.6.2.
- Apple M2 Max, 64 GiB RAM; **x86_64 Python under Rosetta**, not native ARM Python.
- OMP, OpenBLAS, MKL and Accelerate thread limits all set to one. Benchmark
  jobs were run sequentially without concurrent tests or other agent benchmarks.
- Each workload/revision runs in a fresh child process. The first pair is a
  warmup; reported medians use three additional pairs with alternating order.
- Wall time excludes imports, input construction and result comparison. CSTMM
  includes actual output tables and plots. Peak RSS includes imports and input
  construction and native arrays; it is not incremental function memory.

## Measurements

The adjacent JSON files retain individual observations, workload parameters,
environment, ranges, output comparisons and corrected-count hashes.

### Standard workloads

| Component | Old seconds | Current seconds | Time ratio | Old peak MiB | Current peak MiB | RSS ratio |
|---|---:|---:|---:|---:|---:|---:|
| quant_validation | 0.357 | 0.391 | 1.10× | 162.4 | 175.6 | 1.08× |
| fastq_scan | 2.138 | 2.084 | 0.97× | 98.3 | 98.4 | 1.00× |
| wsfilter | 1.729 | 2.440 | 1.41× | 380.9 | 410.0 | 1.08× |
| csfilter | 1.175 | 1.582 | 1.35× | 96.4 | 97.1 | 1.01× |
| tau | 0.311 | 1.032 | 3.31× | 196.6 | 274.8 | 1.40× |
| sva | 2.486 | 2.958 | 1.19× | 227.6 | 242.5 | 1.07× |
| cstmm | 6.872 | 32.680 | 4.76× | 335.9 | 355.5 | 1.06× |

### Larger workloads (400 samples)

| Component | Old seconds | Current seconds | Time ratio | Old peak MiB | Current peak MiB | RSS ratio |
|---|---:|---:|---:|---:|---:|---:|
| wsfilter | 5.130 | 6.345 | 1.24× | 589.6 | 689.3 | 1.17× |
| csfilter | 1.648 | 2.330 | 1.41× | 108.4 | 107.2 | 0.99× |
| tau | 0.345 | 1.190 | 3.45× | 386.9 | 535.3 | 1.38× |
| cstmm | 15.574 | 93.711 | 6.02× | 602.5 | 654.5 | 1.09× |

### Disabled pair diagnostics and separate imputation

| Component | Old seconds | Current seconds | Time ratio | Old peak MiB | Current peak MiB | RSS ratio |
|---|---:|---:|---:|---:|---:|---:|
| cstmm | 5.849 | 5.710 | 0.98× | 335.5 | 332.5 | 0.99× |
| imputation | 0.089 | 0.090 | 1.01× | 157.1 | 163.8 | 1.04× |

CSTMM uses the standard 200-sample workload here. All measured comparisons
passed numeric equivalence; the ten corrected count files also match between
the baseline and both current diagnostic settings. Modest timing changes in
these short runs overlap measurement variation.

### New optional random sampling (current revision only)

| Input reads | Selected reads (10%) | Median seconds | Range seconds | Median peak MiB |
|---:|---:|---:|---:|---:|
| 1,000,000 | 100,000 | 4.529 | 4.491–4.667 | 110.8 |
| 5,000,000 | 500,000 | 22.761 | 22.336–22.905 | 201.9 |

The baseline has no equivalent random sampler, so these are absolute costs,
not regression ratios. Each size has one excluded warmup and three measured
fresh-process runs, with the same single-thread environment and 100-base
single-end gzip fixture. The harness checks the selected read and base counts.
All candidate records are validated and hashed even when only 10% are selected.
Partial Fisher–Yates selection keeps a swap map and selected-position set;
its memory grows with the ending selection rank (the selected count in this
first-round fixture). Larger selections and paired-end inputs need separate
capacity measurements. This feature is optional; the default remains contiguous.

## Why costs increased

1. **CSTMM:** [`write_observed_pair_diagnostics`](../../../amalgkit/cstmm_diagnostics.py) runs a TMM fit
   for every retained sample pair: 19,900 pairs at 200 samples and 79,800 at 400.
   Each fit processes the shared orthologs, including trimming/sorting. The TSV
   writer streams rows, but the plotter subsequently reads the entire pair table.
   This adds quadratic growth in sample count to runtime and diagnostic storage;
   plot memory also grows with the pair count. The existing
   `--tmm_reference_diagnostics no` option omits these diagnostic files/plots.
   It does not change the normalization factors. The disabled-diagnostic
   comparison above measures the remaining work explicitly.
2. **Tau:** [`linear_sample_group_summary`](../../../amalgkit/per_species_common.py) inverse-transforms the whole run
   matrix before averaging. It then builds per-unit means, concatenates them,
   and constructs coverage/weight records. This implements the new
   arithmetic-mean semantics and unit/project controls, but retaining whole-run
   arrays and general per-unit processing costs more than the previous
   transformed group-mean path. A future optimization could aggregate a group
   at a time and specialize the default run-unit case while preserving the new
   numerical and missing-value rules.
3. **Filters:** wsfilter now constructs finite-support counts with an int64
   matrix product, even when the minimum-common-genes gate is disabled.
   csfilter additionally counts finite pairs for every reference group.
   These diagnostic calculations explain additional work in otherwise matching
   default correlation scores; the measurements do not isolate each operation's
   individual contribution.
4. **SVA:** the comparison fixes three surrogate variables and five IRW
   iterations on both revisions. The new implementation uses supported SVD
   directions and additional design/output validation. Automatic variable
   selection and its changed fallback/iteration policy are outside this timing.

No scientific defaults or diagnostic coverage were changed during this
investigation. The correctness fixes are in `93d5f52`; the benchmark and report
form a separate commit.

## Workload and equivalence limits

- The deterministic seed is 719. Expression/filter workloads have ten tissue
  groups. CSTMM and csfilter have ten species.
- Standard expression matrices contain 20,000 genes and 200 samples; the larger
  matrices contain 30,000 genes and 400 samples. csfilter and imputation use
  2,000 ortholog rows; CSTMM uses those 2,000 rows to normalize all genes.
- Tau uses identical within-tissue replicate profiles so the old geometric
  and new arithmetic summary definitions yield equivalent tau values. This
  isolates computation cost; it does **not** imply general output equivalence
  between the old and new biological definitions.
- Correlation tests compare common numeric score columns, not new diagnostics
  absent from the baseline. csfilter has missing values every eleventh row and
  seventh column; wsfilter uses complete values and default run exclusion.
- CSTMM uses no missing values to isolate the reference audit. Its normalization
  factors are compared numerically and all ten species' corrected count files
  are compared by SHA-256. New diagnostic files have no baseline counterpart.
  The workload JSON's `observed_pair_diagnostics` field records the harness's
  requested setting; the baseline has no such audit regardless of that field.
- Imputation is a separate 20%-block-missing fixture with identical complete
  profiles and equal library depths, chosen so raw and library-scaled imputation
  are comparable. It is an easy low-rank case, not a worst-case convergence test.
- FASTQ scanning uses one million synthetic 100-base single-end reads in gzip,
  with low-entropy sequence/quality. Quant validation uses 100,000 transcript
  rows. Neither measures actual kallisto quantification, remote SRA/GSA transfer,
  fasterq-dump, or realistic compression entropy.
- Numeric equivalence uses `allclose(rtol=1e-8, atol=1e-10, equal_nan=True)`.
  This is not a full end-to-end real-data or cross-platform benchmark. Short
  subsecond timings and modest differences are more sensitive to host noise.

## Reproduction

Run from the repository with a supported Python environment containing both
revisions' dependencies. Keep the interpreter and package versions identical.

```bash
mkdir -p /tmp/amalgkit-perf-baseline
git archive 7ac5250 | tar -x -C /tmp/amalgkit-perf-baseline
python benchmarks/benchmark_recent_changes.py \
  --baseline /tmp/amalgkit-perf-baseline --output /tmp/amalgkit-perf-standard
python benchmarks/benchmark_recent_changes.py \
  --baseline /tmp/amalgkit-perf-baseline --cases wsfilter csfilter tau cstmm \
  --genes 30000 --samples 400 --orthologs 2000 --output /tmp/amalgkit-perf-large
python benchmarks/benchmark_recent_changes.py \
  --baseline /tmp/amalgkit-perf-baseline --cases cstmm imputation \
  --without-pairs --output /tmp/amalgkit-perf-no-pairs
for reads in 1000000 5000000; do
  for repeat in 0 1 2 3; do
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
      python benchmarks/benchmark_recent_changes.py \
      --repo "$PWD" --case random_sampling --reads "$reads" \
      --output "/tmp/amalgkit-perf-sampling/$reads/current-$repeat.json"
  done
done
```

Use `--current` to point to an archived `93d5f52` tree if testing from a later
checkout. Raw temporary numeric matrices and generated count/plot files are
not committed; rerunning the harness regenerates them.
