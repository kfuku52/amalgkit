# Benchmarks

`benchmark_core.py` tracks the kernels most likely to affect interactive and
batch throughput: streaming quant-output validation, cross-species correlation,
and getfastq directory inventory construction.

Run a short local sample with:

```bash
python benchmarks/benchmark_core.py --quick --output benchmark-results.json
```

The scheduled end-to-end workflow runs the full sizes and uploads the JSON
result. Compare medians across several runs rather than treating one runner as
a hard absolute baseline; hosted-runner hardware varies. Peak Python allocation
is reported alongside time so an apparent speedup cannot silently replace a
bounded-memory algorithm with a full-file load.

## Statistical filter evaluation

`evaluate_filter_accuracy.py` compares optional wsfilter/csfilter policies using
labelled synthetic interventions. It measures detection and biological retention,
not runtime. Run `PYTHONPATH=. python benchmarks/evaluate_filter_accuracy.py
--seeds 12 --output filter-accuracy.json` with a supported Python environment.
See [filter validation](../.wiki/Filter-validation.md) for interpretation and
why these simulations alone do not justify changing scientific defaults.

## Comparing recent revisions

`benchmark_recent_changes.py` runs the same deterministic workloads in isolated
processes against two source trees. It records wall time, process CPU time,
whole-process peak RSS (including native allocations), repeated measurements,
exact output equality and numerical equivalence. CSTMM also checks hashes of
the exported corrected count tables.
Input generation and imports are outside the timed interval but included in RSS.

```bash
mkdir -p /tmp/amalgkit-baseline
git archive 1cdbd7322a32ecdbe68e07ce803082bc8991b2bb | tar -x -C /tmp/amalgkit-baseline
python benchmarks/benchmark_recent_changes.py \
  --baseline /tmp/amalgkit-baseline --output /tmp/amalgkit-comparison
```

Use one Python environment with both revisions' dependencies installed. Do not
run benchmark jobs concurrently. The runner fixes native thread counts to one,
alternates revision order, excludes the first pair as warmup, and measures three
pairs by default. Inspect `summary.json` equivalence fields before interpreting
timings; scientific definition changes need carefully matched fixtures.

Additional cases include `select_normalization` (repeated plant annotations)
and `gsa_selection` (sparse pending rules in a wide metadata table). String
outputs are compared exactly; numeric outputs also report the maximum absolute
difference. CPU time and recorded load averages help identify scheduling noise,
but throughput claims still require comparable, otherwise idle runs.

The [September 2026 comparison](results/2026-09-10/README.md) documents workload
limits, reproduction commands, results, and the sources of observed increases.

## October 2026 efficiency audit

The audit inspected CLI startup, metadata/taxonomy, selection, FASTQ retrieval
and sampling, quantification/reuse, merge/validation, normalization, batch models,
filters/plots, and test/CI delivery. Profiles identified repeated work in four
kernels; the comparison above uses their pre-change revision as its baseline.

| Kernel | Work removed | Compatibility checks |
| --- | --- | --- |
| wsfilter | repeated numeric column conversion; integer products for all-finite support counts | pairwise correlations/counts, all methods, project exclusion, zero-margin rounding |
| csfilter | repeated metadata scans and sample finite masks | exact correlations/counts, run/species exclusion, lexical groups, missing samples |
| selection normalization | regex scans of resolved rows; repeated organ-text validation | exact annotations, mixed tissues, changed rules, row order |
| deferred GSA selection | row construction/status scans for records without recipes | exact statuses/exclusions, incomplete groups, malformed recipes, array identifiers |

```bash
python benchmarks/benchmark_recent_changes.py \
  --baseline /tmp/amalgkit-baseline \
  --cases wsfilter csfilter select_normalization gsa_selection \
  --genes 20000 --samples 200 --orthologs 2000 --repeats 5 \
  --output /tmp/amalgkit-efficiency
```

Run the core/table benchmarks from CONTRIBUTING.md as additional checks.
Streaming validation and single-pass quant reads already bound table memory;
FASTQ integrity scans and cached content fingerprints protect reuse. These
contracts remain intact. Scientific algorithms/defaults and CI coverage were
retained. Network transfer, external-tool throughput, optional Combat-seq and
platform/dependency CI lanes require their own environments; static inspection
does not establish their performance. Keep generated measurements outside Git.

## Additional efficiency kernels

`benchmark_efficiency.py` compares bounded-memory QC summaries, default tau
aggregation, measured GSA overlays, shared correlation preparation, random
sampling positions, SVA right vectors, classical MDS, all-pair CSTMM diagnostics,
and a planted-factor SVA correction. Use an archived source tree as the baseline:

```bash
mkdir -p /tmp/amalgkit-efficiency-baseline
git archive 6539068259cfd3cb9bbb915fe06a28e88da1d581 | tar -x -C /tmp/amalgkit-efficiency-baseline
python benchmarks/benchmark_efficiency.py \
  --baseline /tmp/amalgkit-efficiency-baseline \
  --output /tmp/amalgkit-efficiency-comparison
```

Each case uses seed 719, one warmup and three repeats in separate processes
with native thread counts fixed to one. JSON records workload sizes, wall/CPU
medians, a separate Python allocation scan, whole-process RSS and load averages.
Input setup/imports are outside timing; RSS includes setup and repeats. Run
without concurrent test/benchmark jobs. System load can still distort wall time;
report that limitation rather than inferring throughput from CPU time alone.

Outputs are compared exactly for selection positions, integer QC summaries,
tau summaries, GSA overlays, correlation statistics and diagnostic TSVs. SVA
subspaces and MDS geometry are compared through their projection/Gram matrices;
the SVA pipeline compares corrected expression with relative tolerance `1e-8`
and absolute tolerance `1e-10`. Contract tests additionally cover deficient
ranks, tiny/large scales, tied MDS spectra, lexical IDs and fractional QC counts.
The latter can differ at rounding precision after chunked summation. Sampling
benchmarks time position preparation, not FASTQ I/O; correlation benchmarks
exclude PDF rendering. All-pair diagnostics still scale quadratically in the
sample count; preserve their full output when comparing alternatives.

## Validating the October candidate implementations

`benchmark_recent_changes.py` also supports `tmm_rounds`, `averaged_inputs`,
`select_exclusion`, `select_exclusion_unique`, and `provenance_hash`. These use
the production functions in both source trees, including the actual base
selection rules, TMM reference choices and averaged-table labels. Comparisons
exit nonzero for unequal output or corrected-count hashes; the warmup is also
checked. The `--allocation-scan` option measures Python allocations in a second
call after timing and process peak RSS, so tracing overhead does not enter the
reported wall/CPU interval.

```bash
python benchmarks/benchmark_recent_changes.py \
  --baseline /tmp/amalgkit-baseline --repeats 5 --allocation-scan \
  --cases tmm_rounds averaged_inputs select_exclusion select_exclusion_unique provenance_hash \
  --genes 20000 --samples 90 --output /tmp/amalgkit-candidates
```

Repeated-text selection is expected to benefit most. The bounded cardinality
probes are a heuristic: inputs engineered to fool both probes still incur one
full factorization before rejection. Hashing retains SHA-256 and the same
mutation/cache checks; its main expected benefit is bounded buffer allocation,
not faster disk I/O. Keep host-load effects and whole-pipeline contribution
separate from these kernel measurements.

QC streaming bounds memory rather than promising faster parsing. Canonical TSV
inputs need a bounded row-width pass before the C parser, which can otherwise
ignore extra fields at chunk boundaries. Other formats use the strict Python
parser. Report that CPU/I/O cost alongside memory savings; remove the extra
pass only after malformed-boundary regressions pass with the underlying parser.
