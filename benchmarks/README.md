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
processes against two source trees. It records wall time, whole-process peak RSS
(including native allocations), repeated measurements, and numerical output
equivalence. CSTMM also checks hashes of the exported corrected count tables.
Input generation and imports are outside the timed interval but included in RSS.

```bash
mkdir -p /tmp/amalgkit-baseline
git archive 7ac5250 | tar -x -C /tmp/amalgkit-baseline
python benchmarks/benchmark_recent_changes.py \
  --baseline /tmp/amalgkit-baseline --output /tmp/amalgkit-comparison
```

Use one Python environment with both revisions' dependencies installed. Do not
run benchmark jobs concurrently. The runner fixes native thread counts to one,
alternates revision order, excludes the first pair as warmup, and measures three
pairs by default. Inspect `summary.json` equivalence fields before interpreting
timings; scientific definition changes need carefully matched fixtures.

The [September 2026 comparison](results/2026-09-10/README.md) documents workload
limits, reproduction commands, results, and the sources of observed increases.
