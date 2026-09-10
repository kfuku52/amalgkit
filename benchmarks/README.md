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
