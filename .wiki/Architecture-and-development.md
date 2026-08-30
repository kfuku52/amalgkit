# Architecture and development

AMALGKIT keeps command modules as stable public facades and moves reusable work
into focused boundaries. The main boundaries are:

- `output_contracts`: streaming quant/BUSCO validation
- `parallel_utils`: deterministic concurrency and preserved task tracebacks
- `subprocess_utils`: timeouts and cached dependency probes
- `cross_species_computation`: numerical kernels
- `getfastq_file_state` and `getfastq_probes`: file inventory and probe parsing
- `logging_utils`: optional structured JSONL events
- `redaction`: credential-safe diagnostic messages and exception chains
- `table_io` and `identifier_validation`: lexical IDs and bounded exact duplicate checks
- `getfastq_resume`: fingerprints and state/file identities
- `fastq_cleanup`: reversible cleanup of managed files and private links

The full maintainer map and migration rules are in
[`ARCHITECTURE.md`](https://github.com/kfuku52/amalgkit/blob/master/ARCHITECTURE.md).

## Diagnostic logging

Normal terminal output remains concise. To record structured events, put the
global options before the command:

```bash
amalgkit --log_file ./logs/quant.jsonl quant --out_dir ./work
amalgkit --log_level DEBUG quant --out_dir ./work
```

Each line is JSON and can include command, run/species context, duration,
return code, or an exception traceback.

URL userinfo, query strings and fragments are removed from JSON fields,
exception chains, timeout stderr and debug traceback output. Raw subprocess
results are retained for programmatic consumers; only displayed diagnostics are
redacted. Avoid putting credentials in URL paths, which remain visible for
diagnosis.

## Fast feedback

```bash
python -m pytest -q -n 2 -m "not integration and not slow and not optional_dependency"
python .github/scripts/check_quality.py
python benchmarks/benchmark_core.py --quick
```

CI additionally runs the complete suite with branch coverage, a dependency-floor
environment, distribution checks, and a scheduled real-tool end-to-end workflow.
