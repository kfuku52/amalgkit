# Changelog

All notable user-visible and maintainer-facing changes are recorded here. GitHub
Releases provide the commit-level generated notes for each version.

## Unreleased

## 0.16.59 - 2026-08-10

- Restored CI compatibility with pandas 3 string dtypes, current NumPy stubs,
  and the minimum Matplotlib dependency set.

## 0.16.58 - 2026-08-10

- Reduced startup dependency-probe latency and added successful probe caching.
- Unified bounded-memory output validation and shared BUSCO/quant schemas.
- Preserved original tracebacks for failures from parallel workflows.
- Added structured JSONL logging, typed computational boundaries, coverage and
  complexity gates, nightly real-tool end-to-end testing, and core benchmarks.
- Changed release automation to publish patch releases with wheel and source
  distributions.
