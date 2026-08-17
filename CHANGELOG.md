# Changelog

All notable user-visible and maintainer-facing changes are recorded here. GitHub
Releases provide the commit-level generated notes for each version.

## Unreleased

- Estimated CSTMM TMM factors from pairwise-observed single-copy counts instead of EM-PCA–imputed ortholog holes.

- Rejected ambiguous metadata merges with empty or duplicate run IDs and
  prevented missing, duplicate, or whitespace-padded quantification target IDs
  from corrupting merged abundance tables.
- Made SVA input scales explicit: raw and fractional estimated counts are
  corrected on the log(count + 1) scale, pre-transformed inputs retain valid
  negative values, and collinear surrogate variables cannot remove protected
  sample-group effects.
- Preserved undefined cross-species correlations without fabricating zeroes,
  retained valid absent/multicopy ortholog cells, and aligned finite-only
  embeddings and dendrogram labels.
- Hardened SRA/ENA download handling by validating each HTTPS redirect before
  connecting, using exact endpoint allowlists and private temporary directories,
  and redacting credentials from URL-like command arguments.
- Restored major/minor-only release automation: patch-only versions no longer
  create Git tags or GitHub Releases and therefore do not trigger Bioconda
  autobump pull requests.

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
