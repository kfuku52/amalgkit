# Changelog

All notable user-visible and maintainer-facing changes are recorded here. GitHub
Releases provide generated notes only for tagged releases. Patch-only updates
remain on the default branch and are recorded below; consult this file rather
than the Releases page for those changes.

## Unreleased

### 0.16.78 - 2026-09-03

- Recorded actual input spot and base counts for private FASTQ runs so
  `getfastq_stats.tsv` no longer reports `num_written=0` for valid outputs.
  Paired private inputs now fail early when mate record counts differ.
- Kept rRNA search sensitivity on MMseqs2 `easy-search` while omitting it from
  nucleotide `createindex`, avoiding invalid k-mer threshold derivation for
  the nucleotide index without patching MMseqs2.

### 0.16.76 - 2026-08-31

- Verify public original FASTQ downloads against ENA's per-file MD5 and byte
  count, then read every gzip member and trailer before publishing or trimming.
  A corrupt transfer is retried once in private staging; persistent corruption
  remains a failure and cannot become a completed run. Legacy sources without
  published checksums still require a complete non-empty gzip stream.

### 0.16.75 - 2026-08-31

- Made `cstmm --metadata` effective, validated species/run alignment, and excluded
  unselected runs from factor estimation and corrected count tables. Explicit
  metadata must retain a row for every input count column; use exclusion flags
  to remove runs. Previously ignored metadata overrides can now change results
  or produce an actionable error. Recompute CSTMM and downstream outputs from
  the matching merge inputs when adopting corrected selection.
- Rejected FPKM from CSTMM counts without valid original `tmm_library_size`
  values, preventing accidental cancellation of TMM when old metadata is used.
- Migrated bundled yeast rules to the current parameter schema, preserved
  reviewed biological group assignments instead of applying plant rules, and
  limited tutorial sampling to three runs per species/group. Standard text
  normalization still applies to group labels.
- Fixed generated species-workspace instructions and CLI examples to pass the
  correct metadata directories and rule file without replacing yeast rules.
- Corrected private FASTQ metadata handoffs, reference/index prerequisites,
  metadata inference, Oarfish abundance scales and compatible normalization,
  run/species batch units, external dependencies, output names, and defaults.
  Made the tutorial metadata edit non-destructive and preserved lexical IDs
  and annotations. Shortened the README and restored a current workflow diagram.
- Added network-free documentation syntax/default/link checks to existing CI,
  executable small-fixture workflow examples, and installed-wheel yeast selection.
- Added a guarded Wiki synchronization tool with version/source-commit footers
  and documented the separate publication step. Aligned security-update and
  changelog guidance with the default-branch patch release policy.

### 0.16.74 - 2026-08-31

- Preserved previous output directories after failed commits, failed rollback
  and interrupts; recovery paths are reported without deleting the last copy.
- Redacted URL credentials in structured log fields, exception chains, timeout
  diagnostics and debug tracebacks while keeping raw subprocess data usable.
- Rejected empty, non-finite and negative quant tables before merge publication.
  Preserved lexical gene/run identifiers through quant, orthology, normalization
  and finalization, including leading zeroes and literal `NA` identifiers.
- Allowed default quant cleanup of managed private-gzip symlinks without
  deleting their source files, and made marker/input rollback interruption-safe.
- Added mate-conversion settings to getfastq resume schema 3. Older state files
  require reprocessing because they cannot prove the setting; `.sra` downloads
  are retained.
- Recognized the getfastq completion manifest during sanity checks instead of
  reporting it as an orphan run output.
- Stabilized TMM trimming ranks and reference ties, verified with independent
  60-digit and edgeR fixtures. Boundary cases can differ from older results;
  regenerate CSTMM and downstream outputs together to adopt the new convention.
- Fixed integer DataFrame correction under pandas 3 in TMM and Combat-seq
  without truncating fractional values or modifying input matrices.
- Bounded duplicate-validation memory with a temporary SQLite inventory and
  retained one canonical target vector while merging quant tables. Added
  reproducible wall-time, allocation/RSS and output-hash benchmarks.
- Extracted resume, FASTQ cleanup, table I/O, identifier validation and redaction
  boundaries with shared formatting/type checks. Reused cached CI installs,
  removed the duplicate Python 3.14 fast lane, gated macOS by relevant paths,
  and extended real-tool E2E through default cleanup, TMM and finalization.
- Verified built wheels outside the checkout by importing their installed
  modules and extracting the bundled FASTA dataset.

### Earlier default-branch updates

- Allowed `quant` to reuse a canonical, existing `getfastq` directory through
  a symbolic-link input root while retaining path-component containment and
  symbolic-link rejection for run directories and output roots.

- Added `tmm_effective_library_size` to `cstmm/metadata.tsv` as the product of
  `tmm_library_size` and `tmm_normalization_factor`, and documented safe edgeR
  handoff without recomputing or double-applying TMM normalization.

- Defined Oarfish abundance as unlength-normalized counts per million, stored annotated length separately from a unit effective-length term, propagated that model through merge and CSTMM outputs, and rejected FPKM for those runs because the transformation is undefined without an effective length.

- Preserved undefined within-species sample correlations instead of fabricating r=0 (including a zero self-correlation) in finalize/overview embeddings.

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
