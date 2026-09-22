# Test pruning review — 2026-09-22

Reviewed the suite by the realistic failure each test can detect, including
its setup, assertions, mocks, and overlap with other tests. Test counts and
coverage were not retention criteria. The normal delivery coverage check is
still run without changing its threshold.

## Removed or consolidated

| Area | Decision and remaining failure detection |
| --- | --- |
| CLI | Remove help-text inventories and duplicate legacy-command rejection. Keep real process entrypoints, topic routing, root help, backend-specific index guidance, and dataset output. Test scientific choices and timeout help in-process. Numeric validators retain their boundary cases on one consumer plus a valid/invalid wiring case on every other consumer. |
| FASTQ extraction | Fold range arguments and avoiding redundant trimming into the existing full/partial-range statistics tests. Drop three dependency-probe variants whose only assertion was the same first executable. Keep tool failures, unsupported versions, paired counts, fallback trimming, filter ordering, and resume regressions. |
| MMseqs | Check configured sensitivity, maximum hits, and memory together at the external command boundary. Keep the separate automatic-option path and fatal-child-error test. Use a non-default memory value so the override check cannot pass by ignoring the option. |
| Text, metadata, selection | Combine blanks, NaN, duplicates, whitespace, and scalar conversion in mixed-input cases. Fold generator and identifier-delimiter checks into existing tests. Remove type-only pivot tests, trivial factory/reorder smoke tests, schema-position assertions, and an encoding-keyword spy. Keep copying, literal identifiers, annotation preservation, invalid metadata, and public selection workflows. |
| Species files and BUSCO | Combine ordinary files, hidden/temporary directories, and multiple species in directory tests. Fold basic BUSCO normalization into the comment/header case and share compressed-input setup. Keep ambiguous matches, missing files, unsafe paths, failed publication, and real plotting. |
| Scientific kernels | Remove shape/nonnegative-only RUV/TMM smoke checks, the weak Leek range assertion, default-dictionary snapshots, and duplicated zero-factor/constant-input tests. Keep independent R/Decimal/hand-calculated references, protected biological signals, labels, failed fits, zero correction, scale restoration, and missing-data behavior. |
| Linear algebra | Remove the internal fallback predicate's branch-by-branch tests, including a copied epsilon formula. Keep actual degenerate-boundary reconstruction tests in imputation and RUV and zero-rank SVA tests. The precise internal fallback decision is not itself a scientific output contract. |
| Imputation | Remove a lower-bound test whose imputed mean never reached the floor; the CSTMM negative-value regression exercises actual clipping. Fold observed-value preservation into the successful iterative-fit test. |
| Correlation | Combine NaN and infinity filtering; remove the direct constant-block case already exercised by the finite-correlation resolver. Retain undefined-pair rejection, extreme scales, and all three correlation methods. |
| Fragment distributions | Replace the integer/default/fractional Cartesian product with non-default integer metadata and fractional file inputs. Keep precedence, invalid inputs, provenance, backend/layout changes, reuse safety, and real kallisto comparisons. |
| Test infrastructure | Remove assertions that rechecked bytes produced by the test-only PDF stub on every emission. Real PDF renderers retain their dedicated tests. |

The removed cross-species exclusion test performed `.eq('no')` in the test
itself: it could detect a pandas filtering failure but not a regression in the
production filtering caller. Exclusion-reason preservation and workflow
filtering tests remain. The removed help-handler mock test likewise repeated
routing already exercised through the real CLI.

## Retention decisions across the rest of the suite

- Download, GSA, SRA, taxonomy, runtime paths, and subprocess tests guard actual
  transfer deadlines, retry limits, hostile paths/XML/URLs, checksum failures,
  truncation, paired-read identity, locks, secret redaction, and cache freshness.
  These are distinct failure modes even when setup looks similar.
- Merge, quantification, integration, sanity, rerun, and output-contract tests
  protect cross-file handoffs, lexical IDs, stale outputs, atomic replacement,
  rollback, and crash recovery. Keep both component and workflow checks where
  the workflow can omit an otherwise correctly tested validator.
- Selection rules, sampling, within/cross-species filtering, normalization,
  finalization, and batch correction retain biologically distinct cases and
  independent reference data. An external tool mock is retained when it checks
  a real command/provenance boundary rather than merely its own return value.
- Cache, directory-scan, worker-budget, and parallel scheduling tests remain
  where they prevent repeated large I/O, retaining matrices, oversubscribing
  cores, losing task results, or corrupting concurrent publications. Internal
  calls are justified here by the concrete resource or correctness failure.
- Documentation workflows, documentation tooling, datasets, dependency-floor
  consistency, logging, and artifact export tests guard installed/public
  behavior or independently maintained files. They are not substitutes for
  testing Python, pandas, argparse, or other libraries themselves.

Deliberately no longer frozen: exhaustive help wording, incidental metadata
column positions, internal default dictionaries, and the exact numerical
fallback predicate. Restoring such snapshots would add maintenance cost
without demonstrating a new user-visible failure.
