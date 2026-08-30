# Architecture

AMALGKIT is a Python CLI whose commands exchange files in a workspace. The
public command modules remain stable import surfaces; reusable kernels and file
contracts live in smaller modules so they can be tested without initializing a
whole workflow.

## Runtime flow

```text
cli_entry / __main__
        |
        v
main.py -> cli_parser.py -> lazy command handler
                              |
                              v
       metadata/select/getfastq/quant/merge/... command module
                              |
                 +------------+------------+
                 |                         |
          shared contracts            external tools
     metadata, paths, outputs,      subprocess timeouts,
       parallel execution              cached probes
```

`main.py` owns process-level error handling and optional JSONL logging.
`cli_parser.py` owns argument compatibility. Command modules own orchestration,
not low-level reusable algorithms.

## Stable boundaries

| Module | Responsibility |
| --- | --- |
| `output_contracts.py` | streaming validation and required schemas for quant/BUSCO files |
| `table_io.py` | lexical identifier/annotation reads, separate from numeric missing values |
| `identifier_validation.py` | exact duplicate detection with bounded RAM and a temporary SQLite inventory |
| `fastq_cleanup.py` | reversible retirement of managed FASTQs and private-input links |
| `getfastq_resume.py` | versioned fingerprints, atomic state serialization and file identities |
| `parallel_utils.py` | CPU allocation, deterministic task results, traceback-preserving failures |
| `subprocess_utils.py` | command execution, timeouts, output decoding, successful probe cache |
| `logging_utils.py` | namespaced JSONL event logging |
| `redaction.py` | URL credential redaction at diagnostic output boundaries |
| `cross_species_computation.py` | correlation, imputation-cache, and embedding kernels |
| `getfastq_file_state.py` | bounded directory inventories for a getfastq run |
| `getfastq_probes.py` | pure fasterq-dump output/version parsing |
| `linalg_utils.py` | shared numerical-stability decisions |
| `text_utils.py` | first-seen-order text normalization |

Historical functions imported directly from command modules are re-exported or
wrapped there. New code should import the stable boundary directly; compatibility
facades can be removed only in a documented breaking release.

## Data and failure contracts

Large tables must be validated or transformed in chunks unless the algorithm
requires the full matrix. `output_contracts.py` is the source of truth for
cross-command schemas. Avoid adding a second validator to a producer or to
`sanity`.

Quant consumers validate their already-loaded columns through the same frame
contract as the streaming validator. Identifier columns are read as lexical
tokens before type inference: `0001`, `1`, and `NA` are different identifiers.
Only empty annotation cells represent missing IDs. Numeric columns retain their
own NA/finite/nonnegative checks. Merge retains one canonical target-ID vector
and at most the worker count of transient input vectors.

The streaming validator retains at most 50,000 prior IDs in Python. Larger
inventories spill to a private SQLite file with a 2 MiB page cache, preserving
exact cross-chunk duplicate detection. This scratch file is removed on success
or error; enough temporary disk space is required for the ID inventory.

Output directory commits keep the previous directory if rollback fails and
attach its recovery path to the original exception. Cleanup must never delete
an un-restored backup. Private FASTQ links are retired as directory entries;
their original source files are not owned or removed by AMALGKIT.

Parallel callers receive failures in input order. When a workflow must stop,
use `raise_task_failures`: the outer exception retains the established concise
message and type, while its `ExceptionGroup` cause holds each original traceback.

External-tool discovery in the startup banner is PATH-only. Commands perform
the authoritative version/help probe with a 30-second default and reuse a
successful probe within the process.

## Incremental modernization

The current codebase predates uniform formatting and typing. CI therefore
enforces Ruff formatting and mypy on migrated boundary modules, while Ruff lint
and the lowered complexity ceiling apply repository-wide. When extracting or
substantially editing another module, add it to the formatted/typed boundary
only after its checks pass; this avoids a repository-wide mechanical diff.

The shared list is `[tool.mypy].files` in `pyproject.toml`; run
`python .github/scripts/check_quality.py` for both formatting and typing checks.
Getfastq orchestration still owns filter execution and restart decisions, while
state serialization and fingerprint changes can be reviewed in isolation.

Performance-sensitive changes should update or run `benchmarks/benchmark_core.py`.
The scheduled end-to-end workflow retains benchmark JSON and exercises real
`sra-tools`, `seqkit`, and `kallisto` binaries.
