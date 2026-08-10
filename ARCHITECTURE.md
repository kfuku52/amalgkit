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
| `parallel_utils.py` | CPU allocation, deterministic task results, traceback-preserving failures |
| `subprocess_utils.py` | command execution, timeouts, output decoding, successful probe cache |
| `logging_utils.py` | namespaced JSONL event logging |
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

Performance-sensitive changes should update or run `benchmarks/benchmark_core.py`.
The scheduled end-to-end workflow retains benchmark JSON and exercises real
`sra-tools`, `seqkit`, and `kallisto` binaries.
