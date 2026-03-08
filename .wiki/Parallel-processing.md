## Overview

AMALGKIT uses a shared CPU-budget model across commands. In most cases, you only need to set `--threads`.

## Key options

| Option | Meaning |
| --- | --- |
| `--threads` | total CPU budget for one AMALGKIT process |
| `--internal_jobs` | override the number of internal workers derived from `--threads` |
| `--internal_cpu_budget` | cap the CPU budget used for automatic worker selection |
| `--batch` | one-based metadata row selector for array-job style execution |

`--threads` is a global budget, not a per-worker thread count. AMALGKIT splits it internally across workers and per-worker CPU use.

When `--batch` is set, `--internal_jobs` is forced to `1`.

## Practical guidance

- Start with `--threads auto` unless you are on a shared node.
- Use `--internal_jobs auto` unless you need to cap concurrent workers manually.
- Prefer `--threads N` over trying to tune every command separately.

## Example: single process with a fixed CPU budget

```bash
amalgkit getfastq --out_dir ./ --threads 8
amalgkit quant --out_dir ./ --threads 8
amalgkit wsfilter --out_dir ./ --threads 8
```

## Example: array jobs

```bash
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --array=1-20

amalgkit quant \
    --out_dir ./ \
    --threads "$SLURM_CPUS_PER_TASK" \
    --batch "$SLURM_ARRAY_TASK_ID"
```

## Notes by command family

- `getfastq` and `quant` parallelize over runs.
- `merge`, `wsfilter`, `csfilter`, and `finalize` parallelize internally over larger work units while respecting the same CPU budget.
- `cstmm` is usually lightweight compared with `getfastq` and `quant`, but still follows the same budget model where applicable.
