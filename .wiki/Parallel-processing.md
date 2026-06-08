## Overview

AMALGKIT uses a shared CPU-budget model. In most runs, set only `--threads` and leave the advanced worker options on `auto`.

`--threads` is the total CPU budget for one AMALGKIT process. It is not a per-worker thread count. Commands that can parallelize divide that budget across internal workers and per-worker work.

## Key Options

| Option | Meaning |
| --- | --- |
| `--threads INT or auto` | total CPU budget for one AMALGKIT process |
| `--internal_jobs INT or auto` | advanced override for internal worker count |
| `--internal_cpu_budget INT or auto` | advanced cap used when automatic worker counts are chosen |
| `--batch INT` | one-based metadata row selector for array-job execution |

When `--batch` is set, AMALGKIT processes one metadata row and forces `--internal_jobs` to `1`. Total CPU demand then scales with the number of concurrent array tasks.

## Practical Defaults

- Use `--threads auto` for local exploratory runs.
- Use `--threads N` on a shared node or scheduler allocation.
- Keep `--internal_jobs auto` unless you need to limit AMALGKIT's internal concurrency.
- Keep `--internal_cpu_budget auto` unless a cluster policy requires a hard cap below the visible CPU count.

## Single-Process Example

```bash
amalgkit getfastq --out_dir ./ --threads 8
amalgkit quant --out_dir ./ --threads 8
amalgkit merge --out_dir ./ --threads 8
amalgkit wsfilter --out_dir ./ --threads 8
amalgkit finalize --out_dir ./ --threads 8
```

## Array-Job Example

`getfastq`, `quant`, `wsfilter`, `csfilter`, and `finalize` accept `--batch`.

```bash
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --array=1-20

amalgkit quant \
    --out_dir ./ \
    --threads "$SLURM_CPUS_PER_TASK" \
    --batch "$SLURM_ARRAY_TASK_ID"
```

If 20 array tasks run at the same time with `--threads 4`, the job may use up to about 80 cores across the cluster.

## Shared Downloads and Locks

Large projects often start many AMALGKIT processes against the same workspace. Use shared cache and lock paths so downloads and database/index builds are not duplicated.

| Option | Default | Role |
| --- | --- | --- |
| `--download_dir` | `out_dir/downloads` | shared download/cache directory |
| `--download_lock_dir` | `download_dir/locks` | filesystem locks and provider semaphores |
| `--ncbi_download_max_concurrency` | `auto` | optional cap for NCBI cloud-object downloads |
| `--aws_download_max_concurrency` | `auto` | optional cap for AWS downloads |
| `--gcp_download_max_concurrency` | `auto` | optional cap for GCP downloads |
| `--ena_download_max_concurrency` | `auto` | optional cap for ENA downloads |
| `--ddbj_download_max_concurrency` | `auto` | optional cap for DDBJ downloads |

For `quant`, shared index-build locks are created under the selected `--index_dir`. This prevents multiple array tasks from building the same species/backend index at the same time.

## Command Families

- `metadata` can parallelize species-wise batch queries and throttle NCBI metadata requests through `--ncbi_metadata_max_concurrency`.
- `getfastq` parallelizes run processing and uses provider fallbacks plus shared download locks.
- `quant` parallelizes run processing and uses shared index locks.
- `merge`, `wsfilter`, `csfilter`, and `finalize` parallelize internal work units while respecting the same CPU budget.
- `cstmm` is usually lighter than download or quantification steps, but still follows its declared resource options.
