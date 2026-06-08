## Overview

`amalgkit rerun` reads `sanity_report.json`, resolves failed pipeline targets, and reruns selected commands.

It also writes a manifest so the rerun plan can be inspected or archived.

## Basic Flow

Run `sanity` first:

```bash
amalgkit sanity --out_dir ./ --all
```

Preview rerun targets:

```bash
amalgkit rerun --out_dir ./ --dry_run yes
```

Execute the rerun:

```bash
amalgkit rerun --out_dir ./
```

## Default Paths

```text
Report:   out_dir/sanity/sanity_report.json
Manifest: out_dir/sanity/rerun_manifest.json
```

`--metadata report` is the default. It uses the metadata path recorded in the sanity report and is usually the safest choice.

## Target Selection

If `--check` is omitted, `rerun` uses the checks recorded in `sanity_report.json`.

Supported targets:

- `getfastq`
- `index`
- `quant`
- `merge`
- `busco`
- `finalize`
- `all`

Rerun selected target families:

```bash
amalgkit rerun --out_dir ./ --check getfastq,quant
```

Limit by run:

```bash
amalgkit rerun \
    --out_dir ./ \
    --check getfastq,quant \
    --run SRR000001,SRR000002
```

Limit by species:

```bash
amalgkit rerun \
    --out_dir ./ \
    --check merge,busco,finalize \
    --species "Homo sapiens"
```

Include warning-level issues, such as missing global summary PDFs:

```bash
amalgkit rerun --out_dir ./ --include_warnings yes
```

## Manifest Control

Write the default manifest:

```bash
amalgkit rerun --out_dir ./ --manifest inferred
```

Write to an explicit path:

```bash
amalgkit rerun --out_dir ./ --manifest ./rerun_manifest.json
```

Disable manifest output:

```bash
amalgkit rerun --out_dir ./ --manifest none
```

## Operational Notes

- `--redo yes` is the default so broken outputs can be overwritten.
- `--download_dir` and `--download_lock_dir` are shared with download-heavy targets.
- Use `--dry_run yes` before large cluster reruns.
- Use `--include_warnings yes` when the goal is to regenerate warning-level missing summaries, not only failed core outputs.
