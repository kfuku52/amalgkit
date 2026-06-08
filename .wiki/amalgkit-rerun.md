## Overview

`amalgkit rerun` reads `sanity_report.json`, resolves failed pipeline targets, and reruns the selected commands. It also writes a manifest so the rerun plan can be inspected or archived.

Run `amalgkit sanity` first:

```bash
amalgkit sanity --out_dir ./ --all
```

Then preview rerun targets:

```bash
amalgkit rerun --out_dir ./ --dry_run yes
```

## Default paths

```text
Report:   out_dir/sanity/sanity_report.json
Manifest: out_dir/sanity/rerun_manifest.json
```

`--metadata report` uses the metadata path stored in the sanity report. This is the default and is usually the safest choice.

## Examples

Rerun all failed targets recorded by the report:

```bash
amalgkit rerun --out_dir ./
```

Rerun only selected target families:

```bash
amalgkit rerun --out_dir ./ --check getfastq,quant
```

Limit by run:

```bash
amalgkit rerun --out_dir ./ --check getfastq,quant --run SRR000001,SRR000002
```

Limit by species:

```bash
amalgkit rerun --out_dir ./ --check merge,busco,finalize --species "Homo sapiens"
```

Include warning-level issues, such as missing global summary PDFs:

```bash
amalgkit rerun --out_dir ./ --include_warnings yes
```

Disable manifest output:

```bash
amalgkit rerun --out_dir ./ --manifest none
```

## Supported targets

- `getfastq`
- `index`
- `quant`
- `merge`
- `busco`
- `finalize`
- `all`

If `--check` is omitted, rerun uses the checks recorded in `sanity_report.json`.

## Notes

- `--redo yes` is the default for rerun so broken outputs can be overwritten.
- `--download_dir` and `--download_lock_dir` are shared with download-heavy targets.
- Use `--dry_run yes` before large cluster reruns.
