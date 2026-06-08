## Overview

`amalgkit sanity` scans a workspace and reports missing or ambiguous outputs for selected pipeline steps.

It writes `sanity_report.json`, which can be consumed by `amalgkit rerun`.

## Basic Use

Check all supported targets:

```bash
amalgkit sanity --out_dir ./ --all
```

Check selected targets:

```bash
amalgkit sanity --out_dir ./ --check getfastq,quant,merge
```

If no target flag or `--check` value is specified, `--all` is assumed.

## Targets

Supported check targets are:

- `getfastq`
- `index`
- `quant`
- `merge`
- `busco`
- `finalize`
- `all`

The older boolean flags are still available:

```bash
amalgkit sanity --out_dir ./ --getfastq --quant --merge
```

`--check` is more convenient for scripts:

```bash
amalgkit sanity --out_dir ./ --check getfastq,index,quant
```

## Limit Scope

Limit run-based checks:

```bash
amalgkit sanity \
    --out_dir ./ \
    --run SRR000001,SRR000002 \
    --check getfastq,quant
```

Limit species-based checks:

```bash
amalgkit sanity \
    --out_dir ./ \
    --species "Homo sapiens,Mus musculus" \
    --check merge,busco,finalize
```

Use explicit directories when outputs live outside the default workspace layout:

```bash
amalgkit sanity \
    --out_dir ./ \
    --check quant,merge \
    --quant_dir ./custom_quant \
    --merge_dir ./custom_merge
```

## Strict Mode

Use strict mode in CI or cluster scripts:

```bash
amalgkit sanity --out_dir ./ --all --strict yes --strict_level error
amalgkit sanity --out_dir ./ --all --strict yes --strict_level warning
```

`--strict_level error` fails only on errors.

`--strict_level warning` fails on warnings or errors.

## Report and Rerun

Default report path:

```text
out_dir/sanity/sanity_report.json
```

Preview reruns from that report:

```bash
amalgkit rerun --out_dir ./ --dry_run yes
```
