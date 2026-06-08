## Overview

`amalgkit sanity` scans a workspace and reports missing or ambiguous outputs for selected pipeline steps. It writes a machine-readable report that can be consumed by `amalgkit rerun`.

## Usage

Check all supported targets:

```bash
amalgkit sanity --out_dir ./ --all
```

Check selected targets:

```bash
amalgkit sanity --out_dir ./ --check getfastq,quant,merge
```

Limit checks by run or species:

```bash
amalgkit sanity --out_dir ./ --run SRR000001,SRR000002 --check getfastq,quant
amalgkit sanity --out_dir ./ --species "Homo sapiens,Mus musculus" --check merge,busco,finalize
```

## Targets

Supported check targets are:

- `getfastq`
- `index`
- `quant`
- `merge`
- `busco`
- `finalize`
- `all`

If no target is specified, `--all` is assumed.

## Strict mode

Use strict mode in CI or cluster scripts:

```bash
amalgkit sanity --out_dir ./ --all --strict yes --strict_level error
amalgkit sanity --out_dir ./ --all --strict_level warning
```

## Report and rerun

The default report path is:

```text
out_dir/sanity/sanity_report.json
```

To preview reruns from that report:

```bash
amalgkit rerun --out_dir ./ --dry_run yes
```
