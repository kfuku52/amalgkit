## Overview

`amalgkit busco` prepares per-species BUSCO tables used by ortholog-based downstream steps.

Typical downstream uses:

- `amalgkit cstmm`
- `amalgkit csfilter`

## Example

Using BUSCO:

```bash
amalgkit busco --out_dir ./ --tool busco
```

Using compleasm:

```bash
amalgkit busco --out_dir ./ --tool compleasm
```

## Main outputs

- `busco/<Species>.tsv` or `busco/<Species>_busco.tsv` depending on workflow context
- summary files needed by `cstmm` and `csfilter`

## Notes

- `busco` itself is a wrapper step; the required external executable depends on `--tool`
- if you already have compatible BUSCO full tables, you can skip this command and pass `--dir_busco` directly to downstream commands
