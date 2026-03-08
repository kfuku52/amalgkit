## Overview

`amalgkit merge` consolidates `amalgkit quant` outputs into per-species tables. These merged tables are used directly by:

- `amalgkit cstmm`
- `amalgkit wsfilter`
- `amalgkit csfilter`
- `amalgkit finalize`

`merge` is Python-only in current releases.

## Example

```bash
amalgkit merge --out_dir ./
```

## Main outputs

For each species:

- `<Species>_est_counts.tsv`
- `<Species>_eff_length.tsv`
- `<Species>_tpm.tsv`

Top-level summary PDFs are also produced, including mapping, library-layout, duplication, insert-size, and expression summary plots.
