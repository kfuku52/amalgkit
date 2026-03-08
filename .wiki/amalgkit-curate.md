## Legacy page

`amalgkit curate` is not exposed by the current CLI. The old monolithic workflow has been replaced by Python implementations of:

- [`amalgkit wsfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-wsfilter)
- [`amalgkit csfilter`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csfilter)
- [`amalgkit finalize`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-finalize)

## Migration guide

Old workflow:

```bash
amalgkit merge
amalgkit cstmm      # optional
amalgkit curate
```

Current workflow:

```bash
amalgkit merge
amalgkit cstmm      # optional
amalgkit wsfilter
amalgkit csfilter   # optional
amalgkit finalize
```

## Batch correction in the current workflow

`amalgkit finalize` now provides Python backends for:

- `no`
- `sva`
- `ruvseq`
- `combatseq`
- `latent_glm`

## Metadata column name

Older documents sometimes refer to `curate_group`. Current releases use `sample_group`.
