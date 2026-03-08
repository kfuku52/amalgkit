## Overview
The `amalgkit select` command assists you in selecting data that is suitable for downstream analysis. Using this command, you can mark data that is unsuitable for evolutionary gene expression analysis (such as data from miRNA-seq) with a `no` value in the `is_qualified` column. Some samples, like those that have been intensively sequenced (e.g., livers of *Bos taurus*), can be subsampled using the `--max_sample` option while preserving the BioProject diversity as much as possible. The excluded data from this operation are marked with a `no` value in the `is_sampled` column, and as a result, these samples will not be included in downstream analysis.

The `amalgkit select` command can also take a directory containing a series of config files as an optional input, specified by the `--config_dir` option. This provides further customization and control over the data selection process. For detailed instructions and further information, consult the [amalgkit config guide](https://github.com/kfuku52/amalgkit/wiki/amalgkit-config).

## Example command
```
amalgkit select --out_dir "./" --config_dir "./config_base"
```
