## Overview

`amalgkit integrate` scans local FASTQ files and writes AMALGKIT-compatible metadata. Use it when a project includes private FASTQ files or non-SRA files that should enter the same downstream pipeline as public runs.

## Common Uses

Create metadata from private FASTQ files only:

```bash
amalgkit integrate --fastq_dir ./private_fastq --out_dir ./ --output_metadata ./metadata/with_private.tsv
```

Append private FASTQ entries to existing metadata:

```bash
amalgkit integrate \
    --fastq_dir ./private_fastq \
    --out_dir ./ \
    --metadata ./metadata/metadata.tsv \
    --output_metadata ./metadata/with_private.tsv
```

Use the first example in a new workspace without `metadata/metadata.tsv`.
If that file already exists, inferred input uses merge mode; the second example
makes the input explicit. Both examples write a separate output for review and
preserve the original public metadata.

## Output Paths

Without `--output_metadata`:

- standalone mode writes `out_dir/metadata_private_fastq.tsv`
- merge mode writes `out_dir/metadata/metadata_updated_for_private_fastq.tsv`

Neither path is the default input of `getfastq`, `quant`, or `merge`.
Without an explicit handoff, private samples can be omitted or the next command
can fail because `metadata/metadata.tsv` does not exist.

## FASTQ Discovery

Files under `--fastq_dir` are scanned recursively. The first subdirectory below `--fastq_dir` is parsed as `scientific_name`, with underscores converted to spaces.

Example:

```text
private_fastq/Homo_sapiens/brain/Sample1_1.fq.gz
private_fastq/Homo_sapiens/brain/Sample1_2.fq.gz
```

This assigns `scientific_name` to `Homo sapiens`.

If the same FASTQ basename appears under multiple species directories, AMALGKIT prefixes the generated run ID with the species directory token.

## FASTQ Naming

Supported extensions:

- `.fastq`
- `.fastq.gz`
- `.fq`
- `.fq.gz`

Typical paired-end names:

- `Sample1_1.fq.gz`
- `Sample1_2.fq.gz`

Typical single-end name:

- `Sample1.fq.gz`

## Useful Options

| Option | Use |
| --- | --- |
| `--fastq_dir` | root directory to scan recursively |
| `--metadata` | existing AMALGKIT metadata to update |
| `--output_metadata` | explicit output metadata path |
| `--seqkit_exe` | path to `seqkit` |
| `--accurate_size yes/no` | exact or sampled size estimation for gzip FASTQ files |
| `--remove_tmp yes/no` | remove temporary files |

## Fields to Review

`integrate` cannot infer every biological label. Before continuing, review the
generated `metadata/with_private.tsv`, including:

- `scientific_name`
- `sample_group`
- `exclusion`
- `is_sampled` (use `yes` for runs intended for processing)
- `lib_layout`
- `read1_path`
- `read2_path`

`sample_group` is used by `select`, `wsfilter`, `csfilter`, and `finalize`.

## Example Rows

| scientific_name | sample_group | run | read1_path | read2_path | exclusion | lib_layout | private_file |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Mus musculus | brain | MMB1 | /path/MMB1_1.fq.gz | /path/MMB1_2.fq.gz | no | paired | yes |
| Mus musculus | liver | MML1 | /path/MML1.fq.gz |  | no | single | yes |

## Next Steps

Install the [external dependencies](https://github.com/kfuku52/amalgkit/wiki/Installation-and-dependencies).
Even private-only `getfastq` currently checks both `fasterq-dump` and SeqKit at
startup. Supply a reference transcriptome for every species under `fasta/` as
described in [quant](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant).

After reviewing the generated metadata, pass the same path through all three
commands. Then use the metadata exported by `merge` for finalization:

```bash
amalgkit getfastq --out_dir ./ --metadata ./metadata/with_private.tsv
amalgkit quant --out_dir ./ --metadata ./metadata/with_private.tsv --build_index yes
amalgkit merge --out_dir ./ --metadata ./metadata/with_private.tsv
amalgkit finalize --out_dir ./ --input_dir ./merge --metadata ./merge/metadata.tsv --batch_effect_alg no
```

This finalization example is for short reads. For Oarfish/long-read data, use
`--norm log2p1-none`; see
[metadata and normalization](https://github.com/kfuku52/amalgkit/wiki/Metadata-and-normalization).
