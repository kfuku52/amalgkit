## Overview

The `amalgkit getfastq` command serves three primary purposes:

1. Retrieving SRA data and extracting FASTQ files.
2. Quality control of FASTQ files using [`fastp`](https://github.com/OpenGene/fastp).
3. Optional filtering with MMseqs2:
   - rRNA removal (`--rrna_filter yes`)
   - taxonomy-aware contaminant removal (`--contam_filter yes`)

To operate, `amalgkit getfastq` requires a metadata table (`--metadata`) produced from one or more of the following commands: [`amalgkit metadata`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata), [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate), or [`amalgkit select`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-select). Alternatively, you can input a specific BioProject, BioSample, or SRA ID directly using `--id`, which will generate the corresponding RNA-seq fastq files.

The resulting fastq files are suitable for tasks like expression level quantification with [`amalgkit quant`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-quant) or for transcriptome assembly using tools such as [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq). 

## Example command

### Generating FASTQ files for all listed data in `metadata.tsv`
If a `metadata/metadata.tsv` file exists in the same working directory where you are running `amalgkit getfastq` (this is the standard output of `both amalgkit metadata` and `amalgkit integrate`), AMALGKIT automatically recognizes it, eliminating the need to specify the `--metadata` option. This command will sequentially process each sample in the `metadata.tsv`. If you have access to a computer cluster, please refer to the relevant chapter.
```
amalgkit getfastq --metadata /PATH/TO/metadata.tsv
```

### Generating FASTQ files for one specified SRA data
```
amalgkit getfastq --id DRR461654
```

## Two-step FASTQ extraction
With `amalgkit getfastq`, when a target FASTQ size is set using `--max_bp`, a unique two-step extraction process is initiated:

1. **First-round extraction:** During this phase, AMALGKIT attempts to extract a FASTQ file that exactly matches the target size specified.
2. **Second-round compensatory extraction:** Since there are inherent losses during processing by [fasterq-dump](https://github.com/ncbi/sra-tools), [fastp](https://github.com/OpenGene/fastp), and optional downstream filters (rRNA/contaminant removal), AMALGKIT factors in this loss. It then initiates a second round of extraction to supplement the FASTQ reads. This ensures the cumulative data size closely aligns with the target size.

This two-fold approach ensures precision while maximizing data utility. This feature is available only for data derived from the SRA and not for private FASTQ files.

## Interpreting `getfastq_stats.tsv`
`getfastq_stats.tsv` contains both count and base metrics. For paired-end libraries, count columns can use different units depending on the stage:

- `num_dumped`, `num_written`, `num_rrna_in/out`, and `num_contam_in/out`: spot counts (paired-end spots = read pairs)
- `num_fastp_in/out`: read counts reported by `fastp` (mates counted separately in paired-end data)
- `bp_*`: total bases

For stage-by-stage removal fractions, `bp_*` columns are the safest values to compare across all stages.

## Cloud options
At present, `amalgkit getfastq` can retrieve samples from three cloud storage services: Amazon (AWS), Google (GCP), and NCBI. While AWS and GCP access might be subject to location constraints and other factors, by default, `amalgkit getfastq` will attempt to download a sample from one of these services. If the sample is not available from the initially chosen service, AMALGKIT will try the next one. If download fails from all configured cloud services, `amalgkit getfastq` exits with an error.

## Local FASTQ files
If you intend to use local FASTQ files that are not available on the SRA in conjunction with AMALGKIT, you must first utilize [`amalgkit integrate`](https://github.com/kfuku52/amalgkit/wiki/amalgkit-integrate) to generate a `metadata.tsv` file containing all pertinent details. Once this is done, `amalgkit getfastq` will be able to process the local files using [fastp](https://github.com/OpenGene/fastp) and carry out other subsequent steps in the same manner as it would for FASTQ files derived from the SRA.

## Parallel processing
The `--batch` option in various AMALGKIT functions is designed for straightforward parallel processing. This argument accepts integer values, where each integer corresponds to the row number of an entry in the `metadata.tsv`. In essence, `--batch 1` processes the first sample, `--batch 2` the second, and so on. For instance, if you wish to process the third sample from the `metadata.tsv`, the command would be:

```
amalgkit getfastq --out_dir "./" --batch 3
```

This option is particularly beneficial when paired with array jobs on computer clusters, such as those managed by SLURM. An example SLURM `sbatch` command to utilize this feature might look like:

```
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --array=1-3

amalgkit getfastq --out_dir "./" --threads $SLURM_CPUS_PER_TASK --batch $SLURM_ARRAY_TASK_ID
```

## Tips for subsequent transcriptome assembly
When it comes to assembly, including more RNA-seq libraries generally yields more transcripts. However, handling an excessive amount of data can be computationally demanding. To address this, `amalgkit getfastq` offers an automatic subsampling feature for RNA-seq reads from various libraries. The volume of data necessary, denoted by the --max_bp parameter, can vary based on several factors, including the specific assembly program in use. For a deeper dive into this topic, you can refer to [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146062).
