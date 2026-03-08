## Overview
`amalgkit quant` provides a seamless interface with the underlying [`kallisto`](https://github.com/pachterlab/kallisto) tool, an RNA-Seq quantification software. AMALGKIT automatically detects and handles single-end or paired-end reads and runs `kallisto` accordingly. `amalgkit quant` needs a kallisto index file for the reference transcriptome. Users can either provide this file or allow AMALGKIT to generate it.

## Example command

Assuming you already have a folder containing the necessary index files for quantification:
```
amalgkit quant --out_dir "./" --index_dir "./index"
```

Assuming you need `amalgkit` to create index files from cds sequences (fasta format):
```
amalgkit quant --out_dir "./" --fasta_dir "./fasta" --build_index yes
```

## Quantification from index files
`kallisto` requires index files for quantification, which are created from CDS sequences in FASTA format for each species. For instance, if your `metadata.tsv` contains a sample with the scientific_name *Mus musculus*, AMALGKIT will search for a file prefixed with `Mus_musculus` in the specified index folder (e.g., `Mus_musculus_version1.0.fasta`). With `--build_index yes`, AMALGKIT generates index files with the `.idx` extension. All CDS sequences for species listed in `metadata.tsv` must be located in the same directory and be in the FASTA format.

## Parallel processing
The `--batch` option in various AMALGKIT functions is designed for straightforward parallel processing. This argument accepts integer values, where each integer corresponds to the row number of an entry in the `metadata.tsv`. In essence, `--batch 1` processes the first sample, `--batch 2` the second, and so on. For instance, if you wish to process the third sample from the `metadata.tsv`, the command would be:

```
amalgkit quant --out_dir "./" --batch 3
```

This option is particularly beneficial when paired with array jobs on computer clusters, such as those managed by SLURM. An example SLURM `sbatch` command to utilize this feature might look like:

```
#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --array=1-3

amalgkit quant --out_dir "./" --threads $SLURM_CPUS_PER_TASK --batch $SLURM_ARRAY_TASK_ID
```

## Output files
* **SRR8819967_abundance.h5**: bootstrap results in `h5dump` format
* **SRR8819967_run_info.json**: contains run info
* **SRR8819967_abundance.tsv**: contains target_id, lentgh, eff_length, est_counts and tpm in human readable .tsv