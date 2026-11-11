# Oarfish placement models on direct RNA SIRV libraries

Using Oarfish 0.11.0 (rammap) through AMALGKIT `quant`, the truncation model
improved isoform quantification on two [Garalde et al. direct RNA
libraries](https://doi.org/10.1038/nmeth.4577). The tested libraries contain only
SIRV controls. E0 is equimolar; E2 contains four concentration levels.

| Library | Input reads | Model | Mean absolute log2 error | RMSE (log2) | Total variation | Spearman |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| E0 | 184,135 | none | 1.526 | 3.533 | 0.213 | — |
| E0 | 184,135 | logistic | 1.383 | 3.253 | 0.211 | — |
| E0 | 184,135 | truncation | 0.522 | 0.834 | 0.149 | — |
| E2 | 302,729 | none | 1.557 | 3.378 | 0.170 | 0.759 |
| E2 | 302,729 | logistic | 1.531 | 3.337 | 0.162 | 0.782 |
| E2 | 302,729 | truncation | 0.746 | 1.615 | 0.123 | 0.885 |

Each row compares all 69 isoforms. Divide each set of estimated counts by its
sum and each set of known molarities by its sum. For the log error, use
`log2(max(estimated_fraction, 0.5 / total_estimated_count) / true_fraction)`.
Total variation is half the sum of absolute differences between the unfloored
fractions. Spearman is undefined for the constant E0 truth.

## Reproduction

1. Obtain the original ENA FASTQs listed below and verify their checksums.
   Keep the full reads with no additional filtering.
2. Download Lexogen's [SIRV Set 1/2 reference and Set 1 concentration
   sheet](https://www.lexogen.com/sirvs/download). Use
   `SIRV_Set1_2_Norm_Sequences_20250402` and
   `SIRV_Set1_Norm_sequence-design-overview_20210507a.xlsx`. From the ZIP's
   seven-region FASTA and GTF, concatenate the exon intervals for each
   `transcript_id` in coordinate order and reverse-complement negative-strand
   transcripts. Use the IDs as FASTA names and no appended poly(A) tails.
   In the workbook, match transcript IDs in column G to E0 in K and E2 in M;
   use only the 69 annotated transcripts.
3. Put each mix in a separate quant workspace per model, with the same reference,
   a single-end ONT run in metadata and the FASTQ under `getfastq/<run>/`.
   Run `quant --threads 1 --internal_jobs 1 --quant_backend oarfish --oarfish_seq_tech ont-drna
   --build_index yes --clean_fastq no` and pass `--oarfish_options` with
   `"--coverage-model none"`, `"--coverage-model logistic"` or
   `"--coverage-model truncation"`. Read the estimated counts from
   `<run>_abundance.tsv`.

| Mix | FASTQ | MD5 |
| --- | --- | --- |
| E0 | [SRR6058584](https://www.ebi.ac.uk/ena/browser/view/SRR6058584) | `ff98d2f899dcf52774b55310fc70b268` |
| E2 | [SRR6058583](https://www.ebi.ac.uk/ena/browser/view/SRR6058583) | `1146f3fbd028d1b6502cb4dd38b8453e` |

Reference ZIP SHA-256:
`225be9142268ba69a0f74767fd32bca26f85783a4dbeae43e8822224fdd58f6b`.
Molarity workbook SHA-256:
`da7c04ffe8e820e7f03c9cf1a963080c7eebb58cd13b194de156926c28e76276`.

Wall times on macOS arm64 with Oarfish invoked as `-j 1` were 60–64 seconds
for E0 and 108–116 seconds for E2. These are single runs per configuration,
not evidence of a speed difference between the models. These two samples
cover one protocol from one study.

## Independent validation on SG-NEx mixed libraries

We also evaluated the first 250,000 complete reads from each of three HEYA8
libraries in the [SG-NEx resource](https://doi.org/10.1038/s41592-025-02623-4).
All contain SIRV-1 E2 and were quantified with Oarfish 0.11.0 using the same
200,648-transcript human, sequin, SIRV, ERCC and long SIRV reference. Error
metrics compare the 69 SIRV-1 isoforms and normalize counts within that set.

| SG-NEx HEYA8 library | No placement model: MAE | Logistic: MAE | Truncation: MAE | Spearman, none → truncation |
| --- | ---: | ---: | ---: | ---: |
| directRNA_replicate1_run1 | 1.127 | 1.109 | 0.747 | 0.785 → 0.931 |
| directRNA_replicate2_run1 | 1.078 | 1.069 | 0.672 | 0.793 → 0.941 |
| directcDNA_replicate1_run3 | 0.993 | 1.030 | 0.903 | 0.807 → 0.842 |

The SG-NEx data was accessed on 2026-11-11 at the [AWS Open Data Registry](https://registry.opendata.aws/sg-nex-data/).
The [sample table](https://github.com/GoekeLab/sg-nex-data/blob/master/docs/samples.tsv)
contains the public FASTQ links. The
[reference FASTA](https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/annotations/transcriptome_fasta/hg38_sequins_SIRV_ERCCs_longSIRVs_cdna.fa)
has SHA-256 `6acae821e99faefb9b8654da7694e3174f32fb82a9a3affd41ea2c45bc4842ba`.
Use the same Oarfish model settings as above, with `--seq-tech ont-drna`
for direct RNA and `--seq-tech ont-cdna` for direct cDNA. Use the E2 molarity
column for all three comparisons.

These are two direct RNA replicates and one direct cDNA library from one cell
line. The prefix subsamples are fixed, not random. Lower error on the controls
does not validate endogenous transcript abundances or other protocols.

### Mapper speed and accuracy on the same direct RNA input

We compared the first library above with the same 250,000 reads, reference,
Oarfish 0.11.0, and truncation model, using `--threads 4 --internal_jobs 1`
on an Apple M2 Max with macOS 27.0.1. The rammap runs read the existing
Oarfish `.mmi` from `--index_dir`; the Spearfish runs built an index in
memory from the workspace `fasta/Test_species.fa`. Timings include the
whole `amalgkit quant` invocation and its normal output validation. The
table gives the median of three runs per mapper after one warm-up each.

| Mapper | Wall time (s) | Peak RSS (GB) | Mapped reads | SIRV MAE (log2) | SIRV Spearman |
| --- | ---: | ---: | ---: | ---: | ---: |
| rammap | 230.53 | 1.67 | 200,420 | 0.7472 | 0.9309 |
| Spearfish | 28.40 | 2.67 | 202,839 | 0.8631 | 0.9065 |

Pass `--oarfish_options="--mapper rammap --coverage-model truncation"` or
`--oarfish_options="--mapper spearfish --coverage-model truncation"` in
separate workspaces with the same `ont-drna` preset. Wall time and peak RSS
are from `/usr/bin/time -l`. Estimated counts differed by up to 0.39 reads
per transcript between repeats with the same mapper; the SIRV metrics were
stable to the shown precision. On this sample Spearfish traded higher RAM
and error for shorter runtime. Its error was still lower than rammap with
no placement model (1.1271). This is a one-library mapper comparison.
