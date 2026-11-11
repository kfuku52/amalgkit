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
   Run `quant --threads 4 --quant_backend oarfish --oarfish_seq_tech ont-drna
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

Four-thread wall times on macOS arm64 were 60–64 seconds for E0 and 108–116
seconds for E2. These are single runs per configuration, not evidence of a
speed difference between the models. Only one protocol from one study was
tested here; the results support an opt-in comparison on similar libraries.
