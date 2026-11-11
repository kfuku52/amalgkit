# RNA-seq quantification benchmarks

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

We also compared the estimated log2 fold change from E0 to E2 for each
isoform with the known fold change in relative molarity. Use the same
normalized fractions and half-read floor independently in each mix.

| Model | Fold-change mean absolute log2 error | Fold-change RMSE (log2) | Spearman |
| --- | ---: | ---: | ---: |
| none | 0.705 | 1.203 | 0.893 |
| logistic | 0.797 | 1.455 | 0.895 |
| truncation | 0.554 | 1.338 | 0.924 |

Truncation improved typical absolute error and rank agreement but increased
the squared error on this fold-change comparison.

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
The direct RNA replicates contain the same E2 mix, so their expected relative
fold change is zero. Their mean absolute estimated log2 fold change was 0.537
with no model, 0.532 with logistic, and 0.650 with truncation. Better estimates
of absolute composition did not translate to closer agreement between these
two replicates.

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

We also compared the two Spearfish reference modes on the same library
using Oarfish directly with `-j 4`. After one warm-up, the median of three
runs was 39.50 s and 2.67 GB peak RSS with `--annotated`, versus 29.83 s and
0.79 GB with `--index`. All runs mapped 202,839 reads and gave the same SIRV
metrics to the displayed precision. Between modes, individual transcript
counts differed by at most 0.57 reads. The prebuilt file was 593 MiB and
took 5.32 s to create with `spearfish index` 0.1.0 and four threads.

To use it through AMALGKIT, save the index as `Test_species.sfi` in the
workspace `index` directory and pass the same Spearfish options. A separate
whole-`quant` run with each reference mode confirmed 2.62 GB versus 0.80 GB
peak RSS. The prebuilt run used only the `.sfi` file, with no FASTA or
adjacent index manifest.

With the first 80,000 reads from the direct cDNA library above and the
`ont-cdna` preset, both reference modes again worked through AMALGKIT.
Both mapped 28,649 reads and returned the same target IDs and lengths.
The sum of absolute per-target count differences was 0.015% of the total
estimated counts. This run checked agreement, not throughput.

### CPU allocation when only one run is eligible

With the first 80,000 reads from the same library, prebuilt rammap index,
and truncation model, we ran `quant --threads 4 --internal_cpu_budget 4`
with `--internal_jobs auto`. The previous allocation passed `-j 1` to
Oarfish; the corrected allocation passed `-j 4`. After one warm-up,
median wall time was 134.2 s before and 54.4 s after; the three paired
improvements ranged from 39% to 59% under varying background load. Peak
RSS was about 1.67 GB in both cases. All runs mapped 33,399 reads, and
the maximum count difference was 0.034 reads per transcript.

## Counterexamples on LRGASP simulated reads

The [LRGASP human cDNA and mouse direct RNA simulations](https://cgl.gi.ucsc.edu/data/LRGASP/)
provide the transcript of origin for each read. Oarfish 0.11.0 (rammap) with
no placement or presence model gave a lower mean log error than either
placement or presence modelling on these inputs. Within the truncation
model, presence modelling assigned less mass to transcripts with no
sampled reads, but detected fewer of the genuinely low-count transcripts.

| Simulation | Placement | Presence | Mean absolute log2 error | Count mass on zero truth (%) | Recall at 1–3 true reads (%) | Total variation |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Human cDNA | none | none | 0.839 | 1.81 | 59.5 | 0.153 |
| Human cDNA | none | spike-slab | 0.862 | 1.74 | 58.8 | 0.155 |
| Human cDNA | none | bayes | 0.881 | 1.73 | 58.3 | 0.158 |
| Human cDNA | truncation | none | 0.934 | 3.69 | 58.4 | 0.151 |
| Human cDNA | truncation | spike-slab | 0.978 | 2.74 | 56.3 | 0.149 |
| Human cDNA | truncation | bayes | 1.013 | 2.24 | 55.5 | 0.152 |
| Mouse direct RNA | none | none | 0.566 | 0.84 | 67.5 | 0.076 |
| Mouse direct RNA | truncation | none | 0.630 | 1.91 | 66.8 | 0.081 |
| Mouse direct RNA | truncation | spike-slab | 0.657 | 1.45 | 65.4 | 0.079 |
| Mouse direct RNA | truncation | bayes | 0.673 | 1.12 | 64.7 | 0.077 |

For each simulation, retain the first 250,000 reads whose exact origin ID
in `ground_truth/ONT.simulated.read_to_isoform.tsv.gz` occurs in the
[GENCODE human v38](https://www.gencodegenes.org/human/release_38.html) or
[GENCODE mouse vM27](https://www.gencodegenes.org/mouse/release_M27.html)
transcript FASTA. Use `human.ONT.simulated.fq.gz` or
`mouse.ONT.simulated.fq.gz` from the matching LRGASP simulation directory.
Use only the versioned transcript ID before the first `|` in each FASTA
header. We skipped 64,402 human and 156,791 mouse records with out-of-reference
origins before reaching 250,000 retained records, and converted `U` to `T`
in the mouse reads. The references contain 237,012 and 142,375 transcripts.

Quantify through `amalgkit quant` with `--threads 4 --internal_jobs 1`,
`--oarfish_seq_tech ont-cdna` for human or `ont-drna` for mouse, and pass
`--coverage-model` and `--presence-model` in `--oarfish_options`.
For mean absolute log2 error, first rescale estimated counts to sum to
250,000, then compare `log2((estimate + 0.5) / (truth + 0.5))` over the
transcripts with positive true counts. Count mass and total variation use
all transcripts. Recall uses the unscaled estimate with a threshold of one
read; its denominator contains only transcripts with 1–3 true reads.
A zero true count means zero in this subsample. These are simulated reads,
and the comparison excludes origins outside the supplied transcriptome.

Changing `--score-prob-denom` also failed to improve the real Garalde E2
library with truncation. The log2 errors for denominators 1, 3 (the default),
10 and 30 were 0.764, 0.746, 0.761 and 0.808, respectively.
## Short-read ERCC controls

With kallisto 0.52.0, sequence-bias correction increased concentration error
in four pure ERCC libraries from the [SEQC project](https://www.ncbi.nlm.nih.gov/bioproject/208369).
Using a fragment distribution estimated from the original pairs instead of
assuming 200/20 bp reduced concentration error when quantifying only read 1.

| SEQC site and run | Mix | Fragment mean/SD (bp) | Paired-end MAE: default → `--bias` | Read-1 MAE: assumed → estimated | Read-1 fragments mapped: assumed → estimated |
| --- | ---: | ---: | ---: | ---: | ---: |
| BGI [SRR896989](https://www.ebi.ac.uk/ena/browser/view/SRR896989) | 1 | 162.6/50.0 | 0.317 → 0.449 | 0.332 → 0.314 | 906,226 → 950,326 |
| BGI [SRR897022](https://www.ebi.ac.uk/ena/browser/view/SRR897022) | 2 | 166.3/49.7 | 0.311 → 0.653 | 0.327 → 0.318 | 898,525 → 945,005 |
| NYG [SRR1216112](https://www.ebi.ac.uk/ena/browser/view/SRR1216112) | 1 | 169.5/49.4 | 0.943 → 1.382 | 0.945 → 0.927 | 915,665 → 952,839 |
| NYG [SRR1216126](https://www.ebi.ac.uk/ena/browser/view/SRR1216126) | 2 | 168.7/46.0 | 1.030 → 1.703 | 1.024 → 1.013 | 905,744 → 948,450 |

For each library, MAE is mean absolute `log2((observed TPM + 0.5) /
(expected TPM + 0.5))` among controls whose expected TPM is at least 20.
Expected TPM is the molarity divided by total molarity within its mix, times
one million. Bias correction did not change paired-end mapping counts.
For the 54 controls meeting the threshold in both mixes, the paired-end
Mix 1/Mix 2 fold-change MAE increased from 0.153 to 0.361 at BGI and 0.252 to
0.454 at NYG. Estimated fragment lengths changed the read-1 fold-change MAE
from 0.171 to 0.154 at BGI and 0.250 to 0.252 at NYG.

To reproduce, take the first 1,000,000 synchronized pairs of each unfiltered
ENA run. Index [ERCC92.fa](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip)
with kallisto's default k-mer size and use the
[vendor's Mix 1 and Mix 2 concentrations](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt).
Run `kallisto quant -t 2` with the two FASTQs, once without extra options and
once with `--bias`. The tested binary accepted `--bias` and logged bias
learning; the [upstream manual](https://pachterlab.github.io/kallisto/manual)
lists it as a legacy option. For the read-1 comparisons, add `--single -l 200
-s 20` and then replace those two numbers with the weighted mean and
population SD of the default paired run's `aux/fld` histogram in
`abundance.h5`. The histogram contained 20,000 fragment lengths in each run.
Compare `abundance.tsv` TPM in all cases.

The inputs are fixed prefixes of pure spike-ins from two sites, with no
endogenous transcripts. The estimated distributions used information from
the original mates in the same data. They do not demonstrate recovery of
fragment lengths from a genuinely single-end library or establish that
sequence-bias correction is ineffective in other protocols.

### 3-prime tag sequencing and length correction

On the first 2 million untrimmed reads of a QuantSeq REV library
([SRR10134316](https://www.ebi.ac.uk/ena/browser/view/SRR10134316), batch 5 of
[Wu et al.](https://doi.org/10.1016/j.celrep.2020.01.068)), estimated counts
were closer to the known ERCC Mix 1 composition than length-corrected TPM.
`--single-overhang` mapped more reads to the controls but increased both
errors. We used kallisto 0.52.0, the same ERCC-only reference, and
`--single --rf-stranded`. Within each row, the only change is the use of
`--single-overhang` in the second run.

| Fragment mean/SD (bp) | Reads mapped to ERCC | TPM MAE | Relative count MAE |
| --- | ---: | ---: | ---: |
| 160/50 | 90,462 → 102,386 | 1.464 → 1.506 | 1.237 → 1.315 |
| 200/20 | 90,400 → 102,386 | 1.508 → 1.552 | 1.236 → 1.315 |
| 250/50 | 89,686 → 102,386 | 1.522 → 1.607 | 1.241 → 1.315 |

Metrics use the same log2 error as above, now for expected TPM of at least
100. Relative counts are `est_counts` divided by their sum over ERCC and
multiplied by one million. The 3′ library is a mixed human and spike-in
sample; this ERCC-only index does not quantify human transcripts.

The direction held at 200/20 bp on one million reads from each of two
batches in the same study. The batch 5 prefix came from ENA; the batch 2
prefix was extracted from the public SRA object with `fastq-dump -X 1000000`.

| Batch and run | Reads mapped to ERCC | TPM MAE | Relative count MAE |
| --- | ---: | ---: | ---: |
| 2 [SRR10134313](https://www.ebi.ac.uk/ena/browser/view/SRR10134313) | 37,308 → 40,330 | 1.684 → 1.707 | 1.521 → 1.550 |
| 5 [SRR10134316](https://www.ebi.ac.uk/ena/browser/view/SRR10134316) | 45,175 → 51,245 | 1.664 → 1.703 | 1.377 → 1.441 |

For the default placement runs, we also ran AMALGKIT `quant` with fragment
mean 200, SD 20, and `--kallisto_options="--rf-stranded"`, then `merge` and
`finalize --norm log2p1-cpm --batch_effect_alg no`. The merged counts matched
the direct kallisto outputs. All 184 expression cells matched independent
`log2(1 + 1e6 * count / library_size)` calculations to within 3.6e-15.

In the four standard SEQC read-1 runs above, `--single-overhang` at 200/20 bp
also mapped more reads, by 6.8 to 7.6 percentage points. Its fold-change
MAE improved at BGI (0.171 → 0.158) and worsened at NYG (0.250 → 0.258).
Mapping rate alone did not indicate which options recovered molarity best.
