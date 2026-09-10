# CSTMM imputation and validation

## Production normalization

CSTMM keeps orthogroups that are single-copy in at least the selected percentage
of eligible species (`--single_copy_threshold`, default 50%). It does not require
any orthogroup to be present in every species. Only single-copy entries that
match a target in the count table supply observations. Zero-copy, multiple-copy
and unmatched entries are missing reference values; an observed count of zero
remains zero.

By default, `--tmm_imputation_scale library_size` divides each reference column
by its original all-target library size and multiplies by one million before
EM-PCA. Missing values are reconstructed in this CPM space and converted back
to the sample's count scale. Observed counts are restored exactly. The existing
two-round TMM implementation receives that completed count matrix and the
original library sizes. Actual correction always uses the resulting factors.
The output is still original per-species counts divided by those factors;
imputed orthogroup expression is not exported as measured expression.

This separates sequencing depth from imputation. It does not prove that
missing ortholog expression is identifiable or that species-specific expression
and annotation differences can be removed by one scalar. A genuinely absent
gene has no hidden expression value to recover: its imputation is an auxiliary
reference estimate, not evidence that the gene is expressed. Single-copy status
alone does not establish the non-DE reference assumption of
[TMM](https://doi.org/10.1186/gb-2010-11-3-r25).

## Settings and diagnostics

| Option | Default | Meaning |
| --- | --- | --- |
| `--tmm_imputation_scale` | `library_size` | CPM-space imputation; `raw` reproduces the former count-space method for comparison |
| `--tmm_imputation_rank` | `4` | Maximum PCA rank, capped by matrix dimensions |
| `--tmm_imputation_max_iter` | `50` | Maximum outer EM-PCA iterations |
| `--tmm_imputation_tol` | `1e-6` | Positive absolute change tolerance in the chosen imputation scale |
| `--tmm_allow_unconverged` | `no` | Explicitly accept nonconvergence or row-mean fallback, with warnings |
| `--tmm_reference_diagnostics` | `yes` | Write observed-only comparisons for every pair of retained samples |

Convergence means the maximum change in missing entries met the tolerance; it
is not a measure of predictive accuracy. A saturated rank can preserve a poor
initial estimate while satisfying that numerical criterion. By default an
unconverged/fallback reconstruction stops the staged run. Increase the iteration
limit or explicitly accept the estimates with `--tmm_allow_unconverged yes`.
No alternative factor estimator is automatically substituted.

`cstmm_normalization.json` records requested/resolved rank, scale, iteration
count, final change, convergence, fallback, number of clipped imputations,
reference rows, complete rows, reference samples, single-copy threshold,
software version, input count directory and count units. `cstmm_missingness.tsv`
gives sample-level missing and observed reference counts and missing fractions,
including columns not used for factor estimation.
`cstmm_orthology_status.tsv` reports candidate orthogroup/species entries as
`observed`, `zero_copy`, `multiple_copy` or `target_not_found`. These are input
facts: `zero_copy` does not distinguish biological loss from incomplete annotation.
Malformed raw counts and empty or duplicate target IDs are rejected rather than
turned into reference missingness. Literal identifiers such as `NA` and `001`
are preserved.

## Observed-only comparisons are reference only

By default, each sample pair uses its jointly
observed entries, the original library sizes and the same TMM kernel. It does
not require a globally complete orthogroup. The table records shared observed,
positive and retained entries, estimation status, the observed-only factor
ratio, the production factor ratio and their log2 difference. The second sample
is the numerator. Unestimable comparisons have blank estimates, not factor 1.

These comparisons **never replace, blend with, or otherwise change production
factors**. They are an audit, not a fallback. Pair-specific gene sets and
references differ from the production reference, so a discrepancy is not itself
proof of imputation error. Shared-positive counts describe support; they are
not calibrated confidence scores or p-values. Disconnected or weakly overlapping
species may have little observed evidence for their relative scale.

The audit evaluates all sample pairs, so its work grows quadratically
with sample count. Rows are streamed to disk. Use `--tmm_reference_diagnostics no`
to skip the comparison table and PDF for large compendia; imputation settings
and missingness diagnostics are always saved.

Normal CSTMM execution writes `cstmm_observed_pair_comparison.pdf` from the
input expression counts and orthology table. The left panel
compares each sample with the **actual global round-2 reference**, recomputing
the observed-only TMM directly in that orientation. The right panel shows all
sample pairs. In each panel, x is `log2(global_factor[sample] / global_factor[reference])`
and y is the log2 observed-only pairwise factor ratio. Agreement lies on the
dashed `y = x` line; a vertical gap of 1 means a twofold difference in the ratio.
Both panels share axis limits. Unestimable pairs are counted and not plotted;
overlapping points and the dependence between pairs should not be read as
independent evidence. Shared gene counts remain available in the TSV.

The global factors have a geometric-mean-one constraint, so matching the
reference sample name alone does not justify comparing separately centered
factor values. Comparing **ratios to that reference** removes the arbitrary
common scale. Both methods must also use the same original library sizes and
TMM options. Even then, observed-only and imputed gene sets and trimming may
differ. The all-pairs panel additionally reflects reference dependence, since
pairwise TMM estimates need not be transitive. This is a sensitivity comparison,
not a direct estimate of imputation error. `global_reference_sample` is recorded
in the TSV; pairs containing that reference always put it in `reference_sample`.
The plot loads the pair table in memory, in addition to the streamed estimation.

## Sensitivity and migration

Run comparisons in separate output directories, pointing `--dir_count` and
`--metadata` at the same original merge inputs. Compare thresholds 50/75/100
(100 may have no usable reference), raw/library-size scales and several ranks.
Do not select a setting only because the same samples look better separated in
PCA. Evaluate factor-ratio and observed-gene fold-change errors on held-out
gene-by-species blocks, not just randomly hidden individual cells. Masking
experiments measure recoverable observations; true gene loss needs a separate
simulation. Where available, compare with independent biological measurements.

The regression suite includes 100 simulated species with no universally observed
orthogroup and a tenfold-depth block-missingness counterexample. Those controlled
examples are not a validation of arbitrary species, tissue or protocol mixtures.

Changing imputation scale/rank/reference selection requires regenerating CSTMM
and all dependent filters, finalization, tables and figures together. To reproduce
the old imputation explicitly, use `--tmm_imputation_scale raw` and, if needed,
`--tmm_allow_unconverged yes`. Input validation now rejects malformed counts.
Changing only downstream abundance to `log2p1-cpm` does not require new factors
if the original `tmm_library_size` and matching counts remain available, but
all downstream stages that use the changed abundance must be rerun. Existing
files are never silently migrated.
