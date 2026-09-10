# Batch correction models, failure policy, and diagnostics

Batch correction in `finalize` is optional (`--batch_effect_alg no` by default).
The corrected tables are exploratory expression summaries; successful fitting,
a clean PCA plot, or a stable selected dimension does not demonstrate biological
validity. BioProject can encode biological differences as well as technical ones.

## Shared contract

`--batch_failure_policy skip` retains every run of an affected species without
batch correction when its model cannot be estimated. Requested normalization and
expression transformation still apply. The program reports the species and reason,
and metadata contains `batch_status` and `batch_skip_reason` in addition to the
existing `batch_corrected` and `batch_alg_used` columns. No subset is silently
corrected and mixed with uncorrected runs. `--batch_failure_policy error` stops on
these model failures. Malformed inputs, missing dependencies, and unexpected
programming errors stop under either policy.

Statuses distinguish `corrected`, `skipped`, and `not_needed` (e.g. explicit k=0).
An unresolved dimension estimate is not reported as an estimated zero. The
historical `stable` fields refer to a particular numerical/selection stage, not
biological validity. No failed fit retries after deleting a protected covariate,
substituting another likelihood, or inventing latent directions.
If fewer factors than an explicitly requested positive dimension can be estimated,
SVA, RUV and the log-linear latent backend follow the failure policy, rather than
reporting the fit as an intentionally selected zero or smaller dimension.

The protected design includes an intercept and treatment-coded `sample_group`.
Additional columns have explicit types:

```bash
amalgkit finalize --out_dir ./ --metadata ./csfilter/metadata.tsv \
  --batch_effect_alg latent_loglinear --latent_k 1 \
  --batch_categorical_covariates sex treatment \
  --batch_continuous_covariates age
```

Categorical reference levels are sorted; continuous covariates are centered and
scaled. Their encoding, run order, numeric matrix, rank, residual degrees of
freedom, and singular values are saved. Missing covariates and duplicate run IDs
are errors. Constant/redundant protected columns make a design unidentifiable and
follow the failure policy. Missing biological information is not filled with an
invented category.

The joint biological/batch design is checked separately, with a group-by-project
table. Complete confounding prevents ComBat-seq fitting. Latent methods do not
fit BioProject directly, so their diagnostic of confounding does not automatically
stop estimation of residual variation. Such residual variation must not be
interpreted as an identified project effect. No method guarantees preservation
of unrecorded biological covariates.

## ComBat-seq

ComBat-seq requires finite, nonnegative **integer raw counts**. It does not accept
CSTMM-divided counts or silently round estimated/fractional counts. Use an
appropriate raw-count input; `finalize` rejects the known CSTMM input path for
this backend. Library offsets/normalization for other backends remain a separate
input contract, not inferred from corrected column sums.

The default `--combatseq_group_model protect` retains sample_group. An explicit
`batch_only` choice omits it; additional protected covariates, if supplied, still
remain. A singleton batch, insufficient batches, confounded design, no residual
degrees of freedom, or a model fitting failure skips correction for the entire
species. The original group model is never replaced by a batch-only retry.

The adapter passes a prebuilt numeric patsy design to InMoose, preserving the
documented encoding rather than reinterpreting category names as formula code.
The intercept is retained for InMoose's within-batch dispersion fits. The fixed
balanced raw-count fixture matches R ComBat-seq exactly with InMoose 0.9.1.
See the [original ComBat-seq implementation](https://github.com/zhangyuqing/ComBat-seq)
and [InMoose covariate preprocessing](https://github.com/epigenelabs/inmoose/blob/master/inmoose/pycombat/covariates.py).

## Experimental log-linear latent removal

`latent_loglinear` is the canonical name; `latent_glm` is a compatibility alias.
This method rescales counts by library totals, applies `log(x + 0.5)`, regresses
on the protected design, and estimates latent directions with a residual SVD.
It subtracts the latent fitted effect, reverses the transform, floors negative
values, and restores the input count scale.

`--latent_weighting uniform` (old `--latent_family poisson`) uses uniform weights.
`dispersion` (old `nb`) uses gene weights derived from a moment dispersion
estimate. These are **not** Poisson/NB likelihood fits. The reported objective is
weighted mean squared error on log residuals, evaluated at the final factors.
Subspace stability is the iteration convergence criterion; failure to converge
follows the shared failure policy.

Specify a nonnegative `--latent_k`. Uncalibrated `auto` is skipped by default;
`--latent_k_selection legacy` explicitly enables the historical 10% leading /
60% cumulative spectral-energy heuristic. Its tendency to select noise factors
is known: the review's independent Poisson(100), 1000-gene, 8-run fixture selects
k=4. It is not a test for technical batch variation. No newly calibrated auto
method is claimed by this change.

## SVA

`--sva_nsv_permutations` and the historical alias `--sva_B` control dimension
estimation only. `--sva_irw_iterations` independently controls IRW iterations
(default 5). This changes historical runs where one B value controlled both.
`--sva_estimation_method be` is the default; `leek` selects the other estimator,
and `be_then_leek` explicitly permits switching when BE is unresolved.

Zero or dependent candidate directions are removed, never filled with an
arbitrary orthogonal basis. IRW reports its actual rank and completed iterations.
If requested manual dimensions cannot be estimated, correction is skipped (or
raises under strict policy). Dimension-selection stability and completed IRW
iterations do not establish IRW convergence. A biological contrast is required
for the IRW full/null model comparison unless nsv=0.

For exploratory correction, only the SV space orthogonal to the protected design
is subtracted. The estimated SVs and this removal basis are exported separately.
R numerical fixtures check the projection independently. The Python estimator's
local false-discovery-rate implementation is not numerically identical to R sva:
on the committed 120-gene fixture the factor-projector distance is about 0.04475
and maximum corrected log-expression difference is about 0.03228. The same
nondegenerate input with matched IRW iterations changes by only about 4e-14
between the reviewed implementation and this repair. These observations are
fixture-specific, not a bound for arbitrary datasets.

## RUVr-based implementation

RUV fits observed counts (without adding 1 to the likelihood input), obtains
deviance residuals, and estimates W from selected residual rows. Its moment-based
Poisson/NB dispersion fitting is an AMALGKIT implementation, not edgeR's empirical
Bayes dispersion procedure. NB failure does not revert to Poisson, and a failed
GLM does not revert to OLS/ANOVA. The log transform for exploratory correction
still uses a pseudocount of 1.

Controls are explicit choices:

- `--ruvseq_control_genes all`: all retained genes, an explicit RUVr choice.
- `empirical` (historical alias `auto`): sufficiently expressed genes ranked by
  non-significance, followed by MAD selection. Non-significance does not establish
  expression invariance.
- `file` with `--ruvseq_control_file PATH`: one unique retained gene ID per line,
  no header. The selected IDs are recorded.

Every strategy must meet `--ruvseq_min_controls` (at least 2). Insufficient controls
skip correction; they never trigger selection of all genes. A GLM failure skips
the whole fit. k=0 does not run a GLM or select controls.

Specify `--ruvseq_k`. Uncalibrated auto requires `--ruvseq_k_selection legacy`;
otherwise it is skipped. The legacy PCA score is an exploratory, same-data
criterion. Ties prefer k=0, then smaller positive k. This compatibility setting
enables the old criterion, not the unsafe old failure behavior or guaranteed
reproduction of historical results.

The original W is exported for downstream models. Exploratory corrected counts
remove only its component orthogonal to the protected design, which is an
AMALGKIT-specific protection step. RUVr kernel reference tests use identical
edgeR residuals to isolate the SVD/correction calculation from the different
GLM implementations. They do not establish end-to-end equivalence to RUVSeq.

For differential expression, use original counts with an appropriate offset and
design containing the biological covariates and W; check that this augmented
design is identifiable. Do not use the exploratory pseudo-count output as if it
were original sequencing counts. See the [RUVSeq vignette](https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html).

## Final matrix and output diagnostics

All after-QC panels and downstream aggregation receive the **final saved matrix**,
after clipping, rounding, scale transformation, and input-zero restoration.
`--maintain_zero yes` remains the default for compatibility. Its mask comes from
observed input zeros, not zeros created by a log transform. Observed zeros need
not mean structural nonexpression. Compare against `--maintain_zero no` when
low expression matters. Clipping within a backend's count reconstruction is
distinct from the finalize `--clip_negative` option.

Schema 2 adds detailed JSON and TSV outputs while retaining existing summary
columns and expression filenames. JSON contains parameter and dependency
versions, input/final scales, design encoding, fitted/retained gene IDs, controls,
factor values, removal basis, and postprocessing events. Each event records its
scale, changed cells/genes/runs, summed absolute change, and maximum change.
Summed event counts can count the same cell more than once.

Linear design protection applies on the model scale before nonlinear operations.
It does not guarantee preservation of arithmetic means, correlations, or tau
after back-transformation, clipping, rounding, or zero restoration. The
aggregation/tau definitions and CSTMM/depth normalization are separate contracts;
this change does not redefine them.

## Migration and scientific validation

Recompute finalize outputs and downstream consumers for affected ComBat fallbacks,
SVA degeneracy/IRW settings, latent fits, and RUV fits. Old results require the
original software and parameters for exact reproduction. A QC-order-only change
does not change numerical expression values, except the corrected observed-zero
mask semantics for plain log transforms. No upstream re-quantification is required
solely for these changes. A changed CSTMM input contract may independently require
upstream recomputation.

The fixed R fixtures and their generation script are in
[tests/fixtures/batch_effect_reference](https://github.com/kfuku52/amalgkit/tree/master/tests/fixtures/batch_effect_reference).
They validate numerical components, not population-level scientific performance.
A future default auto-k requires independent null calibration and biological
effect retention tests across sample size, depth, sparsity, partial confounding,
and unrecorded biological variation, with project-level validation. The repair
does not equate an arbitrary residual direction with a technical batch effect.
