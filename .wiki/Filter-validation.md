# Filter validation

The optional reference, support, and iteration policies in wsfilter/csfilter
address different questions. They are not interchangeable quality diagnoses.
The default small-group policy, reference exclusion, pooled robust z, and
iteration behavior remain unchanged pending independent quality-labelled data.
Automatic exclusion remains enabled. Pairwise correlations do not use a shared
complete-case gene set, and PCA imputation does not feed exclusion scores.

## Reproducible synthetic evaluation

From a checkout with the supported Python/test dependencies:

```bash
PYTHONPATH=. python benchmarks/evaluate_filter_accuracy.py --seeds 12 --output filter-accuracy.json
```

The script generates 200-gene log-scale profiles for two tissues, with two or
five replicates, one species for wsfilter and five for csfilter. Twelve fixed
seeds vary common expression, tissue/species effects, and sample noise. It
compares unchanged default policies with small-group retention, a 50-pair support
cutoff, project/species exclusion, species-group z, and a single wsfilter pass.
The value 50 is an evaluated example, not an established minimum for RNA-seq QC.

Quality labels come from the generating intervention: a known individual swap,
a species-wide technical swap, or a swap with only two observed genes. Controls
include clean expression, a biological lineage-specific shift, and sparse
observation without a swap. The lineage shift and species-wide technical swap
have identical observed expression but different ground truths. This explicitly
tests the limits of identification from expression alone.

The JSON reports TP, FP, FN, TN, and removal of biological-shift samples per
scenario and replicate count. Inspect strata rather than averaging their
arbitrarily chosen prevalences into a single accuracy score. Seeds and runs
sharing reference profiles are not independent biological validation cohorts.

## Adoption criteria and limitations

A support cutoff prevents unstable sparse correlations from triggering
exclusion, but also stops detection of genuinely incorrect sparse samples.
Species-group z can preserve lineage shifts and detect isolated swaps when
replication is adequate, but can miss a species-wide technical problem with the
same expression. Small-group retention similarly trades reduced false removal
for missed quality failures. Species/project exclusion can remove self-related
support without resolving the biological/technical ambiguity. One pass prevents
later removals, but can miss anomalies masked on the initial pass.

Consequently, the new policies are available as explicit options; no scientific
default is changed based on these simulations alone. Implementation fixes
are adopted independently: removal-round evidence survives later wsfilter
rounds and retain their scoring thresholds; newly undefined scores clear stale
values on reruns; prior exclusion reasons are preserved by mapping-rate filtering;
missing project placeholders cannot act as known reference projects;
and `--one_outlier_per_iter yes` enforces both its group and project
limits with deterministic tie handling. Constant finite-pair correlations are
reported as undefined without a numerical warning.

Before changing defaults, evaluate independent quality labels and genuine
biological shifts. Split calibration and evaluation by project/donor and by
species or clade, keep normalization and averaging fixed, and report detection
sensitivity, precision, false removal, biological effect retention, and
uncertainty by stratum. Do not judge success only by increased within-tissue
correlation or cleaner PCA separation. No real-data sensitivity, specificity,
or universally optimal threshold is established by this synthetic evaluation.

## Compatibility

Existing CLI defaults preserve the exclusion policy. Added metadata columns are
additive. The one-per-iteration constraint fix can change results when explicitly
enabled; handling undefined constant correlations also avoids numerical-warning
failures. Removed-run score persistence changes diagnostic output, not which
runs are removed. Compare alternative settings from the same pre-filter
metadata. A changed retained set requires downstream finalization, group means,
tau, and comparative analyses to be recomputed; raw quantification need not be
repeated solely for these filter settings.

See [wsfilter](https://github.com/kfuku52/amalgkit/wiki/amalgkit-wsfilter) and
[csfilter](https://github.com/kfuku52/amalgkit/wiki/amalgkit-csfilter) for options
and output columns. Normalization and averaging changes must be evaluated
separately before thresholds calibrated on them are reused.
