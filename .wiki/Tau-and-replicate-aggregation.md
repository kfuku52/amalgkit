# Tau and replicate aggregation

Tau is computed from **linear-scale arithmetic means**. Each run is inverse-transformed using `--norm` before aggregation. The inverse only removes the logarithm: FPKM, TPM or the supplied normalization remains in place. With batch correction, these are means of inverse-transformed corrected values, not unbiased estimates of the original expression. The final tau table uses the exported final expression matrix, including any zero restoration.

For a complete panel of `n` sample groups, with nonnegative representatives `x`, tau is `sum(1 - x / max(x)) / (n - 1)`. A value of zero denotes equal expression and one denotes expression in only one group. The biological interpretation also depends on normalization, tissue definitions and sampling.

Earlier outputs averaged log expression and then inverse-transformed that average. This is a different representative: for log2p1 it is the geometric mean of `(x + 1)`, minus one. For A = [0, 100] and B = [10, 10], the earlier calculation gave tau approximately 0.095 with B highest. Linear means are A = 50 and B = 10, giving tau = 0.8 with A highest. Neither estimator is universally optimal; this workflow explicitly targets arithmetic mean expression.

## Replicate weights

`wsfilter`, `csfilter` and `finalize` accept the following options for their per-species tau summaries:

| Option | Meaning |
| --- | --- |
| `--tau_unit run` | Default. Each retained run has equal weight. |
| `--tau_unit biosample` | Average runs within each BioSample, then give BioSamples equal weight within each group. |
| `--tau_unit donor` | Average runs within each curated donor, then give donors equal weight within each group. |
| `--tau_balance_projects yes` | Average the selected units within each BioProject, then give projects equal weight. Default is `no`. |

The metadata must contain the selected ID column (`biosample` or `donor`) and, for project balancing, `bioproject`. Missing IDs cause an error; there is no implicit fallback to run weights. Donor IDs must be curated and unambiguous across the species, including across projects. The same donor can occur in different tissues. Project balancing rejects a unit shared across projects within the same group, which would otherwise double-count it.

These options average retained runs within a unit; they do not pool reads or reconstruct libraries. Curate technical replicates and distinct conditions/time points before aggregation. A BioSample accession does not guarantee an independent donor, and project balancing does not remove biological confounding. Sex, age and treatment are not automatically balanced. The default select deduplication key is `(bioproject, biosample)`; it does not establish donor independence or resolve cross-project duplicates. Removed runs cannot be recovered downstream.

For example, after adding verified `donor` IDs to metadata:

```bash
amalgkit finalize --out_dir out --metadata curated_metadata.tsv --sample_group brain,liver --tau_unit donor --tau_balance_projects yes
```

For a weighting sensitivity analysis, rerun with the same retained metadata, input data, normalization and tissue panel, changing only these options. Compare tau differences and highest-tissue changes. These are descriptive sensitivity comparisons, not donor-bootstrap confidence intervals; automatic bootstrap and leave-project-out inference are not implemented.

## Comparable tissues and missing expression

Use the same explicit `--sample_group` list and harmonized tissue labels across species to define a comparison panel. Without this option, groups are resolved from the supplied metadata before species-specific aggregation; the resulting list is recorded. A group is not dropped from the tau panel merely because it is absent in one species or all its runs were excluded. A panel identifier records the exact group set; matching identifiers alone do not establish anatomical equivalence or compatible normalization/weights.

Tau and highest tissue are undefined (`NA`) when any required representative is missing or nonfinite, when fewer than two groups are present, or when expression is all zero or negative. Missing gene/run values propagate through aggregation rather than silently reweighting remaining units. Invalid negative input values are not clipped by the tau calculation. The upstream correction settings determine clipping; negative/invalid runs make their group representative missing. Structural absence is not interpreted as measured zero.

The tau table includes `tau_status`, `num_groups_required`, `num_groups_observed`, and `panel_id`. `highest_ties` lists all exactly tied maxima in sorted label order; `highest` retains the first maximum in panel order for compatibility, and `order` lists positive groups in descending order. Undefined rows have no highest-tissue annotation. The coverage table distinguishes `not_collected`, `all_excluded`, and `missing_runs`; per-gene invalid values are reflected in the linear mean and tau status. Weights are fixed across genes; a missing value makes the affected representative undefined instead of redistributing its weight.

## Outputs and migration

[Finalize](amalgkit-finalize.md) exports the linear tau inputs, coverage, run weights and a JSON definition next to the tau table. These artifacts are also generated by the per-species preparation workers. `skip_curation` uses the uncorrected expression matrix and still produces the tau tables.

The ordinary `sample_group_mean` tables remain arithmetic means on the selected transformation scale. Cross-species mean correlation/embedding plots also aggregate run expression on that scale. They are not tau inputs, and the `--tau_*` options do not change correlation filtering or those averaged plots. Tau histogram panels use the linear representative definition. The before/after QC panel's matrix stage is separate from the final tau export; its zero-restoration ordering is not changed here.

Recalculate existing tau tables, their histograms, highest-tissue assignments and downstream tau-based gene sets. Existing final run matrices and curated metadata suffice if normalization, correction and retained runs are unchanged. A log-average table alone cannot recover a linear arithmetic mean. Changing select membership or upstream normalization/correction requires rerunning those stages. Historical tau thresholds must be reassessed for the new representative definition.
