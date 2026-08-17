"""Numerical kernels used by cross-species filtering and plotting."""

from __future__ import annotations

import numpy
import pandas

from amalgkit.imputation import impute_expression


def safe_correlation(left_values, right_values, method: str) -> float:
    """Calculate a finite-aware correlation, returning NaN when undefined."""
    if str(method).lower() == "pearson":
        try:
            left = numpy.asarray(left_values, dtype=float).reshape(-1)
            right = numpy.asarray(right_values, dtype=float).reshape(-1)
        except (TypeError, ValueError):
            left = pandas.to_numeric(pandas.Series(left_values), errors="coerce").to_numpy(dtype=float)
            right = pandas.to_numeric(pandas.Series(right_values), errors="coerce").to_numpy(dtype=float)
        valid = numpy.isfinite(left) & numpy.isfinite(right)
        if int(numpy.sum(valid)) <= 1:
            return numpy.nan
        left = left[valid]
        right = right[valid]
        left = left - numpy.mean(left)
        right = right - numpy.mean(right)
        denominator = float(numpy.sqrt(numpy.dot(left, left) * numpy.dot(right, right)))
        if denominator <= 0.0:
            return numpy.nan
        return float(numpy.dot(left, right) / denominator)
    left = pandas.to_numeric(pandas.Series(left_values), errors="coerce")
    right = pandas.to_numeric(pandas.Series(right_values), errors="coerce")
    valid = left.notna() & right.notna()
    if int(valid.sum()) <= 1 or left.loc[valid].nunique() <= 1 or right.loc[valid].nunique() <= 1:
        return numpy.nan
    try:
        return float(left.loc[valid].corr(right.loc[valid], method=method))
    except ValueError:
        return numpy.nan


def calculate_correlation_within_group(
    df_metadata: pandas.DataFrame,
    ortholog_matrix: pandas.DataFrame,
    correction_label: str,
) -> pandas.DataFrame:
    """Add within-group and best non-group correlation columns to metadata."""
    target_col = f"within_group_cor_{correction_label}"
    nongroup_col = f"max_nongroup_cor_{correction_label}"
    out = df_metadata.copy()
    out.loc[:, target_col] = numpy.nan
    out.loc[:, nongroup_col] = numpy.nan
    if ortholog_matrix.shape[1] == 0:
        return out
    kept = out.loc[out["exclusion"].eq("no"), :].copy()
    kept.loc[:, "sample_id"] = kept["species_tag"].astype(str) + "_" + kept["run"].astype(str)
    kept = kept.loc[kept["sample_id"].isin(ortholog_matrix.columns), :].copy()
    if kept.shape[0] == 0:
        return out
    group_order = list(dict.fromkeys(kept["sample_group"].astype(str).tolist()))
    ortholog_med = pandas.DataFrame(index=ortholog_matrix.index, columns=group_order, dtype=float)
    species_profiles_by_group = {}
    species_order = list(dict.fromkeys(kept["species_tag"].astype(str).tolist()))
    for sample_group in group_order:
        species_profiles = pandas.DataFrame(index=ortholog_matrix.index, dtype=float)
        for species_tag in species_order:
            sample_ids = (
                kept.loc[
                    kept["species_tag"].astype(str).eq(species_tag) & kept["sample_group"].astype(str).eq(sample_group),
                    "sample_id",
                ]
                .astype(str)
                .tolist()
            )
            sample_ids = [sample_id for sample_id in sample_ids if sample_id in ortholog_matrix.columns]
            if len(sample_ids) == 0:
                continue
            tc_species_group = ortholog_matrix.loc[:, sample_ids].apply(
                pandas.to_numeric,
                errors="coerce",
            )
            species_profiles.loc[:, species_tag] = tc_species_group.mean(axis=1, skipna=True)
        if species_profiles.shape[1] == 0:
            continue
        species_profiles_by_group[sample_group] = species_profiles
        ortholog_med.loc[:, sample_group] = species_profiles.median(axis=1, skipna=True)
    for row_idx, row in kept.iterrows():
        sample_id = str(row["sample_id"])
        sample_group = str(row["sample_group"])
        if sample_id not in ortholog_matrix.columns or sample_group not in ortholog_med.columns:
            continue
        sample_values = ortholog_matrix.loc[:, sample_id]
        within_profiles = species_profiles_by_group.get(sample_group)
        if within_profiles is None:
            within_cor = numpy.nan
        else:
            within_profiles = within_profiles.copy()
            species_tag = str(row["species_tag"])
            other_sample_ids = (
                kept.loc[
                    kept["species_tag"].astype(str).eq(species_tag)
                    & kept["sample_group"].astype(str).eq(sample_group)
                    & ~kept["sample_id"].astype(str).eq(sample_id),
                    "sample_id",
                ]
                .astype(str)
                .tolist()
            )
            other_sample_ids = [
                other_sample_id for other_sample_id in other_sample_ids if other_sample_id in ortholog_matrix.columns
            ]
            if len(other_sample_ids) == 0:
                within_profiles = within_profiles.drop(columns=[species_tag], errors="ignore")
            else:
                within_profiles.loc[:, species_tag] = (
                    ortholog_matrix.loc[:, other_sample_ids]
                    .apply(
                        pandas.to_numeric,
                        errors="coerce",
                    )
                    .mean(axis=1, skipna=True)
                )
            if within_profiles.shape[1] == 0:
                within_cor = numpy.nan
            else:
                within_reference = within_profiles.median(axis=1, skipna=True)
                within_cor = safe_correlation(sample_values, within_reference, method="pearson")
        other_groups = [group for group in group_order if group != sample_group]
        nongroup_values = [
            safe_correlation(sample_values, ortholog_med.loc[:, other_group], method="pearson")
            for other_group in other_groups
        ]
        nongroup_values = [value for value in nongroup_values if numpy.isfinite(value)]
        max_nongroup = max(nongroup_values) if len(nongroup_values) > 0 else numpy.nan
        out.loc[row_idx, target_col] = within_cor
        out.loc[row_idx, nongroup_col] = max_nongroup
    return out


def resolve_matrix_for_embedding(
    matrix_df: pandas.DataFrame,
    missing_strategy: str,
    cache: dict | None = None,
) -> pandas.DataFrame:
    """Apply one missing-value strategy, optionally caching the result."""
    cache_key = ("filled", id(matrix_df), str(missing_strategy).lower())
    if cache is not None and cache_key in cache:
        return cache[cache_key]
    if str(missing_strategy).lower() == "strict":
        out = matrix_df.dropna(axis=0, how="any").copy()
    else:
        out = impute_expression(matrix_df=matrix_df, strategy=missing_strategy)
    if cache is not None:
        cache[cache_key] = out
    return out


def resolve_correlation_matrix(
    matrix_df: pandas.DataFrame,
    missing_strategy: str = "row_mean",
    cache: dict | None = None,
) -> pandas.DataFrame:
    """Return a correlation matrix sharing the embedding imputation cache."""
    strategy_key = str(missing_strategy).lower()
    cache_key = ("correlation", id(matrix_df), strategy_key)
    if cache is not None and cache_key in cache:
        return cache[cache_key]
    filled = resolve_matrix_for_embedding(
        matrix_df,
        missing_strategy=missing_strategy,
        cache=cache,
    )
    corr = filled.corr(method="pearson")
    if cache is not None:
        cache[cache_key] = corr
        if strategy_key == "row_mean":
            cache.pop(("filled", id(matrix_df), strategy_key), None)
    return corr


def finite_correlation_block(corr_df: pandas.DataFrame) -> pandas.DataFrame:
    """Drop samples whose self-correlation is undefined and return a finite block.

    A constant sample produces an all-NaN row/column, including a NaN diagonal.
    Fabricating r=0 there corrupts PCA, MDS, and dendrograms.  Correlations
    between otherwise defined samples must remain finite; reject partial NaN
    patterns instead of silently discarding additional usable samples.
    """
    if corr_df.empty:
        return corr_df.copy()
    if corr_df.shape[0] != corr_df.shape[1] or not corr_df.index.equals(corr_df.columns):
        raise ValueError("Correlation matrix must be square with matching index and columns.")
    values = corr_df.to_numpy(dtype=float)
    defined = numpy.isfinite(numpy.diag(values))
    block = corr_df.loc[defined, defined]
    if block.size > 0 and not numpy.isfinite(block.to_numpy(dtype=float)).all():
        raise ValueError("Correlations between defined samples must be finite.")
    return block


def resolve_finite_correlation_matrix(
    matrix_df: pandas.DataFrame,
    missing_strategy: str = "row_mean",
    cache: dict | None = None,
) -> pandas.DataFrame:
    """Return the fully defined sample subset required by eigendecompositions."""
    strategy_key = str(missing_strategy).lower()
    cache_key = ("finite_correlation", id(matrix_df), strategy_key)
    if cache is not None and cache_key in cache:
        return cache[cache_key]
    corr = resolve_correlation_matrix(
        matrix_df,
        missing_strategy=missing_strategy,
        cache=cache,
    )
    finite_corr = finite_correlation_block(corr)
    if finite_corr.size > 0 and not numpy.isfinite(finite_corr.to_numpy(dtype=float)).all():
        raise ValueError("Non-degenerate sample correlations must be finite after imputation.")
    if cache is not None:
        cache[cache_key] = finite_corr
    return finite_corr


def resolve_tsne_perplexity(num_samples: int) -> int | None:
    """Choose a conservative t-SNE perplexity valid for the sample count."""
    if int(num_samples) < 4:
        return None
    max_perplexity = int((int(num_samples) - 1) // 3)
    if max_perplexity < 1:
        return None
    return min(30, max_perplexity)
