from dataclasses import dataclass

import numpy
import pandas

from amalgkit.batch_effect_common import BatchEffectResult, normalize_run_ids


EPS = 1e-8
MAX_EXPONENT = 50.0


@dataclass
class _LatentFit:
    corrected_df: pandas.DataFrame
    latent_df: pandas.DataFrame
    objective: float
    iterations: int
    converged: bool
    actual_k: int


def _build_design_matrix(metadata_df, sample_group_column='sample_group'):
    sample_count = metadata_df.shape[0]
    intercept = numpy.ones((sample_count, 1), dtype=float)
    if sample_group_column not in metadata_df.columns:
        return intercept, ['Intercept']
    sample_groups = metadata_df.loc[:, sample_group_column].fillna('unknown').astype(str).str.strip()
    sample_groups = sample_groups.mask(sample_groups.eq(''), 'unknown')
    dummies = pandas.get_dummies(sample_groups, drop_first=True, dtype=float)
    if dummies.shape[1] == 0:
        return intercept, ['Intercept']
    matrix = numpy.column_stack([intercept, dummies.to_numpy(dtype=float)])
    return matrix, ['Intercept'] + list(dummies.columns)


def _orthonormal_basis(matrix):
    if matrix.size == 0:
        return numpy.zeros((matrix.shape[0], 0), dtype=float)
    q, r = numpy.linalg.qr(matrix)
    if r.ndim != 2 or r.shape[0] == 0 or r.shape[1] == 0:
        return numpy.zeros((matrix.shape[0], 0), dtype=float)
    diag = numpy.abs(numpy.diag(r))
    rank = int(numpy.sum(diag > 1e-8))
    if rank <= 0:
        return numpy.zeros((matrix.shape[0], 0), dtype=float)
    return q[:, :rank]


def _project_orthogonal_to_design(matrix, design_matrix):
    basis = _orthonormal_basis(design_matrix)
    if basis.shape[1] == 0:
        return matrix.copy()
    return matrix - matrix.dot(basis).dot(basis.T)


def _orthonormalize_latent(latent_matrix, design_matrix):
    if latent_matrix.size == 0:
        return numpy.zeros((design_matrix.shape[0], 0), dtype=float)
    basis = _orthonormal_basis(design_matrix)
    if basis.shape[1] == 0:
        projected = latent_matrix.copy()
    else:
        projected = latent_matrix - basis.dot(basis.T).dot(latent_matrix)
    if numpy.allclose(projected, 0.0):
        return numpy.zeros((design_matrix.shape[0], 0), dtype=float)
    q, r = numpy.linalg.qr(projected)
    if r.ndim != 2 or r.shape[0] == 0 or r.shape[1] == 0:
        return numpy.zeros((design_matrix.shape[0], 0), dtype=float)
    diag = numpy.abs(numpy.diag(r))
    rank = int(numpy.sum(diag > 1e-8))
    if rank <= 0:
        return numpy.zeros((design_matrix.shape[0], 0), dtype=float)
    return q[:, :rank]


def _subspace_distance(left, right):
    if (left.shape[1] == 0) and (right.shape[1] == 0):
        return 0.0
    left_proj = left.dot(left.T) if left.shape[1] > 0 else numpy.zeros((left.shape[0], left.shape[0]), dtype=float)
    right_proj = right.dot(right.T) if right.shape[1] > 0 else numpy.zeros((right.shape[0], right.shape[0]), dtype=float)
    return float(numpy.linalg.norm(left_proj - right_proj, ord='fro'))


def _fit_linear_effects(response_matrix, design_matrix):
    coefficients = numpy.linalg.pinv(design_matrix).dot(response_matrix.T)
    fitted = design_matrix.dot(coefficients).T
    return coefficients, fitted


def _estimate_gene_dispersion(normalized_counts):
    if normalized_counts.shape[1] <= 1:
        return numpy.zeros((normalized_counts.shape[0],), dtype=float)
    means = normalized_counts.mean(axis=1)
    variances = normalized_counts.var(axis=1, ddof=1)
    alpha = (variances - means) / numpy.maximum(means ** 2, EPS)
    alpha[~numpy.isfinite(alpha)] = 0.0
    alpha = numpy.maximum(alpha, 0.0)
    return alpha


def _weighted_residuals(residual_matrix, normalized_counts, fitted_matrix, family):
    if family == 'poisson':
        gene_weights = numpy.ones((residual_matrix.shape[0], 1), dtype=float)
        return residual_matrix, gene_weights
    alpha = _estimate_gene_dispersion(normalized_counts)
    fitted_mean = numpy.exp(numpy.clip(fitted_matrix, -MAX_EXPONENT, MAX_EXPONENT)).mean(axis=1)
    gene_weights = 1.0 / numpy.sqrt(1.0 + alpha * numpy.maximum(fitted_mean, EPS))
    gene_weights = gene_weights.reshape(-1, 1)
    return residual_matrix * gene_weights, gene_weights


def _update_latent_factors(weighted_residual_matrix, design_matrix, k):
    sample_count = design_matrix.shape[0]
    if k <= 0:
        return numpy.zeros((sample_count, 0), dtype=float)
    projected = _project_orthogonal_to_design(weighted_residual_matrix, design_matrix)
    if projected.size == 0 or numpy.allclose(projected, 0.0):
        return numpy.zeros((sample_count, 0), dtype=float)
    _u, singular_values, vt = numpy.linalg.svd(projected, full_matrices=False)
    rank = int(numpy.sum(singular_values > 1e-8))
    if rank <= 0:
        return numpy.zeros((sample_count, 0), dtype=float)
    latent = vt[:min(k, rank), :].T
    return _orthonormalize_latent(latent, design_matrix)


def _objective_value(weighted_residual_matrix):
    if weighted_residual_matrix.size == 0:
        return 0.0
    return float(numpy.mean(weighted_residual_matrix ** 2))


def _prepare_input_arrays(counts_df):
    observed_runs = normalize_run_ids(counts_df.columns)
    counts = counts_df.loc[:, observed_runs].to_numpy(dtype=float)
    counts = numpy.maximum(counts, 0.0)
    library_sizes = counts.sum(axis=0)
    if numpy.allclose(library_sizes, 0.0):
        library_sizes = numpy.ones_like(library_sizes)
    positive_library_sizes = library_sizes[library_sizes > 0]
    scale = float(numpy.median(positive_library_sizes)) if positive_library_sizes.size > 0 else 1.0
    if scale <= 0:
        scale = 1.0
    offsets = numpy.log(numpy.maximum(library_sizes / scale, EPS))
    normalized_counts = counts / numpy.exp(offsets).reshape(1, -1)
    response = numpy.log(normalized_counts + 0.5)
    return observed_runs, counts, normalized_counts, offsets, response


def _fit_latent_model(
    counts_df,
    design_matrix,
    response_matrix,
    normalized_counts,
    offsets,
    k,
    family,
    max_iter,
    tol,
):
    run_ids = list(counts_df.columns)
    if k <= 0:
        coefficients, fitted = _fit_linear_effects(response_matrix, design_matrix)
        corrected_log_counts = design_matrix.dot(coefficients).T + offsets.reshape(1, -1)
        corrected = numpy.exp(numpy.clip(corrected_log_counts, -MAX_EXPONENT, MAX_EXPONENT))
        corrected_df = pandas.DataFrame(corrected, index=counts_df.index, columns=run_ids)
        latent_df = pandas.DataFrame(index=run_ids)
        objective = _objective_value(response_matrix - fitted)
        return _LatentFit(corrected_df, latent_df, objective, 0, True, 0)

    weighted_residuals, _gene_weights = _weighted_residuals(
        _project_orthogonal_to_design(response_matrix, design_matrix),
        normalized_counts,
        response_matrix,
        family,
    )
    latent = _update_latent_factors(weighted_residuals, design_matrix, k)
    if latent.shape[1] == 0:
        return _fit_latent_model(
            counts_df=counts_df,
            design_matrix=design_matrix,
            response_matrix=response_matrix,
            normalized_counts=normalized_counts,
            offsets=offsets,
            k=0,
            family=family,
            max_iter=max_iter,
            tol=tol,
        )

    converged = False
    iterations = 0
    objective = numpy.nan
    for iteration in range(int(max_iter)):
        design_augmented = numpy.column_stack([design_matrix, latent])
        coefficients, fitted = _fit_linear_effects(response_matrix, design_augmented)
        residuals = response_matrix - fitted
        weighted_residuals, _gene_weights = _weighted_residuals(
            residuals,
            normalized_counts=normalized_counts,
            fitted_matrix=fitted,
            family=family,
        )
        updated_latent = _update_latent_factors(weighted_residuals, design_matrix, latent.shape[1])
        iterations = iteration + 1
        objective = _objective_value(weighted_residuals)
        if _subspace_distance(latent, updated_latent) <= float(tol):
            latent = updated_latent
            converged = True
            break
        latent = updated_latent
        if latent.shape[1] == 0:
            break

    actual_k = int(latent.shape[1])
    if actual_k <= 0:
        return _fit_latent_model(
            counts_df=counts_df,
            design_matrix=design_matrix,
            response_matrix=response_matrix,
            normalized_counts=normalized_counts,
            offsets=offsets,
            k=0,
            family=family,
            max_iter=max_iter,
            tol=tol,
        )

    design_augmented = numpy.column_stack([design_matrix, latent])
    coefficients, fitted = _fit_linear_effects(response_matrix, design_augmented)
    design_only = design_matrix.dot(coefficients[:design_matrix.shape[1], :]).T
    corrected_log_counts = design_only + offsets.reshape(1, -1)
    corrected = numpy.exp(numpy.clip(corrected_log_counts, -MAX_EXPONENT, MAX_EXPONENT))
    corrected_df = pandas.DataFrame(corrected, index=counts_df.index, columns=run_ids)
    latent_df = pandas.DataFrame(
        latent,
        index=run_ids,
        columns=['latent_{}'.format(idx + 1) for idx in range(actual_k)],
    )
    if numpy.isnan(objective):
        weighted_residuals, _gene_weights = _weighted_residuals(
            response_matrix - fitted,
            normalized_counts=normalized_counts,
            fitted_matrix=fitted,
            family=family,
        )
        objective = _objective_value(weighted_residuals)
    return _LatentFit(corrected_df, latent_df, objective, iterations, converged, actual_k)


def _resolve_manual_k(k_setting, max_k):
    if str(k_setting).lower() == 'auto':
        return None
    try:
        requested = int(k_setting)
    except (TypeError, ValueError):
        requested = 0
    return max(0, min(requested, max_k))


def _resolve_auto_k(response_matrix, normalized_counts, design_matrix, family, max_allowed_k):
    if max_allowed_k <= 0:
        return 0
    residual_matrix = _project_orthogonal_to_design(response_matrix, design_matrix)
    weighted_residuals, _gene_weights = _weighted_residuals(
        residual_matrix,
        normalized_counts=normalized_counts,
        fitted_matrix=response_matrix,
        family=family,
    )
    if weighted_residuals.size == 0 or numpy.allclose(weighted_residuals, 0.0):
        return 0
    _u, singular_values, _vt = numpy.linalg.svd(weighted_residuals, full_matrices=False)
    energy = singular_values ** 2
    total_energy = float(energy.sum())
    if total_energy <= EPS:
        return 0
    leading_share = float(energy[0] / total_energy)
    if leading_share < 0.10:
        return 0
    usable_energy = energy[:max_allowed_k]
    cumulative = numpy.cumsum(usable_energy) / total_energy
    resolved = int(numpy.searchsorted(cumulative, 0.60) + 1)
    return max(1, min(resolved, max_allowed_k))


def run_latent_glm_backend(
    counts_df,
    metadata_df,
    family='nb',
    k_setting='auto',
    k_max=5,
    sample_group_column='sample_group',
    max_iter=200,
    tol=1e-5,
):
    if counts_df.shape[1] == 0:
        raise ValueError('latent_glm backend requires at least one sample.')
    if family not in {'poisson', 'nb'}:
        raise ValueError('Unsupported latent_glm family: {}'.format(family))
    counts = counts_df.copy()
    counts.index = counts.index.map(str)
    counts.columns = counts.columns.map(str)
    metadata = metadata_df.copy()
    if 'run' not in metadata.columns:
        raise ValueError('latent_glm backend requires a "run" column in metadata.')
    metadata.loc[:, 'run'] = metadata.loc[:, 'run'].astype(str)
    metadata = metadata.set_index('run', drop=False).loc[list(counts.columns), :].reset_index(drop=True)

    run_ids, _counts_array, normalized_counts, offsets, response_matrix = _prepare_input_arrays(counts)
    design_matrix, _design_columns = _build_design_matrix(metadata, sample_group_column=sample_group_column)
    max_allowed_k = min(int(k_max), max(0, counts.shape[1] - _orthonormal_basis(design_matrix).shape[1]))
    manual_k = _resolve_manual_k(k_setting=k_setting, max_k=max_allowed_k)

    if manual_k is not None:
        fit = _fit_latent_model(
            counts_df=counts,
            design_matrix=design_matrix,
            response_matrix=response_matrix,
            normalized_counts=normalized_counts,
            offsets=offsets,
            k=manual_k,
            family=family,
            max_iter=max_iter,
            tol=tol,
        )
        resolved_k = int(fit.actual_k)
        method = 'manual'
    else:
        resolved_auto_k = _resolve_auto_k(
            response_matrix=response_matrix,
            normalized_counts=normalized_counts,
            design_matrix=design_matrix,
            family=family,
            max_allowed_k=max_allowed_k,
        )
        fit = _fit_latent_model(
            counts_df=counts,
            design_matrix=design_matrix,
            response_matrix=response_matrix,
            normalized_counts=normalized_counts,
            offsets=offsets,
            k=resolved_auto_k,
            family=family,
            max_iter=max_iter,
            tol=tol,
        )
        resolved_k = int(fit.actual_k)
        method = 'auto'

    if resolved_k <= 0:
        summary = BatchEffectResult(
            backend='latent_glm',
            method=method,
            skip_reason='latent_k_zero',
            stable=True,
            corrected_run_ids=[],
            uncorrected_run_ids=list(run_ids),
            resolved_latent_k=0,
            latent_family=family,
            latent_iterations=int(fit.iterations),
            latent_objective=float(fit.objective),
            negative_values_before_clip=0,
            negative_values_after_clip=0,
            extra={'latent_converged': True},
        ).to_jsonable()
        return counts.copy(), pandas.DataFrame(index=run_ids), summary

    corrected_df = fit.corrected_df.loc[:, run_ids]
    corrected_df = corrected_df.clip(lower=0.0)
    summary = BatchEffectResult(
        backend='latent_glm',
        method=method,
        skip_reason='',
        stable=bool(fit.converged),
        corrected_run_ids=list(run_ids),
        uncorrected_run_ids=[],
        resolved_latent_k=int(resolved_k),
        latent_family=family,
        latent_iterations=int(fit.iterations),
        latent_objective=float(fit.objective),
        negative_values_before_clip=0,
        negative_values_after_clip=0,
        extra={'latent_converged': bool(fit.converged)},
    ).to_jsonable()
    return corrected_df, fit.latent_df, summary


__all__ = [
    'run_latent_glm_backend',
]
