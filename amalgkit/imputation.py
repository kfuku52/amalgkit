import warnings

import numpy
import pandas

from amalgkit.linalg_utils import gram_components_require_svd_fallback as _gram_components_require_svd_fallback


IMPUTATION_STRATEGIES = ('em_pca', 'nipals', 'row_mean')


def _warn_row_mean_fallback(strategy, reason):
    warnings.warn(
        '{} imputation fell back to row-mean imputation ({}). '
        'Imputation quality may be degraded.'.format(strategy, reason)
    )


def _row_mean_impute(values):
    out = numpy.asarray(values, dtype=float).copy()
    finite = numpy.isfinite(out)
    row_sums = numpy.where(finite, out, 0.0).sum(axis=1)
    row_counts = finite.sum(axis=1)
    row_means = numpy.divide(
        row_sums,
        row_counts,
        out=numpy.zeros_like(row_sums, dtype=float),
        where=row_counts > 0,
    )
    missing_rows, missing_columns = numpy.where(~finite)
    out[missing_rows, missing_columns] = row_means[missing_rows]
    return out


def _truncated_svd_reconstruction(centered, num_pc):
    """Return the rank-``num_pc`` SVD reconstruction without a full SVD.

    Expression matrices are normally much taller than they are wide.  The
    non-zero singular vectors can therefore be recovered exactly from the
    smaller Gram matrix, avoiding computation of singular vectors that the
    imputation step immediately discards.
    """
    values = numpy.asarray(centered, dtype=float)
    if values.ndim != 2:
        raise ValueError('centered must be a two-dimensional matrix.')
    resolved_pc = min(max(0, int(num_pc)), min(values.shape))
    if resolved_pc < 1:
        raise ValueError('PCA could not resolve a principal component.')
    if values.shape[0] >= values.shape[1]:
        eigenvalues, eigenvectors = numpy.linalg.eigh(values.T @ values)
        order = numpy.argsort(eigenvalues)[::-1]
        descending = eigenvalues[order]
        if _gram_components_require_svd_fallback(
            descending,
            resolved_pc,
            values.shape,
        ):
            return _full_svd_reconstruction(values, resolved_pc)
        top = eigenvectors[:, order[:resolved_pc]]
        return (values @ top) @ top.T
    eigenvalues, eigenvectors = numpy.linalg.eigh(values @ values.T)
    order = numpy.argsort(eigenvalues)[::-1]
    descending = eigenvalues[order]
    if _gram_components_require_svd_fallback(
        descending,
        resolved_pc,
        values.shape,
    ):
        return _full_svd_reconstruction(values, resolved_pc)
    top = eigenvectors[:, order[:resolved_pc]]
    return top @ (top.T @ values)


def _full_svd_reconstruction(values, num_pc):
    left, singular_values, right = numpy.linalg.svd(values, full_matrices=False)
    return (
        left[:, :num_pc] * singular_values[:num_pc].reshape(1, -1)
    ).dot(right[:num_pc, :])


def _fit_nipals(values, num_pc, max_iter, tol):
    residual = numpy.asarray(values, dtype=float).copy()
    scores = []
    loadings = []
    max_components = min(int(num_pc), min(residual.shape))
    for _component in range(max_components):
        column_variance = numpy.var(residual, axis=0, ddof=1)
        column_variance[~numpy.isfinite(column_variance)] = -numpy.inf
        initial_column = int(numpy.argmax(column_variance))
        if not numpy.isfinite(column_variance[initial_column]):
            break
        score = residual[:, initial_column].copy()
        if (not numpy.isfinite(score).all()) or float(score @ score) <= numpy.finfo(float).eps:
            column_energy = numpy.sum(residual ** 2, axis=0)
            candidates = numpy.flatnonzero(column_energy > numpy.finfo(float).eps)
            if candidates.size == 0:
                break
            score = residual[:, int(candidates[0])].copy()
        loading = numpy.zeros((residual.shape[1],), dtype=float)
        for _iteration in range(max(20, int(max_iter))):
            denominator = float(score @ score)
            if (not numpy.isfinite(denominator)) or denominator <= numpy.finfo(float).eps:
                break
            loading = residual.T.dot(score) / denominator
            loading_norm = float(numpy.linalg.norm(loading))
            if (not numpy.isfinite(loading_norm)) or loading_norm <= numpy.finfo(float).eps:
                break
            loading = loading / loading_norm
            updated_score = residual.dot(loading)
            if not numpy.isfinite(updated_score).all():
                raise ValueError('NIPALS produced non-finite scores.')
            delta = float(numpy.max(numpy.abs(updated_score - score)))
            score = updated_score
            if delta < float(tol):
                break
        if (
            (not numpy.isfinite(score).all())
            or (not numpy.isfinite(loading).all())
            or float(numpy.sum(numpy.abs(score))) <= numpy.finfo(float).eps
            or float(numpy.sum(numpy.abs(loading))) <= numpy.finfo(float).eps
        ):
            break
        scores.append(score)
        loadings.append(loading)
        residual = residual - numpy.outer(score, loading)
    if not scores:
        raise ValueError('NIPALS could not resolve a principal component.')
    return numpy.column_stack(scores), numpy.column_stack(loadings)


def _iterative_pca_impute(values, missing_mask, num_pc, max_iter, tol, strategy, diagnostics=None):
    imputed = _row_mean_impute(values)
    if diagnostics is not None:
        diagnostics.update(converged=False, iterations=0, final_delta=None)
    for _iteration in range(int(max_iter)):
        column_means = numpy.mean(imputed, axis=0)
        centered = imputed - column_means.reshape(1, -1)
        if strategy == 'nipals':
            scores, loadings = _fit_nipals(
                values=centered,
                num_pc=num_pc,
                max_iter=max_iter,
                tol=tol,
            )
            reconstructed = scores.dot(loadings.T)
        else:
            reconstructed = _truncated_svd_reconstruction(
                centered=centered,
                num_pc=num_pc,
            )
        reconstructed = reconstructed + column_means.reshape(1, -1)
        old_values = imputed[missing_mask].copy()
        new_values = reconstructed[missing_mask]
        if not numpy.isfinite(new_values).all():
            raise ValueError('PCA imputation produced non-finite values.')
        imputed[missing_mask] = new_values
        if old_values.size == 0:
            break
        delta = float(numpy.max(numpy.abs(old_values - new_values)))
        if diagnostics is not None:
            diagnostics.update(iterations=_iteration + 1, final_delta=delta, converged=delta < float(tol))
        if (not numpy.isfinite(delta)) or delta < float(tol):
            break
    return imputed


def impute_expression(
    matrix_df,
    strategy='em_pca',
    num_pc=4,
    max_iter=50,
    tol=1e-6,
    minimum_imputed_value=None,
    return_diagnostics=False,
):
    strategy = str(strategy).strip().lower()
    if strategy not in IMPUTATION_STRATEGIES:
        raise ValueError('Unknown missing-value strategy: {}'.format(strategy))
    if int(num_pc) < 1 or int(max_iter) < 1 or not numpy.isfinite(tol) or float(tol) <= 0:
        raise ValueError('num_pc, max_iter and tol must be positive and tol must be finite.')
    diagnostics = dict(strategy=strategy, requested_rank=int(num_pc), resolved_rank=0,
                       max_iter=int(max_iter), tolerance=float(tol), converged=True,
                       iterations=0, final_delta=None, fallback=False, clipped_cells=0,
                       missing_cells=0)
    if minimum_imputed_value is not None:
        minimum_imputed_value = float(minimum_imputed_value)
        if not numpy.isfinite(minimum_imputed_value):
            raise ValueError('minimum_imputed_value must be finite.')
    numeric = matrix_df.apply(pandas.to_numeric, errors='coerce')
    if numeric.shape[0] == 0 or numeric.shape[1] == 0:
        return (numeric.copy(), diagnostics) if return_diagnostics else numeric.copy()
    values = numeric.to_numpy(dtype=float)
    missing_mask = ~numpy.isfinite(values)
    diagnostics['missing_cells'] = int(missing_mask.sum())
    if not missing_mask.any():
        return (numeric.copy(), diagnostics) if return_diagnostics else numeric.copy()
    values_for_imputation = values.copy()
    values_for_imputation[missing_mask] = numpy.nan
    if strategy == 'row_mean':
        imputed = _row_mean_impute(values_for_imputation)
    else:
        max_pc = min(values.shape[0] - 1, values.shape[1] - 1)
        if max_pc < 1:
            diagnostics.update(fallback=True, converged=False)
            _warn_row_mean_fallback(
                strategy,
                'a {}x{} matrix is too small to resolve a principal component'.format(
                    values.shape[0], values.shape[1]
                ),
            )
            imputed = _row_mean_impute(values_for_imputation)
        else:
            resolved_pc = min(max(1, int(num_pc)), max_pc)
            diagnostics['resolved_rank'] = resolved_pc
            try:
                imputed = _iterative_pca_impute(
                    values=values_for_imputation,
                    missing_mask=missing_mask,
                    num_pc=resolved_pc,
                    max_iter=max_iter,
                    tol=tol,
                    strategy=strategy,
                    diagnostics=diagnostics,
                )
            except (ValueError, numpy.linalg.LinAlgError) as exc:
                diagnostics.update(fallback=True, converged=False)
                _warn_row_mean_fallback(strategy, 'failed to converge/resolve: {}'.format(exc))
                imputed = _row_mean_impute(values_for_imputation)
    if minimum_imputed_value is not None:
        diagnostics['clipped_cells'] = int((imputed[missing_mask] < minimum_imputed_value).sum())
        imputed[missing_mask] = numpy.maximum(
            imputed[missing_mask],
            minimum_imputed_value,
        )
    imputed[~missing_mask] = values[~missing_mask]
    result = pandas.DataFrame(imputed, index=numeric.index, columns=numeric.columns)
    return (result, diagnostics) if return_diagnostics else result


__all__ = [
    'IMPUTATION_STRATEGIES',
    'impute_expression',
]
