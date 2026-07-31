import numpy
import pandas


IMPUTATION_STRATEGIES = ('em_pca', 'nipals', 'row_mean')


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


def _iterative_pca_impute(values, missing_mask, num_pc, max_iter, tol, strategy):
    imputed = _row_mean_impute(values)
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
            left, singular_values, right = numpy.linalg.svd(centered, full_matrices=False)
            resolved_pc = min(int(num_pc), singular_values.size)
            if resolved_pc < 1:
                raise ValueError('PCA could not resolve a principal component.')
            reconstructed = (
                left[:, :resolved_pc] * singular_values[:resolved_pc].reshape(1, -1)
            ).dot(right[:resolved_pc, :])
        reconstructed = reconstructed + column_means.reshape(1, -1)
        old_values = imputed[missing_mask].copy()
        new_values = reconstructed[missing_mask]
        if not numpy.isfinite(new_values).all():
            raise ValueError('PCA imputation produced non-finite values.')
        imputed[missing_mask] = new_values
        if old_values.size == 0:
            break
        delta = float(numpy.max(numpy.abs(old_values - new_values)))
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
):
    strategy = str(strategy).strip().lower()
    if strategy not in IMPUTATION_STRATEGIES:
        raise ValueError('Unknown missing-value strategy: {}'.format(strategy))
    if minimum_imputed_value is not None:
        minimum_imputed_value = float(minimum_imputed_value)
        if not numpy.isfinite(minimum_imputed_value):
            raise ValueError('minimum_imputed_value must be finite.')
    numeric = matrix_df.apply(pandas.to_numeric, errors='coerce')
    if numeric.shape[0] == 0 or numeric.shape[1] == 0:
        return numeric.copy()
    values = numeric.to_numpy(dtype=float)
    missing_mask = ~numpy.isfinite(values)
    if not missing_mask.any():
        return numeric.copy()
    values_for_imputation = values.copy()
    values_for_imputation[missing_mask] = numpy.nan
    if strategy == 'row_mean':
        imputed = _row_mean_impute(values_for_imputation)
    else:
        max_pc = min(values.shape[0] - 1, values.shape[1] - 1)
        if max_pc < 1:
            imputed = _row_mean_impute(values_for_imputation)
        else:
            resolved_pc = min(max(1, int(num_pc)), max_pc)
            try:
                imputed = _iterative_pca_impute(
                    values=values_for_imputation,
                    missing_mask=missing_mask,
                    num_pc=resolved_pc,
                    max_iter=max_iter,
                    tol=tol,
                    strategy=strategy,
                )
            except (ValueError, numpy.linalg.LinAlgError):
                imputed = _row_mean_impute(values_for_imputation)
    if minimum_imputed_value is not None:
        imputed[missing_mask] = numpy.maximum(
            imputed[missing_mask],
            minimum_imputed_value,
        )
    imputed[~missing_mask] = values[~missing_mask]
    return pandas.DataFrame(imputed, index=numeric.index, columns=numeric.columns)


__all__ = [
    'IMPUTATION_STRATEGIES',
    'impute_expression',
]
