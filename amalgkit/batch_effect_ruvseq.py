import math

import numpy
import pandas
from scipy.special import gammaln, xlogy
from scipy.stats import chi2, f

from amalgkit.batch_effect_common import align_metadata_to_counts
from amalgkit.batch_effect_contract import BatchModelError, batch_backend, removal_basis, matrix_change
from amalgkit.normalization_tmm import calc_factor_quantile


RUVSEQ_SCORE_TOLERANCE = 1e-12
RUVSEQ_POISSON_ALPHA_THRESHOLD = 1e-6
# Absolute cutoff in GLM deviance-residual space, used to recognise numerical
# no-effect fits.
#
# The control residuals come from statsmodels' fit.resid_deviance. They are not
# raw-count quantities, so the cutoff must NOT be scaled by the raw count
# magnitude: that mixes units and lets the threshold grow without bound with
# library depth, eventually classifying nontrivial residual structure as
# numerical noise.
#
# Identical runs leave residuals at ~1.5e-8, the convergence-scale artifact of
# the statsmodels IRLS fit. The 1e-6 floor sits safely above that measured
# artifact. A regression test separately pins that a residual magnitude of 0.5
# remains actionable even when normalized counts are in the millions.
RUVSEQ_RESIDUAL_NOISE_FLOOR = 1e-6


def _align_metadata_to_counts(counts_df, metadata_df):
    return align_metadata_to_counts(counts_df=counts_df, metadata_df=metadata_df)


def _build_sample_group_design(aligned_metadata, sample_group_column='sample_group'):
    if sample_group_column not in aligned_metadata.columns:
        raise ValueError('Missing required metadata column: {}'.format(sample_group_column))
    sample_groups = aligned_metadata.loc[:, sample_group_column].fillna('').astype(str).str.strip()
    if (sample_groups == '').any():
        raise ValueError('sample_group contains empty values.')
    levels = sorted(sample_groups.unique().tolist())
    design = pandas.DataFrame({'Intercept': numpy.ones((len(sample_groups),), dtype=float)}, index=aligned_metadata['run'])
    for level in levels[1:]:
        design['sample_group_{}'.format(level)] = sample_groups.eq(level).astype(float).to_numpy()
    return design, sample_groups.tolist()


def _build_ruvseq_skip_output(counts_df, method, skip_reason):
    return (
        counts_df.copy(),
        pandas.DataFrame(index=counts_df.columns),
        {
            'backend': 'ruvseq',
            'method': method,
            'skip_reason': str(skip_reason),
            'stable': None,
            'corrected_run_ids': [],
            'uncorrected_run_ids': [str(run_id) for run_id in counts_df.columns],
            'resolved_ruv_k': None,
            'resolved_ruv_controls': None,
            'ruv_baseline_score': math.nan,
            'ruv_selected_score': math.nan,
            'ruv_selected_penalized_score': math.nan,
            'ruv_penalty': math.nan,
            'ruv_residual_method': 'not_run',
            'ruv_pvalue_method': 'not_run',
            'ruv_fallback_used': False,
            'ruv_fallback_reason': '',
            'ruv_nb_fallback_genes': 0,
            'ruv_anova_failure_genes': 0,
        },
    )


def compute_factor_r2(values, factor_values):
    y = pandas.to_numeric(pandas.Series(values), errors='coerce').to_numpy(dtype=float)
    factors = pandas.Series(factor_values).astype(str).str.strip().to_numpy(dtype=object)
    valid = numpy.isfinite(y) & (factors != '') & (factors != 'nan')
    if valid.sum() < 2:
        return math.nan
    y = y[valid]
    factors = factors[valid]
    levels = [level for level in sorted(numpy.unique(factors).tolist()) if level != '']
    if len(levels) <= 1:
        return math.nan
    design = numpy.ones((len(y), len(levels)), dtype=float)
    for idx, level in enumerate(levels[1:], start=1):
        design[:, idx] = (factors == level).astype(float)
    beta, _, _, _ = numpy.linalg.lstsq(design, y, rcond=None)
    fitted = design @ beta
    rss = float(numpy.sum((y - fitted) ** 2))
    sst = float(numpy.sum((y - numpy.mean(y)) ** 2))
    if (not numpy.isfinite(sst)) or (sst <= 0):
        return math.nan
    r2 = 1.0 - (rss / sst)
    return float(r2)


def _leading_pca_scores(centered, n_pc):
    values = numpy.asarray(centered, dtype=float)
    resolved_pc = min(max(0, int(n_pc)), min(values.shape))
    if resolved_pc == 0:
        return numpy.zeros((values.shape[0], 0), dtype=float)
    eigenvalues, eigenvectors = numpy.linalg.eigh(values @ values.T)
    order = numpy.argsort(eigenvalues)[::-1]
    descending = numpy.maximum(eigenvalues[order], 0.0)
    scale = float(descending[0]) if descending.size > 0 else 0.0
    tolerance = numpy.finfo(float).eps * max(values.shape) * scale
    boundary_index = resolved_pc - 1
    requires_svd = (
        (not numpy.isfinite(scale))
        or scale <= 0.0
        or descending[boundary_index] <= tolerance
    )
    if resolved_pc < descending.size:
        boundary_gap = descending[boundary_index] - descending[resolved_pc]
        requires_svd = requires_svd or boundary_gap <= tolerance
    if requires_svd:
        left, singular_values, _right = numpy.linalg.svd(values, full_matrices=False)
        return (
            left[:, :resolved_pc]
            * singular_values[:resolved_pc].reshape(1, -1)
        )
    singular_values = numpy.sqrt(descending[:resolved_pc])
    return eigenvectors[:, order[:resolved_pc]] * singular_values.reshape(1, -1)


def score_ruvseq_components(mat_df, metadata_df, n_pc=3, batch_column='bioproject', sample_group_column='sample_group'):
    if (mat_df.shape[1] < 3) or (mat_df.shape[0] < 2):
        return {'score': math.nan, 'group_score': math.nan, 'batch_score': math.nan}
    aligned_metadata = _align_metadata_to_counts(counts_df=mat_df, metadata_df=metadata_df)
    samples_by_genes = mat_df.transpose().to_numpy(dtype=float)
    centered = samples_by_genes - numpy.mean(samples_by_genes, axis=0, keepdims=True)
    try:
        pcs = _leading_pca_scores(centered, n_pc=n_pc)
    except numpy.linalg.LinAlgError:
        return {'score': math.nan, 'group_score': math.nan, 'batch_score': math.nan}
    if pcs.shape[1] == 0:
        return {'score': math.nan, 'group_score': math.nan, 'batch_score': math.nan}
    batch_values = (
        aligned_metadata.loc[:, batch_column].fillna('not_provided').astype(str).tolist()
        if batch_column in aligned_metadata.columns else
        ['not_provided'] * aligned_metadata.shape[0]
    )
    group_values = aligned_metadata.loc[:, sample_group_column].astype(str).tolist()
    batch_r2 = [compute_factor_r2(pcs[:, idx], batch_values) for idx in range(pcs.shape[1])]
    group_r2 = [compute_factor_r2(pcs[:, idx], group_values) for idx in range(pcs.shape[1])]
    batch_score = 0.0 if all(not numpy.isfinite(v) for v in batch_r2) else float(numpy.nanmean(batch_r2))
    group_score = 0.0 if all(not numpy.isfinite(v) for v in group_r2) else float(numpy.nanmean(group_r2))
    return {
        'score': group_score - batch_score,
        'group_score': group_score,
        'batch_score': batch_score,
    }


def score_ruvseq_matrix(mat_df, metadata_df, n_pc=3, batch_column='bioproject', sample_group_column='sample_group'):
    return score_ruvseq_components(
        mat_df=mat_df,
        metadata_df=metadata_df,
        n_pc=n_pc,
        batch_column=batch_column,
        sample_group_column=sample_group_column,
    )['score']


def _upperquartile_normalize(counts_df, round_counts=True):
    lib_sizes = counts_df.sum(axis=0).astype(float)
    raw_factors = calc_factor_quantile(counts_df, lib_sizes, p=0.75).replace(0, numpy.nan)
    if raw_factors.isna().any():
        raw_factors = raw_factors.fillna(1.0)
    norm_factors = raw_factors / math.exp(float(numpy.mean(numpy.log(raw_factors.to_numpy(dtype=float)))))
    normalized = counts_df.astype(float).copy()
    normalized.loc[:, :] = normalized.to_numpy(dtype=float) / norm_factors.reindex(normalized.columns).to_numpy(dtype=float)
    if round_counts:
        normalized.loc[:, :] = numpy.round(normalized.to_numpy(dtype=float))
    effective_lib_sizes = lib_sizes * norm_factors
    return normalized, norm_factors, effective_lib_sizes


def _between_lane_normalize_upper(counts_df, round_counts=True):
    quantiles = counts_df.apply(lambda col: float(numpy.quantile(col.to_numpy(dtype=float), 0.75, method='linear')), axis=0)
    mean_quantile = float(numpy.mean(quantiles.to_numpy(dtype=float)))
    if (not numpy.isfinite(mean_quantile)) or (mean_quantile == 0):
        scales = pandas.Series(numpy.ones((counts_df.shape[1],), dtype=float), index=counts_df.columns, dtype=float)
    else:
        scales = quantiles / mean_quantile
    # A sparse lane can have a 75th percentile of zero while the mean quantile
    # is positive, which would divide that column by zero and produce inf/NaN.
    # Such a lane carries no usable scale information, so it is left unscaled -
    # the same policy _upperquartile_normalize already applies to its own
    # zero/non-finite factors.
    scales = scales.astype(float)
    unusable = ~numpy.isfinite(scales.to_numpy(dtype=float)) | (scales.to_numpy(dtype=float) <= 0)
    if unusable.any():
        scales.loc[unusable] = 1.0
    normalized = counts_df.astype(float).copy()
    normalized.loc[:, :] = normalized.to_numpy(dtype=float) / scales.reindex(normalized.columns).to_numpy(dtype=float)
    if round_counts:
        normalized.loc[:, :] = numpy.round(normalized.to_numpy(dtype=float))
    return normalized, scales


def _counts_per_million(counts_df, effective_lib_sizes):
    denom = pandas.Series(effective_lib_sizes, index=counts_df.columns, dtype=float).replace(0, numpy.nan)
    cpm = counts_df.astype(float).copy()
    cpm.loc[:, :] = (cpm.to_numpy(dtype=float) / denom.reindex(cpm.columns).to_numpy(dtype=float)) * 1e6
    return cpm


def _load_statsmodels():
    try:
        import statsmodels.api as sm
    except ImportError:
        return None
    return sm


def _estimate_nb_alpha_from_poisson_fit(y, mu):
    y = numpy.asarray(y, dtype=float).reshape(-1)
    mu = numpy.asarray(mu, dtype=float).reshape(-1)
    valid = numpy.isfinite(y) & numpy.isfinite(mu) & (mu > 0)
    if valid.sum() < 2:
        return 0.0
    alpha_terms = ((y[valid] - mu[valid]) ** 2 - mu[valid]) / numpy.maximum(mu[valid] ** 2, 1e-12)
    alpha = float(numpy.mean(alpha_terms))
    if (not numpy.isfinite(alpha)) or (alpha <= 0):
        return 0.0
    return alpha


def _record_ruv_fallback(diagnostics, reason):
    diagnostics['ruv_fallback_used'] = True
    reasons = diagnostics.setdefault('_ruv_fallback_reasons', [])
    reason = str(reason).strip()
    if reason != '' and reason not in reasons:
        reasons.append(reason)


def _compute_glm_pvalues_and_residuals(counts_df, design_df, effective_lib_sizes, diagnostics=None):
    if diagnostics is None:
        diagnostics = {}
    sm = _load_statsmodels()
    if sm is None:
        raise ImportError('RUVSeq requires statsmodels.')
    from statsmodels.tools.sm_exceptions import PerfectSeparationError
    x_full = design_df.to_numpy(dtype=float)
    offset = numpy.log(pandas.Series(effective_lib_sizes, index=counts_df.columns, dtype=float).to_numpy(dtype=float))
    offset_scale = numpy.exp(offset)
    offset_scale_sum = float(numpy.sum(offset_scale))
    residuals = numpy.empty(counts_df.shape, dtype=float)
    pvalues = numpy.ones((counts_df.shape[0],), dtype=float)
    rank_full = int(numpy.linalg.matrix_rank(x_full))
    df_diff = max(1, rank_full - 1)
    for row_idx in range(counts_df.shape[0]):
        y = counts_df.iloc[row_idx, :].to_numpy(dtype=float)
        if not y.any():
            residuals[row_idx, :] = 0.0
            pvalues[row_idx] = numpy.nan
            continue
        try:
            poisson_full = sm.GLM(y, x_full, family=sm.families.Poisson(), offset=offset).fit(maxiter=100, disp=0)
        except (ValueError, FloatingPointError, numpy.linalg.LinAlgError, PerfectSeparationError) as exc:
            raise BatchModelError('ruvseq_glm_failed', str(exc), {'failed_gene': str(counts_df.index[row_idx])}) from exc
        if not poisson_full.converged:
            raise BatchModelError('ruvseq_glm_not_converged', 'Poisson fit did not converge.', {'failed_gene': str(counts_df.index[row_idx])})
        alpha = _estimate_nb_alpha_from_poisson_fit(y=y, mu=poisson_full.fittedvalues)
        fit_full = poisson_full
        y_sum = float(numpy.sum(y))
        null_mean = offset_scale * (y_sum / offset_scale_sum)
        poisson_null_llf = float(
            numpy.sum(xlogy(y, null_mean) - null_mean - gammaln(y + 1.0))
        )
        null_llf = poisson_null_llf
        if alpha > RUVSEQ_POISSON_ALPHA_THRESHOLD:
            try:
                nb_family = sm.families.NegativeBinomial(alpha=max(alpha, 1e-8))
                fit_full = sm.GLM(y, x_full, family=nb_family, offset=offset).fit(maxiter=100, disp=0)
                fit_null = sm.GLM(
                    y,
                    numpy.ones((design_df.shape[0], 1), dtype=float),
                    family=nb_family,
                    offset=offset,
                ).fit(maxiter=100, disp=0)
                null_llf = float(fit_null.llf)
            except (ValueError, FloatingPointError, numpy.linalg.LinAlgError, PerfectSeparationError) as exc:
                raise BatchModelError('ruvseq_nb_glm_failed', str(exc), {'failed_gene': str(counts_df.index[row_idx])}) from exc
            if not fit_full.converged or not fit_null.converged:
                raise BatchModelError('ruvseq_glm_not_converged', 'NB fit did not converge.', {'failed_gene': str(counts_df.index[row_idx])})
        residuals[row_idx, :] = numpy.asarray(fit_full.resid_deviance, dtype=float).reshape(-1)
        llf_stat = max(0.0, 2.0 * float(fit_full.llf - null_llf))
        pvalues[row_idx] = float(chi2.sf(llf_stat, df_diff))
        if not numpy.isfinite(pvalues[row_idx]) or not numpy.isfinite(residuals[row_idx]).all():
            raise BatchModelError('ruvseq_nonfinite_fit', 'GLM produced nonfinite diagnostics.', {'failed_gene': str(counts_df.index[row_idx])})
    residuals_df = pandas.DataFrame(residuals, index=counts_df.index, columns=counts_df.columns)
    pvalues_series = pandas.Series(pvalues, index=counts_df.index, dtype=float)
    return pvalues_series, residuals_df


def _compute_group_pvalues(seq_uq_df, sample_groups, diagnostics=None):
    if diagnostics is None:
        diagnostics = {}
    groups = pandas.Series(sample_groups, index=seq_uq_df.columns).astype(str)
    levels = [level for level in sorted(groups.unique().tolist()) if level != '']
    if len(levels) <= 1:
        return pandas.Series(numpy.ones((seq_uq_df.shape[0],), dtype=float), index=seq_uq_df.index)
    with numpy.errstate(divide='ignore', invalid='ignore'):
        log_mat = numpy.log(seq_uq_df.to_numpy(dtype=float))
    finite_log = numpy.isfinite(log_mat)
    finite_counts = numpy.sum(finite_log, axis=1)
    row_centers = numpy.divide(
        numpy.sum(numpy.where(finite_log, log_mat, 0.0), axis=1),
        finite_counts,
        out=numpy.zeros((seq_uq_df.shape[0],), dtype=float),
        where=finite_counts > 0,
    )
    centered_log_mat = numpy.where(
        finite_log,
        log_mat - row_centers.reshape(-1, 1),
        0.0,
    )
    sums = numpy.zeros((seq_uq_df.shape[0], len(levels)), dtype=float)
    sum_squares = numpy.zeros_like(sums)
    counts = numpy.zeros_like(sums)
    for level_idx, level in enumerate(levels):
        level_mask = groups.eq(level).to_numpy()
        values = centered_log_mat[:, level_mask]
        finite = finite_log[:, level_mask]
        sums[:, level_idx] = numpy.sum(values, axis=1)
        sum_squares[:, level_idx] = numpy.sum(values ** 2, axis=1)
        counts[:, level_idx] = numpy.sum(finite, axis=1)
    active = counts > 0
    active_groups = numpy.sum(active, axis=1)
    total_counts = numpy.sum(counts, axis=1)
    safe_counts = numpy.where(active, counts, 1.0)
    group_correction = numpy.where(active, (sums ** 2) / safe_counts, 0.0)
    total_sums = numpy.sum(sums, axis=1)
    total_correction = numpy.divide(
        total_sums ** 2,
        total_counts,
        out=numpy.zeros_like(total_sums),
        where=total_counts > 0,
    )
    ss_between = numpy.maximum(numpy.sum(group_correction, axis=1) - total_correction, 0.0)
    ss_within = numpy.maximum(numpy.sum(sum_squares - group_correction, axis=1), 0.0)
    df_between = active_groups - 1
    df_within = total_counts - active_groups
    valid = (active_groups > 1) & (df_within > 0)
    f_stat = numpy.full((seq_uq_df.shape[0],), numpy.nan, dtype=float)
    regular = valid & (ss_within > 0)
    f_stat[regular] = (
        ss_between[regular] / df_between[regular]
    ) / (
        ss_within[regular] / df_within[regular]
    )
    f_stat[valid & (ss_within == 0) & (ss_between > 0)] = numpy.inf
    pvalues = numpy.full((seq_uq_df.shape[0],), numpy.nan, dtype=float)
    pvalues[valid] = f.sf(f_stat[valid], df_between[valid], df_within[valid])
    return pandas.Series(pvalues, index=seq_uq_df.index, dtype=float)


def select_ruvseq_controls(
    counts_df,
    seq_uq_df,
    pvalues,
    design_df,
    mode='auto',
    top_n=1000,
    min_controls=100,
    effective_lib_sizes=None,
):
    num_genes = seq_uq_df.shape[0]
    controls = numpy.ones((num_genes,), dtype=bool)
    if str(mode).lower() == 'all':
        return controls
    if design_df.shape[1] <= 1:
        raise BatchModelError('ruvseq_controls_unidentifiable', 'Empirical controls require a biological contrast; explicitly select all or file.')
    cpm_mat = _counts_per_million(
        counts_df=counts_df,
        effective_lib_sizes=effective_lib_sizes if effective_lib_sizes is not None else counts_df.sum(axis=0),
    )
    min_samples = max(2, int(math.floor(seq_uq_df.shape[1] / 4.0)))
    is_expressed = (cpm_mat > 1).sum(axis=1).to_numpy(dtype=int) >= min_samples
    pvalues_array = pandas.to_numeric(pandas.Series(pvalues, index=seq_uq_df.index), errors='coerce').to_numpy(dtype=float)
    eligible = is_expressed & numpy.isfinite(pvalues_array)
    num_eligible = int(eligible.sum())
    if num_eligible < int(min_controls):
        raise BatchModelError('ruvseq_insufficient_controls', 'Too few eligible empirical controls.', {'ruv_eligible_controls': num_eligible})
    n_select = min(int(top_n), num_eligible)
    ord_idx = numpy.argsort(pvalues_array[eligible])[::-1]
    idx_stage1 = numpy.where(eligible)[0][ord_idx[:n_select]]
    mad_vals = []
    for row_idx in idx_stage1:
        values = seq_uq_df.iloc[row_idx, :].to_numpy(dtype=float)
        median = numpy.nanmedian(values)
        mad_vals.append(float(numpy.nanmedian(numpy.abs(values - median)) * 1.4826))
    mad_vals = numpy.asarray(mad_vals, dtype=float)
    ord_mad = numpy.argsort(mad_vals)
    keep_n = max(int(min_controls), int(math.floor(len(ord_mad) * 0.5)))
    keep_n = min(keep_n, len(ord_mad))
    if keep_n < int(min_controls):
        raise BatchModelError('ruvseq_insufficient_controls', 'Control selection did not meet min_controls.', {'ruv_eligible_controls': num_eligible})
    chosen = idx_stage1[ord_mad[:keep_n]]
    controls = numpy.zeros((num_genes,), dtype=bool)
    controls[chosen] = True
    return controls


def compute_design_residuals(seq_uq_df, design_df):
    with numpy.errstate(divide='ignore', invalid='ignore'):
        samples_by_genes = numpy.log1p(
            seq_uq_df.to_numpy(dtype=float)
        ).transpose()
    x = design_df.to_numpy(dtype=float)
    beta, _, _, _ = numpy.linalg.lstsq(x, samples_by_genes, rcond=None)
    fitted = x @ beta
    residuals = samples_by_genes - fitted
    return pandas.DataFrame(residuals.transpose(), index=seq_uq_df.index, columns=seq_uq_df.columns)


def _compute_ruvr_basis(residuals_df, controls, center=True, tolerance=1e-8):
    residuals = residuals_df.to_numpy(dtype=float).transpose()
    if center:
        residuals = residuals - numpy.mean(residuals, axis=0, keepdims=True)
    controls = numpy.asarray(controls, dtype=bool).reshape(-1)
    if controls.size != residuals.shape[1]:
        raise ValueError('controls must have one element per gene.')
    e_controls = residuals[:, controls]
    residual_scale = float(numpy.max(numpy.abs(e_controls))) if e_controls.size > 0 else 0.0
    if residual_scale <= RUVSEQ_RESIDUAL_NOISE_FLOOR:
        return numpy.zeros((residuals.shape[0], 0), dtype=float)
    eigenvectors, singular_values, _right_vectors = numpy.linalg.svd(
        e_controls,
        full_matrices=False,
    )
    positive = singular_values > float(tolerance)
    return eigenvectors[:, positive]


def ruvr_correct_counts(seq_uq_df, controls, k, residuals_df, center=True, round_counts=True, epsilon=1.0, tolerance=1e-8, is_log=False, residual_basis=None, design_matrix=None):
    x = seq_uq_df.to_numpy(dtype=float)
    if (not is_log) and numpy.any(numpy.abs(x - numpy.round(x)) > 1e-8):
        pass
    y = x.transpose() if is_log else numpy.log(x + float(epsilon)).transpose()
    controls = numpy.asarray(controls, dtype=bool).reshape(-1)
    if controls.size != x.shape[0]:
        raise ValueError('controls must have one element per gene.')
    if int(k) <= 0:
        return seq_uq_df.copy(), pandas.DataFrame(index=seq_uq_df.columns)
    basis = residual_basis
    if basis is None:
        basis = _compute_ruvr_basis(
            residuals_df=residuals_df,
            controls=controls,
            center=center,
            tolerance=tolerance,
        )
    basis = numpy.asarray(basis, dtype=float)
    if basis.ndim != 2 or basis.shape[0] != x.shape[1]:
        raise ValueError('residual_basis must have one row per sample.')
    if basis.shape[1] == 0:
        return seq_uq_df.copy(), pandas.DataFrame(index=seq_uq_df.columns)
    resolved_k = min(int(k), basis.shape[1])
    if resolved_k <= 0:
        return seq_uq_df.copy(), pandas.DataFrame(index=seq_uq_df.columns)
    w = basis[:, :resolved_k]
    # W is retained for downstream GLMs; exploration removes only the part
    # orthogonal to the explicitly protected design.
    remove = w if design_matrix is None else removal_basis(w, design_matrix)
    alpha, _, _, _ = numpy.linalg.lstsq(remove, y, rcond=None)
    corrected_y = y - (remove @ alpha)
    if is_log:
        corrected = corrected_y.transpose()
    else:
        corrected = numpy.exp(corrected_y).transpose() - float(epsilon)
        if round_counts:
            corrected = numpy.round(corrected)
            corrected[corrected < 0] = 0
    corrected_df = pandas.DataFrame(corrected, index=seq_uq_df.index, columns=seq_uq_df.columns)
    if not is_log:
        corrected_df.attrs['postprocessing'] = [matrix_change(numpy.exp(corrected_y).T - float(epsilon), corrected, 'round_and_clip', 'upper_quartile_counts')]
    w_df = pandas.DataFrame(
        w,
        index=seq_uq_df.columns,
        columns=['W_{}'.format(i + 1) for i in range(w.shape[1])],
    )
    return corrected_df, w_df


def resolve_ruvseq_k_and_matrix(
    seq_uq_df,
    controls,
    residuals_df,
    metadata_df,
    k_setting='auto',
    k_max=5,
    batch_column='bioproject',
    sample_group_column='sample_group',
    design_matrix=None,
):
    if str(k_setting) != 'auto':
        selected_k = int(k_setting)
        if selected_k < 0:
            selected_k = 1
        if selected_k == 0:
            corrected_df = seq_uq_df.copy()
            w_df = pandas.DataFrame(index=seq_uq_df.columns)
        else:
            residual_basis = _compute_ruvr_basis(
                residuals_df=residuals_df,
                controls=controls,
            )
            corrected_df, w_df = ruvr_correct_counts(
                seq_uq_df,
                controls,
                selected_k,
                residuals_df,
                residual_basis=residual_basis,
                design_matrix=design_matrix,
            )
        resolved_k = int(w_df.shape[1])
        comp = score_ruvseq_components(
            mat_df=corrected_df,
            metadata_df=metadata_df,
            batch_column=batch_column,
            sample_group_column=sample_group_column,
        )
        return {
            'k': resolved_k,
            'matrix': corrected_df,
            'w': w_df,
            'score': comp['score'],
            'group_score': comp['group_score'],
            'batch_score': comp['batch_score'],
            'penalized_score': comp['score'],
            'baseline_score': math.nan,
            'baseline_group_score': math.nan,
            'penalty': 0.0,
        }
    residual_basis = _compute_ruvr_basis(
        residuals_df=residuals_df,
        controls=controls,
    )
    max_k = max(1, int(k_max))
    design_rank = 1 if design_matrix is None else numpy.linalg.matrix_rank(design_matrix)
    max_allowed = max(0, seq_uq_df.shape[1] - design_rank - 1)
    max_k = min(max_k, max_allowed)
    baseline_comp = score_ruvseq_components(
        mat_df=seq_uq_df,
        metadata_df=metadata_df,
        batch_column=batch_column,
        sample_group_column=sample_group_column,
    )
    baseline_score = baseline_comp['score']
    if not numpy.isfinite(baseline_score):
        raise BatchModelError('ruvseq_auto_score_unavailable', 'Cannot compare automatic k candidates to a finite baseline.')
    baseline_group_score = baseline_comp['group_score']
    best_k = 0
    best_score = baseline_score
    best_group_score = baseline_group_score
    best_batch_score = baseline_comp['batch_score']
    best_penalty = 0.0
    best_penalized_score = baseline_score
    best_matrix = seq_uq_df.copy()
    best_w = pandas.DataFrame(index=seq_uq_df.columns)
    for k in range(1, max_k + 1):
        corrected_df, w_df = ruvr_correct_counts(
            seq_uq_df,
            controls,
            k,
            residuals_df,
            residual_basis=residual_basis,
            design_matrix=design_matrix,
        )
        resolved_k = int(w_df.shape[1])
        if resolved_k <= 0:
            continue
        comp = score_ruvseq_components(
            mat_df=corrected_df,
            metadata_df=metadata_df,
            batch_column=batch_column,
            sample_group_column=sample_group_column,
        )
        score = comp['score']
        group_score = comp['group_score']
        if numpy.isfinite(baseline_group_score) and numpy.isfinite(group_score):
            penalty = max(0.0, float(baseline_group_score - group_score)) * 2.0
        elif numpy.isfinite(baseline_group_score) and (not numpy.isfinite(group_score)):
            penalty = 1.0
        else:
            penalty = 0.0
        penalized_score = score - penalty if numpy.isfinite(score) else math.nan
        if not numpy.isfinite(best_penalized_score):
            best_k = resolved_k
            best_score = score
            best_group_score = group_score
            best_batch_score = comp['batch_score']
            best_penalty = penalty
            best_penalized_score = penalized_score
            best_matrix = corrected_df
            best_w = w_df
            continue
        if not numpy.isfinite(penalized_score):
            continue
        if penalized_score > (best_penalized_score + RUVSEQ_SCORE_TOLERANCE):
            best_k = resolved_k
            best_score = score
            best_group_score = group_score
            best_batch_score = comp['batch_score']
            best_penalty = penalty
            best_penalized_score = penalized_score
            best_matrix = corrected_df
            best_w = w_df
            continue
        if abs(float(penalized_score - best_penalized_score)) <= RUVSEQ_SCORE_TOLERANCE:
            should_replace = False
            if (
                (int(best_k) > 0)
                and (int(resolved_k) > 0)
                and (int(resolved_k) < int(best_k))
            ):
                should_replace = True
            if not should_replace:
                continue
            best_k = resolved_k
            best_score = score
            best_group_score = group_score
            best_batch_score = comp['batch_score']
            best_penalty = penalty
            best_penalized_score = penalized_score
            best_matrix = corrected_df
            best_w = w_df
    return {
        'k': best_k,
        'matrix': best_matrix,
        'w': best_w,
        'score': best_score,
        'group_score': best_group_score,
        'batch_score': best_batch_score,
        'penalized_score': best_penalized_score,
        'baseline_score': baseline_score,
        'baseline_group_score': baseline_group_score,
        'penalty': best_penalty,
    }


@batch_backend('ruvseq')
def run_ruvseq_backend(
    counts_df,
    metadata_df,
    control_mode='auto',
    k_setting='auto',
    k_max=5,
    top_n=1000,
    min_controls=100,
    batch_column='bioproject',
    sample_group_column='sample_group',
    k_selection='manual',
    control_gene_ids=None,
    protected_design=None,
):
    method = 'manual' if str(k_setting) != 'auto' else 'auto'
    if counts_df.shape[0] == 0:
        return _build_ruvseq_skip_output(
            counts_df=counts_df,
            method=method,
            skip_reason='no_expressed_genes',
        )
    if control_mode not in {'auto', 'empirical', 'all', 'file'}:
        raise ValueError('Unknown RUV control mode.')
    if control_gene_ids is not None and control_mode != 'file':
        raise ValueError('A control gene file requires RUV control mode file.')
    if k_selection not in {'manual', 'legacy'} or int(k_max) < 0 or int(min_controls) < 2 or int(top_n) < 1:
        raise ValueError('Invalid RUV selection options.')
    if str(k_setting) != 'auto' and int(k_setting) < 0:
        raise ValueError('RUV k must be nonnegative.')
    if str(k_setting) == '0':
        output = _build_ruvseq_skip_output(counts_df, method, 'ruvseq_k_zero')
        output[-1]['resolved_ruv_k'] = 0
        return output
    if str(k_setting) == 'auto' and k_selection != 'legacy':
        raise BatchModelError('ruvseq_auto_not_calibrated', 'Specify --ruvseq_k INT; uncalibrated auto requires --ruvseq_k_selection legacy.')
    aligned_metadata = metadata_df
    design_df = protected_design.matrix
    max_k = min(int(k_max), max(0, counts_df.shape[1] - design_df.shape[1] - 1))
    if max_k == 0 or (str(k_setting) != 'auto' and int(k_setting) > max_k):
        raise BatchModelError('ruvseq_insufficient_df', 'Requested k leaves no residual degrees of freedom.')
    if (counts_df.sum(axis=0) <= 0).any():
        raise BatchModelError('ruvseq_zero_library', 'Cannot estimate offsets for a zero-count library.')
    # The count GLM uses the observed counts. The pseudo-count belongs only
    # to the log transformation used for the exploratory corrected matrix.
    _edge_uq_df, _uq_factors, effective_lib_sizes = _upperquartile_normalize(counts_df, round_counts=False)
    seq_uq_df, seq_uq_scales = _between_lane_normalize_upper(counts_df, round_counts=True)
    diagnostics = {
        'ruv_residual_method': 'glm_deviance',
        'ruv_pvalue_method': 'glm_lrt',
        'ruv_fallback_used': False,
        'ruv_nb_fallback_genes': 0,
        'ruv_anova_failure_genes': 0,
    }
    pvalues, residuals_df = _compute_glm_pvalues_and_residuals(
        counts_df=counts_df,
        design_df=design_df,
        effective_lib_sizes=effective_lib_sizes,
        diagnostics=diagnostics,
    )
    if (pvalues is None) or (residuals_df is None):
        raise BatchModelError('ruvseq_glm_failed', 'GLM did not produce usable residuals and p-values.', diagnostics)
    if control_mode == 'file':
        if not control_gene_ids:
            raise ValueError('RUV file mode requires control gene IDs.')
        controls = counts_df.index.isin(control_gene_ids)
    else:
        controls = select_ruvseq_controls(
            counts_df=counts_df, seq_uq_df=seq_uq_df, pvalues=pvalues,
            design_df=design_df, mode=control_mode, top_n=top_n,
            min_controls=min_controls, effective_lib_sizes=effective_lib_sizes,
        )
    if int(controls.sum()) < int(min_controls):
        raise BatchModelError('ruvseq_insufficient_controls', 'Selected controls do not meet min_controls.',
                              {'resolved_ruv_controls': int(controls.sum()), 'ruv_control_mode': control_mode})
    resolved = resolve_ruvseq_k_and_matrix(
        seq_uq_df=seq_uq_df,
        controls=controls,
        residuals_df=residuals_df,
        metadata_df=aligned_metadata,
        k_setting=k_setting,
        k_max=max_k,
        batch_column=batch_column,
        sample_group_column=sample_group_column,
        design_matrix=design_df.to_numpy(dtype=float),
    )
    postprocessing = list(resolved['matrix'].attrs.get('postprocessing', []))
    if str(k_setting) != 'auto' and int(resolved['k']) < int(k_setting):
        raise BatchModelError('ruvseq_degenerate_factors', 'Requested RUV dimensions could not be estimated.',
                              {'requested_ruv_k': int(k_setting), 'resolved_ruv_k': int(resolved['k'])})
    corrected_df = counts_df.copy()
    corrected_run_ids = []
    uncorrected_run_ids = [str(run_id) for run_id in counts_df.columns]
    skip_reason = 'ruvseq_k_zero'
    if int(resolved['k']) > 0:
        # The correction ran on the between-lane upper-quartile normalized
        # matrix; multiply back by the per-lane scale factors so the returned
        # matrix is on the original count scale (consistent with the k=0 path
        # and with the raw-count nonexpressed genes concatenated downstream in
        # per_species_finalize_python). Previously the corrected matrix was
        # returned on the normalized scale, silently rescaling every corrected
        # count and mixing two scales in the final table.
        normalized_corrected = resolved['matrix'].reindex(index=counts_df.index, columns=counts_df.columns)
        scale_values = seq_uq_scales.reindex(counts_df.columns).to_numpy(dtype=float)
        corrected_values = normalized_corrected.to_numpy(dtype=float) * scale_values[numpy.newaxis, :]
        before_rounding = corrected_values.copy()
        corrected_values = numpy.round(corrected_values)
        corrected_values[corrected_values < 0] = 0
        postprocessing.append(matrix_change(before_rounding, corrected_values, 'round_and_clip', 'original_counts'))
        corrected_df = pandas.DataFrame(
            corrected_values,
            index=counts_df.index,
            columns=counts_df.columns,
        )
        corrected_run_ids = [str(run_id) for run_id in counts_df.columns]
        uncorrected_run_ids = []
        skip_reason = ''
    summary = {
        'backend': 'ruvseq',
        'method': method,
        'skip_reason': skip_reason,
        'stable': None,
        'corrected_run_ids': corrected_run_ids,
        'uncorrected_run_ids': uncorrected_run_ids,
        'resolved_ruv_k': int(resolved['k']),
        'resolved_ruv_controls': int(controls.sum()),
        'ruv_control_mode': control_mode,
        'ruv_control_gene_ids': counts_df.index[controls].tolist(),
        'ruv_missing_control_gene_ids': [] if control_gene_ids is None else [
            gene for gene in control_gene_ids if gene not in counts_df.index],
        'ruv_glm_model': 'poisson_or_moment_nb',
        'ruv_input_pseudocount': 0.0,
        'ruv_log_pseudocount': 1.0,
        'ruv_effective_library_sizes': effective_lib_sizes.to_dict(),
        'ruv_between_lane_scales': seq_uq_scales.to_dict(),
        'postprocessing': postprocessing,
        'k_selection': k_selection,
        'ruv_baseline_score': resolved['baseline_score'],
        'ruv_selected_score': resolved['score'],
        'ruv_selected_penalized_score': resolved['penalized_score'],
        'ruv_penalty': resolved['penalty'],
        'ruv_residual_method': diagnostics['ruv_residual_method'],
        'ruv_pvalue_method': diagnostics['ruv_pvalue_method'],
        'ruv_fallback_used': bool(diagnostics['ruv_fallback_used']),
        'ruv_fallback_reason': '|'.join(diagnostics.get('_ruv_fallback_reasons', [])),
        'ruv_nb_fallback_genes': int(diagnostics['ruv_nb_fallback_genes']),
        'ruv_anova_failure_genes': int(diagnostics['ruv_anova_failure_genes']),
    }
    return corrected_df, resolved['w'], summary


__all__ = [
    'compute_factor_r2',
    'compute_design_residuals',
    'resolve_ruvseq_k_and_matrix',
    'run_ruvseq_backend',
    'ruvr_correct_counts',
    'score_ruvseq_components',
    'score_ruvseq_matrix',
    'select_ruvseq_controls',
]
