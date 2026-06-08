import math

import numpy
import pandas
from scipy.stats import chi2, f_oneway

from amalgkit.normalization_tmm import calc_factor_quantile


RUVSEQ_SCORE_TOLERANCE = 1e-12
RUVSEQ_POISSON_ALPHA_THRESHOLD = 1e-6


def _align_metadata_to_counts(counts_df, metadata_df):
    if 'run' not in metadata_df.columns:
        raise ValueError('Missing required metadata column: run')
    run_indexed = metadata_df.copy()
    run_indexed['run'] = run_indexed['run'].astype(str)
    missing_runs = [run_id for run_id in counts_df.columns if run_id not in set(run_indexed['run'])]
    if missing_runs:
        raise ValueError('Metadata is missing rows for runs: {}'.format(', '.join(missing_runs)))
    aligned = run_indexed.drop_duplicates(subset=['run'], keep='first').set_index('run')
    return aligned.loc[list(counts_df.columns), :].reset_index()


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


def score_ruvseq_components(mat_df, metadata_df, n_pc=3, batch_column='bioproject', sample_group_column='sample_group'):
    if (mat_df.shape[1] < 3) or (mat_df.shape[0] < 2):
        return {'score': math.nan, 'group_score': math.nan, 'batch_score': math.nan}
    aligned_metadata = _align_metadata_to_counts(counts_df=mat_df, metadata_df=metadata_df)
    samples_by_genes = mat_df.transpose().to_numpy(dtype=float)
    centered = samples_by_genes - numpy.mean(samples_by_genes, axis=0, keepdims=True)
    try:
        u, s, _vh = numpy.linalg.svd(centered, full_matrices=False)
    except numpy.linalg.LinAlgError:
        return {'score': math.nan, 'group_score': math.nan, 'batch_score': math.nan}
    if s.size == 0:
        return {'score': math.nan, 'group_score': math.nan, 'batch_score': math.nan}
    pcs = u[:, : min(int(n_pc), u.shape[1])] * s[: min(int(n_pc), s.size)].reshape(1, -1)
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


def _compute_glm_pvalues_and_residuals(counts_df, design_df, effective_lib_sizes):
    sm = _load_statsmodels()
    if sm is None:
        return None, None
    x_full = design_df.to_numpy(dtype=float)
    x_null = numpy.ones((design_df.shape[0], 1), dtype=float)
    offset = numpy.log(pandas.Series(effective_lib_sizes, index=counts_df.columns, dtype=float).to_numpy(dtype=float))
    residuals = numpy.empty(counts_df.shape, dtype=float)
    pvalues = numpy.ones((counts_df.shape[0],), dtype=float)
    rank_full = int(numpy.linalg.matrix_rank(x_full))
    rank_null = int(numpy.linalg.matrix_rank(x_null))
    df_diff = max(1, rank_full - rank_null)
    for row_idx in range(counts_df.shape[0]):
        y = counts_df.iloc[row_idx, :].to_numpy(dtype=float)
        try:
            poisson_full = sm.GLM(y, x_full, family=sm.families.Poisson(), offset=offset).fit(maxiter=100, disp=0)
            poisson_null = sm.GLM(y, x_null, family=sm.families.Poisson(), offset=offset).fit(maxiter=100, disp=0)
        except Exception:
            return None, None
        alpha = _estimate_nb_alpha_from_poisson_fit(y=y, mu=poisson_full.fittedvalues)
        fit_full = poisson_full
        fit_null = poisson_null
        if alpha > RUVSEQ_POISSON_ALPHA_THRESHOLD:
            try:
                nb_family = sm.families.NegativeBinomial(alpha=max(alpha, 1e-8))
                fit_full = sm.GLM(y, x_full, family=nb_family, offset=offset).fit(maxiter=100, disp=0)
                fit_null = sm.GLM(y, x_null, family=nb_family, offset=offset).fit(maxiter=100, disp=0)
            except Exception:
                fit_full = poisson_full
                fit_null = poisson_null
        residuals[row_idx, :] = numpy.asarray(fit_full.resid_deviance, dtype=float).reshape(-1)
        llf_stat = max(0.0, 2.0 * float(fit_full.llf - fit_null.llf))
        pvalues[row_idx] = float(chi2.sf(llf_stat, df_diff))
    residuals_df = pandas.DataFrame(residuals, index=counts_df.index, columns=counts_df.columns)
    pvalues_series = pandas.Series(pvalues, index=counts_df.index, dtype=float)
    return pvalues_series, residuals_df


def _compute_group_pvalues(seq_uq_df, sample_groups):
    groups = pandas.Series(sample_groups, index=seq_uq_df.columns).astype(str)
    levels = [level for level in sorted(groups.unique().tolist()) if level != '']
    if len(levels) <= 1:
        return pandas.Series(numpy.ones((seq_uq_df.shape[0],), dtype=float), index=seq_uq_df.index)
    pvalues = numpy.ones((seq_uq_df.shape[0],), dtype=float)
    log_mat = numpy.log(seq_uq_df.to_numpy(dtype=float))
    for idx in range(seq_uq_df.shape[0]):
        arrays = []
        for level in levels:
            level_values = log_mat[idx, groups.eq(level).to_numpy()]
            level_values = level_values[numpy.isfinite(level_values)]
            if level_values.size == 0:
                continue
            arrays.append(level_values)
        if len(arrays) <= 1:
            pvalues[idx] = 1.0
            continue
        try:
            pvalues[idx] = float(f_oneway(*arrays).pvalue)
        except Exception:
            pvalues[idx] = 1.0
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
        return controls
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
        return controls
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
        return controls
    chosen = idx_stage1[ord_mad[:keep_n]]
    controls = numpy.zeros((num_genes,), dtype=bool)
    controls[chosen] = True
    return controls


def compute_design_residuals(seq_uq_df, design_df):
    samples_by_genes = numpy.log(seq_uq_df.to_numpy(dtype=float)).transpose()
    x = design_df.to_numpy(dtype=float)
    beta, _, _, _ = numpy.linalg.lstsq(x, samples_by_genes, rcond=None)
    fitted = x @ beta
    residuals = samples_by_genes - fitted
    return pandas.DataFrame(residuals.transpose(), index=seq_uq_df.index, columns=seq_uq_df.columns)


def ruvr_correct_counts(seq_uq_df, controls, k, residuals_df, center=True, round_counts=True, epsilon=1.0, tolerance=1e-8, is_log=False):
    x = seq_uq_df.to_numpy(dtype=float)
    residuals = residuals_df.to_numpy(dtype=float)
    if (not is_log) and numpy.any(numpy.abs(x - numpy.round(x)) > 1e-8):
        pass
    y = x.transpose() if is_log else numpy.log(x + float(epsilon)).transpose()
    e = residuals.transpose()
    if center:
        e = e - numpy.mean(e, axis=0, keepdims=True)
    controls = numpy.asarray(controls, dtype=bool).reshape(-1)
    if controls.size != x.shape[0]:
        raise ValueError('controls must have one element per gene.')
    if int(k) <= 0:
        return seq_uq_df.copy(), pandas.DataFrame(index=seq_uq_df.columns)
    e_controls = e[:, controls]
    u, s, _vh = numpy.linalg.svd(e_controls, full_matrices=False)
    positive = numpy.where(s > float(tolerance))[0]
    if positive.size == 0:
        return seq_uq_df.copy(), pandas.DataFrame(index=seq_uq_df.columns)
    resolved_k = min(int(k), int(positive[-1] + 1))
    if resolved_k <= 0:
        return seq_uq_df.copy(), pandas.DataFrame(index=seq_uq_df.columns)
    w = u[:, :resolved_k]
    alpha, _, _, _ = numpy.linalg.lstsq(w, y, rcond=None)
    corrected_y = y - (w @ alpha)
    if is_log:
        corrected = corrected_y.transpose()
    else:
        corrected = numpy.exp(corrected_y).transpose() - float(epsilon)
        if round_counts:
            corrected = numpy.round(corrected)
            corrected[corrected < 0] = 0
    corrected_df = pandas.DataFrame(corrected, index=seq_uq_df.index, columns=seq_uq_df.columns)
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
):
    if str(k_setting) != 'auto':
        selected_k = int(k_setting)
        if selected_k < 0:
            selected_k = 1
        corrected_df, w_df = ruvr_correct_counts(seq_uq_df, controls, selected_k, residuals_df)
        comp = score_ruvseq_components(
            mat_df=corrected_df,
            metadata_df=metadata_df,
            batch_column=batch_column,
            sample_group_column=sample_group_column,
        )
        return {
            'k': selected_k,
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
    max_k = max(1, int(k_max))
    max_allowed = max(0, seq_uq_df.shape[1] - 1)
    max_k = min(max_k, max_allowed)
    baseline_comp = score_ruvseq_components(
        mat_df=seq_uq_df,
        metadata_df=metadata_df,
        batch_column=batch_column,
        sample_group_column=sample_group_column,
    )
    baseline_score = baseline_comp['score']
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
        corrected_df, w_df = ruvr_correct_counts(seq_uq_df, controls, k, residuals_df)
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
            best_k = k
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
            best_k = k
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
            if (int(best_k) <= 0) and (int(k) > 0):
                should_replace = True
            elif (int(best_k) > 0) and (int(k) > 0) and (int(k) < int(best_k)):
                should_replace = True
            if not should_replace:
                continue
            best_k = k
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
):
    aligned_metadata = _align_metadata_to_counts(counts_df=counts_df, metadata_df=metadata_df)
    method = 'manual' if str(k_setting) != 'auto' else 'auto'
    try:
        design_df, sample_groups = _build_sample_group_design(
            aligned_metadata=aligned_metadata,
            sample_group_column=sample_group_column,
        )
    except ValueError:
        return _build_ruvseq_skip_output(
            counts_df=counts_df,
            method=method,
            skip_reason='ruvseq_design_failed',
        )
    if len(set(sample_groups)) <= 1:
        return _build_ruvseq_skip_output(
            counts_df=counts_df,
            method=method,
            skip_reason='ruvseq_design_failed',
        )
    counts_plus_one = counts_df.astype(float) + 1.0
    _edge_uq_df, _uq_factors, effective_lib_sizes = _upperquartile_normalize(counts_plus_one, round_counts=True)
    seq_uq_df, _seq_uq_scales = _between_lane_normalize_upper(counts_plus_one, round_counts=True)
    pvalues, residuals_df = _compute_glm_pvalues_and_residuals(
        counts_df=counts_plus_one,
        design_df=design_df,
        effective_lib_sizes=effective_lib_sizes,
    )
    if (pvalues is None) or (residuals_df is None):
        residuals_df = compute_design_residuals(seq_uq_df=seq_uq_df, design_df=design_df)
        pvalues = _compute_group_pvalues(seq_uq_df=seq_uq_df, sample_groups=sample_groups)
    controls = select_ruvseq_controls(
        counts_df=counts_plus_one,
        seq_uq_df=seq_uq_df,
        pvalues=pvalues,
        design_df=design_df,
        mode=control_mode,
        top_n=top_n,
        min_controls=min_controls,
        effective_lib_sizes=effective_lib_sizes,
    )
    if int(controls.sum()) < 2:
        controls = numpy.ones((seq_uq_df.shape[0],), dtype=bool)
    resolved = resolve_ruvseq_k_and_matrix(
        seq_uq_df=seq_uq_df,
        controls=controls,
        residuals_df=residuals_df,
        metadata_df=aligned_metadata,
        k_setting=k_setting,
        k_max=k_max,
        batch_column=batch_column,
        sample_group_column=sample_group_column,
    )
    corrected_df = counts_df.copy()
    corrected_run_ids = []
    uncorrected_run_ids = [str(run_id) for run_id in counts_df.columns]
    skip_reason = 'ruvseq_k_zero'
    if int(resolved['k']) > 0:
        corrected_df = resolved['matrix'].reindex(index=counts_df.index, columns=counts_df.columns)
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
        'ruv_baseline_score': resolved['baseline_score'],
        'ruv_selected_score': resolved['score'],
        'ruv_selected_penalized_score': resolved['penalized_score'],
        'ruv_penalty': resolved['penalty'],
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
