import math
from dataclasses import dataclass, field
from typing import Callable, List, Optional

import numpy
import pandas
from scipy import interpolate, stats


@dataclass
class SVAEstimate:
    nsv: Optional[int]
    method: str


@dataclass
class SVAParameterResolution:
    nsv: int
    B: int
    stable: bool
    method: str
    trace_B: List[int] = field(default_factory=list)
    trace_nsv: List[int] = field(default_factory=list)
    trace_method: List[str] = field(default_factory=list)


def build_sample_group_design_matrix(sample_groups):
    if sample_groups is None:
        raise ValueError('sample_groups is required.')
    values = pandas.Series(sample_groups).fillna('').astype(str).str.strip()
    if values.shape[0] == 0:
        raise ValueError('sample_groups is empty.')
    if (values == '').any():
        raise ValueError('sample_groups contains empty values.')
    levels = sorted(values.drop_duplicates().tolist())
    intercept = numpy.ones((values.shape[0], 1), dtype=float)
    column_names = ['(Intercept)']
    if len(levels) == 1:
        return intercept, column_names
    dummy_columns = []
    for level in levels[1:]:
        dummy_columns.append((values == level).astype(float).to_numpy().reshape(-1, 1))
        column_names.append('sample_group{}'.format(level))
    return numpy.hstack([intercept] + dummy_columns), column_names


def build_intercept_only_design_matrix(num_rows):
    num_rows = _coerce_int(num_rows)
    if (num_rows is None) or (num_rows < 0):
        raise ValueError('num_rows must be a non-negative integer.')
    return numpy.ones((num_rows, 1), dtype=float), ['(Intercept)']


def clean_y_matrix(y_matrix, mod_matrix, sv_matrix):
    y = numpy.asarray(y_matrix, dtype=float)
    mod = numpy.asarray(mod_matrix, dtype=float)
    svs = numpy.asarray(sv_matrix, dtype=float)
    if y.ndim != 2:
        raise ValueError('y_matrix must be two-dimensional.')
    if mod.ndim != 2:
        raise ValueError('mod_matrix must be two-dimensional.')
    if svs.ndim != 2:
        raise ValueError('sv_matrix must be two-dimensional.')
    if mod.shape[0] != y.shape[1]:
        raise ValueError('mod_matrix row count must match sample count.')
    if svs.shape[0] != y.shape[1]:
        raise ValueError('sv_matrix row count must match sample count.')
    if svs.shape[1] == 0:
        return y.copy()
    X = numpy.hstack([mod, svs])
    hat = numpy.linalg.solve(X.T @ X, X.T)
    beta = hat @ y.T
    P = mod.shape[1]
    adjusted = y - (X[:, P:] @ beta[P:, :]).T
    return adjusted


def _coerce_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _compute_hat_matrix(design_matrix):
    design = numpy.asarray(design_matrix, dtype=float)
    if design.ndim != 2:
        raise ValueError('design_matrix must be two-dimensional.')
    xtx = design.T @ design
    try:
        solved = numpy.linalg.solve(xtx, design.T)
    except numpy.linalg.LinAlgError:
        solved = numpy.linalg.pinv(xtx) @ design.T
    return design @ solved


def _residualize_by_design(data_matrix, hat_matrix):
    data = numpy.asarray(data_matrix, dtype=float)
    hat = numpy.asarray(hat_matrix, dtype=float)
    if data.ndim != 2:
        raise ValueError('data_matrix must be two-dimensional.')
    if hat.ndim != 2:
        raise ValueError('hat_matrix must be two-dimensional.')
    if data.shape[1] != hat.shape[0]:
        raise ValueError('hat_matrix row count must match sample count.')
    return data - (hat @ data.T).T


def _normalized_singular_energy(singular_values, ndf):
    ndf = _coerce_int(ndf)
    if (ndf is None) or (ndf <= 0):
        return numpy.zeros((0,), dtype=float)
    values = numpy.asarray(singular_values, dtype=float)
    values = values[:ndf]
    if values.shape[0] < ndf:
        values = numpy.pad(values, (0, ndf - values.shape[0]), mode='constant')
    denom = float(numpy.sum(values ** 2))
    if denom <= 0:
        return numpy.zeros((ndf,), dtype=float)
    return (values ** 2) / denom


def _permute_rows_without_replacement(matrix, rng):
    values = numpy.asarray(matrix, dtype=float)
    if values.ndim != 2:
        raise ValueError('matrix must be two-dimensional.')
    if values.size == 0:
        return values.copy()
    order = numpy.argsort(rng.random(values.shape), axis=1)
    return numpy.take_along_axis(values, order, axis=1)


def _row_vars(matrix):
    values = numpy.asarray(matrix, dtype=float)
    if values.ndim != 2:
        raise ValueError('matrix must be two-dimensional.')
    if values.shape[1] <= 1:
        return numpy.zeros((values.shape[0],), dtype=float)
    return numpy.var(values, axis=1, ddof=1)


def _modefunc(values):
    array = numpy.asarray(values).reshape(-1)
    if array.size == 0:
        return 0
    unique_values, counts = numpy.unique(array, return_counts=True)
    order = numpy.lexsort((unique_values, -counts))
    return unique_values[order[0]]


def _orthogonal_complement_basis(design_matrix, tol=1e-10):
    design = numpy.asarray(design_matrix, dtype=float)
    if design.ndim != 2:
        raise ValueError('design_matrix must be two-dimensional.')
    q_matrix, _r_matrix = numpy.linalg.qr(design, mode='complete')
    rank = int(numpy.linalg.matrix_rank(design, tol=tol))
    return q_matrix[:, rank:]


def _normalize_vector(vector, tol=1e-10):
    values = numpy.asarray(vector, dtype=float).reshape(-1)
    norm = float(numpy.linalg.norm(values))
    if norm <= tol:
        return None
    return values / norm


def _augmented_design_has_full_rank(design_matrix, extra_columns, tol=1e-10):
    design = numpy.asarray(design_matrix, dtype=float)
    extras = numpy.asarray(extra_columns, dtype=float)
    if extras.ndim == 1:
        extras = extras.reshape(-1, 1)
    if extras.shape[1] == 0:
        return True
    augmented = numpy.hstack([design, extras])
    return int(numpy.linalg.matrix_rank(augmented, tol=tol)) == augmented.shape[1]


def _stabilize_surrogate_matrix(sv_matrix, design_matrix, tol=1e-10):
    sv = numpy.asarray(sv_matrix, dtype=float)
    design = numpy.asarray(design_matrix, dtype=float)
    if sv.ndim != 2:
        raise ValueError('sv_matrix must be two-dimensional.')
    if design.ndim != 2:
        raise ValueError('design_matrix must be two-dimensional.')
    if sv.shape[0] != design.shape[0]:
        raise ValueError('sv_matrix row count must match design row count.')
    if sv.shape[1] == 0:
        return sv.copy()

    basis = _orthogonal_complement_basis(design, tol=tol)
    stabilized_columns = []
    basis_index = 0

    def accept_or_fallback(candidate_vector):
        nonlocal basis_index
        candidate = _normalize_vector(candidate_vector, tol=tol)
        if candidate is not None:
            extras = stabilized_columns + [candidate]
            if _augmented_design_has_full_rank(design, numpy.column_stack(extras), tol=tol):
                return candidate
        while basis_index < basis.shape[1]:
            fallback = basis[:, basis_index].copy()
            basis_index += 1
            for prev in stabilized_columns:
                fallback = fallback - (prev @ fallback) * prev
            fallback = _normalize_vector(fallback, tol=tol)
            if fallback is None:
                continue
            extras = stabilized_columns + [fallback]
            if _augmented_design_has_full_rank(design, numpy.column_stack(extras), tol=tol):
                return fallback
        return None

    for col_index in range(sv.shape[1]):
        accepted = accept_or_fallback(sv[:, col_index])
        if accepted is None:
            continue
        stabilized_columns.append(accepted)

    if len(stabilized_columns) == 0:
        return numpy.zeros((sv.shape[0], 0), dtype=float)
    return numpy.column_stack(stabilized_columns)


def estimate_num_sv_be(
    data_matrix,
    mod_matrix,
    B_value,
    max_nsv,
    random_seed=None,
):
    data = numpy.asarray(data_matrix, dtype=float)
    mod = numpy.asarray(mod_matrix, dtype=float)
    B_value = _coerce_int(B_value)
    max_nsv = _coerce_int(max_nsv)
    if data.ndim != 2:
        raise ValueError('data_matrix must be two-dimensional.')
    if mod.ndim != 2:
        raise ValueError('mod_matrix must be two-dimensional.')
    if data.shape[1] != mod.shape[0]:
        raise ValueError('mod_matrix row count must match sample count.')
    if (B_value is None) or (B_value < 1):
        B_value = 20
    if (max_nsv is None) or (max_nsv < 0):
        max_nsv = max(0, min(data.shape) - mod.shape[1] - 1)
    if (data.shape[0] == 0) or (data.shape[1] == 0):
        return SVAEstimate(nsv=0, method='be')

    hat = _compute_hat_matrix(mod)
    residual = _residualize_by_design(data, hat)
    ndf = min(data.shape) - int(math.ceil(float(numpy.trace(hat))))
    if ndf <= 0:
        return SVAEstimate(nsv=0, method='be')

    singular_values = numpy.linalg.svd(residual, full_matrices=False, compute_uv=False)
    if float(numpy.sum(singular_values[:ndf] ** 2)) <= 0:
        return SVAEstimate(nsv=None, method='be')
    observed = _normalized_singular_energy(singular_values, ndf=ndf)
    if observed.shape[0] == 0:
        return SVAEstimate(nsv=None, method='be')

    rng = numpy.random.default_rng(_coerce_int(random_seed))
    permuted_stats = numpy.zeros((B_value, ndf), dtype=float)
    for idx in range(B_value):
        residual_perm = _permute_rows_without_replacement(residual, rng)
        residual_perm = _residualize_by_design(residual_perm, hat)
        permuted_singular_values = numpy.linalg.svd(residual_perm, full_matrices=False, compute_uv=False)
        if float(numpy.sum(permuted_singular_values[:ndf] ** 2)) <= 0:
            return SVAEstimate(nsv=None, method='be')
        permuted_stats[idx, :] = _normalized_singular_energy(
            permuted_singular_values,
            ndf=ndf,
        )

    psv = numpy.ones((ndf,), dtype=float)
    for idx in range(ndf):
        psv[idx] = float(numpy.mean(permuted_stats[:, idx] >= observed[idx]))
    if ndf > 1:
        psv = numpy.maximum.accumulate(psv)
    nsv = int(numpy.sum(psv <= 0.10))
    return SVAEstimate(nsv=min(max_nsv, max(0, nsv)), method='be')


def estimate_num_sv_leek(data_matrix, mod_matrix, max_nsv):
    data = numpy.asarray(data_matrix, dtype=float)
    mod = numpy.asarray(mod_matrix, dtype=float)
    max_nsv = _coerce_int(max_nsv)
    if data.ndim != 2:
        raise ValueError('data_matrix must be two-dimensional.')
    if mod.ndim != 2:
        raise ValueError('mod_matrix must be two-dimensional.')
    if data.shape[1] != mod.shape[0]:
        raise ValueError('mod_matrix row count must match sample count.')
    if (max_nsv is None) or (max_nsv < 0):
        max_nsv = max(0, min(data.shape) - mod.shape[1] - 1)
    if data.shape[0] < 10:
        raise ValueError('Leek estimator requires at least 10 rows.')

    dims = data.shape
    a = numpy.linspace(0.0, 2.0, num=100)
    block_size = int(math.floor(dims[0] / 10))
    if block_size < 1:
        raise ValueError('Leek estimator requires at least one row per block.')
    rhat = numpy.zeros((100, 10), dtype=float)
    projection = numpy.eye(dims[1]) - _compute_hat_matrix(mod)

    for j in range(10):
        dats = data[: (j + 1) * block_size, :]
        gram = dats.T @ dats
        gram_values = numpy.linalg.eigvalsh(gram)
        gram_values = numpy.sort(gram_values)[::-1]
        sigbar = float(gram_values[-1] / ((j + 1) * block_size))
        residual = dats @ projection
        wm = (residual.T @ residual) / float((j + 1) * block_size) - (projection * sigbar)
        wm_values = numpy.linalg.eigvalsh(wm)
        wm_values = numpy.sort(wm_values)[::-1]
        combined = numpy.concatenate([
            a * (((j + 1) * block_size) ** (-1.0 / 3.0)) * dims[1],
            wm_values,
        ])
        markers = numpy.concatenate([
            numpy.ones((100,), dtype=bool),
            numpy.zeros((dims[1],), dtype=bool),
        ])
        order = numpy.argsort(combined)[::-1]
        ordered_markers = markers[order]
        positions = numpy.nonzero(ordered_markers)[0] + 1
        rhat[:, j] = (positions - numpy.arange(1, 101))[::-1]

    ss = _row_vars(rhat)
    bump_candidates = numpy.nonzero(ss > (2.0 * ss[0]))[0]
    bumpstart = int(bump_candidates[0] + 1) if bump_candidates.size > 0 else 1
    start_probe = numpy.concatenate([
        numpy.full((bumpstart,), 1e5, dtype=float),
        ss[bumpstart:],
    ])
    start_candidates = numpy.nonzero(start_probe < (0.5 * ss[0]))[0]
    start = int(start_candidates[0] + 1) if start_candidates.size > 0 else 1
    finish_mask = ss * numpy.concatenate([
        numpy.zeros((start,), dtype=float),
        numpy.ones((100 - start,), dtype=float),
    ])
    finish_candidates = numpy.nonzero(finish_mask > ss[0])[0]
    finish = int(finish_candidates[0] + 1) if finish_candidates.size > 0 else 1
    if finish == 1:
        finish = 100
    nsv = int(_modefunc(rhat[(start - 1):finish, 9]))
    nsv = min(max_nsv, max(0, nsv))
    return SVAEstimate(nsv=nsv, method='leek')


def estimate_num_sv_at_B(
    data_matrix,
    mod_matrix,
    B_value,
    max_nsv,
    random_seed=None,
):
    try:
        be_estimate = estimate_num_sv_be(
            data_matrix=data_matrix,
            mod_matrix=mod_matrix,
            B_value=B_value,
            max_nsv=max_nsv,
            random_seed=random_seed,
        )
        if be_estimate.nsv is not None:
            return be_estimate
    except (FloatingPointError, ValueError, numpy.linalg.LinAlgError):
        pass
    try:
        return estimate_num_sv_leek(
            data_matrix=data_matrix,
            mod_matrix=mod_matrix,
            max_nsv=max_nsv,
        )
    except (FloatingPointError, ValueError, numpy.linalg.LinAlgError):
        return SVAEstimate(nsv=None, method='failed')


def f_pvalue(data_matrix, mod_matrix, mod0_matrix):
    data = numpy.asarray(data_matrix, dtype=float)
    mod = numpy.asarray(mod_matrix, dtype=float)
    mod0 = numpy.asarray(mod0_matrix, dtype=float)
    if data.ndim != 2:
        raise ValueError('data_matrix must be two-dimensional.')
    if mod.ndim != 2:
        raise ValueError('mod_matrix must be two-dimensional.')
    if mod0.ndim != 2:
        raise ValueError('mod0_matrix must be two-dimensional.')
    if (data.shape[1] != mod.shape[0]) or (data.shape[1] != mod0.shape[0]):
        raise ValueError('Design matrix row count must match sample count.')

    sample_count = data.shape[1]
    df1 = mod.shape[1]
    df0 = mod0.shape[1]
    identity = numpy.eye(sample_count)

    resid = data @ (identity - _compute_hat_matrix(mod))
    rss1 = numpy.sum(resid * resid, axis=1)

    resid0 = data @ (identity - _compute_hat_matrix(mod0))
    rss0 = numpy.sum(resid0 * resid0, axis=1)

    numerator_df = df1 - df0
    denominator_df = sample_count - df1
    if numerator_df <= 0:
        raise ValueError('mod_matrix must have more columns than mod0_matrix.')
    if denominator_df <= 0:
        raise ValueError('Not enough residual degrees of freedom.')

    with numpy.errstate(divide='ignore', invalid='ignore'):
        fstats = ((rss0 - rss1) / float(numerator_df)) / (rss1 / float(denominator_df))
    fstats = numpy.nan_to_num(fstats, nan=0.0, posinf=0.0, neginf=0.0)
    fstats = numpy.maximum(fstats, 0.0)
    return 1.0 - stats.f.cdf(fstats, dfn=numerator_df, dfd=denominator_df)


def edge_lfdr(
    p_values,
    trunc=True,
    monotone=True,
    transf='probit',
    adj=1.5,
    eps=1e-8,
    lambda_value=0.8,
):
    p = numpy.asarray(p_values, dtype=float).reshape(-1)
    if p.size == 0:
        return p.copy()
    pi0 = float(numpy.mean(p >= lambda_value) / (1.0 - lambda_value))
    pi0 = min(pi0, 1.0)
    if p.size < 2:
        return numpy.full(p.shape, fill_value=pi0, dtype=float)
    transf = str(transf).strip().lower()

    try:
        if transf == 'probit':
            p = numpy.clip(p, eps, 1.0 - eps)
            x = stats.norm.ppf(p)
            if numpy.allclose(x, x[0]):
                x = x + numpy.linspace(-eps, eps, num=x.shape[0])
            kde = stats.gaussian_kde(x, bw_method=lambda obj: obj.scotts_factor() * adj)
            grid = numpy.linspace(float(numpy.min(x)), float(numpy.max(x)), max(512, p.size))
            density_values = kde(grid)
            spline = interpolate.UnivariateSpline(grid, density_values, s=len(grid))
            y = spline(x)
            y = numpy.maximum(y, eps)
            lfdr = pi0 * stats.norm.pdf(x) / y
        elif transf == 'logit':
            x = numpy.log((p + eps) / (1.0 - p + eps))
            if numpy.allclose(x, x[0]):
                x = x + numpy.linspace(-eps, eps, num=x.shape[0])
            kde = stats.gaussian_kde(x, bw_method=lambda obj: obj.scotts_factor() * adj)
            grid = numpy.linspace(float(numpy.min(x)), float(numpy.max(x)), max(512, p.size))
            density_values = kde(grid)
            spline = interpolate.UnivariateSpline(grid, density_values, s=len(grid))
            y = spline(x)
            y = numpy.maximum(y, eps)
            dx = numpy.exp(x) / (1.0 + numpy.exp(x)) ** 2
            lfdr = pi0 * dx / y
        else:
            raise ValueError('Unsupported edge_lfdr transform: {}'.format(transf))
    except (ValueError, numpy.linalg.LinAlgError):
        lfdr = numpy.full(p.shape, fill_value=pi0, dtype=float)

    if trunc:
        lfdr = numpy.minimum(lfdr, 1.0)
    if monotone:
        order = numpy.argsort(p, kind='mergesort')
        ranked = lfdr[order]
        ranked = numpy.maximum.accumulate(ranked)
        inverse = numpy.empty_like(order)
        inverse[order] = numpy.arange(order.shape[0])
        lfdr = ranked[inverse]
    return lfdr


def irwsva_build(data_matrix, mod_matrix, mod0_matrix=None, nsv=1, B_iterations=5):
    data = numpy.asarray(data_matrix, dtype=float)
    mod = numpy.asarray(mod_matrix, dtype=float)
    if mod0_matrix is None:
        mod0 = mod[:, [0]]
    else:
        mod0 = numpy.asarray(mod0_matrix, dtype=float)
    nsv = _coerce_int(nsv)
    B_iterations = _coerce_int(B_iterations)
    if data.ndim != 2:
        raise ValueError('data_matrix must be two-dimensional.')
    if mod.ndim != 2:
        raise ValueError('mod_matrix must be two-dimensional.')
    if mod0.ndim != 2:
        raise ValueError('mod0_matrix must be two-dimensional.')
    if (data.shape[1] != mod.shape[0]) or (data.shape[1] != mod0.shape[0]):
        raise ValueError('Design matrix row count must match sample count.')
    if (nsv is None) or (nsv <= 0):
        raise ValueError('nsv must be a positive integer.')
    if (B_iterations is None) or (B_iterations < 1):
        B_iterations = 5

    sample_count = data.shape[1]
    identity = numpy.eye(sample_count)
    residual = data @ (identity - _compute_hat_matrix(mod))
    eigvals, eigvecs = numpy.linalg.eigh(residual.T @ residual)
    order = numpy.argsort(eigvals)[::-1]
    eigvecs = eigvecs[:, order]
    current_sv = _stabilize_surrogate_matrix(eigvecs[:, :nsv], mod)
    if current_sv.shape[1] == 0:
        current_sv = _orthogonal_complement_basis(mod)[:, :nsv]

    pprob_gam = numpy.zeros((data.shape[0],), dtype=float)
    pprob_b = numpy.zeros((data.shape[0],), dtype=float)
    dats = data.copy()
    for _ in range(B_iterations):
        mod_b = numpy.hstack([mod, current_sv])
        mod0_b = numpy.hstack([mod0, current_sv])
        ptmp_b = f_pvalue(data, mod_b, mod0_b)
        pprob_b = 1.0 - edge_lfdr(ptmp_b)

        mod_gam = numpy.hstack([mod0, current_sv])
        ptmp_gam = f_pvalue(data, mod_gam, mod0)
        pprob_gam = 1.0 - edge_lfdr(ptmp_gam)

        pprob = pprob_gam * (1.0 - pprob_b)
        dats = data * pprob.reshape(-1, 1)
        dats = dats - numpy.mean(dats, axis=1, keepdims=True)
        _svd_u, _svd_s, svd_vh = numpy.linalg.svd(dats, full_matrices=False)
        current_sv = _stabilize_surrogate_matrix(svd_vh.T[:, :nsv], mod)
        if current_sv.shape[1] == 0:
            current_sv = _orthogonal_complement_basis(mod)[:, :nsv]

    _svd_u, _svd_s, svd_vh = numpy.linalg.svd(dats, full_matrices=False)
    sv = _stabilize_surrogate_matrix(svd_vh.T[:, :nsv], mod)
    if sv.shape[1] == 0:
        sv = _orthogonal_complement_basis(mod)[:, :nsv]
    return {
        'sv': sv,
        'pprob_gam': pprob_gam,
        'pprob_b': pprob_b,
        'n_svs': nsv,
    }


def resolve_sva_B_value(B_setting='auto', sample_count=None, auto_max=100):
    auto_max = _coerce_int(auto_max)
    if (auto_max is None) or (auto_max < 5):
        auto_max = 100
    setting_text = str(B_setting).strip().lower()
    if setting_text != 'auto':
        manual = _coerce_int(setting_text)
        if (manual is not None) and (manual >= 1):
            return int(manual)
    sample_count = _coerce_int(sample_count)
    if (sample_count is None) or (sample_count < 2):
        return int(min(auto_max, 20))
    auto_value = int(math.ceil(120 / sample_count))
    auto_value = max(5, min(auto_max, auto_value))
    return int(auto_value)


def resolve_sva_parameters(
    num_samples,
    design_columns,
    nsv_setting='auto',
    B_setting='auto',
    B_auto_max=100,
    estimate_nsv_at_B: Optional[Callable[[int, int], SVAEstimate]] = None,
):
    num_samples = _coerce_int(num_samples)
    design_columns = _coerce_int(design_columns)
    if (num_samples is None) or (num_samples < 0):
        raise ValueError('num_samples must be a non-negative integer.')
    if (design_columns is None) or (design_columns < 0):
        raise ValueError('design_columns must be a non-negative integer.')
    max_nsv = max(0, num_samples - design_columns - 1)
    B_setting_text = str(B_setting).strip().lower()
    nsv_setting_text = str(nsv_setting).strip().lower()
    B_default = resolve_sva_B_value(
        B_setting=B_setting_text,
        sample_count=num_samples,
        auto_max=B_auto_max,
    )
    if nsv_setting_text != 'auto':
        manual_nsv = _coerce_int(nsv_setting_text)
        if (manual_nsv is not None) and (manual_nsv >= 0):
            return SVAParameterResolution(
                nsv=min(max_nsv, manual_nsv),
                B=B_default,
                stable=True,
                method='manual_nsv',
            )
    if max_nsv == 0:
        return SVAParameterResolution(
            nsv=0,
            B=B_default,
            stable=True,
            method='max_nsv_zero',
        )
    if estimate_nsv_at_B is None:
        raise ValueError('estimate_nsv_at_B is required when nsv_setting is auto and max_nsv > 0.')

    def normalize_estimate(raw_estimate):
        if isinstance(raw_estimate, SVAEstimate):
            nsv_value = raw_estimate.nsv
            method = raw_estimate.method
        else:
            nsv_value = getattr(raw_estimate, 'nsv', None)
            method = getattr(raw_estimate, 'method', 'failed')
        nsv_value = _coerce_int(nsv_value)
        if nsv_value is None:
            return SVAEstimate(nsv=None, method=str(method))
        return SVAEstimate(
            nsv=min(max_nsv, max(0, nsv_value)),
            method=str(method),
        )

    if B_setting_text != 'auto':
        estimate = normalize_estimate(estimate_nsv_at_B(B_default, max_nsv))
        return SVAParameterResolution(
            nsv=0 if estimate.nsv is None else int(estimate.nsv),
            B=B_default,
            stable=True,
            method='manual_B_{}'.format(estimate.method),
            trace_B=[int(B_default)],
            trace_nsv=[-1 if estimate.nsv is None else int(estimate.nsv)],
            trace_method=[estimate.method],
        )

    B_auto_max = _coerce_int(B_auto_max)
    if (B_auto_max is None) or (B_auto_max < 5):
        B_auto_max = 100
    base_B = resolve_sva_B_value(B_setting='auto', sample_count=num_samples, auto_max=B_auto_max)
    candidates = sorted({
        int(base_B),
        int(math.ceil(base_B * 1.5)),
        int(base_B * 2),
        int(B_auto_max),
    })
    candidates = [cand for cand in candidates if (cand >= 1) and (cand <= B_auto_max)]
    if len(candidates) == 0:
        candidates = [int(min(B_auto_max, 20))]

    prev_nsv = None
    selected_nsv = None
    selected_B = int(candidates[-1])
    selected_method = 'failed'
    stable = False
    trace_B = []
    trace_nsv = []
    trace_method = []

    for candidate in candidates:
        estimate = normalize_estimate(estimate_nsv_at_B(candidate, max_nsv))
        trace_B.append(int(candidate))
        trace_nsv.append(-1 if estimate.nsv is None else int(estimate.nsv))
        trace_method.append(estimate.method)
        if estimate.nsv is None:
            continue
        selected_nsv = int(estimate.nsv)
        selected_B = int(candidate)
        selected_method = estimate.method
        if (prev_nsv is not None) and (selected_nsv == prev_nsv):
            stable = True
            break
        prev_nsv = selected_nsv

    if selected_nsv is None:
        selected_nsv = 0

    return SVAParameterResolution(
        nsv=int(selected_nsv),
        B=int(selected_B),
        stable=stable,
        method='auto_{}'.format(selected_method),
        trace_B=trace_B,
        trace_nsv=trace_nsv,
        trace_method=trace_method,
    )


def run_sva_backend(
    counts_df,
    metadata_df,
    nsv_setting='auto',
    B_setting='auto',
    B_auto_max=100,
    sample_group_column='sample_group',
    random_seed=None,
):
    if sample_group_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(sample_group_column))
    sample_groups = metadata_df.loc[:, sample_group_column]
    mod_matrix, _mod_names = build_sample_group_design_matrix(sample_groups)
    _mod0_matrix, _mod0_names = build_intercept_only_design_matrix(metadata_df.shape[0])
    counts_matrix = counts_df.to_numpy(dtype=float, copy=False)
    resolved = resolve_sva_parameters(
        num_samples=counts_df.shape[1],
        design_columns=mod_matrix.shape[1],
        nsv_setting=nsv_setting,
        B_setting=B_setting,
        B_auto_max=B_auto_max,
        estimate_nsv_at_B=lambda B_value, max_nsv: estimate_num_sv_at_B(
            data_matrix=counts_matrix,
            mod_matrix=mod_matrix,
            B_value=B_value,
            max_nsv=max_nsv,
            random_seed=random_seed,
        ),
    )
    if resolved.nsv <= 0:
        empty_sv = pandas.DataFrame(index=counts_df.columns)
        summary = {
            'backend': 'sva',
            'method': resolved.method,
            'skip_reason': 'sva_nsv_zero',
            'stable': resolved.stable,
            'corrected_run_ids': [],
            'uncorrected_run_ids': [str(col) for col in counts_df.columns],
            'resolved_sva_nsv': resolved.nsv,
            'resolved_sva_B': resolved.B,
            'sva_estimation_method': resolved.method,
            'sva_stable': resolved.stable,
            'trace_B': resolved.trace_B,
            'trace_nsv': resolved.trace_nsv,
            'trace_method': resolved.trace_method,
        }
        return counts_df.copy(), empty_sv, summary
    irw = irwsva_build(
        data_matrix=counts_matrix,
        mod_matrix=mod_matrix,
        mod0_matrix=_mod0_matrix,
        nsv=resolved.nsv,
        B_iterations=resolved.B,
    )
    sv_matrix = numpy.asarray(irw['sv'], dtype=float)
    corrected_matrix = clean_y_matrix(
        y_matrix=counts_matrix,
        mod_matrix=mod_matrix,
        sv_matrix=sv_matrix,
    )
    sv_columns = ['sv{}'.format(idx + 1) for idx in range(sv_matrix.shape[1])]
    corrected_df = pandas.DataFrame(
        corrected_matrix,
        index=counts_df.index,
        columns=counts_df.columns,
    )
    sv_df = pandas.DataFrame(
        sv_matrix,
        index=counts_df.columns,
        columns=sv_columns,
    )
    summary = {
        'backend': 'sva',
        'method': 'irw',
        'skip_reason': '',
        'stable': resolved.stable,
        'corrected_run_ids': [str(col) for col in counts_df.columns],
        'uncorrected_run_ids': [],
        'resolved_sva_nsv': resolved.nsv,
        'resolved_sva_B': resolved.B,
        'sva_estimation_method': resolved.method,
        'sva_stable': resolved.stable,
        'trace_B': resolved.trace_B,
        'trace_nsv': resolved.trace_nsv,
        'trace_method': resolved.trace_method,
    }
    return corrected_df, sv_df, summary
