import math
from dataclasses import dataclass

import numpy
import pandas
from scipy.stats import rankdata


@dataclass
class TMMRoundTripResult:
    round1_factors: pandas.Series
    round1_reference_column: int
    median_reference_columns: list
    round2_factors: pandas.Series
    library_sizes: pandas.Series


def _as_matrix_with_columns(counts):
    if isinstance(counts, pandas.DataFrame):
        matrix = counts.to_numpy(dtype=float, copy=False)
        columns = [str(col) for col in counts.columns]
    else:
        matrix = numpy.asarray(counts, dtype=float)
        if matrix.ndim != 2:
            raise ValueError('counts must be two-dimensional.')
        columns = ['sample{}'.format(i + 1) for i in range(matrix.shape[1])]
    if matrix.ndim != 2:
        raise ValueError('counts must be two-dimensional.')
    if numpy.isnan(matrix).any():
        raise ValueError('NA counts are not permitted.')
    return matrix, columns


def _coerce_library_sizes(lib_size, matrix, columns):
    nsamples = matrix.shape[1]
    if lib_size is None:
        values = matrix.sum(axis=0, dtype=float)
    elif isinstance(lib_size, pandas.Series):
        values = lib_size.reindex(columns).to_numpy(dtype=float)
    else:
        values = numpy.asarray(lib_size, dtype=float).reshape(-1)
        if values.size != nsamples:
            values = numpy.resize(values, nsamples)
    if numpy.isnan(values).any():
        raise ValueError('NA lib.sizes are not permitted.')
    return pandas.Series(values, index=columns, dtype=float)


def _remove_all_zero_rows(matrix):
    keep = numpy.count_nonzero(matrix > 0, axis=1) != 0
    if keep.all():
        return matrix
    return matrix[keep, :]


def _recycle(values, size):
    array = numpy.asarray(values, dtype=float).reshape(-1)
    if array.size == 0:
        return array
    if array.size == size:
        return array
    return numpy.resize(array, size)


def calc_factor_quantile(data, lib_size, p=0.75):
    matrix, columns = _as_matrix_with_columns(data)
    libs = _coerce_library_sizes(lib_size, matrix, columns).to_numpy(dtype=float)
    factors = numpy.ones((matrix.shape[1],), dtype=float)
    for idx in range(matrix.shape[1]):
        factors[idx] = float(numpy.quantile(matrix[:, idx], q=p, method='linear'))
    return pandas.Series(factors / libs, index=columns, dtype=float)


def calc_factor_tmm(
    obs,
    ref,
    libsize_obs=None,
    libsize_ref=None,
    logratio_trim=0.3,
    sum_trim=0.05,
    do_weighting=True,
    acutoff=-1e10,
):
    obs_array = numpy.asarray(obs, dtype=float)
    ref_array = numpy.asarray(ref, dtype=float)
    if obs_array.ndim > 2 or ref_array.ndim > 2:
        raise ValueError('obs and ref must be one- or two-dimensional.')
    obs_vector = obs_array.reshape(-1, order='F')
    ref_vector = ref_array.reshape(-1, order='F')
    n_obs = float(obs_vector.sum()) if libsize_obs is None else float(numpy.asarray(libsize_obs, dtype=float).reshape(-1)[0])
    if libsize_ref is None:
        n_ref_values = numpy.array([float(ref_vector.sum())], dtype=float)
    else:
        n_ref_values = numpy.asarray(libsize_ref, dtype=float).reshape(-1)
        if n_ref_values.size == 0:
            n_ref_values = numpy.array([float(ref_vector.sum())], dtype=float)
    target_size = max(obs_vector.size, ref_vector.size)
    obs_scaled = _recycle(obs_vector / n_obs, target_size)
    ref_scaled = _recycle(ref_vector, target_size) / _recycle(n_ref_values, target_size)
    log_r = numpy.log2(obs_scaled / ref_scaled)
    abs_e = (numpy.log2(obs_scaled) + numpy.log2(ref_scaled)) / 2.0
    left_var = _recycle((n_obs - obs_vector) / n_obs / obs_vector, target_size)
    right_var = (
        (_recycle(n_ref_values, target_size) - _recycle(ref_vector, target_size))
        / _recycle(n_ref_values, target_size)
        / _recycle(ref_vector, target_size)
    )
    variances = left_var + right_var
    finite = numpy.isfinite(log_r) & numpy.isfinite(abs_e) & (abs_e > acutoff)
    log_r = log_r[finite]
    abs_e = abs_e[finite]
    variances = variances[finite]
    if log_r.size == 0:
        return 1.0
    if float(numpy.max(numpy.abs(log_r))) < 1e-6:
        return 1.0
    n = log_r.size
    lo_l = int(math.floor(n * logratio_trim) + 1)
    hi_l = int(n + 1 - lo_l)
    lo_s = int(math.floor(n * sum_trim) + 1)
    hi_s = int(n + 1 - lo_s)
    rank_log_r = rankdata(log_r, method='average')
    rank_abs_e = rankdata(abs_e, method='average')
    keep = (
        (rank_log_r >= lo_l) & (rank_log_r <= hi_l) &
        (rank_abs_e >= lo_s) & (rank_abs_e <= hi_s)
    )
    if do_weighting:
        numerator = numpy.nansum(log_r[keep] / variances[keep])
        denominator = numpy.nansum(1.0 / variances[keep])
        factor_log = numerator / denominator
    else:
        factor_log = float(numpy.nanmean(log_r[keep]))
    if numpy.isnan(factor_log):
        factor_log = 0.0
    return float(2.0 ** factor_log)


def _resolve_tmm_reference_column(matrix, lib_sizes):
    if isinstance(lib_sizes, pandas.Series):
        lib_size_values = lib_sizes.to_numpy(dtype=float)
    else:
        lib_size_values = numpy.asarray(lib_sizes, dtype=float).reshape(-1)
    f75 = calc_factor_quantile(matrix, lib_size_values).to_numpy(dtype=float)
    if float(numpy.median(f75)) < 1e-20:
        ref_column = int(numpy.argmax(numpy.sum(numpy.sqrt(matrix), axis=0)))
    else:
        ref_column = int(numpy.argmin(numpy.abs(f75 - numpy.mean(f75))))
    return ref_column


def calc_norm_factors_tmm(
    counts,
    lib_size=None,
    ref_column=None,
    logratio_trim=0.3,
    sum_trim=0.05,
    do_weighting=True,
    acutoff=-1e10,
):
    matrix, columns = _as_matrix_with_columns(counts)
    libs = _coerce_library_sizes(lib_size, matrix, columns)
    matrix = _remove_all_zero_rows(matrix)
    nsamples = matrix.shape[1]
    if (matrix.shape[0] == 0) or (nsamples == 1):
        return pandas.Series(numpy.ones((nsamples,), dtype=float), index=columns, dtype=float)

    if ref_column is None:
        resolved_ref = _resolve_tmm_reference_column(matrix, libs)
    else:
        resolved_ref = ref_column
    if numpy.isscalar(resolved_ref):
        ref_indices = [int(resolved_ref)]
        ref_matrix = matrix[:, ref_indices[0]]
        ref_lib_size = float(libs.iloc[ref_indices[0]])
    else:
        ref_indices = [int(idx) for idx in list(resolved_ref)]
        ref_matrix = matrix[:, ref_indices]
        ref_lib_size = libs.iloc[ref_indices].to_numpy(dtype=float)

    factors = numpy.empty((nsamples,), dtype=float)
    for idx in range(nsamples):
        factors[idx] = calc_factor_tmm(
            obs=matrix[:, idx],
            ref=ref_matrix,
            libsize_obs=float(libs.iloc[idx]),
            libsize_ref=ref_lib_size,
            logratio_trim=logratio_trim,
            sum_trim=sum_trim,
            do_weighting=do_weighting,
            acutoff=acutoff,
        )
    factors = factors / math.exp(float(numpy.mean(numpy.log(factors))))
    return pandas.Series(factors, index=columns, dtype=float)


def run_tmm_rounds_for_cstmm(counts, lib_size=None):
    matrix, columns = _as_matrix_with_columns(counts)
    libs = _coerce_library_sizes(lib_size, matrix, columns)
    round1_counts = counts if isinstance(counts, pandas.DataFrame) else pandas.DataFrame(matrix, columns=columns)
    round1 = calc_norm_factors_tmm(round1_counts, lib_size=libs, ref_column=None)
    round1_values = round1.to_numpy(dtype=float)
    round1_reference_column = _resolve_tmm_reference_column(_remove_all_zero_rows(matrix), libs)
    median_value = float(numpy.sort(round1_values)[int(math.ceil(round1_values.size / 2.0) - 1)])
    median_reference_columns = [
        int(idx) for idx, value in enumerate(round1_values) if value == median_value
    ]
    round2 = calc_norm_factors_tmm(round1_counts, lib_size=libs, ref_column=median_reference_columns)
    return TMMRoundTripResult(
        round1_factors=round1,
        round1_reference_column=round1_reference_column,
        median_reference_columns=median_reference_columns,
        round2_factors=round2,
        library_sizes=libs,
    )


def apply_tmm_factors(counts, norm_factors):
    if not isinstance(counts, pandas.DataFrame):
        raise ValueError('counts must be a pandas.DataFrame.')
    if isinstance(norm_factors, pandas.Series):
        factors = norm_factors.reindex(counts.columns).to_numpy(dtype=float)
    else:
        factors = numpy.asarray(norm_factors, dtype=float).reshape(-1)
        if factors.size != counts.shape[1]:
            raise ValueError('norm_factors must match the number of sample columns.')
    corrected = counts.copy()
    corrected.loc[:, :] = counts.to_numpy(dtype=float, copy=False) / factors.reshape(1, -1)
    return corrected


__all__ = [
    'TMMRoundTripResult',
    'apply_tmm_factors',
    'calc_factor_quantile',
    'calc_factor_tmm',
    'calc_norm_factors_tmm',
    'run_tmm_rounds_for_cstmm',
]
