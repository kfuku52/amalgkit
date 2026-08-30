import math
from dataclasses import dataclass
from fractions import Fraction

import numpy
import pandas


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
    if not numpy.isfinite(matrix).all():
        raise ValueError('Counts must contain only finite values.')
    if (matrix < 0).any():
        raise ValueError('Negative counts are not permitted.')
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
            raise ValueError('lib.size must match the number of sample columns.')
    if not numpy.isfinite(values).all():
        raise ValueError('lib.sizes must contain only finite values.')
    if (values <= 0).any():
        raise ValueError('lib.sizes must contain only positive values.')
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


def _rank_count_expression(obs, ref, ref_libs, *, ratio):
    """Average ranks before logs/library scaling can split equal count ratios.

    Positive products/ratios are represented as mantissa + binary exponent to
    avoid overflow/underflow. No decimal rounding or tolerance merges distinct
    keys. A common library factor does not affect rank and is omitted.
    """
    obs_m, obs_e = numpy.frexp(obs)
    ref_m, ref_e = numpy.frexp(ref)
    mantissas = obs_m / ref_m if ratio else obs_m * ref_m
    exponents = obs_e - ref_e if ratio else obs_e + ref_e
    if not numpy.all(ref_libs == ref_libs[0]):
        # Preserve the historical vector-reference API, where library sizes
        # can vary per recycled element rather than being one common factor.
        lib_m, lib_e = numpy.frexp(ref_libs)
        mantissas = mantissas * lib_m if ratio else mantissas / lib_m
        exponents = exponents + lib_e if ratio else exponents - lib_e
    mantissas, shift = numpy.frexp(mantissas)
    exponents = exponents + shift
    order = numpy.lexsort((mantissas, exponents))
    different = (exponents[order][1:] != exponents[order][:-1]) | (mantissas[order][1:] != mantissas[order][:-1])
    starts = numpy.r_[0, numpy.flatnonzero(different) + 1]
    ends = numpy.r_[starts[1:], len(order)]
    ranks = numpy.empty(len(order), dtype=float)
    ranks[order] = numpy.repeat((starts + ends + 1) / 2.0, ends - starts)
    return ranks


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
    if (not numpy.isfinite(obs_array).all()) or (not numpy.isfinite(ref_array).all()):
        raise ValueError('TMM count vectors must contain only finite values.')
    if (obs_array < 0).any() or (ref_array < 0).any():
        raise ValueError('TMM count vectors must not contain negative values.')
    obs_vector = obs_array.reshape(-1, order='F')
    ref_vector = ref_array.reshape(-1, order='F')
    if libsize_obs is None:
        n_obs = float(obs_vector.sum())
    else:
        n_obs_values = numpy.asarray(libsize_obs, dtype=float).reshape(-1)
        if n_obs_values.size == 0:
            raise ValueError('libsize_obs must contain one finite positive value.')
        n_obs = float(n_obs_values[0])
    if libsize_ref is None:
        n_ref_values = numpy.array([float(ref_vector.sum())], dtype=float)
    else:
        n_ref_values = numpy.asarray(libsize_ref, dtype=float).reshape(-1)
        if n_ref_values.size == 0:
            n_ref_values = numpy.array([float(ref_vector.sum())], dtype=float)
    if (not numpy.isfinite(n_obs)) or n_obs <= 0:
        raise ValueError('libsize_obs must be finite and positive.')
    if (not numpy.isfinite(n_ref_values).all()) or (n_ref_values <= 0).any():
        raise ValueError('libsize_ref must contain only finite positive values.')
    target_size = max(obs_vector.size, ref_vector.size)
    if not obs_vector.size or not ref_vector.size:
        return 1.0
    obs_counts = _recycle(obs_vector, target_size)
    ref_counts = _recycle(ref_vector, target_size)
    ref_libs = _recycle(n_ref_values, target_size)
    positive = (obs_counts > 0) & (ref_counts > 0)
    obs_counts, ref_counts, ref_libs = obs_counts[positive], ref_counts[positive], ref_libs[positive]
    log_obs = numpy.log2(obs_counts) - math.log2(n_obs)
    log_ref = numpy.log2(ref_counts) - numpy.log2(ref_libs)
    log_r = log_obs - log_ref
    abs_e = (log_obs + log_ref) / 2.0
    with numpy.errstate(over='ignore', divide='ignore', invalid='ignore'):
        variances = (n_obs - obs_counts) / n_obs / obs_counts + (ref_libs - ref_counts) / ref_libs / ref_counts
    finite = numpy.isfinite(log_r) & numpy.isfinite(abs_e) & (abs_e > acutoff)
    log_r = log_r[finite]
    variances = variances[finite]
    obs_counts, ref_counts, ref_libs = obs_counts[finite], ref_counts[finite], ref_libs[finite]
    if log_r.size == 0:
        return 1.0
    if float(numpy.max(numpy.abs(log_r))) < 1e-6:
        return 1.0
    n = log_r.size
    lo_l = int(math.floor(n * logratio_trim) + 1)
    hi_l = int(n + 1 - lo_l)
    lo_s = int(math.floor(n * sum_trim) + 1)
    hi_s = int(n + 1 - lo_s)
    rank_log_r = _rank_count_expression(obs_counts, ref_counts, ref_libs, ratio=True)
    rank_abs_e = _rank_count_expression(obs_counts, ref_counts, ref_libs, ratio=False)
    keep = (
        (rank_log_r >= lo_l) & (rank_log_r <= hi_l) &
        (rank_abs_e >= lo_s) & (rank_abs_e <= hi_s)
    )
    if not keep.any():
        return 1.0
    with numpy.errstate(divide='ignore', invalid='ignore'):
        if do_weighting:
            numerator = numpy.nansum(log_r[keep] / variances[keep])
            denominator = numpy.nansum(1.0 / variances[keep])
            factor_log = numerator / denominator
        else:
            factor_log = float(numpy.mean(log_r[keep]))
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
        # Compare distances to the exact mean of the represented quartiles;
        # forming a rounded mean can break an otherwise exact reference tie.
        quartiles = [Fraction(float(value)) for value in f75]
        mean_quartile = sum(quartiles) / len(quartiles)
        ref_column = min(range(len(quartiles)), key=lambda index: abs(quartiles[index] - mean_quartile))
    return ref_column


def _resolve_median_reference_column(factors):
    values = numpy.asarray(factors, dtype=float).reshape(-1)
    if values.size == 0:
        raise ValueError('At least one TMM normalization factor is required.')
    if not numpy.isfinite(values).all() or (values <= 0).any():
        raise ValueError('TMM normalization factors must be finite and positive.')
    ordered = numpy.sort(values)
    central = ordered[[(values.size - 1) // 2, values.size // 2]]
    # For even n both central values are equally close to the true median.
    # Choose the first input column, without subtracting a rounded midpoint.
    return int(numpy.flatnonzero(numpy.isin(values, central))[0])


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
    if matrix.shape[1] == 0:
        raise ValueError('TMM normalization requires at least one sample column.')
    if _remove_all_zero_rows(matrix).shape[0] == 0:
        raise ValueError('TMM normalization requires at least one positive count.')
    libs = _coerce_library_sizes(lib_size, matrix, columns)
    round1_counts = counts if isinstance(counts, pandas.DataFrame) else pandas.DataFrame(matrix, columns=columns)
    round1 = calc_norm_factors_tmm(round1_counts, lib_size=libs, ref_column=None)
    round1_values = round1.to_numpy(dtype=float)
    round1_reference_column = _resolve_tmm_reference_column(_remove_all_zero_rows(matrix), libs)
    median_reference_columns = [_resolve_median_reference_column(round1_values)]
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
    _matrix, _columns = _as_matrix_with_columns(counts)
    if isinstance(norm_factors, pandas.Series):
        factors = norm_factors.reindex(counts.columns).to_numpy(dtype=float)
    else:
        factors = numpy.asarray(norm_factors, dtype=float).reshape(-1)
        if factors.size != counts.shape[1]:
            raise ValueError('norm_factors must match the number of sample columns.')
    if (not numpy.isfinite(factors).all()) or (factors <= 0).any():
        raise ValueError('norm_factors must contain only finite positive values.')
    return pandas.DataFrame(_matrix / factors.reshape(1, -1), index=counts.index, columns=counts.columns)


__all__ = [
    'TMMRoundTripResult',
    'apply_tmm_factors',
    'calc_factor_quantile',
    'calc_factor_tmm',
    'calc_norm_factors_tmm',
    'run_tmm_rounds_for_cstmm',
]
