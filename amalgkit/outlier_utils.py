import numpy
import pandas


def compute_groupwise_robust_z(values, groups, mad_scale=1.4826, min_group_size=3):
    values_arr = pandas.to_numeric(pandas.Series(values), errors='coerce').to_numpy(dtype=float)
    groups_arr = pandas.Series(groups).fillna('').astype(str).str.strip().to_numpy(dtype=object)
    out = numpy.full((values_arr.shape[0],), numpy.nan, dtype=float)
    if values_arr.shape[0] == 0:
        return out
    valid_group = groups_arr != ''
    group_levels = list(dict.fromkeys(groups_arr[valid_group].tolist()))
    for group_value in group_levels:
        idx = numpy.where(valid_group & (groups_arr == group_value) & numpy.isfinite(values_arr))[0]
        if idx.size == 0:
            continue
        if idx.size < int(min_group_size):
            out[idx] = numpy.nan
            continue
        x = values_arr[idx]
        median_value = float(numpy.median(x))
        mad_raw = float(numpy.median(numpy.abs(x - median_value)))
        scale = mad_raw * float(mad_scale)
        if (not numpy.isfinite(scale)) or (scale <= 0):
            sd_value = float(numpy.std(x, ddof=1)) if idx.size > 1 else numpy.nan
            if numpy.isfinite(sd_value) and (sd_value > 0):
                scale = sd_value
            else:
                out[idx] = 0.0
                continue
        out[idx] = (x - median_value) / scale
    return out


def flag_margin_outliers(
    df,
    margin_col,
    group_col='sample_group',
    margin_threshold=0.0,
    robust_z_threshold=-2.5,
    min_group_size=3,
    robust_z_col='robust_z',
    outlier_col='outlier_flag',
):
    if margin_col not in df.columns:
        raise ValueError('Margin column not found: {}'.format(margin_col))
    if group_col not in df.columns:
        raise ValueError('Group column not found: {}'.format(group_col))
    out = df.copy()
    margin_values = pandas.to_numeric(out.loc[:, margin_col], errors='coerce').to_numpy(dtype=float)
    robust_z = compute_groupwise_robust_z(
        values=margin_values,
        groups=out.loc[:, group_col],
        min_group_size=min_group_size,
    )
    is_finite_margin = numpy.isfinite(margin_values)
    outlier_flag = (
        is_finite_margin
        & (margin_values < float(margin_threshold))
        & ((robust_z <= float(robust_z_threshold)) | numpy.isnan(robust_z))
    )
    out.loc[:, robust_z_col] = robust_z
    out.loc[:, outlier_col] = outlier_flag
    return out


__all__ = [
    'compute_groupwise_robust_z',
    'flag_margin_outliers',
]
