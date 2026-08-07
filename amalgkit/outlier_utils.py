import warnings

import numpy
import pandas


SMALL_GROUP_POLICIES = ('margin_fallback', 'retain')


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
    small_group_policy='margin_fallback',
    small_group_col='small_group',
    report_small_groups=True,
):
    """Flag low-margin samples, optionally falling back to the absolute margin.

    A sample is flagged when its margin is finite, below ``margin_threshold``,
    and its robust z-score crosses ``robust_z_threshold``.

    ``compute_groupwise_robust_z`` cannot score a group with fewer than
    ``min_group_size`` members and returns NaN for it. Two-replicate groups are
    common in cross-species RNA-seq, so how that NaN is treated decides whether
    those groups get any margin-based screening at all. ``small_group_policy``
    makes the choice explicit:

    ``'margin_fallback'`` (default)
        When the robust z-score is unavailable, fall back to the absolute
        ``margin_threshold`` alone. This is the established behaviour and
        matches the pre-port R implementation, so an explicitly supplied
        ``--margin_threshold`` keeps working for small groups.
    ``'retain'``
        Never flag a sample whose group could not be robustly scored. Margin
        based screening is disabled for small groups, and ``margin_threshold``
        has no effect there.

    Groups that could not be robustly scored are marked in ``small_group_col``
    and, when ``report_small_groups`` is set, named in a warning.
    """
    if margin_col not in df.columns:
        raise ValueError('Margin column not found: {}'.format(margin_col))
    if group_col not in df.columns:
        raise ValueError('Group column not found: {}'.format(group_col))
    policy = str(small_group_policy).strip().lower()
    if policy not in SMALL_GROUP_POLICIES:
        raise ValueError(
            'Unknown small_group_policy: {} (expected one of {})'.format(
                small_group_policy, ', '.join(SMALL_GROUP_POLICIES)
            )
        )
    out = df.copy()
    margin_values = pandas.to_numeric(out.loc[:, margin_col], errors='coerce').to_numpy(dtype=float)
    group_values = pandas.Series(out.loc[:, group_col]).fillna('').astype(str).str.strip()
    robust_z = compute_groupwise_robust_z(
        values=margin_values,
        groups=group_values,
        min_group_size=min_group_size,
    )
    is_finite_margin = numpy.isfinite(margin_values)
    is_finite_z = numpy.isfinite(robust_z)
    # A sample with a usable margin and a named group but no robust z-score is
    # one whose group was too small to score.
    is_small_group = is_finite_margin & (~is_finite_z) & (group_values.to_numpy(dtype=object) != '')

    crosses_margin = is_finite_margin & (margin_values < float(margin_threshold))
    crosses_z = is_finite_z & (robust_z <= float(robust_z_threshold))
    if policy == 'margin_fallback':
        outlier_flag = crosses_margin & (crosses_z | is_small_group)
    else:
        outlier_flag = crosses_margin & crosses_z

    if report_small_groups and bool(is_small_group.any()):
        affected_groups = sorted(set(group_values.to_numpy(dtype=object)[is_small_group].tolist()))
        warnings.warn(
            '{} sample(s) in {} sample_group(s) could not be robustly scored because the group '
            'has fewer than {} usable members: {}. Under small_group_policy="{}", {}.'.format(
                int(is_small_group.sum()),
                len(affected_groups),
                int(min_group_size),
                ', '.join(affected_groups),
                policy,
                'they are screened on the absolute margin threshold alone'
                if policy == 'margin_fallback'
                else 'they are retained without margin-based screening',
            )
        )

    out.loc[:, robust_z_col] = robust_z
    out.loc[:, outlier_col] = outlier_flag
    if small_group_col is not None:
        out.loc[:, small_group_col] = is_small_group
    return out


__all__ = [
    'SMALL_GROUP_POLICIES',
    'compute_groupwise_robust_z',
    'flag_margin_outliers',
]
