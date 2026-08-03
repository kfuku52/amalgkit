import warnings

import numpy
import pandas
import pytest

from amalgkit.outlier_utils import compute_groupwise_robust_z, flag_margin_outliers


def _small_group_frame():
    # Two 2-sample groups: below the default min_group_size of 3, so the robust
    # z-score is NaN for every row.
    return pandas.DataFrame(
        {
            'sample_group': ['g1', 'g1', 'g2', 'g2'],
            'ws_margin': [-0.1, -0.2, 0.0, 0.4],
        }
    )


def test_compute_groupwise_robust_z_small_group_is_nan():
    values = [-0.1, -0.2, 0.3, -0.4, 1.2]
    groups = ['g1', 'g1', 'g2', 'g2', 'g2']
    z = compute_groupwise_robust_z(values, groups, min_group_size=3)
    # g1 has 2 members (< min_group_size) -> not scoreable.
    assert numpy.isnan(z[0])
    assert numpy.isnan(z[1])
    # g2 has 3 members -> finite z-scores.
    assert numpy.isfinite(z[2:]).all()


def test_flag_margin_outliers_default_policy_falls_back_to_absolute_margin():
    # Established behaviour (and the pre-port R behaviour): when the robust
    # z-score is unavailable, the absolute margin threshold decides on its own.
    # At the default margin_threshold=0.0, a negative margin still flags.
    df = _small_group_frame()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = flag_margin_outliers(df, margin_col='ws_margin', group_col='sample_group', min_group_size=3)

    assert out.loc[:, 'outlier_flag'].tolist() == [True, True, False, False]
    assert out.loc[:, 'robust_z'].isna().all()


def test_flag_margin_outliers_default_policy_respects_explicit_threshold():
    # An explicitly supplied --margin_threshold must stay effective for small
    # groups. A margin of exactly 0.0 does not cross the default threshold
    # (0 < 0 is false) but does cross an explicit 0.1.
    df = _small_group_frame()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        default_out = flag_margin_outliers(
            df, margin_col='ws_margin', group_col='sample_group', min_group_size=3
        )
        explicit_out = flag_margin_outliers(
            df, margin_col='ws_margin', group_col='sample_group', min_group_size=3,
            margin_threshold=0.1,
        )

    # g2's first sample has margin 0.0: not flagged by default, flagged at 0.1.
    assert not bool(default_out.loc[2, 'outlier_flag'])
    assert bool(explicit_out.loc[2, 'outlier_flag'])
    # The 0.4 margin crosses neither threshold.
    assert not bool(explicit_out.loc[3, 'outlier_flag'])


def test_flag_margin_outliers_retain_policy_disables_small_group_screening():
    # The opt-in policy: a group that cannot be robustly scored is never
    # flagged, at the default threshold or an explicit one.
    df = _small_group_frame()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        default_out = flag_margin_outliers(
            df, margin_col='ws_margin', group_col='sample_group', min_group_size=3,
            small_group_policy='retain',
        )
        explicit_out = flag_margin_outliers(
            df, margin_col='ws_margin', group_col='sample_group', min_group_size=3,
            margin_threshold=0.1, small_group_policy='retain',
        )

    assert not default_out.loc[:, 'outlier_flag'].any()
    assert not explicit_out.loc[:, 'outlier_flag'].any()


def test_flag_margin_outliers_reports_unscoreable_groups():
    # Issue #171 asks for the affected groups to be reported rather than
    # handled silently, whichever policy is in force.
    df = _small_group_frame()

    with pytest.warns(UserWarning, match='could not be robustly scored') as record:
        out = flag_margin_outliers(
            df, margin_col='ws_margin', group_col='sample_group', min_group_size=3
        )

    message = str(record[0].message)
    assert 'g1' in message
    assert 'g2' in message
    assert 'margin_fallback' in message
    assert out.loc[:, 'small_group'].tolist() == [True, True, True, True]


def test_flag_margin_outliers_does_not_report_when_every_group_is_scoreable():
    # Negative control: a report that fires on every run carries no signal.
    df = pandas.DataFrame(
        {
            'sample_group': ['g1'] * 5,
            'ws_margin': [-3.0, -0.1, -0.2, 0.5, 1.0],
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        out = flag_margin_outliers(df, margin_col='ws_margin', group_col='sample_group', min_group_size=3)

    assert not out.loc[:, 'small_group'].any()


def test_flag_margin_outliers_large_group_flags_only_crossing_samples():
    df = pandas.DataFrame(
        {
            'sample_group': ['g1'] * 5,
            'ws_margin': [-3.0, -0.1, -0.2, 0.5, 1.0],
        }
    )
    out = flag_margin_outliers(df, margin_col='ws_margin', group_col='sample_group', min_group_size=3)
    # Scoreable group: both the margin and the z-score must cross. The policy
    # is irrelevant here, so it must not change the result.
    retained = flag_margin_outliers(
        df, margin_col='ws_margin', group_col='sample_group', min_group_size=3,
        small_group_policy='retain',
    )
    assert out.loc[:, 'outlier_flag'].tolist() == [True, False, False, False, False]
    assert retained.loc[:, 'outlier_flag'].tolist() == out.loc[:, 'outlier_flag'].tolist()


def test_flag_margin_outliers_nonfinite_margin_never_flags():
    df = pandas.DataFrame(
        {
            'sample_group': ['g1', 'g1', 'g1', 'g1'],
            'ws_margin': [numpy.nan, numpy.inf, -numpy.inf, -0.3],
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = flag_margin_outliers(df, margin_col='ws_margin', group_col='sample_group', min_group_size=3)
    # Non-finite margins are never outliers and are never counted as an
    # unscoreable small group.
    assert not out.loc[:, 'outlier_flag'].iloc[0:3].any()
    assert not out.loc[:, 'small_group'].iloc[0:3].any()


def test_flag_margin_outliers_rejects_unknown_policy():
    df = _small_group_frame()
    with pytest.raises(ValueError, match='Unknown small_group_policy'):
        flag_margin_outliers(
            df, margin_col='ws_margin', group_col='sample_group',
            small_group_policy='keep_everything',
        )


def test_flag_margin_outliers_returns_copies():
    df = pandas.DataFrame(
        {
            'sample_group': ['g1', 'g1', 'g1'],
            'ws_margin': [-0.1, 0.2, 0.3],
        }
    )
    out = flag_margin_outliers(df, margin_col='ws_margin', group_col='sample_group')
    # The input frame is not mutated.
    assert list(df.columns) == ['sample_group', 'ws_margin']
    assert 'outlier_flag' not in df.columns
    assert 'outlier_flag' in out.columns
