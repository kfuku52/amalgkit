from types import SimpleNamespace

import numpy
import pandas
import pytest

from amalgkit.per_species_python import (
    _apply_within_group_filter,
    _compute_sample_group_correlation_metrics,
    _should_stop_within_group_filter,
)


def test_within_group_reference_leaves_evaluated_sample_out():
    counts = pandas.DataFrame(
        {
            'A_good': [5.0, 4.0, 3.0, 2.0, 1.0],
            'A_bad': [1.0, 2.0, 3.0, 4.0, 5.0],
            'B1': [1.0, 2.0, 3.0, 4.0, 5.0],
            'B2': [1.0, 2.0, 3.0, 4.0, 5.0],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5'],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['A_good', 'A_bad', 'B1', 'B2'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['P1', 'P1', 'P2', 'P2'],
            'exclusion': ['no', 'no', 'no', 'no'],
        }
    )

    observed = _compute_sample_group_correlation_metrics(
        tc=counts,
        sra=metadata,
        selected_sample_groups=['A', 'B'],
        dist_method='pearson',
    ).set_index('run')

    assert observed.loc['A_bad', 'ws_within_group_cor'] == pytest.approx(-1.0)
    assert observed.loc['A_bad', 'ws_max_nongroup_cor'] == pytest.approx(1.0)
    assert observed.loc['A_bad', 'ws_margin'] == pytest.approx(-2.0)


def test_within_group_correlation_rejects_unknown_method():
    counts = pandas.DataFrame({'A1': [1.0, 2.0], 'B1': [2.0, 1.0]})
    metadata = pandas.DataFrame(
        {'run': ['A1', 'B1'], 'sample_group': ['A', 'B']}
    )

    with pytest.raises(ValueError, match='Unsupported correlation method'):
        _compute_sample_group_correlation_metrics(
            tc=counts,
            sra=metadata,
            selected_sample_groups=['A', 'B'],
            dist_method='typo',
        )


@pytest.mark.parametrize('method', ['pearson', 'spearman', 'kendall'])
@pytest.mark.parametrize('reference_exclusion', ['run', 'bioproject'])
@pytest.mark.parametrize('missing', [False, True])
def test_correlation_metrics_match_direct_pairwise_references(method, reference_exclusion, missing):
    rng = numpy.random.default_rng(719)
    counts = pandas.DataFrame(rng.normal(size=(521, 6)), columns=['a1', 'a2', 'a3', 'b1', 'b2', 'b3'])
    counts['b3'] = 1.0
    if missing:
        counts.loc[::7, 'a1'] = numpy.nan
        counts.loc[::11, 'b1'] = numpy.inf
    metadata = pandas.DataFrame(dict(
        run=counts.columns, sample_group=['A'] * 3 + ['B'] * 3,
        bioproject=['P1', 'P1', 'P2', 'P3', 'P3', 'P4'], exclusion='no',
    ))
    original = counts.copy(deep=True)
    result = _compute_sample_group_correlation_metrics(
        counts, metadata, ['A', 'B'], method, reference_exclusion=reference_exclusion,
    ).set_index('run')
    for row in metadata.itertuples(index=False):
        correlations, common = {}, {}
        for group in ['A', 'B']:
            members = metadata.loc[metadata.sample_group.eq(group)]
            if reference_exclusion == 'bioproject':
                members = members.loc[members.bioproject.ne(row.bioproject)]
            else:
                members = members.loc[members.run.ne(row.run)]
            reference = counts[members.run].mean(axis=1)
            valid = numpy.isfinite(counts[row.run]) & numpy.isfinite(reference)
            common[group] = int(valid.sum())
            left, right = counts.loc[valid, row.run], reference.loc[valid]
            correlations[group] = (
                left.corr(right, method=method) if left.nunique() > 1 and right.nunique() > 1 else numpy.nan
            )
        other = 'B' if row.sample_group == 'A' else 'A'
        numpy.testing.assert_allclose(
            result.loc[row.run, ['ws_within_group_cor', 'ws_max_nongroup_cor']].to_numpy(dtype=float),
            [correlations[row.sample_group], correlations[other]], rtol=1e-14, atol=1e-14,
        )
        assert result.loc[row.run, 'ws_within_common_genes'] == common[row.sample_group]
        assert result.loc[row.run, 'ws_min_nongroup_common_genes'] == common[other]
    pandas.testing.assert_frame_equal(counts, original)


def test_within_group_filter_stops_when_an_iteration_makes_no_progress():
    current = pandas.DataFrame({'RUN1': [1.0], 'RUN2': [2.0]}, index=['G1'])
    unchanged = current.copy()

    assert _should_stop_within_group_filter(
        current_tc=current,
        next_tc=unchanged,
        excluded_runs=['RUN1'],
    )


def test_within_group_filter_clears_stale_candidates_after_all_runs_removed():
    counts = pandas.DataFrame(
        {
            'A1': [1.0, 2.0, 3.0, 4.0],
            'A2': [1.0, 2.0, 3.0, 4.0],
            'B1': [1.0, 2.0, 3.0, 4.0],
            'B2': [1.0, 2.0, 3.0, 4.0],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['A1', 'A2', 'B1', 'B2'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['P1', 'P1', 'P2', 'P2'],
            'exclusion': ['no', 'no', 'no', 'no'],
        }
    )
    args = SimpleNamespace(
        dist_method='pearson',
        margin_threshold=0.1,
        robust_z_threshold=-2.5,
        one_outlier_per_iter=False,
    )

    empty_counts, excluded_metadata, excluded_runs = _apply_within_group_filter(
        tc=counts,
        sra=metadata,
        args=args,
        selected_sample_groups=['A', 'B'],
    )
    assert empty_counts.shape[1] == 0
    assert excluded_runs == ['A1', 'A2', 'B1', 'B2']
    assert excluded_metadata['ws_small_group'].tolist() == [True, True, True, True]

    next_counts, next_metadata, next_excluded_runs = _apply_within_group_filter(
        tc=empty_counts,
        sra=excluded_metadata,
        args=args,
        selected_sample_groups=['A', 'B'],
    )

    assert next_counts.shape[1] == 0
    assert next_excluded_runs == []
    # No new removals, but the evidence from the removal round survives.
    assert next_metadata['ws_outlier_candidate'].fillna(False).all()
    assert next_metadata['ws_small_group'].fillna(False).all()
