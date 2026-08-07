from types import SimpleNamespace

import pandas

from amalgkit.per_species_python import (
    _apply_within_group_filter,
    _should_stop_within_group_filter,
)


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
    assert not next_metadata['ws_outlier_candidate'].fillna(False).any()
    assert not next_metadata['ws_small_group'].fillna(False).any()
