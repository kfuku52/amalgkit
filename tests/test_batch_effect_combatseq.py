import pytest
import pandas

from amalgkit.batch_effect_combatseq import run_combatseq_backend


pytest.importorskip('inmoose.pycombat')


def _balanced_counts():
    return pandas.DataFrame(
        {
            'RUN1': [100, 101, 102, 103, 10, 11, 12, 13, 50, 51],
            'RUN2': [95, 96, 97, 98, 15, 16, 17, 18, 55, 56],
            'RUN3': [90, 91, 92, 93, 20, 21, 22, 23, 40, 41],
            'RUN4': [85, 86, 87, 88, 25, 26, 27, 28, 45, 46],
        },
        index=[f'G{i}' for i in range(10)],
    )


def test_run_combatseq_backend_matches_expected_balanced_group_case():
    counts = _balanced_counts()
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'B', 'A', 'B'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2'],
        }
    )
    corrected, summary = run_combatseq_backend(counts_df=counts, metadata_df=metadata)
    expected = pandas.DataFrame(
        {
            'RUN1': [97, 98, 99, 100, 14, 15, 16, 18, 46, 47],
            'RUN2': [92, 93, 94, 95, 21, 22, 23, 24, 51, 52],
            'RUN3': [94, 95, 96, 97, 15, 16, 17, 18, 44, 45],
            'RUN4': [89, 90, 91, 92, 19, 20, 21, 22, 50, 51],
        },
        index=counts.index,
    )
    pandas.testing.assert_frame_equal(corrected, expected)
    assert summary['method'] == 'group'
    assert summary['skip_reason'] == ''
    assert summary['corrected_run_ids'] == ['RUN1', 'RUN2', 'RUN3', 'RUN4']
    assert summary['uncorrected_run_ids'] == []
    assert summary['group_model_used'] is True
    assert summary['group_fallback_used'] is False


def test_run_combatseq_backend_falls_back_without_group_when_confounded():
    counts = pandas.DataFrame(
        {
            'RUN1': [100, 101, 102, 103, 10, 11, 12, 13, 50, 51],
            'RUN2': [100, 101, 102, 103, 10, 11, 12, 13, 50, 51],
            'RUN3': [5, 6, 7, 8, 100, 101, 102, 103, 40, 41],
            'RUN4': [5, 6, 7, 8, 100, 101, 102, 103, 40, 41],
        },
        index=[f'G{i}' for i in range(10)],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2'],
        }
    )
    corrected, summary = run_combatseq_backend(counts_df=counts, metadata_df=metadata)
    expected = pandas.DataFrame(
        {
            'RUN1': [24, 26, 28, 30, 33, 35, 37, 38, 47, 48],
            'RUN2': [24, 26, 28, 30, 33, 35, 37, 38, 47, 48],
            'RUN3': [22, 24, 26, 28, 31, 33, 34, 36, 44, 45],
            'RUN4': [22, 24, 26, 28, 31, 33, 34, 36, 44, 45],
        },
        index=counts.index,
    )
    pandas.testing.assert_frame_equal(corrected, expected)
    assert summary['method'] == 'no_group'
    assert summary['group_model_used'] is False
    assert summary['group_fallback_used'] is True
    assert summary['skip_reason'] == ''


def test_run_combatseq_backend_keeps_all_singleton_batches_uncorrected():
    counts = pandas.DataFrame(
        {
            'RUN1': [10, 20],
            'RUN2': [30, 40],
            'RUN3': [50, 60],
        },
        index=['G1', 'G2'],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3'],
            'sample_group': ['A', 'A', 'B'],
            'bioproject': ['BP1', 'BP2', 'BP3'],
        }
    )
    corrected, summary = run_combatseq_backend(counts_df=counts, metadata_df=metadata)
    pandas.testing.assert_frame_equal(corrected, counts)
    assert summary['method'] == 'all_singleton'
    assert summary['skip_reason'] == 'combatseq_all_singleton'
    assert summary['corrected_run_ids'] == []
    assert summary['uncorrected_run_ids'] == ['RUN1', 'RUN2', 'RUN3']


def test_run_combatseq_backend_keeps_singleton_batch_uncorrected():
    counts = pandas.DataFrame(
        {
            'RUN1': [100, 120, 10, 12],
            'RUN2': [90, 110, 20, 22],
            'RUN3': [95, 115, 15, 17],
            'RUN4': [85, 105, 25, 27],
            'RUN5': [60, 80, 30, 32],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4', 'RUN5'],
            'sample_group': ['A', 'B', 'A', 'B', 'A'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2', 'BP3'],
        }
    )
    corrected, summary = run_combatseq_backend(counts_df=counts, metadata_df=metadata)
    assert corrected.loc[:, 'RUN5'].tolist() == counts.loc[:, 'RUN5'].tolist()
    assert summary['skip_reason'] == 'combatseq_singleton_kept'
    assert summary['uncorrected_run_ids'] == ['RUN5']
    assert summary['corrected_run_ids'] == ['RUN1', 'RUN2', 'RUN3', 'RUN4']
