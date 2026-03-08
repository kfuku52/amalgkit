import numpy
import pandas

from amalgkit.batch_effect_ruvseq import run_ruvseq_backend, ruvr_correct_counts


def test_ruvr_correct_counts_preserves_shape_and_nonnegative_values():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [101.0, 95.0, 11.0, 10.0],
            'RUN2': [98.0, 94.0, 13.0, 9.0],
            'RUN3': [14.0, 15.0, 103.0, 100.0],
            'RUN4': [12.0, 10.0, 98.0, 97.0],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    residuals_df = pandas.DataFrame(
        {
            'RUN1': [0.2, -0.1, 0.3, -0.2],
            'RUN2': [0.1, 0.0, 0.4, -0.1],
            'RUN3': [-0.3, 0.2, -0.2, 0.1],
            'RUN4': [-0.2, 0.1, -0.1, 0.2],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    controls = numpy.array([True, True, False, True], dtype=bool)

    corrected_df, w_df = ruvr_correct_counts(
        seq_uq_df=counts_df,
        controls=controls,
        k=1,
        residuals_df=residuals_df,
    )

    assert corrected_df.shape == counts_df.shape
    assert list(corrected_df.columns) == list(counts_df.columns)
    assert numpy.all(corrected_df.to_numpy(dtype=float) >= 0)
    assert w_df.shape == (counts_df.shape[1], 1)
    assert list(w_df.index) == list(counts_df.columns)


def test_run_ruvseq_backend_returns_nonnegative_corrected_matrix_and_summary():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [100.0, 110.0, 12.0, 10.0, 9.0],
            'RUN2': [95.0, 105.0, 15.0, 11.0, 8.0],
            'RUN3': [11.0, 13.0, 102.0, 98.0, 7.0],
            'RUN4': [10.0, 12.0, 99.0, 96.0, 6.0],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2'],
        }
    )

    corrected_df, w_df, summary = run_ruvseq_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        control_mode='auto',
        k_setting='1',
        k_max=3,
        top_n=5,
        min_controls=2,
    )

    assert corrected_df.shape == counts_df.shape
    assert list(corrected_df.columns) == list(counts_df.columns)
    assert numpy.all(corrected_df.to_numpy(dtype=float) >= 0)
    assert int(summary['resolved_ruv_k']) == 1
    assert int(summary['resolved_ruv_controls']) >= 2
    assert summary['skip_reason'] == ''
    assert list(w_df.index) == list(counts_df.columns)


def test_run_ruvseq_backend_auto_k_selects_positive_k_for_balanced_fixture():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [100.0, 110.0, 12.0, 10.0, 9.0, 8.0],
            'RUN2': [95.0, 105.0, 15.0, 11.0, 8.0, 7.0],
            'RUN3': [11.0, 13.0, 102.0, 98.0, 7.0, 6.0],
            'RUN4': [10.0, 12.0, 99.0, 96.0, 6.0, 5.0],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5', 'G6'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP2', 'BP1', 'BP2'],
        }
    )

    corrected_df, _w_df, summary = run_ruvseq_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        control_mode='auto',
        k_setting='auto',
        k_max=3,
        top_n=6,
        min_controls=2,
    )

    assert corrected_df.shape == counts_df.shape
    assert int(summary['resolved_ruv_k']) >= 1
    assert int(summary['resolved_ruv_k']) <= 3
    assert int(summary['resolved_ruv_controls']) >= 2


def test_run_ruvseq_backend_single_group_design_skip():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [100.0, 110.0, 12.0, 10.0],
            'RUN2': [95.0, 105.0, 15.0, 11.0],
            'RUN3': [11.0, 13.0, 102.0, 98.0],
            'RUN4': [10.0, 12.0, 99.0, 96.0],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'A', 'A', 'A'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2'],
        }
    )

    corrected_df, w_df, summary = run_ruvseq_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        control_mode='auto',
        k_setting='auto',
        k_max=3,
        top_n=4,
        min_controls=2,
    )

    pandas.testing.assert_frame_equal(corrected_df, counts_df)
    assert summary['skip_reason'] == 'ruvseq_design_failed'
    assert pandas.isna(summary['resolved_ruv_k'])
    assert pandas.isna(summary['resolved_ruv_controls'])
    assert w_df.shape == (counts_df.shape[1], 0)
