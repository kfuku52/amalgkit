import numpy
import pandas

from amalgkit.batch_effect_ruvseq import (
    _between_lane_normalize_upper,
    run_ruvseq_backend,
    ruvr_correct_counts,
)


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


def test_run_ruvseq_backend_rank_zero_returns_raw_counts_and_reports_k_zero():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [10.0, 20.0, 30.0],
            'RUN2': [10.0, 20.0, 30.0],
            'RUN3': [10.0, 20.0, 30.0],
            'RUN4': [10.0, 20.0, 30.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': list(counts_df.columns),
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2'],
        }
    )

    for k_setting in ('1', 'auto'):
        corrected_df, w_df, summary = run_ruvseq_backend(
            counts_df=counts_df,
            metadata_df=metadata_df,
            k_setting=k_setting,
            k_max=2,
            min_controls=2,
        )

        pandas.testing.assert_frame_equal(corrected_df, counts_df)
        assert w_df.shape == (counts_df.shape[1], 0)
        assert summary['resolved_ruv_k'] == 0
        assert summary['skip_reason'] == 'ruvseq_k_zero'
        assert summary['corrected_run_ids'] == []


def test_run_ruvseq_backend_empty_gene_table_returns_no_expressed_genes_skip():
    counts_df = pandas.DataFrame(columns=['RUN1', 'RUN2'], dtype=float)
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2'],
            'sample_group': ['A', 'B'],
            'bioproject': ['BP1', 'BP2'],
        }
    )

    corrected_df, w_df, summary = run_ruvseq_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
    )

    pandas.testing.assert_frame_equal(corrected_df, counts_df)
    assert list(w_df.index) == ['RUN1', 'RUN2']
    assert w_df.shape == (2, 0)
    assert summary['skip_reason'] == 'no_expressed_genes'
    assert summary['corrected_run_ids'] == []
    assert summary['uncorrected_run_ids'] == ['RUN1', 'RUN2']


def test_run_ruvseq_backend_no_op_correction_does_not_inflate_counts():
    # Regression test: the +1 pseudo-count used for the GLM fit must not leak
    # into the corrected matrix. When the correction is a no-op (W@alpha ~ 0),
    # the corrected counts must equal the raw counts, not raw + 1.
    counts_df = pandas.DataFrame(
        {
            'RUN1': [10.0, 20.0, 30.0],
            'RUN2': [10.0, 20.0, 30.0],
            'RUN3': [10.0, 20.0, 30.0],
            'RUN4': [10.0, 20.0, 30.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': list(counts_df.columns),
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2'],
        }
    )

    # A forced k>0 with numerically-zero residuals must still return the raw
    # counts (the noise-floor guard resolves k to 0 instead of fitting on float
    # noise from the GLM fit).
    corrected_df, w_df, summary = run_ruvseq_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        k_setting='1',
        k_max=2,
        min_controls=2,
    )

    pandas.testing.assert_frame_equal(corrected_df, counts_df)
    assert w_df.shape == (counts_df.shape[1], 0)
    assert summary['resolved_ruv_k'] == 0
    assert summary['skip_reason'] == 'ruvseq_k_zero'


def test_ruvr_correct_counts_returns_raw_counts_on_noise_floor_residuals():
    # Direct unit test for the noise-floor guard: residuals at float-noise
    # magnitude must be treated as "no batch effect" (k resolves to 0) rather
    # than producing a spurious correction.
    counts_df = pandas.DataFrame(
        {
            'RUN1': [10.0, 20.0, 30.0],
            'RUN2': [10.0, 20.0, 30.0],
            'RUN3': [10.0, 20.0, 30.0],
            'RUN4': [10.0, 20.0, 30.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    residuals_df = pandas.DataFrame(
        {
            'RUN1': [1e-8, -1e-8, 1e-8],
            'RUN2': [-1e-8, 1e-8, -1e-8],
            'RUN3': [1e-8, 1e-8, -1e-8],
            'RUN4': [-1e-8, -1e-8, 1e-8],
        },
        index=counts_df.index,
    )
    controls = numpy.array([True, True, True], dtype=bool)

    corrected_df, w_df = ruvr_correct_counts(
        seq_uq_df=counts_df,
        controls=controls,
        k=1,
        residuals_df=residuals_df,
    )

    pandas.testing.assert_frame_equal(corrected_df, counts_df)
    assert w_df.shape == (counts_df.shape[1], 0)


def test_run_ruvseq_backend_restores_raw_count_scale_with_unequal_library_sizes():
    # Regression for the between-lane-scale finding: the RUVr correction runs
    # on the upper-quartile normalized matrix, but run_ruvseq_backend must
    # multiply back by the per-lane scale factors so the returned "corrected
    # counts" stay on the original count scale. Previously the corrected
    # matrix was returned on the normalized scale (silently rescaling every
    # count ~2x on high-depth runs and ~30x on low-depth runs) and was
    # concatenated with raw-count nonexpressed genes downstream, mixing two
    # scales in the final table.
    counts_df = pandas.DataFrame(
        {
            'RUN1': [1100.0, 200.0, 800.0, 900.0, 300.0],
            'RUN2': [1050.0, 190.0, 780.0, 880.0, 280.0],
            'RUN3': [13.0, 4.0, 12.0, 11.0, 5.0],
            'RUN4': [12.0, 3.0, 11.0, 10.0, 4.0],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5'],
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
        k_setting='1',
        k_max=3,
        top_n=5,
        min_controls=2,
    )

    assert corrected_df.shape == counts_df.shape
    assert numpy.all(corrected_df.to_numpy(dtype=float) >= 0)
    assert int(summary['resolved_ruv_k']) == 1
    raw_sums = counts_df.sum(axis=1).to_numpy(dtype=float)
    corrected_sums = corrected_df.sum(axis=1).to_numpy(dtype=float)
    # Corrected counts must stay on the raw-count scale: per-gene totals
    # within ~2x of the raw totals (the correction removes batch signal, it
    # must not rescale library depth).
    assert numpy.all(corrected_sums >= raw_sums * 0.5)
    assert numpy.all(corrected_sums <= raw_sums * 2.0)


def test_between_lane_normalize_upper_survives_a_zero_quantile_lane():
    # A sparse lane can have a 75th percentile of zero while the mean quantile
    # is positive. Dividing by that per-lane scale produced inf/NaN, which then
    # propagated into the design residuals and the corrected matrix.
    counts_df = pandas.DataFrame(
        {
            'RUN1': [0.0, 0.0, 0.0, 0.0, 100.0],
            'RUN2': [10.0, 12.0, 11.0, 13.0, 14.0],
            'RUN3': [11.0, 13.0, 10.0, 12.0, 15.0],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5'],
    )

    normalized, scales = _between_lane_normalize_upper(counts_df)

    values = normalized.to_numpy(dtype=float)
    assert numpy.isfinite(values).all()
    # The unusable lane is left unscaled rather than divided by zero.
    assert float(scales.loc['RUN1']) == 1.0
    assert normalized.loc[:, 'RUN1'].tolist() == counts_df.loc[:, 'RUN1'].tolist()
    # The lanes that do carry a usable quantile are still scaled.
    assert float(scales.loc['RUN2']) > 0.0
    assert float(scales.loc['RUN3']) > 0.0


def test_run_ruvseq_backend_survives_a_zero_quantile_lane_end_to_end():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [0.0, 0.0, 0.0, 0.0, 100.0],
            'RUN2': [120.0, 140.0, 110.0, 130.0, 150.0],
            'RUN3': [118.0, 138.0, 112.0, 128.0, 148.0],
            'RUN4': [50.0, 60.0, 45.0, 55.0, 65.0],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP2', 'BP1', 'BP2'],
        }
    )

    corrected_df, _w_df, _summary = run_ruvseq_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        k_setting='1',
        k_max=2,
        top_n=5,
        min_controls=2,
    )

    values = corrected_df.to_numpy(dtype=float)
    assert numpy.isfinite(values).all()
    assert numpy.all(values >= 0)


def test_run_ruvseq_backend_detects_batch_effect_at_high_counts():
    # The residual noise floor must not scale with count magnitude. The control
    # residuals are log-scale, so comparing them against a cutoff derived from
    # the maximum raw count classified a genuine batch effect as numerical
    # noise once library depth was large enough: with counts in the millions
    # the cutoff reached ~2.0 while the real residual was ~1.3, so k collapsed
    # to 0 and the correction was silently skipped.
    rng = numpy.random.default_rng(0)
    values = rng.gamma(shape=4.0, scale=2e6 / 4.0, size=(40, 6))
    # A genuine 2.5x batch effect on the second half of the runs.
    values[:, 3:] *= 2.5
    counts_df = pandas.DataFrame(
        numpy.round(values),
        index=['G{:02d}'.format(i) for i in range(40)],
        columns=['RUN{}'.format(i + 1) for i in range(6)],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': list(counts_df.columns),
            'sample_group': ['A', 'A', 'A', 'B', 'B', 'B'],
            'bioproject': ['BP1', 'BP1', 'BP1', 'BP2', 'BP2', 'BP2'],
        }
    )
    assert counts_df.to_numpy(dtype=float).max() > 1e6

    _corrected_df, w_df, summary = run_ruvseq_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        k_setting='1',
        k_max=2,
        min_controls=5,
    )

    assert int(summary['resolved_ruv_k']) == 1
    assert summary['skip_reason'] == ''
    assert w_df.shape[1] == 1


def test_run_ruvseq_backend_reconstructs_raw_counts_when_correction_is_a_noop():
    # Stronger form of the raw-scale contract than a 0.5x-2x band: with no
    # batch effect present the correction is a no-op, so the returned matrix
    # must reproduce the raw counts exactly, including the unequal per-lane
    # upper-quartile scales that the RUVr step divides out and must multiply
    # back.
    counts_df = pandas.DataFrame(
        {
            'RUN1': [1100.0, 200.0, 800.0, 900.0, 300.0],
            'RUN2': [1100.0, 200.0, 800.0, 900.0, 300.0],
            'RUN3': [13.0, 4.0, 12.0, 11.0, 5.0],
            'RUN4': [13.0, 4.0, 12.0, 11.0, 5.0],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G5'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP2', 'BP1', 'BP2'],
        }
    )
    # The lanes really do sit on very different depths, so a missing
    # multiply-back would be visible rather than masked by equal scales.
    lane_totals = counts_df.sum(axis=0).to_numpy(dtype=float)
    assert lane_totals.max() / lane_totals.min() > 50

    corrected_df, _w_df, summary = run_ruvseq_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        k_setting='1',
        k_max=2,
        top_n=5,
        min_controls=2,
    )

    assert int(summary['resolved_ruv_k']) == 0
    pandas.testing.assert_frame_equal(corrected_df, counts_df)
