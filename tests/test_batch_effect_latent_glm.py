import numpy
import pandas

import amalgkit.batch_effect_latent_glm as latent_glm
from amalgkit.batch_effect_latent_glm import run_latent_glm_backend


def _latent_fixture():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [120, 118, 122, 119, 12, 14, 15, 13, 80, 82, 78, 81],
            'RUN2': [121, 117, 123, 120, 13, 15, 14, 12, 32, 30, 34, 31],
            'RUN3': [14, 13, 12, 11, 131, 129, 133, 130, 79, 83, 77, 80],
            'RUN4': [13, 12, 11, 10, 132, 130, 134, 131, 33, 31, 35, 32],
        },
        index=['G{:02d}'.format(i + 1) for i in range(12)],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP2', 'BP1', 'BP2'],
        }
    )
    return counts_df, metadata_df


def test_run_latent_glm_backend_manual_k_returns_nonnegative_counts_and_latent_scores():
    counts_df, metadata_df = _latent_fixture()
    corrected_df, latent_df, summary = run_latent_glm_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        family='nb',
        k_setting='1',
        k_max=3,
        max_iter=50,
        tol=1e-6,
    )

    assert corrected_df.shape == counts_df.shape
    assert latent_df.shape == (counts_df.shape[1], 1)
    assert numpy.all(corrected_df.to_numpy(dtype=float) >= 0.0)
    assert summary['backend'] == 'latent_glm'
    assert summary['method'] == 'manual'
    assert summary['skip_reason'] == ''
    assert summary['resolved_latent_k'] == 1
    assert summary['latent_family'] == 'nb'
    assert summary['corrected_run_ids'] == ['RUN1', 'RUN2', 'RUN3', 'RUN4']
    assert corrected_df.loc['G09', 'RUN1'] > corrected_df.loc['G09', 'RUN2']


def test_run_latent_glm_backend_auto_k_selects_positive_latent_dimension_for_strong_batch_signal():
    counts_df, metadata_df = _latent_fixture()
    corrected_df, latent_df, summary = run_latent_glm_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        family='poisson',
        k_setting='auto',
        k_max=3,
        max_iter=50,
        tol=1e-6,
    )

    assert corrected_df.shape == counts_df.shape
    assert int(summary['resolved_latent_k']) >= 1
    assert latent_df.shape[1] == int(summary['resolved_latent_k'])
    assert summary['skip_reason'] == ''
    assert numpy.all(corrected_df.to_numpy(dtype=float) >= 0.0)


def test_run_latent_glm_backend_k_zero_returns_noop_summary():
    counts_df, metadata_df = _latent_fixture()
    corrected_df, latent_df, summary = run_latent_glm_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        family='nb',
        k_setting='0',
        k_max=3,
        max_iter=10,
        tol=1e-5,
    )

    pandas.testing.assert_frame_equal(corrected_df, counts_df)
    assert latent_df.shape == (counts_df.shape[1], 0)
    assert summary['skip_reason'] == 'latent_k_zero'
    assert summary['resolved_latent_k'] == 0
    assert summary['corrected_run_ids'] == []


def test_latent_glm_removes_only_latent_effect_and_inverts_pseudocount():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [40.0, 20.0, 30.0, 30.0, 10.0],
            'RUN2': [30.0, 30.0, 40.0, 20.0, 10.0],
            'RUN3': [20.0, 40.0, 30.0, 30.0, 10.0],
            'RUN4': [20.0, 40.0, 30.0, 30.0, 10.0],
            'RUN5': [30.0, 30.0, 20.0, 40.0, 10.0],
            'RUN6': [40.0, 20.0, 30.0, 30.0, 10.0],
        },
        index=['G1', 'G2', 'G3', 'G4', 'G_CONSTANT'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': list(counts_df.columns),
            'sample_group': ['A', 'A', 'A', 'B', 'B', 'B'],
        }
    )

    corrected_df, _latent_df, summary = run_latent_glm_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        family='poisson',
        k_setting='1',
        k_max=1,
        max_iter=20,
        tol=1e-8,
    )

    assert summary['latent_converged'] is True
    numpy.testing.assert_allclose(
        corrected_df.loc['G_CONSTANT', :].to_numpy(dtype=float),
        counts_df.loc['G_CONSTANT', :].to_numpy(dtype=float),
        atol=1e-8,
    )
    assert corrected_df.loc[:, ['RUN1', 'RUN2', 'RUN3']].nunique(axis=1).max() > 1


def test_latent_glm_nonconvergence_returns_original_counts():
    counts_df, metadata_df = _latent_fixture()

    corrected_df, latent_df, summary = run_latent_glm_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        family='nb',
        k_setting='1',
        k_max=1,
        max_iter=0,
        tol=1e-12,
    )

    pandas.testing.assert_frame_equal(corrected_df, counts_df)
    assert latent_df.shape == (counts_df.shape[1], 1)
    assert summary['skip_reason'] == 'latent_not_converged'
    assert summary['latent_converged'] is False
    assert summary['corrected_run_ids'] == []


def _low_count_batch_fixture():
    # Low counts with a strong multiplicative batch effect. Removing the latent
    # effect drives exp(corrected_response) below the 0.5 pseudo-count for some
    # low-expression cells, which is what produces pre-clip negatives.
    rng = numpy.random.default_rng(0)
    values = rng.poisson(lam=1.0, size=(20, 6)).astype(float)
    batch_multiplier = numpy.where(numpy.array([0, 0, 0, 1, 1, 1]) == 1, 3.0, 1.0)
    values *= batch_multiplier.reshape(1, -1)
    counts_df = pandas.DataFrame(
        values,
        index=['G{:02d}'.format(i + 1) for i in range(20)],
        columns=['RUN{}'.format(i + 1) for i in range(6)],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': list(counts_df.columns),
            'sample_group': ['A', 'B'] * 3,
        }
    )
    return counts_df, metadata_df


def test_run_latent_glm_backend_reports_real_pre_clip_negative_counts():
    # Regression pin for the negative-value status counters.
    #
    # The back-transform in _fit_latent_model is exp(corrected_response) - 0.5,
    # which goes negative for low-count genes, and is floored at zero *there*.
    # Measuring negatives on the returned corrected matrix therefore always
    # yields 0 regardless of how much clipping actually happened, so the counter
    # has to be carried out of the fit. This test fails (0 != > 0) if the count
    # is ever re-derived from the post-floor matrix again.
    counts_df, metadata_df = _low_count_batch_fixture()

    corrected_df, _latent_df, summary = run_latent_glm_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        k_setting=1,
        family='nb',
    )

    # The correction must actually have been applied, not skipped.
    assert summary['skip_reason'] == ''
    assert summary['resolved_latent_k'] == 1

    assert summary['negative_values_before_clip'] > 0
    assert summary['negative_values_after_clip'] == 0
    assert (corrected_df.to_numpy(dtype=float) >= 0).all()


def test_run_latent_glm_backend_reports_zero_negatives_when_no_clipping_occurs():
    # Negative control: a counter that is nonzero for every input would carry no
    # information. High, well-separated counts never cross the pseudo-count
    # boundary, so both counters must read zero here.
    counts_df, metadata_df = _latent_fixture()

    corrected_df, _latent_df, summary = run_latent_glm_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        k_setting=1,
        family='nb',
    )

    assert summary['skip_reason'] == ''
    assert summary['negative_values_before_clip'] == 0
    assert summary['negative_values_after_clip'] == 0
    assert (corrected_df.to_numpy(dtype=float) >= 0).all()


def test_latent_updates_use_current_iteration_residuals(monkeypatch):
    # The latent pattern is deliberately mixed with the sample-group design
    # effect, so each latent update must use the residual after fitting the
    # current design-plus-latent model rather than the initial design residual.
    counts_df = pandas.DataFrame(
        {
            'R1': [120, 24, 18, 30, 42, 21, 15, 33],
            'R2': [100, 30, 22, 27, 38, 20, 17, 29],
            'R3': [80, 36, 26, 24, 35, 18, 19, 25],
            'R4': [35, 90, 48, 22, 20, 40, 28, 16],
            'R5': [30, 75, 42, 25, 18, 34, 24, 19],
            'R6': [25, 60, 36, 28, 16, 30, 21, 22],
        },
        index=['G{}'.format(i + 1) for i in range(8)],
        dtype=float,
    )
    metadata_df = pandas.DataFrame(
        {
            'run': list(counts_df.columns),
            'sample_group': ['A', 'A', 'A', 'B', 'B', 'B'],
        }
    )
    calls = []
    original_weighted_residuals = latent_glm._weighted_residuals

    def recording_weighted_residuals(residual_matrix, normalized_counts, fitted_matrix, family):
        calls.append((residual_matrix.copy(), fitted_matrix.copy()))
        return original_weighted_residuals(residual_matrix, normalized_counts, fitted_matrix, family)

    monkeypatch.setattr(latent_glm, '_weighted_residuals', recording_weighted_residuals)
    corrected_df, latent_df, summary = run_latent_glm_backend(
        counts_df=counts_df,
        metadata_df=metadata_df,
        family='poisson',
        k_setting='1',
        k_max=1,
        max_iter=4,
        tol=1e-6,
    )

    assert summary['resolved_latent_k'] == 1
    assert latent_df.shape == (6, 1)
    assert numpy.isfinite(corrected_df.to_numpy(dtype=float)).all()
    assert len(calls) >= 2
    _initial_residuals, _initial_fitted = calls[0]
    for residuals, fitted in calls[1:]:
        expected = latent_glm._prepare_input_arrays(counts_df)[4] - fitted
        numpy.testing.assert_allclose(residuals, expected)
