import numpy
import pandas

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
