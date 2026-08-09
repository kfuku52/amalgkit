import warnings

import numpy
import pandas
import pytest

from amalgkit.imputation import _truncated_svd_reconstruction, impute_expression


def test_truncated_svd_reconstruction_matches_full_svd_at_degenerate_boundary():
    values = numpy.zeros((8, 4), dtype=float)
    values[:4, :4] = numpy.eye(4)
    left, singular_values, right = numpy.linalg.svd(values, full_matrices=False)
    expected = (left[:, :1] * singular_values[:1]).dot(right[:1, :])

    observed = _truncated_svd_reconstruction(values, num_pc=1)

    numpy.testing.assert_allclose(observed, expected, atol=0.0, rtol=0.0)


def test_impute_expression_row_mean_fills_missing_values():
    df = pandas.DataFrame(
        {
            'RUN1': [1.0, 2.0, 3.0],
            'RUN2': [4.0, numpy.nan, 6.0],
            'RUN3': [7.0, 8.0, numpy.nan],
        },
        index=['G1', 'G2', 'G3'],
    )
    imputed = impute_expression(df, strategy='row_mean')

    assert numpy.isfinite(imputed.to_numpy(dtype=float)).all()
    assert imputed.loc['G2', 'RUN2'] == pytest.approx(5.0)
    # G3 row mean: (3 + 6) / 2 = 4.5
    assert imputed.loc['G3', 'RUN3'] == pytest.approx(4.5)
    # Observed values are never altered.
    assert imputed.loc['G1', 'RUN1'] == 1.0


def test_impute_expression_preserves_observed_values():
    df = pandas.DataFrame(
        {
            'RUN1': [1.0, 2.0, 3.0],
            'RUN2': [4.0, numpy.nan, 6.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    imputed = impute_expression(df, strategy='em_pca')

    assert imputed.loc['G1', 'RUN1'] == 1.0
    assert imputed.loc['G1', 'RUN2'] == 4.0
    assert imputed.loc['G3', 'RUN1'] == 3.0
    assert imputed.loc['G3', 'RUN2'] == 6.0


def test_impute_expression_em_pca_fallback_warns_on_degenerate_input():
    # A zero-variance matrix cannot resolve a NIPALS component; the fallback
    # to row-mean imputation must be surfaced with a warning, not silently
    # presented as EM-PCA output.
    df = pandas.DataFrame(
        {
            'RUN1': [5.0, 5.0, 5.0],
            'RUN2': [5.0, numpy.nan, 5.0],
            'RUN3': [5.0, 5.0, numpy.nan],
        },
        index=['G1', 'G2', 'G3'],
    )

    with pytest.warns(UserWarning, match='row-mean'):
        imputed = impute_expression(df, strategy='nipals', num_pc=1)

    assert numpy.isfinite(imputed.to_numpy(dtype=float)).all()
    assert imputed.loc['G2', 'RUN2'] == pytest.approx(5.0)
    assert imputed.loc['G3', 'RUN3'] == pytest.approx(5.0)


def test_impute_expression_single_row_warns_on_row_mean_fallback():
    # max_pc = min(nrow - 1, ncol - 1) = 0 for a one-row matrix, so the
    # EM-PCA path cannot run at all. That degradation must be surfaced rather
    # than returning row-mean values labelled as EM-PCA output.
    df = pandas.DataFrame(
        {
            'RUN1': [1.0],
            'RUN2': [numpy.nan],
            'RUN3': [3.0],
        },
        index=['G1'],
    )

    with pytest.warns(UserWarning, match='row-mean'):
        imputed = impute_expression(df, strategy='em_pca')

    assert numpy.isfinite(imputed.to_numpy(dtype=float)).all()
    # Row mean of the observed values: (1 + 3) / 2 = 2.0
    assert imputed.loc['G1', 'RUN2'] == pytest.approx(2.0)
    assert imputed.loc['G1', 'RUN1'] == 1.0
    assert imputed.loc['G1', 'RUN3'] == 3.0


def test_impute_expression_single_column_warns_on_row_mean_fallback():
    # The mirror case: one column gives max_pc = 0 as well.
    df = pandas.DataFrame(
        {'RUN1': [1.0, numpy.nan, 3.0]},
        index=['G1', 'G2', 'G3'],
    )

    with pytest.warns(UserWarning, match='row-mean'):
        imputed = impute_expression(df, strategy='nipals')

    assert numpy.isfinite(imputed.to_numpy(dtype=float)).all()
    # G2 has no observed value in any column, so its row mean is 0.
    assert imputed.loc['G2', 'RUN1'] == pytest.approx(0.0)
    assert imputed.loc['G1', 'RUN1'] == 1.0


def test_impute_expression_does_not_warn_when_em_pca_succeeds():
    # Negative control: the fallback warning must not fire on the happy path,
    # otherwise the warning carries no information.
    df = pandas.DataFrame(
        {
            'RUN1': [1.0, 2.0, 3.0, 4.0],
            'RUN2': [2.0, 4.1, 6.0, 8.2],
            'RUN3': [3.0, numpy.nan, 9.1, 12.0],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )

    with warnings.catch_warnings():
        warnings.simplefilter('error')
        imputed = impute_expression(df, strategy='em_pca', num_pc=1)

    assert numpy.isfinite(imputed.to_numpy(dtype=float)).all()


def test_impute_expression_minimum_imputed_value_floor():
    df = pandas.DataFrame(
        {
            'RUN1': [1.0, 2.0, 3.0],
            'RUN2': [4.0, numpy.nan, 6.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    imputed = impute_expression(df, strategy='row_mean', minimum_imputed_value=0.5)

    assert imputed.loc['G2', 'RUN2'] >= 0.5
    assert imputed.loc['G1', 'RUN1'] == 1.0
