import numpy
import pandas
import pytest

from amalgkit.cross_species_computation import (
    finite_correlation_block,
    resolve_correlation_matrix,
    resolve_finite_correlation_matrix,
    resolve_matrix_for_embedding,
    resolve_tsne_perplexity,
    safe_correlation,
)


def test_safe_correlation_returns_nan_for_constant_vectors():
    assert numpy.isnan(safe_correlation([1.0, 1.0, 1.0], [1.0, 2.0, 3.0], 'pearson'))
    assert numpy.isnan(safe_correlation([1.0, 2.0, 3.0], [4.0, 4.0, 4.0], 'spearman'))


def test_safe_correlation_uses_only_finite_pairs():
    observed = safe_correlation([1.0, numpy.nan, 3.0, 4.0], [2.0, 9.0, 6.0, 8.0], 'pearson')
    expected = pandas.Series([1.0, 3.0, 4.0]).corr(pandas.Series([2.0, 6.0, 8.0]))
    assert numpy.isclose(observed, expected)


def test_resolve_matrix_for_embedding_strict_and_cache():
    matrix = pandas.DataFrame([[1.0, 2.0], [numpy.nan, 4.0]], columns=['A', 'B'])
    cache = {}

    strict = resolve_matrix_for_embedding(matrix, 'strict', cache=cache)
    imputed = resolve_matrix_for_embedding(matrix, 'row_mean', cache=cache)

    assert strict.index.tolist() == [0]
    assert numpy.isfinite(imputed.to_numpy()).all()
    assert resolve_matrix_for_embedding(matrix, 'row_mean', cache=cache) is imputed


def test_resolve_correlation_matrix_preserves_undefined_constant_sample():
    matrix = pandas.DataFrame(
        {
            'sample_a': [1.0, 2.0, 3.0, 4.0],
            'sample_b': [2.0, 4.0, 6.0, 8.0],
            'constant_sample': [5.0, 5.0, 5.0, 5.0],
        }
    )

    correlation = resolve_correlation_matrix(matrix)

    assert correlation.index.tolist() == ['sample_a', 'sample_b', 'constant_sample']
    assert correlation.loc['constant_sample'].isna().all()
    assert correlation.loc[:, 'constant_sample'].isna().all()


def test_resolve_finite_correlation_matrix_drops_undefined_constant_sample():
    matrix = pandas.DataFrame(
        {
            'sample_a': [1.0, 2.0, 3.0, 4.0],
            'sample_b': [2.0, 4.0, 6.0, 8.0],
            'constant_sample': [5.0, 5.0, 5.0, 5.0],
        }
    )

    correlation = resolve_finite_correlation_matrix(matrix)

    assert correlation.index.tolist() == ['sample_a', 'sample_b']
    assert numpy.isfinite(correlation.to_numpy()).all()


def test_resolve_tsne_perplexity_is_valid_for_sample_count():
    assert resolve_tsne_perplexity(3) is None
    assert resolve_tsne_perplexity(4) == 1
    assert resolve_tsne_perplexity(10) == 3
    assert resolve_tsne_perplexity(100) == 30


def test_finite_correlation_block_drops_constant_sample():
    corr = pandas.DataFrame(
        {
            'sample_a': [1.0, 1.0, numpy.nan],
            'sample_b': [1.0, 1.0, numpy.nan],
            'constant_sample': [numpy.nan, numpy.nan, numpy.nan],
        },
        index=['sample_a', 'sample_b', 'constant_sample'],
    )

    block = finite_correlation_block(corr)

    assert block.index.tolist() == ['sample_a', 'sample_b']
    assert numpy.isfinite(block.to_numpy()).all()


def test_finite_correlation_block_rejects_partial_undefined_correlations():
    corr = pandas.DataFrame(
        [
            [1.0, 0.4, 0.5],
            [0.4, 1.0, numpy.nan],
            [0.5, numpy.nan, 1.0],
        ],
        index=['sample_a', 'sample_b', 'sample_c'],
        columns=['sample_a', 'sample_b', 'sample_c'],
    )

    with pytest.raises(ValueError, match='Correlations between defined samples must be finite'):
        finite_correlation_block(corr)
