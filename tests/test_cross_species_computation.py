import numpy
import pandas
import pytest

from amalgkit.cross_species_computation import (
    calculate_correlation_within_group,
    finite_correlation_block,
    resolve_correlation_matrix,
    resolve_finite_correlation_matrix,
    resolve_matrix_for_embedding,
    resolve_tsne_perplexity,
    safe_correlation,
)


@pytest.mark.parametrize('reference_exclusion', ['run', 'species'])
def test_cross_group_lookup_preserves_lexical_groups_and_pairwise_support(reference_exclusion):
    matrix = pandas.DataFrame({
        'S_0001': [1., 2., numpy.nan, 4., 8.],
        'S_2': [2., 3., 6., 7., 9.],
        'T_1': [1., 3., 5., numpy.nan, 9.],
        'U_1': [3., 2., 8., 6., 7.],
        'S_3': [9., 7., 6., 4., 1.],
        'T_2': [8., 6., numpy.nan, 3., 2.],
        'U_2': [7., 8., 5., 2., 1.],
    }, index=['0001', 'NA', 'g3', 'g4', 'g5'])
    metadata = pandas.DataFrame({
        'species_tag': ['S', 'T', 'S', 'U', 'S', 'T', 'U', 'Z'],
        'run': ['0001', '1', '2', '1', '3', '2', '2', 'absent'],
        'sample_group': ['NA'] * 4 + ['0001'] * 4,
        'exclusion': 'no',
    }, index=[9, 1, 8, 2, 7, 3, 6, 4])
    original = metadata.copy(deep=True)
    within_columns = ['S_2', 'T_1', 'U_1'] if reference_exclusion == 'run' else ['T_1', 'U_1']
    other_columns = ['S_3', 'T_2', 'U_2'] if reference_exclusion == 'run' else ['T_2', 'U_2']
    within = matrix[within_columns].median(axis=1)
    other = matrix[other_columns].median(axis=1)
    result = calculate_correlation_within_group(metadata, matrix, 'corrected', reference_exclusion, 3)
    assert result.index.tolist() == metadata.index.tolist()
    assert result.at[9, 'within_group_cor_corrected'] == pytest.approx(matrix.S_0001.corr(within))
    assert result.at[9, 'max_nongroup_cor_corrected'] == pytest.approx(matrix.S_0001.corr(other))
    assert result.at[9, 'within_common_genes_corrected'] == 4
    assert result.at[9, 'min_nongroup_common_genes_corrected'] == 4
    assert numpy.isnan(result.at[4, 'within_group_cor_corrected'])
    pandas.testing.assert_frame_equal(metadata, original)


def test_safe_correlation_returns_nan_for_constant_vectors():
    assert numpy.isnan(safe_correlation([1.0, 1.0, 1.0], [1.0, 2.0, 3.0], 'pearson'))
    assert numpy.isnan(safe_correlation([1.0, 2.0, 3.0], [4.0, 4.0, 4.0], 'spearman'))


@pytest.mark.parametrize('method', ['pearson', 'spearman', 'kendall'])
def test_safe_correlation_excludes_nonfinite_pairs(method):
    observed = safe_correlation([1, numpy.inf, 3, 4, numpy.nan], [1, 0, 3, 4, 9], method)
    assert observed == pytest.approx(1.0)


@pytest.mark.parametrize('scale', [1e-100, 1e100, 1e307])
@pytest.mark.parametrize('direction', [-1, 1])
def test_safe_correlation_is_invariant_to_scale(scale, direction):
    left = numpy.array([-3., -1., 1., 3.]) * scale
    right = direction * left
    assert safe_correlation(left, right, 'pearson') == pytest.approx(direction)


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
