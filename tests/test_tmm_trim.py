"""TMM partial selection must preserve exact count-key ranks and tie boundaries."""

import runpy
import weakref
from decimal import localcontext
from fractions import Fraction
from pathlib import Path

import numpy
import pandas
import pytest

from amalgkit import normalization_tmm as tmm


@pytest.mark.parametrize('size', [31, 1023, 1024, 1025, 4096])
@pytest.mark.parametrize('kind', ['integer_ties', 'fractional', 'extreme', 'all_tied'])
@pytest.mark.parametrize('ratio', [False, True])
@pytest.mark.parametrize('variable_library', [False, True])
def test_trim_mask_matches_full_ranks(size, kind, ratio, variable_library):
    rng = numpy.random.default_rng(31273 + size)
    if kind == 'integer_ties':
        obs, ref = rng.integers(1, 30, size=(2, size)).astype(float)
    elif kind == 'fractional':
        obs, ref = numpy.exp(rng.normal(size=(2, size)))
    elif kind == 'extreme':
        obs = numpy.ldexp(rng.uniform(.5, 1., size), rng.integers(-1070, 1024, size))
        ref = numpy.ldexp(rng.uniform(.5, 1., size), rng.integers(-1070, 1024, size))
    else:
        obs, ref = numpy.full(size, 4.), numpy.full(size, 3.)
    libraries = numpy.exp(rng.uniform(1., 10., size)) if variable_library else numpy.full(size, 1e8)
    if kind == 'extreme' and variable_library:
        libraries = numpy.ldexp(rng.uniform(.5, 1., size), rng.integers(-1070, 1024, size))
    before = [array.copy() for array in (obs, ref, libraries)]
    ranks = tmm._rank_count_expression(obs, ref, libraries, ratio=ratio)

    for trim in [-.1, 0., .05, .3, .499, .5, .6]:
        low = int(numpy.floor(size * trim) + 1)
        high = size + 1 - low
        actual = tmm._trim_count_expression(obs, ref, libraries, ratio=ratio, low=low, high=high)
        numpy.testing.assert_array_equal(actual, (ranks >= low) & (ranks <= high))
    for array, original in zip((obs, ref, libraries), before):
        numpy.testing.assert_array_equal(array, original)


def _fraction_average_ranks(values):
    # Independent exact rational ordering; no production key or rank helper.
    occurrences = {}
    for rank, index in enumerate(sorted(range(len(values)), key=values.__getitem__), 1):
        occurrences.setdefault(values[index], []).append(rank)
    means = {value: Fraction(sum(ranks), len(ranks)) for value, ranks in occurrences.items()}
    return [means[value] for value in values]


@pytest.mark.parametrize('ratio', [False, True])
@pytest.mark.parametrize('variable_library', [False, True])
def test_partition_boundary_ties_match_independent_fraction_ranks(ratio, variable_library, monkeypatch):
    obs = numpy.repeat([1., 2., 3., 4.], 256)
    ref = numpy.ones(1024)
    libraries = numpy.tile([8., 16.], 512) if variable_library else numpy.full(1024, 8.)
    keys = []
    for observed, reference, library in zip(obs, ref, libraries):
        observed, reference, library = map(Fraction, (float(observed), float(reference), float(library)))
        key = observed / reference if ratio else observed * reference
        if variable_library:
            key = key * library if ratio else key / library
        keys.append(key)
    ranks = _fraction_average_ranks(keys)

    def reject_full_rank(*args, **kwargs):
        raise AssertionError('The partition path must be exercised directly.')

    monkeypatch.setattr(tmm, '_rank_count_expression', reject_full_rank)
    # Exercise whole-tie inclusion/exclusion on either side of half ranks,
    # and independently disable the lower or upper trim.
    for low, high in [(384, 641), (385, 641), (384, 640), (1, 640), (385, 1024), (1, 1024)]:
        actual = tmm._trim_count_expression(obs, ref, libraries, ratio=ratio, low=low, high=high)
        numpy.testing.assert_array_equal(actual, [low <= rank <= high for rank in ranks])


def test_partition_empty_bounds_need_no_count_keys():
    empty = numpy.array([], dtype=float)
    actual = tmm._trim_count_expression(empty, empty, empty, ratio=True, low=1, high=0)
    numpy.testing.assert_array_equal(actual, numpy.array([], dtype=bool))


@pytest.mark.parametrize('ratio', [False, True])
def test_extreme_variable_libraries_match_fraction_oracle(ratio, monkeypatch):
    rng = numpy.random.default_rng(83147)
    exponents = [-1074, -1022, -100, -1, 0, 100, 1000, 1023]
    obs, ref, libraries = [numpy.ldexp(numpy.ones(1024), rng.choice(exponents, 1024)) for _ in range(3)]
    keys = []
    for observed, reference, library in zip(obs, ref, libraries):
        observed, reference, library = map(Fraction, (float(observed), float(reference), float(library)))
        keys.append(observed / reference * library if ratio else observed * reference / library)
    ranks = _fraction_average_ranks(keys)

    def reject_full_rank(*args, **kwargs):
        raise AssertionError('Extreme exponents must also use partial selection.')

    monkeypatch.setattr(tmm, '_rank_count_expression', reject_full_rank)
    for low, high in [(52, 973), (308, 717)]:
        actual = tmm._trim_count_expression(obs, ref, libraries, ratio=ratio, low=low, high=high)
        numpy.testing.assert_array_equal(actual, [low <= rank <= high for rank in ranks])


@pytest.fixture(scope='module')
def independent_pair_oracle():
    path = Path(__file__).parent / 'reference' / 'generate_tmm_fixtures.py'
    return runpy.run_path(str(path))['pair_factor']


@pytest.mark.parametrize('kind', ['integer_ties', 'fractional'])
def test_large_pair_matches_independent_decimal_oracle(kind, independent_pair_oracle, monkeypatch):
    rng = numpy.random.default_rng(20261023)
    if kind == 'integer_ties':
        obs, ref = rng.integers(1, 40, size=(2, 2048)).astype(float)
    else:
        obs, ref = numpy.exp(rng.normal(size=(2, 2048)))
    libraries = [float(obs.sum() * 3), float(ref.sum() * 2)]
    with localcontext() as context:
        context.prec = 60
        expected = float(independent_pair_oracle(obs, ref, *libraries))

    def reject_full_rank(*args, **kwargs):
        raise AssertionError('The independent oracle must cover partial selection.')

    monkeypatch.setattr(tmm, '_rank_count_expression', reject_full_rank)
    observed = tmm.calc_factor_tmm(obs, ref, *libraries)
    numpy.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize('kind', ['dense', 'sparse', 'recycled_vector', 'matrix_reference'])
@pytest.mark.parametrize('weighted', [False, True])
def test_factor_and_diagnostics_preserve_full_rank_behavior(kind, weighted, monkeypatch):
    rng = numpy.random.default_rng(41273)
    obs, ref = numpy.exp(rng.normal(3., 2., size=(2, 4096)))
    libraries = [float(obs.sum() * 3), float(ref.sum() * 2)]
    if kind == 'sparse':
        obs[::3], ref[::5] = 0., 0.
    elif kind == 'recycled_vector':
        ref = ref[:3]
        libraries[1] = numpy.array([200_000., 400_000., 300_000.])
    elif kind == 'matrix_reference':
        ref = ref.reshape((2048, 2))
        libraries[1] = numpy.array([200_000., 400_000.])
    observed_diagnostics = {}
    observed = tmm.calc_factor_tmm(obs, ref, *libraries, do_weighting=weighted, diagnostics=observed_diagnostics)

    def full_rank_trim(observed, reference, reference_libraries, *, ratio, low, high):
        ranks = tmm._rank_count_expression(observed, reference, reference_libraries, ratio=ratio)
        return (ranks >= low) & (ranks <= high)

    monkeypatch.setattr(tmm, '_trim_count_expression', full_rank_trim)
    expected_diagnostics = {}
    expected = tmm.calc_factor_tmm(obs, ref, *libraries, do_weighting=weighted, diagnostics=expected_diagnostics)
    assert observed == expected
    assert observed_diagnostics == expected_diagnostics


@pytest.mark.parametrize('value', [0., -0., .2, numpy.nan, numpy.inf, -numpy.inf, 5e-324, 1e308])
@pytest.mark.parametrize('size', [0, 1, 7, 1024])
def test_scalar_recycling_retains_float_bits_and_writable_storage(value, size):
    source = numpy.array([value], dtype=float)
    expected = numpy.resize(source, size)
    actual = tmm._recycle(source, size)
    assert actual.tobytes() == expected.tobytes()
    assert actual.dtype == expected.dtype
    assert actual.shape == (size,)
    assert actual.flags.writeable
    if size > 1:
        actual[0] = 42.
        assert source.tobytes() == numpy.array([value], dtype=float).tobytes()


@pytest.mark.parametrize('values', [[], [3., 7.], [[1., 2.], [3., 4.]]])
@pytest.mark.parametrize('size', [0, 1, 9])
def test_nonscalar_recycling_preserves_resize_and_flattening(values, size):
    source = numpy.asarray(values, dtype=float).reshape(-1)
    expected = source if source.size == 0 else numpy.resize(source, size)
    actual = tmm._recycle(values, size)
    numpy.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('single_column', [False, True])
def test_rounds_resolve_reference_once_and_release_filtered_copy(single_column, monkeypatch):
    columns = ['0001'] if single_column else [1, '1', '0001', 'NA']
    values = numpy.array([[0., 0., 0., 0.], [10., 20., 30., 40.], [30., 10., 40., 20.]])
    counts = pandas.DataFrame(values[:, :len(columns)], columns=columns)
    original_counts = counts.copy(deep=True)
    libraries = counts.sum().iloc[::-1] * 3
    references = []
    filtered_references = []
    original_resolve = tmm._resolve_tmm_reference_column
    original_filter = tmm._remove_all_zero_rows
    original_normalize = tmm.calc_norm_factors_tmm

    def record_reference(matrix, lib_sizes):
        reference = original_resolve(matrix, lib_sizes)
        references.append(reference)
        return reference

    def record_filter(matrix):
        filtered = original_filter(matrix)
        if filtered is not matrix:
            filtered_references.append(weakref.ref(filtered))
        return filtered

    def check_released_copy(*args, **kwargs):
        assert all(reference() is None for reference in filtered_references)
        return original_normalize(*args, **kwargs)

    monkeypatch.setattr(tmm, '_resolve_tmm_reference_column', record_reference)
    monkeypatch.setattr(tmm, '_remove_all_zero_rows', record_filter)
    monkeypatch.setattr(tmm, 'calc_norm_factors_tmm', check_released_copy)
    result = tmm.run_tmm_rounds_for_cstmm(counts, lib_size=libraries)

    assert references == [result.round1_reference_column]
    assert filtered_references and all(reference() is None for reference in filtered_references)
    assert result.round1_factors.index.tolist() == columns
    assert result.round2_factors.index.tolist() == columns
    pandas.testing.assert_series_equal(result.library_sizes, libraries.reindex(columns))
    pandas.testing.assert_frame_equal(counts, original_counts)
    if single_column:
        assert result.round1_reference_column == 0
        assert result.median_reference_columns == [0]
        numpy.testing.assert_array_equal(result.round1_factors, [1.])
        numpy.testing.assert_array_equal(result.round2_factors, [1.])
