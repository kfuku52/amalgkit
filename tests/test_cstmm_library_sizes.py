import numpy
import pandas

from amalgkit.cstmm_python import _get_library_sizes


def test_library_sizes_preserve_selection_labels_and_last_species_values():
    selected = pandas.DataFrame(
        columns=pandas.Index(['NA', '0001', 'missing', 'NA'], name='sample_id'),
    )
    counts = {
        'first': pandas.DataFrame({'0001': [1.25, numpy.nan, 2.5], 'NA': [100., 200., 300.],
                                   'unselected': [900., 900., 900.]}),
        'last': pandas.DataFrame({'NA': [10., numpy.nan, 20.]}),
    }
    observed = _get_library_sizes(selected, counts)
    expected = pandas.Series(
        [30., 3.75, numpy.nan, 30.], index=selected.columns, dtype=float,
    )
    pandas.testing.assert_series_equal(observed, expected, check_exact=True)


def test_library_sizes_preserve_column_sum_rounding():
    rng = numpy.random.default_rng(1082308)
    values = rng.lognormal(3.5, 2.0, (20000, 5))
    values[::37, ::2] = numpy.nan
    counts = pandas.DataFrame(values, columns=['a', 'b', 'c', 'd', 'e'])
    selected = counts.loc[:, ['e', 'a', 'c']]
    # Whole-frame column reduction can use a different floating summation order.
    expected = pandas.Series(
        [float(counts['e'].sum()), float(counts['a'].sum()), float(counts['c'].sum())],
        index=selected.columns,
    )
    observed = _get_library_sizes(selected, {'species': counts})
    pandas.testing.assert_series_equal(observed, expected, check_exact=True)


def test_library_sizes_preserve_partial_multiindex_selection():
    selected = pandas.DataFrame(columns=pandas.MultiIndex.from_tuples(
        [('a', '1'), ('b', '2'), ('a', '3')], names=['species', 'run'],
    ))
    observed = _get_library_sizes(selected, {'species': pandas.DataFrame({'a': [1., 2.], 'b': [4., 5.]})})
    expected = pandas.Series([3., 9., 3.], index=selected.columns)
    pandas.testing.assert_series_equal(observed, expected, check_exact=True)
