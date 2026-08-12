"""Direct tests for amalgkit.text_utils.normalize_unique_text."""

import pandas

from amalgkit.text_utils import normalize_unique_text


def test_deduplicates_preserving_first_seen_order():
    assert normalize_unique_text(['b', 'a', 'b', 'c', 'a']) == ['b', 'a', 'c']


def test_strips_and_drops_blank_values():
    assert normalize_unique_text(['  x ', '', '  ', 'y']) == ['x', 'y']


def test_returns_empty_list_for_empty_input():
    assert normalize_unique_text([]) == []


def test_drops_none_and_nan():
    assert normalize_unique_text(['a', None, 'b', pandas.NA]) == ['a', 'b']


def test_coerces_non_string_scalars():
    assert normalize_unique_text([1, 2.5, 'x']) == ['1', '2.5', 'x']


def test_nan_marker_dropped_not_stringified():
    # pandas.isna must catch the NaN before str() turns it into 'nan'.
    assert normalize_unique_text(['ok', float('nan'), 'ok']) == ['ok']
