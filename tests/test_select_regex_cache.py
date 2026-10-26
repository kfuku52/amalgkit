"""Selection regex compression must preserve lexical text and rule behavior."""

import re
from types import SimpleNamespace

import numpy
import pandas
import pytest

from amalgkit import select


def _rule(pattern, columns=('sample_title',), **kwargs):
    return {
        'pattern': pattern,
        'regex': re.compile(pattern, re.IGNORECASE),
        'columns': list(columns),
        **kwargs,
    }


def _pandas_match(frame, rule):
    expected = pandas.Series(False, index=frame.index)
    for column in rule['columns']:
        if column in frame:
            text = frame[column].fillna('').astype(str)
            expected |= text.str.contains(rule['regex'], regex=True).fillna(False)
    return expected


@pytest.mark.parametrize('rows', [0, 32, 4095, 4096, 8192])
@pytest.mark.parametrize('pattern', [r'^$', r'\b(?:leaf|single[\s_-]cell)\b', r'^NA$|^0001$',
                                    r'葉|é|e\u0301', r'(?s)single.cell', r'(?<=a)[ _-]b'])
def test_selection_regex_matches_pandas_for_lexical_unicode_and_duplicate_indexes(rows, pattern):
    values = [None, numpy.nan, '', 'NA', '0001', 'Leaf', '葉', 'é', 'e\u0301',
              'single\ncell', 'a_b', 'a b', 'a-b', 12, True]
    frame = pandas.DataFrame({
        'sample_title': [values[i % len(values)] for i in range(rows)],
        'other': [values[(i + 3) % len(values)] for i in range(rows)],
    }, index=pandas.Index([i % 7 for i in range(rows)], name='row'))
    original = frame.copy(deep=True)
    rule = _rule(pattern, ('sample_title', 'absent', 'other'))
    for cache in (None, {}):
        actual = select.build_select_rule_match_mask_with_cache(frame, rule, cache)
        pandas.testing.assert_series_equal(actual, _pandas_match(frame, rule), check_exact=True)
    pandas.testing.assert_frame_equal(frame, original, check_exact=True)


def test_high_cardinality_text_matches_pandas():
    frame = pandas.DataFrame({'sample_title': [f'leaf accession_{i}' for i in range(8192)]})
    rule = _rule(r'accession_0$|accession_8191$')
    actual = select.build_select_rule_match_mask_with_cache(frame, rule, {})
    pandas.testing.assert_series_equal(actual, _pandas_match(frame, rule), check_exact=True)


@pytest.mark.parametrize('unique_count', [256, 257])
def test_compact_codes_preserve_identifiers_across_integer_width_boundary(unique_count):
    values = ['0000'] * 8192
    positions = [i for i in range(256, len(values)) if i % 32]
    for number, position in enumerate(positions):
        values[position] = f'{number % unique_count:04d}'
    frame = pandas.DataFrame({'sample_title': values})
    rule = _rule(f'^{unique_count - 1:04d}$')
    cache = {}
    actual = select.build_select_rule_match_mask_with_cache(frame, rule, cache)
    pandas.testing.assert_series_equal(actual, _pandas_match(frame, rule), check_exact=True)
    assert actual.any()


def test_new_invocation_after_text_mutation_uses_current_values():
    frame = pandas.DataFrame({'sample_title': ['leaf'] * 4096})
    cache = {}
    rule = _rule(r'leaf')
    assert select.build_select_rule_match_mask_with_cache(frame, rule, cache).all()
    # Callers clear the invocation-local snapshot before editing source text.
    frame.loc[:, 'sample_title'] = 'root'
    cache.clear()
    assert not select.build_select_rule_match_mask_with_cache(frame, rule, cache).any()
    frame.loc[:, 'sample_title'] = 'leaf'
    assert select.build_select_rule_match_mask_with_cache(frame, rule, {}).all()


def test_factor_cache_is_independent_of_user_column_names_and_order():
    frame = pandas.DataFrame({
        'sample_title': ['leaf', 'root'] * 4096,
        '_SELECT_REGEX_FACTOR_CACHE_KEY': ['0001', 'NA'] * 4096,
    })
    cache = {}
    rules = [_rule(r'leaf', ('sample_title',)),
             _rule(r'^NA$', ('_SELECT_REGEX_FACTOR_CACHE_KEY',)),
             _rule(r'^0001$|root', ('_SELECT_REGEX_FACTOR_CACHE_KEY', 'sample_title'))]
    for rule in rules:
        actual = select.build_select_rule_match_mask_with_cache(frame, rule, cache)
        pandas.testing.assert_series_equal(actual, _pandas_match(frame, rule), check_exact=True)


def test_exclude_rule_priority_stop_and_existing_evidence_are_preserved(monkeypatch):
    frame = pandas.DataFrame({
        'sample_title': ['leaf single cell', 'leaf', 'root', 'NA'] * 2048,
        'exclusion': ['no', 'no', 'existing_evidence', 'no'] * 2048,
        'run': [f'{i:06}' for i in range(8192)],
    }, index=pandas.Index(numpy.arange(8192) * 2 + 1, name='row'))
    rules = [
        _rule(r'leaf', stage='exclude', target_column='exclusion', outcome='first', stop_on_match=True),
        _rule(r'single cell|root', stage='exclude', target_column='exclusion', outcome='second', stop_on_match=True),
        _rule(r'^NA$', stage='exclude', target_column='exclusion', outcome='temporary', stop_on_match=False),
        _rule(r'^NA$', stage='exclude', target_column='exclusion', outcome='last', stop_on_match=True),
    ]
    actual = select.apply_select_exclude_rules(SimpleNamespace(df=frame.copy()), rules).df
    monkeypatch.setattr(select, 'build_select_rule_match_mask_with_cache',
                        lambda df, rule, text_cache=None: _pandas_match(df, rule))
    expected = select.apply_select_exclude_rules(SimpleNamespace(df=frame.copy()), rules).df
    pandas.testing.assert_frame_equal(actual, expected, check_exact=True)
    assert actual.exclusion.tolist() == ['first', 'first', 'existing_evidence', 'last'] * 2048
    assert actual['run'].iloc[0] == '000000'
