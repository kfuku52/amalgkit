import json
from types import SimpleNamespace

import pandas
import pytest

from amalgkit.gsa_select import DEFERRED_COLUMN, resolve_deferred_selection


def test_sparse_deferred_rules_preserve_order_and_unrelated_rows():
    threshold = dict(kind='min_spots', threshold=3, target='exclusion', outcome='low_spots')
    dedup = dict(kind='dedup', columns=['biosample'], outcome='redundant')
    recipes = json.dumps([threshold, dedup], sort_keys=True)
    frame = pandas.DataFrame({
        'run': ['0001', 'CRR1', 'CRR2', 'CRR3', 'CRR4', 'NA'],
        'data_source': ['sra', 'gsa', 'gsa', 'gsa', 'gsa', 'sra'],
        'biosample': ['NA', 'B1', 'B1', 'B2', 'B3', '0001'],
        'total_spots': [8, 10, 20, None, 1, 4],
        'read_count_status': ['', 'measured', 'measured', 'unknown', 'measured', ''],
        'exclusion': 'no', 'is_sampled': 'yes',
        DEFERRED_COLUMN: ['', recipes, recipes, recipes, recipes, '[]'],
        'gsa_selection_status': ['', 'pending_input_counts', 'pending_input_counts',
                                 'pending_input_counts', 'pending_input_counts', 'resolved'],
    }, index=[17, 2, 40, 8, 11, 3])
    untouched = frame.loc[[17, 3]].copy()
    result = resolve_deferred_selection(SimpleNamespace(df=frame)).df
    assert result.index.tolist() == [17, 2, 40, 8, 11, 3]
    assert result.exclusion.tolist() == ['no', 'redundant', 'no', 'no', 'low_spots', 'no']
    assert result.gsa_selection_status.tolist() == ['', 'resolved', 'resolved', 'pending_input_counts', 'resolved', 'resolved']
    assert json.loads(result.at[8, DEFERRED_COLUMN]) == [threshold, dedup]
    pandas.testing.assert_frame_equal(result.loc[[17, 3]], untouched)
    resolved = result.copy(deep=True)
    resolve_deferred_selection(SimpleNamespace(df=result))
    pandas.testing.assert_frame_equal(result, resolved)


@pytest.mark.parametrize('value', ['{', '{}', '[1]'])
def test_invalid_deferred_rules_are_rejected_among_empty_rows(value):
    metadata = SimpleNamespace(df=pandas.DataFrame({DEFERRED_COLUMN: ['', '[]', value, '']}))
    with pytest.raises(ValueError):
        resolve_deferred_selection(metadata)


def test_resolved_selection_leaves_empty_recipes_and_statuses_unchanged():
    frame = pandas.DataFrame({DEFERRED_COLUMN: ['', '[]', None],
                              'gsa_selection_status': ['', 'resolved', '']}, index=[8, 2, 9])
    original = frame.copy(deep=True)
    resolve_deferred_selection(SimpleNamespace(df=frame))
    pandas.testing.assert_frame_equal(frame, original)
