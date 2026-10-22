"""Exact identifier inventories across memory and SQL batch boundaries."""

import os

import pytest

from amalgkit.identifier_validation import TargetIdTracker


@pytest.mark.parametrize('limit', [0, 2, 50_000])
def test_inventory_retains_every_lexical_identifier_across_batches(limit):
    ids = ['0001', '1', 'NA', 'a\x00b', 'é', 'e\u0301', "x'); DROP TABLE ids; --"]
    ids += [f'gene_{i:08d}' for i in range(1601)]
    with TargetIdTracker(memory_limit=limit) as tracker:
        assert tracker.add(ids) == set()
        for start in range(0, len(ids), 137):
            batch = ids[start:start + 137]
            assert tracker.add(batch) == set(batch)
        assert tracker.add(['new']) == set()
        assert tracker.add(['new']) == {'new'}


def test_spilled_duplicate_batch_does_not_insert_any_new_ids():
    with TargetIdTracker(memory_limit=2) as tracker:
        assert tracker.add(['first', 'last', 'NA']) == set()
        values = [f'new_{i}' for i in range(1601)]
        values[501], values[1500] = 'first', 'last'
        assert tracker.add(values) == {'first', 'last'}
        novel = [value for value in values if value not in {'first', 'last'}]
        assert tracker.add(novel) == set()
        assert tracker.add(novel) == set(novel)


@pytest.mark.parametrize('raise_error', [False, True])
def test_spilled_inventory_is_removed_after_success_or_failure(raise_error):
    scratch = None
    try:
        with TargetIdTracker(memory_limit=1) as tracker:
            tracker.add([f'ID{i}' for i in range(1001)])
            assert tracker._db.execute('PRAGMA cache_size').fetchone()[0] == -2048
            assert len(tracker._seen) <= tracker.memory_limit
            scratch = tracker._scratch.name
            assert os.path.isdir(scratch)
            if raise_error:
                raise ValueError('controlled failure')
    except ValueError as exc:
        assert raise_error and str(exc) == 'controlled failure'
    assert scratch is not None
    assert not os.path.exists(scratch)
