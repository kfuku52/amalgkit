import numpy as np
import pytest
from amalgkit import merge_plots


def test_explicit_histogram_edges_are_finite_unique_and_sorted():
    np.testing.assert_array_equal(merge_plots._normalize_hist_breaks([], [3, np.nan, 1, 3, np.inf]), [1, 3])
    assert merge_plots._normalize_hist_breaks([np.nan, np.inf]) is None


@pytest.mark.parametrize('mode', ['normal', 'exception', 'degenerate'])
def test_constant_histogram_has_finite_edges_enclosing_data(monkeypatch, mode):
    if mode == 'exception':
        def fail(*args, **kwargs):
            raise ValueError('unavailable edges')
        monkeypatch.setattr(np, 'histogram_bin_edges', fail)
    elif mode == 'degenerate':
        monkeypatch.setattr(np, 'histogram_bin_edges', lambda *a, **kw: np.array([5, 5, np.nan]))
    edges = merge_plots._normalize_hist_breaks([5, 5, np.nan], bin_breaks=[np.nan], bins=3)
    assert np.isfinite(edges).all()
    assert (np.diff(edges) > 0).all()
    assert edges[0] < 5 < edges[-1]


@pytest.mark.parametrize('values,step', [([100], 10), ([300], 25), ([10, 110], 25), ([10, 410], 50), ([10, 810], 100), ([10, 2010], 200)])
def test_insert_axis_ticks_cover_data(values, step):
    ticks = merge_plots._insert_axis_breaks(values + [np.nan])
    assert len(ticks) >= 3
    assert ticks[0] <= min(values) <= max(values) <= ticks[-1]
    assert (np.diff(ticks) == step).all()


def test_insert_axis_ticks_skip_nonfinite_data():
    assert merge_plots._insert_axis_breaks([np.nan, np.inf]) is None
