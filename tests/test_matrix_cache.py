"""Identity reuse is simulated so these regressions do not depend on allocation luck."""
import gc
import weakref

import numpy
import pandas
import pytest

from amalgkit import cross_species_computation as computation
from amalgkit import cross_species_filter as filtering


@pytest.fixture(params=['filled', 'correlation', 'finite_correlation', 'tsne'])
def resolver(request):
    return {
        'filled': computation.resolve_matrix_for_embedding,
        'correlation': computation.resolve_correlation_matrix,
        'finite_correlation': computation.resolve_finite_correlation_matrix,
        'tsne': filtering._compute_tsne_coordinates,
    }[request.param]


def matrix(offset=0):
    return pandas.DataFrame(
        numpy.random.default_rng(42 + offset).normal(size=(6, 4)),
        columns=['A', 'B', 'C', 'D'],
    )


def test_cache_reuses_only_same_live_source(resolver, monkeypatch):
    monkeypatch.setattr(computation, 'id', lambda _: 7, raising=False)
    monkeypatch.setattr(filtering, 'id', lambda _: 7, raising=False)
    cache = {}
    original = matrix()
    replacement = matrix(1)
    first = resolver(original, 'row_mean', cache=cache)
    assert resolver(original, 'row_mean', cache=cache) is first
    old_refs = [entry[0] for entry in cache.values()]
    expected = resolver(replacement, 'row_mean')
    actual = resolver(replacement, 'row_mean', cache=cache)
    assert actual is not first
    pandas.testing.assert_frame_equal(actual, expected)
    del original
    gc.collect()
    assert all(ref() is None for ref in old_refs)
    # Old-source callbacks must not evict replacement entries sharing the ID.
    assert resolver(replacement, 'row_mean', cache=cache) is actual


def test_cache_evicts_results_when_input_is_collected(resolver):
    cache = {}
    source = matrix()
    source_ref = weakref.ref(source)
    result = resolver(source, 'row_mean', cache=cache)
    result_ref = weakref.ref(result)
    assert cache
    del source, result
    gc.collect()
    assert source_ref() is None
    assert cache == {}
    assert result_ref() is None


def test_dead_source_entry_is_not_reused():
    source = matrix()
    ref = weakref.ref(source)
    del source
    gc.collect()
    replacement = matrix(1)
    key = ('filled', id(replacement), 'row_mean')
    cache = {key: (ref, pandas.DataFrame())}
    assert computation.get_cached_matrix(cache, key, replacement) is None
    assert cache == {}


def test_discarded_cache_does_not_retain_results():
    source = matrix()
    cache = {}
    result = computation.resolve_matrix_for_embedding(source, 'row_mean', cache=cache)
    result_ref = weakref.ref(result)
    del cache, result
    gc.collect()
    assert result_ref() is None


def test_explicit_eviction_before_input_mutation_recomputes():
    source = matrix()
    cache = {}
    before = computation.resolve_correlation_matrix(source, cache=cache)
    filtering._evict_embedding_intermediates(cache, source)
    source['B'] = source['A']
    after = computation.resolve_correlation_matrix(source, cache=cache)
    assert after is not before
    assert after.loc['A', 'B'] == pytest.approx(1)
    pandas.testing.assert_frame_equal(after, source.corr())
