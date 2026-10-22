"""Contracts for the averaged expression tables consumed by cross-species plots."""

import numpy
import pandas
import pytest

from amalgkit.cross_species_filter import _build_averaged_cross_species_inputs


@pytest.mark.parametrize('duplicate_index', [False, True])
@pytest.mark.parametrize('with_project', [False, True])
def test_averages_keep_lexical_groups_weighted_runs_and_inputs(duplicate_index, with_project):
    metadata = pandas.DataFrame({
        'species_tag': ['A', 'A', 'A', 'B', 'B', 'B'],
        'scientific_name': ['Species A'] * 3 + ['Species B'] * 3,
        'run': ['a', 'a', 'b', 'c', 'd', 'missing'],
        'sample_group': ['0001'] * 3 + ['NA', 'NA', 'NA'],
        'exclusion': ['no'] * 5 + ['removed'],
    }, index=[3] * 6 if duplicate_index else [8, 5, 4, 7, 2, 1])
    if with_project:
        metadata['bioproject'] = ['P1', 'P1', 'P2', 'P3', 'P3', 'P4']
    matrix = pandas.DataFrame({
        'A_a': ['3', 'bad', '6'], 'A_b': ['12', '8', '0'],
        'B_c': ['4', '2', 'bad'], 'B_d': ['8', 'bad', 'bad'],
    }, index=['OG1', 'NA', '0001'])
    original_metadata, original_matrix = metadata.copy(deep=True), matrix.copy(deep=True)
    observed = _build_averaged_cross_species_inputs(metadata, {'uncorrected': matrix, 'corrected': matrix})
    expected = pandas.DataFrame({'A_0001': [6., 8., 4.], 'B_NA': [6., 2., numpy.nan]}, index=matrix.index)
    for correction in ['uncorrected', 'corrected']:
        pandas.testing.assert_frame_equal(observed[correction], expected, check_exact=True)
    labels = observed['labels']
    assert labels['averaged_id'].tolist() == ['A_0001', 'B_NA']
    assert labels['num_run'].tolist() == [2, 2]
    assert labels['num_bp'].tolist() == ([2, 1] if with_project else [0, 0])
    observed['corrected'].iloc[0, 0] = -99.
    pandas.testing.assert_frame_equal(metadata, original_metadata)
    pandas.testing.assert_frame_equal(matrix, original_matrix)
    pandas.testing.assert_frame_equal(observed['uncorrected'], expected, check_exact=True)


@pytest.mark.parametrize('dtype', ['int64', 'float64'])
def test_singleton_averages_preserve_dtype_and_do_not_alias_source(dtype):
    metadata = pandas.DataFrame(dict(species_tag=['A'], scientific_name=['Species A'],
                                    run=['a'], sample_group=['leaf'], exclusion=['no']))
    matrix = pandas.DataFrame({'A_a': [2, 4]}, dtype=dtype)
    original = matrix.copy(deep=True)
    observed = _build_averaged_cross_species_inputs(metadata, dict(uncorrected=matrix, corrected=matrix))
    pandas.testing.assert_frame_equal(observed['corrected'], matrix.rename(columns={'A_a': 'A_leaf'}))
    observed['corrected'].iloc[0, 0] = -1
    pandas.testing.assert_frame_equal(matrix, original)


def test_averages_keep_shared_tag_groups_and_repeated_labels():
    # Distinct scientific names sharing one tag historically share the same
    # expression reference and retain both label rows in the exported plot data.
    metadata = pandas.DataFrame(dict(species_tag=['A', 'A'], scientific_name=['Alpha', 'Beta'],
                                    run=['a', 'b'], sample_group=['leaf', 'leaf'], exclusion=['no', 'no']))
    matrix = pandas.DataFrame({'A_a': [2.], 'A_b': [6.]})
    observed = _build_averaged_cross_species_inputs(metadata, dict(uncorrected=matrix, corrected=matrix))
    assert observed['labels']['scientific_name'].tolist() == ['Alpha', 'Beta']
    assert observed['uncorrected'].columns.tolist() == ['A_leaf', 'A_leaf', 'A_leaf', 'A_leaf']
    numpy.testing.assert_array_equal(observed['uncorrected'].to_numpy(), [[4., 4., 4., 4.]])


def test_averages_handle_many_groups_without_fragmentation():
    rows = [dict(species_tag=f'S{i}', scientific_name=f'Species {i}', run=str(g),
                 sample_group=f'T{g}', exclusion='no') for i in range(20) for g in range(6)]
    metadata = pandas.DataFrame(rows)
    matrix = pandas.DataFrame({f'S{i}_{g}': [i + g] for i in range(20) for g in range(6)})
    observed = _build_averaged_cross_species_inputs(metadata, dict(uncorrected=matrix, corrected=matrix))
    assert len(observed['labels']) == 120
    for correction in ['uncorrected', 'corrected']:
        assert observed[correction].shape == (1, 120)
        for i, g in [(0, 0), (7, 4), (19, 5)]:
            assert observed[correction].loc[0, f'S{i}_T{g}'] == i + g
