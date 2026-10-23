from types import SimpleNamespace

import numpy
import pandas
import pytest

from amalgkit.cstmm_python import save_corrected_output_files_python


@pytest.mark.parametrize('kind', ['full', 'partial', 'absent', 'reordered', 'repeated'])
@pytest.mark.parametrize('shape', ['mixed', 'empty_rows', 'empty_columns', 'duplicate_columns'])
def test_export_scaling_exact_and_owns_input(tmp_path, kind, shape):
    species = 'Species_A'
    counts = pandas.DataFrame({
        'Species_A_0001': pandas.Series([3, 7, 9], dtype='int64'),
        'Species_A_NA': pandas.Series([0.1, 1.7, 5.2], dtype='float32'),
        'unprefixed': pandas.Series(['4.4', '5.5', '6.6'], dtype=object),
    })
    counts.index = pandas.Index(['0001', 'NA', 'g3'], name='original')
    if shape == 'empty_rows':
        counts = counts.iloc[:0]
    elif shape == 'empty_columns':
        counts = counts.iloc[:, :0]
    elif shape == 'duplicate_columns':
        counts.columns = ['Species_A_0001', 'Species_A_0001', 'unprefixed']
    factors = pandas.Series([1.7, 0.3, 9.9], index=['Species_A_0001', 'Species_A_NA', 'unprefixed'])
    if kind == 'partial':
        factors = factors.iloc[:1]
    elif kind == 'absent':
        factors = pandas.Series([None], index=['missing'])
    elif kind == 'reordered':
        factors = factors.iloc[::-1]
    elif kind == 'repeated':
        factors = pandas.Series([1.7, 0.3, 9.9], index=['Species_A_0001'] * 3)
    counts.attrs = {'provenance': {'models': ['effective']}}
    before = counts.copy(deep=True)
    expected = counts.astype(float).copy()
    for sample_id, factor in factors.items():
        if sample_id in expected.columns:
            expected.loc[:, sample_id] = expected.loc[:, sample_id].to_numpy(dtype=float) / float(factor)
    expected.insert(0, 'target_id', expected.index)
    expected.columns = [col.removeprefix('Species_A_') for col in expected.columns]
    source = tmp_path / 'merge' / species
    source.mkdir(parents=True)
    for name in ['eff_length', 'quant_model']:
        (source / f'{species}_{name}.tsv').write_bytes(b'target_id\t0001\tNA\n0001\t1\t2\nNA\t3\t4\n')
    actual = save_corrected_output_files_python(
        {species: counts}, SimpleNamespace(round2_factors=factors), str(tmp_path / 'merge'), str(tmp_path / 'out'))
    pandas.testing.assert_frame_equal(actual[species], expected.set_index('target_id'), check_exact=True)
    assert (tmp_path / 'out' / species / f'{species}_cstmm_counts.tsv').read_text() == expected.to_csv(sep='\t', index=False)
    pandas.testing.assert_frame_equal(counts, before, check_exact=True)
    for name in ['eff_length', 'quant_model']:
        assert (tmp_path / 'out' / species / f'{species}_{name}.tsv').read_bytes() == (source / f'{species}_{name}.tsv').read_bytes()
    actual[species].attrs['provenance']['models'].append('changed')
    assert counts.attrs == {'provenance': {'models': ['effective']}}
    if actual[species].size:
        actual[species].iloc[0, 0] = numpy.pi
        pandas.testing.assert_frame_equal(counts, before, check_exact=True)


def test_export_scaling_preserves_matched_factor_conversion_error(tmp_path):
    counts = pandas.DataFrame({'Species_A_0001': [1, 2]})
    with pytest.raises(TypeError):
        save_corrected_output_files_python(
            {'Species_A': counts}, SimpleNamespace(round2_factors=pandas.Series([None], index=counts.columns)),
            str(tmp_path), str(tmp_path / 'out'))


def test_export_scaling_preserves_division_warning(tmp_path):
    counts = pandas.DataFrame({'Species_A_0001': [1.0, 2.0]})
    with pytest.warns(RuntimeWarning, match='divide by zero encountered'):
        result = save_corrected_output_files_python(
            {'Species_A': counts}, SimpleNamespace(round2_factors=pandas.Series([0.0], index=counts.columns)),
            str(tmp_path), str(tmp_path / 'out'))
    assert numpy.isinf(result['Species_A'].to_numpy()).all()
