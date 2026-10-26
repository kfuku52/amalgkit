from types import SimpleNamespace

import pandas
import pytest

from amalgkit.cstmm_python import save_corrected_output_files_python


@pytest.mark.parametrize('kind', ['full', 'partial', 'absent', 'reordered'])
def test_export_scaling_exact_and_owns_input(tmp_path, kind):
    species = 'Species_A'
    counts = pandas.DataFrame({
        'Species_A_0001': pandas.Series([3, 7, 9], dtype='int64'),
        'Species_A_NA': pandas.Series([0.1, 1.7, 5.2], dtype='float32'),
        'Species_A_unscaled': [4.4, 5.5, 6.6],
    })
    counts.index = pandas.Index(['0001', 'NA', 'g3'], name='original')
    factors = pandas.Series([1.7, 0.3, 9.9], index=['Species_A_0001', 'Species_A_NA', 'Species_B_0001'])
    if kind == 'partial':
        factors = factors.iloc[:1]
    elif kind == 'absent':
        factors = factors.iloc[2:]
    elif kind == 'reordered':
        factors = factors.iloc[::-1]
    counts.attrs = {'provenance': {'models': ['effective']}}
    before = counts.copy(deep=True)
    # Samples with no reference factor retain their original counts.
    expected = counts.astype(float).div(factors.reindex(counts.columns).fillna(1.0), axis='columns')
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
    actual[species].iloc[0, 0] = 0.0
    pandas.testing.assert_frame_equal(counts, before, check_exact=True)
