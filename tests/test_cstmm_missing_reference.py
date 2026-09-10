import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from amalgkit.cstmm import cstmm_main
from amalgkit.cstmm_diagnostics import write_observed_pair_diagnostics, plot_observed_pair_diagnostics
from amalgkit.cstmm_python import _get_df_nonzero, _get_df_exp_single_copy_ortholog, _read_est_counts
from amalgkit.imputation import impute_expression
from amalgkit.normalization_tmm import run_tmm_rounds_for_cstmm


def test_species_block_imputation_respects_tenfold_depth_and_preserves_observations():
    genes = np.arange(1, 21, dtype=float) * 100
    complete = pd.DataFrame(dict(A1=genes, A2=genes, B1=genes * 10, B2=genes * 10))
    masked = complete.copy()
    masked.iloc[4:12, 2:] = np.nan
    libraries = complete.sum()
    raw = _get_df_nonzero(masked, libraries, scale='raw')
    scaled, diagnostics = _get_df_nonzero(masked, libraries, scale='library_size', return_diagnostics=True)
    np.testing.assert_array_equal(scaled.to_numpy()[masked.notna()], complete.to_numpy()[masked.notna()])
    np.testing.assert_allclose(scaled, complete, rtol=1e-12)
    np.testing.assert_allclose(run_tmm_rounds_for_cstmm(scaled, libraries).round2_factors, 1, rtol=1e-12)
    np.testing.assert_allclose(run_tmm_rounds_for_cstmm(raw, libraries).round2_factors,
                               [np.sqrt(10), np.sqrt(10), 1 / np.sqrt(10), 1 / np.sqrt(10)])
    assert diagnostics['converged'] and diagnostics['missing_cells'] == 16


def test_hundred_species_with_no_universally_observed_ortholog():
    complete = pd.DataFrame(np.arange(1, 201)[:, None] * np.geomspace(1, 100, 100)[None, :],
                            columns=[f'Species_{i}_R1' for i in range(100)])
    masked = complete.copy()
    # Half of the species are missing each orthogroup; every pair still overlaps.
    missing = (np.arange(200)[:, None] + np.arange(100)[None, :]) % 101 < 50
    masked[missing] = np.nan
    assert not masked.notna().all(axis=1).any()
    imputed, diagnostics = _get_df_nonzero(masked, complete.sum(), scale='library_size', return_diagnostics=True)
    assert diagnostics['complete_rows'] == 0
    np.testing.assert_allclose(imputed, complete, rtol=1e-10)
    np.testing.assert_allclose(run_tmm_rounds_for_cstmm(imputed, complete.sum()).round2_factors, 1, rtol=1e-10)


def test_observed_pair_comparisons_keep_unestimable_pairs_missing(tmp_path):
    counts = pd.DataFrame({'a': [10., 20, np.nan], 'b': [100., 200, np.nan],
                           'c': [np.nan, np.nan, 30.]})
    factors = pd.Series({'a': 2., 'b': .5, 'c': 1.})
    original = factors.copy()
    output = tmp_path / 'pairs.tsv'
    write_observed_pair_diagnostics(counts, pd.Series({'a': 100., 'b': 1000., 'c': 100.}), factors, output)
    pairs = pd.read_csv(output, sep='\t')
    assert set(pairs.purpose) == {'reference_only'}
    assert pairs.loc[0, 'observed_factor_ratio'] == pytest.approx(1)
    assert pairs.loc[0, 'applied_factor_ratio'] == pytest.approx(.25)
    assert pairs.loc[0, 'log2_ratio_difference'] == pytest.approx(2)
    assert pairs.loc[1:, 'observed_factor_ratio'].isna().all()
    pd.testing.assert_series_equal(factors, original)


def test_imputation_convergence_diagnostics_do_not_change_shared_api():
    matrix = pd.DataFrame([[1., 2., np.nan], [2., np.nan, 7.], [3., 8., 2.], [5., 3., 4.]])
    legacy = impute_expression(matrix, num_pc=1, max_iter=1)
    result, diagnostic = impute_expression(matrix, num_pc=1, max_iter=1, return_diagnostics=True)
    pd.testing.assert_frame_equal(result, legacy)
    assert not diagnostic['converged'] and diagnostic['iterations'] == 1
    assert diagnostic['final_delta'] > diagnostic['tolerance']


def test_orthology_missingness_distinguishes_copy_number_and_identifier(tmp_path):
    gc = tmp_path / 'gc.tsv'
    og = tmp_path / 'og.tsv'
    pd.DataFrame({'orthogroup_id': ['og0', 'og1', 'og2', 'og3'],
                  'A': [1, 1, 1, 1], 'B': [0, 2, 1, 1]}).to_csv(gc, sep='\t', index=False)
    pd.DataFrame({'orthogroup_id': ['og0', 'og1', 'og2', 'og3'],
                  'A': ['a0', 'a1', 'a2', 'a3'], 'B': ['', 'b1,b2', 'missing', 'b3']}).to_csv(og, sep='\t', index=False)
    counts = {'A': pd.DataFrame({'A_run': [1., 2, 3, 4]}, index=['a0', 'a1', 'a2', 'a3']),
              'B': pd.DataFrame({'B_run': [1., 2, 0]}, index=['b1', 'b2', 'b3'])}
    result = _get_df_exp_single_copy_ortholog(gc, og, tmp_path, counts)
    audit = pd.DataFrame(result.attrs['orthology_audit'])
    assert audit.loc[audit.species.eq('B'), 'status'].tolist() == [
        'zero_copy', 'multiple_copy', 'target_not_found', 'observed']
    assert result.B_run.iloc[:3].isna().all()
    assert result.B_run.iloc[3] == 0


@pytest.mark.integration
def test_diagnostics_never_replace_applied_imputation_factors(tmp_path, stub_pdf_rendering):
    # A/B have the review's block missingness. All original targets remain in
    # merge, so the library sizes stay fixed while orthology annotation varies.
    merge = tmp_path / 'merge'
    genes = np.arange(1, 21) * 100.
    for species, depth in [('A', 1), ('B', 10)]:
        directory = merge / f'Species_{species}'
        directory.mkdir(parents=True)
        pd.DataFrame({'target_id': [f'g{i}' for i in range(20)], 'R1': genes * depth,
                      'R2': genes * depth}).to_csv(directory / f'Species_{species}_est_counts.tsv', sep='\t', index=False)
    pd.DataFrame({'scientific_name': ['Species A', 'Species A', 'Species B', 'Species B'],
                  'run': ['R1', 'R2', 'R1', 'R2'], 'exclusion': 'no', 'sample_group': 'tissue'}).to_csv(
        merge / 'metadata.tsv', sep='\t', index=False)
    og = tmp_path / 'og.tsv'
    pd.DataFrame({'busco_id': [f'og{i}' for i in range(20)], 'Species_A': [f'g{i}' for i in range(20)],
                  'Species_B': ['' if 4 <= i < 12 else f'g{i}' for i in range(20)]}).to_csv(og, sep='\t', index=False)
    args = SimpleNamespace(out_dir=str(tmp_path), dir_count='inferred', orthogroup_table=str(og), dir_busco=None,
                           redo=True, tmm_imputation_scale='raw', tmm_reference_diagnostics=False)
    cstmm_main(args)
    count_path = tmp_path / 'cstmm' / 'Species_B' / 'Species_B_cstmm_counts.tsv'
    before = count_path.read_bytes()
    metadata_before = (tmp_path / 'cstmm' / 'metadata.tsv').read_bytes()
    assert not (tmp_path / 'cstmm' / 'cstmm_observed_pair_comparison.pdf').exists()
    del args.tmm_reference_diagnostics  # Ordinary execution enables the real-input comparison.
    cstmm_main(args)
    assert (tmp_path / 'cstmm' / 'cstmm_observed_pair_comparison.pdf').is_file()
    assert count_path.read_bytes() == before
    assert (tmp_path / 'cstmm' / 'metadata.tsv').read_bytes() == metadata_before
    manifest = json.loads((tmp_path / 'cstmm' / 'cstmm_normalization.json').read_text())
    assert manifest['factor_source'] == 'imputed_reference'
    assert manifest['single_copy_threshold'] == 50
    pairs = pd.read_csv(tmp_path / 'cstmm' / 'cstmm_observed_pair_diagnostics.tsv', sep='\t')
    ab = pairs.loc[pairs.reference_sample.eq('Species_A_R1') & pairs['sample'].eq('Species_B_R1')].iloc[0]
    assert ab.observed_factor_ratio == pytest.approx(1)
    assert ab.applied_factor_ratio == pytest.approx(.1)
    corrected = pd.read_csv(count_path, sep='\t')
    # Only original targets are saved, not the internal orthogroup imputation.
    np.testing.assert_allclose(corrected.R1, genes * 10 * np.sqrt(10))


@pytest.mark.integration
def test_production_accepts_real_absences_without_any_universal_ortholog(tmp_path, stub_pdf_rendering):
    orthology = {'busco_id': [f'og{i}' for i in range(8)]}
    metadata = []
    for species in range(8):
        token = f'Species_{species}'
        genes = [f'g{i}' for i in range(8) if (i + species) % 8 < 4]
        directory = tmp_path / 'merge' / token
        directory.mkdir(parents=True)
        pd.DataFrame({'target_id': genes, 'R1': 10. ** (species % 3)}).to_csv(
            directory / f'{token}_est_counts.tsv', sep='\t', index=False)
        orthology[token] = [f'g{i}' if f'g{i}' in genes else '' for i in range(8)]
        metadata.append(dict(scientific_name=f'Species {species}', run='R1', exclusion='no', sample_group='tissue'))
    og = tmp_path / 'og.tsv'
    pd.DataFrame(orthology).to_csv(og, sep='\t', index=False)
    pd.DataFrame(metadata).to_csv(tmp_path / 'merge' / 'metadata.tsv', sep='\t', index=False)
    cstmm_main(SimpleNamespace(out_dir=str(tmp_path), dir_count='inferred', orthogroup_table=str(og), dir_busco=None))
    manifest = json.loads((tmp_path / 'cstmm' / 'cstmm_normalization.json').read_text())
    assert manifest['complete_rows'] == 0 and manifest['missing_cells'] == 32
    assert manifest['scale'] == 'library_size' and manifest['converged']
    output = pd.read_csv(tmp_path / 'cstmm' / 'metadata.tsv', sep='\t')
    assert output.exclusion.eq('no').all()
    np.testing.assert_allclose(output.tmm_normalization_factor, 1, rtol=1e-12)
    for species in range(8):
        token = f'Species_{species}'
        corrected = pd.read_csv(tmp_path / 'cstmm' / token / f'{token}_cstmm_counts.tsv', sep='\t')
        assert len(corrected) == 4
        np.testing.assert_allclose(corrected.R1, 10. ** (species % 3), rtol=1e-12)


def test_global_reference_is_used_directly_even_when_not_first(tmp_path):
    counts = pd.DataFrame({'A': [10., 20., 30.], 'NA': [20., 25., 10.], 'C': [30., 15., 20.]})
    libraries = pd.Series([100., 120., 200.], index=counts.columns)
    from amalgkit.normalization_tmm import calc_norm_factors_tmm
    factors = calc_norm_factors_tmm(counts, libraries, ref_column=1)
    path = tmp_path / 'pairs.tsv'
    write_observed_pair_diagnostics(counts, libraries, factors, path, global_reference='NA')
    table = pd.read_csv(path, sep='\t', keep_default_na=False)
    fixed = table.loc[table.reference_sample.eq('NA')]
    assert set(fixed['sample']) == {'A', 'C'}
    np.testing.assert_allclose(fixed.observed_factor_ratio, fixed.applied_factor_ratio, rtol=1e-12)
    fig = plot_observed_pair_diagnostics(path, tmp_path / 'comparison.pdf', 'NA')
    xy = fig.axes[0].collections[0].get_offsets()
    np.testing.assert_allclose(xy[:, 0], xy[:, 1], atol=1e-12)
    assert (tmp_path / 'comparison.pdf').read_bytes().startswith(b'%PDF')
    # Global centering is arbitrary; multiplying all f by a constant changes nothing.
    write_observed_pair_diagnostics(counts, libraries, factors * 7, path, global_reference='NA')
    rescaled = pd.read_csv(path, sep='\t', keep_default_na=False)
    np.testing.assert_allclose(rescaled.applied_factor_ratio, table.applied_factor_ratio, rtol=1e-12)


@pytest.mark.parametrize('singleton', [False, True])
def test_plot_counts_unestimable_comparisons_without_inventing_points(tmp_path, singleton):
    counts = pd.DataFrame({'A': [10., np.nan], 'B': [np.nan, 20.]})
    if singleton:
        counts = counts[['A']]
    path = tmp_path / 'pairs.tsv'
    write_observed_pair_diagnostics(counts, counts.sum(), pd.Series(1., index=counts.columns), path, global_reference='A')
    fig = plot_observed_pair_diagnostics(path, tmp_path / 'empty.png', 'A')
    assert not fig.axes[0].collections
    assert f'{0 if singleton else 1} not estimable' in fig.axes[0].get_title()


@pytest.mark.parametrize('target_id', ['', '   '])
def test_cstmm_rejects_empty_target_ids_without_treating_literal_na_as_missing(tmp_path, target_id):
    directory = tmp_path / 'Species_A'
    directory.mkdir()
    path = directory / 'Species_A_est_counts.tsv'
    path.write_text(f'target_id\tR1\n{target_id}\t10\nNA\t20\n')
    with pytest.raises(ValueError, match='nonempty target IDs'):
        _read_est_counts(tmp_path, 'Species_A')
    path.write_text('target_id\tR1\n001\t10\nNA\t20\n')
    counts = _read_est_counts(tmp_path, 'Species_A')
    assert counts.index.tolist() == ['001', 'NA']
    np.testing.assert_array_equal(counts.Species_A_R1, [10, 20])
