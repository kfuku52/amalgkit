"""Numerical and saved-data contracts affected by the efficiency kernels."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from amalgkit.batch_effect_sva import _supported_right_vectors
from amalgkit.getfastq_sampling import _sorted_selected_indices, selected_indices, random_sampling, random_sampling_mask
from amalgkit.gsa_fastq import has_gsa_rows, all_gsa_rows, is_gsa_row
from amalgkit.linalg_utils import classical_mds_coordinates
from amalgkit.merge_plots import _read_est_count_summary


@pytest.mark.parametrize('rank', [0, 1, 2, 3])
def test_sva_retains_direct_svd_rank_and_subspace(rank):
    rng = np.random.default_rng(719)
    values = rng.normal(size=(400, rank)) @ rng.normal(size=(rank, 30))
    if rank:
        values += rng.normal(size=values.shape) * 1e-18
    _, singular, right = np.linalg.svd(values, full_matrices=False)
    tolerance = np.finfo(float).eps * max(values.shape) * max(1., np.linalg.norm(values))
    expected = right[:min(3, int((singular > tolerance).sum()))].T
    observed = _supported_right_vectors(values, 3)
    assert observed.shape == expected.shape
    np.testing.assert_allclose(observed @ observed.T, expected @ expected.T, atol=1e-12)


@pytest.mark.parametrize('scale', [1e-80, 1., 1e80])
def test_sva_separated_full_rank_vectors_match_direct_svd(scale):
    values = np.random.default_rng(719).normal(size=(1000, 40)) * scale
    _, singular, right = np.linalg.svd(values, full_matrices=False)
    tolerance = np.finfo(float).eps * max(values.shape) * max(1., np.linalg.norm(values))
    expected = right[:min(3, int((singular > tolerance).sum()))].T
    observed = _supported_right_vectors(values, 3)
    assert observed.shape == expected.shape
    np.testing.assert_allclose(observed @ observed.T, expected @ expected.T, atol=1e-11)


@pytest.mark.parametrize('total,start,end', [(101, 1, 50), (101, 21, 101), (100000, 11, 70),
                                           (31, 1, 31), (2**100, 1, 10)])
def test_compact_sampling_preserves_permutation_ranks(total, start, end):
    for seed in [0, 719, 42]:
        observed = list(_sorted_selected_indices(total, start, end, seed, '0001'))
        assert observed == sorted(selected_indices(total, start, end, seed, '0001'))


@pytest.mark.parametrize('method,private', [('contiguous', False), ('random', False), ('random', True)])
def test_vectorized_source_and_sampling_flags_keep_row_semantics(method, private):
    frame = pd.DataFrame({'private_file': ['yes', ' YES ', '', None, pd.NA, 1],
                          'data_source': ['GSA', ' gsa ', 'sra', None, pd.NA, 1]}, index=[8, 2, 9, 4, 3, 6])
    args = SimpleNamespace(sampling_method=method, sampling_private=private)
    assert random_sampling_mask(args, frame).tolist() == [random_sampling(args, row) for _, row in frame.iterrows()]
    assert has_gsa_rows(frame) == any(is_gsa_row(row) for _, row in frame.iterrows())
    assert all_gsa_rows(frame) == all(is_gsa_row(row) for _, row in frame.iterrows())
    assert not has_gsa_rows(frame.drop(columns='data_source'))
    assert all_gsa_rows(frame.iloc[:0])


def test_streaming_plot_summary_keeps_missing_and_fractional_counts(tmp_path):
    rng = np.random.default_rng(719)
    frame = pd.DataFrame(rng.gamma(2, 20, size=(12001, 4)), columns=['0001', 'NA', 'mixed', 'empty'])
    frame['empty'] = np.nan
    frame.loc[::11, 'NA'] = np.nan
    frame['mixed'] = frame['mixed'].astype(object)
    frame.loc[11000, 'mixed'] = 'unparseable'
    frame.insert(0, 'target_id', [f'g{i}' for i in range(len(frame))])
    frame.loc[5000, 'target_id'] = 'g\twith\nquoted "identifier"'
    path = tmp_path / 'Species_A' / 'Species_A_est_counts.tsv'
    path.parent.mkdir()
    frame.to_csv(path, sep='\t', index=False)
    numeric = pd.read_csv(path, sep='\t').drop(columns='target_id').apply(pd.to_numeric, errors='coerce')
    observed = _read_est_count_summary(path)
    assert observed.run.tolist() == ['0001', 'NA', 'mixed', 'empty']
    np.testing.assert_allclose(observed.lib_size, numeric.sum(), rtol=1e-14, atol=1e-10)
    np.testing.assert_allclose(observed.mean_before, numeric.mean(), rtol=1e-14, equal_nan=True)


@pytest.mark.parametrize('compressed', [False, True])
def test_streaming_plot_summary_does_not_hide_malformed_rows(tmp_path, compressed):
    import gzip
    path = tmp_path / ('bad.tsv.gz' if compressed else 'bad.tsv')
    text = 'target_id\tR1\tR2\n' + 'g1\t1\t2\n' * 5000 + 'bad\t3\t4\textra\n'
    if compressed:
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
    else:
        path.write_text(text)
    with pytest.warns(UserWarning, match='Failed to read est_counts file'):
        assert _read_est_count_summary(path) is None


@pytest.mark.parametrize('size', [0, 1, 20, 100])
def test_mds_coordinates_preserve_distance_geometry(size):
    rng = np.random.default_rng(719)
    points = rng.normal(size=(size, 2))
    distances = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
    coords = classical_mds_coordinates(distances)
    assert coords.shape == (size, 2)
    observed = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    if size:
        np.testing.assert_allclose(observed, distances, atol=1e-11)


def test_mds_tied_spectrum_keeps_legacy_embedding():
    distances = np.ones((80, 80)) - np.eye(80)
    center = np.eye(80) - np.ones((80, 80)) / 80
    values, vectors = np.linalg.eigh(-.5 * center.dot(distances**2).dot(center))
    order = np.argsort(values)[::-1][:2]
    expected = vectors[:, order] * np.sqrt(np.maximum(values[order], 0))
    np.testing.assert_array_equal(classical_mds_coordinates(distances), expected)


def test_mds_partial_solver_failure_retains_full_solver(monkeypatch):
    import scipy.linalg
    points = np.random.default_rng(719).normal(size=(80, 2))
    distances = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
    def failed(*args, **kwargs):
        raise np.linalg.LinAlgError('partial solver failed')
    monkeypatch.setattr(scipy.linalg, 'eigh', failed)
    coords = classical_mds_coordinates(distances)
    observed = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    np.testing.assert_allclose(observed, distances, atol=1e-11)


@pytest.mark.parametrize('returned', [0, 1, 2])
def test_mds_incomplete_partial_result_retains_legacy_solver(monkeypatch, returned):
    import scipy.linalg
    distances = np.ones((80, 80)) - np.eye(80)
    center = np.eye(80) - np.ones((80, 80)) / 80
    values, vectors = np.linalg.eigh(-.5 * center.dot(distances**2).dot(center))
    order = np.argsort(values)[::-1][:2]
    expected = vectors[:, order] * np.sqrt(np.maximum(values[order], 0))
    monkeypatch.setattr(scipy.linalg, 'eigh', lambda *args, **kwargs: (values[:returned], vectors[:, :returned]))
    np.testing.assert_array_equal(classical_mds_coordinates(distances), expected)


def test_overview_reuses_correlation_in_metadata_order(tmp_path, monkeypatch, stub_pdf_encoding):
    from amalgkit import per_species_outputs as outputs
    counts = pd.DataFrame({'0001': [1., 2., 4.], 'NA': [4., 3., 1.], 'later': [2., 4., 3.]})
    metadata = pd.DataFrame({'run': ['NA', 'later', '0001'], 'sample_group': ['A', 'B', 'A'],
                             'bioproject': ['P1', 'P2', 'P2'], 'exclusion': 'no'})
    correlation = counts.corr()
    expected = correlation.loc[metadata.run, metadata.run]
    def forbidden(*args, **kwargs):
        raise AssertionError('correlation was recomputed')
    def capture(frame):
        pd.testing.assert_frame_equal(frame, expected)
        return np.zeros((len(frame), 2))
    monkeypatch.setattr(outputs, '_compute_corr_matrix', forbidden)
    monkeypatch.setattr(outputs, '_compute_pca_coords', capture)
    monkeypatch.setattr(outputs, '_compute_mds_coords', capture)
    monkeypatch.setattr(outputs, '_compute_tsne_coords', lambda frame: np.zeros((frame.shape[1], 2)))
    outputs.save_state_overview_pdf(counts, metadata, ['A', 'B'], str(tmp_path / 'plot.pdf'),
                                    transform_method='none-none', precomputed_tc_dist_matrix=correlation)


def test_observed_pair_diagnostics_keep_all_pairs_and_reference_direction(tmp_path):
    from amalgkit.cstmm_diagnostics import write_observed_pair_diagnostics
    counts = pd.DataFrame({name: np.arange(1., 17.) for name in ['0001', 'NA', 'later', 'last']})
    libraries = counts.sum().to_dict()
    factors = pd.Series([1., 2., 3., 4.], index=counts.columns)
    path = tmp_path / 'pairs.tsv'
    write_observed_pair_diagnostics(counts, libraries, factors, path, global_reference='later')
    result = pd.read_csv(path, sep='\t', dtype={'reference_sample':str, 'sample':str}, keep_default_na=False)
    assert list(zip(result.reference_sample, result['sample'])) == [
        ('0001','NA'), ('later','0001'), ('0001','last'),
        ('later','NA'), ('NA','last'), ('later','last'),
    ]
    assert result.shared_observed.eq(16).all()
    np.testing.assert_allclose(result.observed_factor_ratio.astype(float), 1.)
    for row in result.itertuples(index=False):
        assert row.applied_factor_ratio == pytest.approx(factors[row.sample] / factors[row.reference_sample])
