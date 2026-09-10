"""Numerical kernels tested against fixed, independently generated R outputs."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from amalgkit.batch_effect_sva import clean_y_matrix, run_sva_backend
from amalgkit.batch_effect_ruvseq import ruvr_correct_counts
from amalgkit.batch_effect_combatseq import run_combatseq_backend


ROOT = Path(__file__).parent / 'fixtures' / 'batch_effect_reference'


def read(name):
    return pd.read_csv(ROOT / name, sep='\t', index_col=0)


def projector(matrix):
    q = np.linalg.qr(np.asarray(matrix))[0]
    return q @ q.T


def test_sva_protected_projection_matches_independent_r_calculation():
    counts = read('counts.tsv')
    metadata = pd.read_csv(ROOT / 'metadata.tsv', sep='\t')
    design = np.column_stack([np.ones(len(metadata)), metadata['sample_group']])
    corrected = clean_y_matrix(np.log1p(counts), design, read('sva_factors.tsv'))
    np.testing.assert_allclose(corrected, read('sva_cleaned.tsv'), rtol=1e-12, atol=1e-12)


def test_sva_manual_fit_keeps_data_supported_signal_and_separates_irw_from_permutations():
    counts = np.log1p(read('counts.tsv'))
    metadata = pd.read_csv(ROOT / 'metadata.tsv', sep='\t')
    first, w, summary = run_sva_backend(counts, metadata, nsv_setting=1, B_setting=20,
                                       irw_iterations=5, failure_policy='error')
    second, w2, _ = run_sva_backend(counts, metadata, nsv_setting=1, B_setting=100,
                                   irw_iterations=5, failure_policy='error')
    np.testing.assert_allclose(first, second, rtol=0, atol=1e-12)
    np.testing.assert_allclose(projector(w), projector(w2), atol=1e-12)
    assert summary['sva_irw_iterations_completed'] == 5
    assert summary['sva_irw_converged'] is None
    # Independent generating signal, not the score used to choose k.
    assert abs(np.corrcoef(w.iloc[:, 0], metadata['bioproject'])[0, 1]) > 0.95


def test_ruvr_kernel_matches_upstream_with_identical_edgeR_residuals():
    counts = read('counts.tsv')
    corrected, w = ruvr_correct_counts(counts, np.ones(len(counts), dtype=bool), 1, read('ruv_residuals.tsv'))
    np.testing.assert_allclose(projector(w), projector(read('ruv_factors.tsv')), atol=1e-12)
    np.testing.assert_array_equal(corrected.to_numpy(), read('ruv_counts.tsv').to_numpy())


@pytest.mark.optional_dependency
def test_combat_matches_r_reference_on_balanced_raw_counts():
    pytest.importorskip('inmoose.pycombat')
    counts = read('counts.tsv')
    metadata = pd.read_csv(ROOT / 'metadata.tsv', sep='\t')
    corrected, _ = run_combatseq_backend(counts, metadata, failure_policy='error')
    np.testing.assert_array_equal(corrected.to_numpy(), read('combat_counts.tsv').to_numpy())
