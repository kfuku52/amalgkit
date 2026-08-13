"""Direct tests for amalgkit.linalg_utils.gram_components_require_svd_fallback.

These exercise the numerical-stability boundaries that decide whether a
decomposition should fall back to an SVD: empty/zero/non-finite eigenvalue
spectra, an eigenvalue that sits on the tolerance boundary, and a degenerate
spectral gap between the kept and dropped components.
"""

import numpy

from amalgkit.linalg_utils import gram_components_require_svd_fallback

EPS = numpy.finfo(float).eps


def _tolerance(scale, shape):
    return EPS * max(shape) * scale


def test_empty_spectrum_requires_fallback():
    assert gram_components_require_svd_fallback(numpy.array([]), num_components=1, matrix_shape=(4, 4)) is True


def test_zero_leading_eigenvalue_requires_fallback():
    # numpy.maximum clamps negatives/zeros to 0 -> scale is 0 -> fallback.
    assert gram_components_require_svd_fallback([0.0, 0.0, 0.0], num_components=1, matrix_shape=(4, 4)) is True


def test_negative_leading_eigenvalue_requires_fallback():
    assert gram_components_require_svd_fallback([-5.0, 1.0], num_components=1, matrix_shape=(4, 4)) is True


def test_non_finite_spectrum_requires_fallback():
    assert gram_components_require_svd_fallback([numpy.nan, 1.0], num_components=1, matrix_shape=(4, 4)) is True


def test_boundary_eigenvalue_at_or_below_tolerance_requires_fallback():
    scale = 1.0
    shape = (3, 3)
    tolerance = _tolerance(scale, shape)
    # Eigenvalue[1] (the component boundary) sits below the tolerance.
    eigenvalues = [scale, tolerance / 2.0]
    assert gram_components_require_svd_fallback(eigenvalues, num_components=2, matrix_shape=shape) is True


def test_degenerate_spectral_gap_requires_fallback():
    # Kept and dropped components are numerically identical -> zero gap.
    eigenvalues = [1.0, 1.0, 1.0, 1.0]
    assert gram_components_require_svd_fallback(eigenvalues, num_components=3, matrix_shape=(3, 3)) is True


def test_stable_spectrum_does_not_require_fallback():
    eigenvalues = [100.0, 50.0, 10.0]
    assert gram_components_require_svd_fallback(eigenvalues, num_components=2, matrix_shape=(3, 3)) is False
