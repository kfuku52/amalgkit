"""Small numerical stability checks shared by decomposition backends."""

from __future__ import annotations

from collections.abc import Sequence

import numpy


def gram_components_require_svd_fallback(
    eigenvalues_descending,
    num_components: int,
    matrix_shape: Sequence[int],
) -> bool:
    """Return whether an eigendecomposition boundary is numerically unstable."""
    eigenvalues = numpy.maximum(
        numpy.asarray(eigenvalues_descending, dtype=float),
        0.0,
    )
    if eigenvalues.size == 0:
        return True
    scale = float(eigenvalues[0])
    if not numpy.isfinite(scale) or scale <= 0.0:
        return True
    tolerance = numpy.finfo(float).eps * max(matrix_shape) * scale
    boundary_index = int(num_components) - 1
    if eigenvalues[boundary_index] <= tolerance:
        return True
    if int(num_components) < eigenvalues.size:
        boundary_gap = eigenvalues[boundary_index] - eigenvalues[int(num_components)]
        if boundary_gap <= tolerance:
            return True
    return False
