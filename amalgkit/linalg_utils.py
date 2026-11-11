"""Small numerical stability checks shared by decomposition backends."""

from __future__ import annotations

from collections.abc import Sequence

import numpy


def classical_mds_coordinates(dist):
    """Return the leading two classical MDS coordinates with bounded copies.

    Large, separated spectra use analytic double centering and a partial
    eigensolve. A tied/undefined boundary retains the established full solver.
    """
    n = dist.shape[0]
    coords = numpy.zeros((n, 2), dtype=float)
    if n <= 1:
        return coords
    squared = dist**2
    if n >= 64 and numpy.isfinite(squared).all():
        from scipy.linalg import eigh

        gram = -0.5 * (squared - squared.mean(axis=1)[:, None] - squared.mean(axis=0)[None, :] + squared.mean())
        # Accept only the requested complete subspace. Some backends return
        # fewer eigenpairs at an exactly tied spectral boundary.
        try:
            eigvals, eigvecs = eigh(gram, subset_by_index=(n - 3, n - 1), driver="evx")
        except numpy.linalg.LinAlgError:
            pass  # The full solver below retains its original result/error.
        else:
            if (
                eigvals.shape == (3,)
                and eigvecs.shape == (n, 3)
                and numpy.isfinite(eigvals).all()
                and numpy.isfinite(eigvecs).all()
            ):
                scale = max(1.0, float(numpy.max(numpy.abs(eigvals))))
                tolerance = numpy.finfo(float).eps * n * scale
                if eigvals[1] > tolerance and eigvals[1] - eigvals[0] > tolerance:
                    return eigvecs[:, ::-1][:, :2] * numpy.sqrt(numpy.maximum(eigvals[::-1][:2], 0.0))
    centering = numpy.eye(n) - numpy.ones((n, n), dtype=float) / float(n)
    gram = -0.5 * centering.dot(squared).dot(centering)
    eigvals, eigvecs = numpy.linalg.eigh(gram)
    order = numpy.argsort(eigvals)[::-1][:2]
    coords[:, : len(order)] = eigvecs[:, order] * numpy.sqrt(numpy.maximum(eigvals[order], 0.0))
    return coords


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


def leading_gram_components(
    gram: numpy.ndarray, num_components: int, matrix_shape: Sequence[int], *, partial: bool
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Return descending eigenvalues and the leading vectors of a Gram matrix.

    For a small leading subspace of a large Gram matrix, keep one additional
    eigenvalue to verify the gap to the excluded directions. Use the full
    spectrum when the partial solver cannot establish that gap.
    """
    size = gram.shape[0]
    if partial and size >= 512 and num_components <= size // 8:
        from scipy.linalg import eigh

        width = num_components + 1
        try:
            values, vectors = eigh(gram, subset_by_index=(size - width, size - 1), driver="evx")
        except (numpy.linalg.LinAlgError, ValueError):
            pass  # The full solver retains its result or authoritative error.
        else:
            if (
                values.shape == (width,)
                and vectors.shape == (size, width)
                and numpy.isfinite(values).all()
                and numpy.isfinite(vectors).all()
                and not gram_components_require_svd_fallback(values[::-1], num_components, matrix_shape)
            ):
                return values[::-1], vectors[:, ::-1][:, :num_components]
    values, vectors = numpy.linalg.eigh(gram)
    order = numpy.argsort(values)[::-1]
    return values[order], vectors[:, order[:num_components]]
