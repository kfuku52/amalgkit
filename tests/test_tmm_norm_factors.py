"""Direct numeric tests for amalgkit's TMM normalization factors.

The production TMM implementation (amalgkit/normalization_tmm.py) is a
faithful re-implementation of edgeR's ``calcNormFactors(method="TMM")``
(Robinson & Oshlack 2010).  The tests here pin it to real numeric values
that were verified independently:

* ``calc_factor_tmm`` on a fixed pair of count vectors, hand-derived to
  ``352/2047`` (the trimmed mean of M-values collapses to a single gene).
* ``calc_norm_factors_tmm`` on the classic edgeR toy example (50 genes,
  four libraries, controls 10 vs cases 20) whose published factors are
  ``1/sqrt(2)`` for the controls and ``sqrt(2)`` for the cases.

They also retain the original structural checks (positive factors,
reference-column bounds, symmetric identical-column behavior).
"""

import numpy
import pandas

from amalgkit.normalization_tmm import (
    _resolve_median_reference_column,
    apply_tmm_factors,
    calc_factor_tmm,
    calc_norm_factors_tmm,
    run_tmm_rounds_for_cstmm,
)


def test_run_tmm_rounds_for_cstmm_returns_positive_factors_and_reference_columns():
    counts = pandas.DataFrame(
        {
            'SP1_RUN1': [100, 110, 120, 130, 10, 11, 12, 13],
            'SP1_RUN2': [90, 100, 110, 120, 20, 21, 22, 23],
            'SP2_RUN1': [80, 85, 90, 95, 40, 41, 42, 43],
            'SP2_RUN2': [60, 65, 70, 75, 50, 51, 52, 53],
        },
        index=[f'OG{i}' for i in range(1, 9)],
    )
    library_sizes = pandas.Series(
        [2000.0, 2200.0, 1800.0, 1900.0],
        index=list(counts.columns),
        dtype=float,
    )

    observed = run_tmm_rounds_for_cstmm(counts=counts, lib_size=library_sizes)

    assert observed.round1_factors.index.tolist() == list(counts.columns)
    assert observed.round2_factors.index.tolist() == list(counts.columns)
    assert numpy.all(observed.round1_factors.to_numpy(dtype=float) > 0)
    assert numpy.all(observed.round2_factors.to_numpy(dtype=float) > 0)
    assert 0 <= observed.round1_reference_column < counts.shape[1]
    assert len(observed.median_reference_columns) == 1
    assert all(0 <= idx < counts.shape[1] for idx in observed.median_reference_columns)


def test_run_tmm_rounds_for_cstmm_keeps_identical_columns_balanced():
    counts = pandas.DataFrame(
        {
            'RUN1': [10, 20, 30, 40],
            'RUN2': [10, 20, 30, 40],
            'RUN3': [10, 20, 30, 40],
            'RUN4': [10, 20, 30, 40],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    library_sizes = pandas.Series([100.0, 100.0, 100.0, 100.0], index=list(counts.columns), dtype=float)

    observed = run_tmm_rounds_for_cstmm(counts=counts, lib_size=library_sizes)

    numpy.testing.assert_allclose(observed.round1_factors.to_numpy(dtype=float), numpy.ones(4), rtol=0.0, atol=1e-12)
    numpy.testing.assert_allclose(observed.round2_factors.to_numpy(dtype=float), numpy.ones(4), rtol=0.0, atol=1e-12)
    assert observed.median_reference_columns == [0]


def test_median_tmm_reference_uses_true_even_median_and_one_column():
    assert _resolve_median_reference_column([101.0, 1.0, 99.0, 2.0]) == 2


def test_apply_tmm_factors_divides_each_sample_column():
    counts = pandas.DataFrame(
        {
            'RUN1': [10.0, 20.0],
            'RUN2': [30.0, 40.0],
        },
        index=['G1', 'G2'],
    )
    factors = pandas.Series([2.0, 4.0], index=['RUN1', 'RUN2'])

    corrected = apply_tmm_factors(counts, factors)

    expected = pandas.DataFrame(
        {
            'RUN1': [5.0, 10.0],
            'RUN2': [7.5, 10.0],
        },
        index=['G1', 'G2'],
    )
    pandas.testing.assert_frame_equal(corrected, expected)


def test_calc_factor_tmm_matches_hand_computed_trimmed_mean_of_m_values():
    """Hand-derived oracle for calc_factor_tmm.

    obs = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024] (lib.size 2047)
    ref = [1] * 11                                      (lib.size 11)

    With n = 11, the default trims give log-ratio rank bounds of [4, 8] and
    absolute-level rank bounds of [1, 11], so exactly the five middle genes
    (counts 8, 16, 32, 64, 128) survive.  Their log-ratios are
    k + log2(11/2047) for k = 3..7, whose (unweighted) mean is
    5 + log2(11/2047).  The TMM factor is therefore

        2 ** (5 + log2(11/2047)) = 2 ** 5 * 11 / 2047 = 352 / 2047.
    """
    obs = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    ref = [1] * 11

    factor = calc_factor_tmm(
        obs=obs,
        ref=ref,
        libsize_obs=2047.0,
        libsize_ref=11.0,
        do_weighting=False,
    )

    expected = 352.0 / 2047.0
    assert abs(factor - expected) < 1e-12


def test_calc_norm_factors_tmm_matches_edger_toy_example():
    """Classic edgeR toy example: 50 genes, controls 10 vs cases 20.

    Every gene has the same fold change between the two groups, so the raw
    pairwise TMM factors are 2 for each case relative to each control and 1
    within a group.  The symmetry (geometric-mean) adjustment then divides by
    (2 * 2 * 1 * 1) ** 0.25 = sqrt(2), yielding published edgeR factors of
    1/sqrt(2) for the controls and sqrt(2) for the cases.
    """
    rows = [[10, 10, 20, 20]] * 50
    counts = pandas.DataFrame(rows, columns=['cont1', 'cont2', 'case1', 'case2'])
    library_sizes = pandas.Series([500.0, 500.0, 500.0, 500.0], index=counts.columns, dtype=float)

    factors = calc_norm_factors_tmm(counts, lib_size=library_sizes)

    expected = pandas.Series(
        [1.0 / numpy.sqrt(2.0)] * 2 + [numpy.sqrt(2.0)] * 2,
        index=counts.columns,
        dtype=float,
    )
    pandas.testing.assert_series_equal(
        factors,
        expected,
        check_exact=False,
        rtol=0.0,
        atol=1e-9,
    )
    # Symmetry constraint: normalized factors multiply to 1.
    numpy.testing.assert_allclose(factors.prod(), 1.0, rtol=0.0, atol=1e-9)
