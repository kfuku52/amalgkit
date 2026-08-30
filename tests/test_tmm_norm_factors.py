"""Direct numeric tests for amalgkit's TMM normalization factors.

The production TMM implementation (amalgkit/normalization_tmm.py) follows
Robinson & Oshlack (2010), with deterministic ties at trimming/reference
boundaries that can differ from edgeR's floating-point ordering. The tests
pin it to numeric values that were verified independently:

* ``calc_factor_tmm`` on a fixed pair of count vectors, hand-derived to
  ``352/2047`` (the trimmed mean of M-values collapses to a single gene).
* ``calc_norm_factors_tmm`` on the classic edgeR toy example (50 genes,
  four libraries, controls 10 vs cases 20) whose published factors are
  ``1/sqrt(2)`` for the controls and ``sqrt(2)`` for the cases.

They also retain the original structural checks (positive factors,
reference-column bounds, symmetric identical-column behavior).
"""

import json
from pathlib import Path

import numpy
import pandas
import pytest

from amalgkit.normalization_tmm import (
    _resolve_median_reference_column,
    apply_tmm_factors,
    calc_factor_tmm,
    calc_norm_factors_tmm,
    run_tmm_rounds_for_cstmm,
)

FIXTURES = Path(__file__).parent / 'fixtures'
ORACLE_CASES = json.loads((FIXTURES / 'tmm.json').read_text())['cases']


@pytest.mark.parametrize('case', ORACLE_CASES, ids=lambda case: 'oracle-' + str(case['case']))
def test_tmm_rounds_match_independent_high_precision_oracle(case):
    counts = numpy.array(case['counts'])
    result = run_tmm_rounds_for_cstmm(counts, lib_size=case['library_sizes'])
    assert result.round1_reference_column == case['round1_reference']
    assert result.median_reference_columns == [case['round2_reference']]
    numpy.testing.assert_allclose(result.round1_factors, case['round1'], rtol=1e-12, atol=1e-12)
    numpy.testing.assert_allclose(result.round2_factors, case['round2'], rtol=1e-12, atol=1e-12)
    fixed = calc_norm_factors_tmm(counts, lib_size=case['library_sizes'], ref_column=0)
    numpy.testing.assert_allclose(fixed, case['fixed_first'], rtol=1e-12, atol=1e-12)
    # Joint power-of-two scaling preserves rates and relative binomial weights.
    scaled = calc_norm_factors_tmm(counts[::-1] * 8, lib_size=numpy.array(case['library_sizes']) * 8, ref_column=0)
    numpy.testing.assert_allclose(scaled, fixed, rtol=1e-12, atol=1e-12)


def test_fixed_reference_non_boundary_cases_match_edger():
    reference = json.loads((FIXTURES / 'tmm-edger.json').read_text())
    assert reference['edgeR_version']
    for case in ORACLE_CASES:
        expected = reference['fixed_first'].get(str(case['case']))
        if expected is None:
            continue
        observed = calc_norm_factors_tmm(case['counts'], lib_size=case['library_sizes'], ref_column=0)
        numpy.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)


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


@pytest.mark.parametrize('dtype', ['int64', 'float64'])
def test_apply_tmm_factors_divides_each_sample_column(dtype):
    counts = pandas.DataFrame(
        {
            'RUN1': [10.0, 20.0],
            'RUN2': [30.0, 40.0],
        },
        index=['G1', 'G2'],
    )
    counts = counts.astype(dtype)
    before = counts.copy(deep=True)
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
    pandas.testing.assert_frame_equal(counts, before)


def test_even_median_reference_does_not_depend_on_rounded_midpoint():
    # A rounded midpoint makes the second value one ULP closer by subtraction.
    values = [0.7, numpy.nextafter(0.8, 1.0)]
    assert _resolve_median_reference_column(values) == 0
    assert _resolve_median_reference_column(list(reversed(values))) == 0


def test_tmm_discards_zero_overlap_without_runtime_warnings():
    assert calc_factor_tmm([0, 2, 0], [3, 0, 0], 10, 20) == 1.0


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
