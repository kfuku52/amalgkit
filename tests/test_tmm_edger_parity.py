import numpy
import pandas

from amalgkit.normalization_tmm import (
    _resolve_tmm_reference_column,
    apply_tmm_factors,
    calc_norm_factors_tmm,
    run_tmm_rounds_for_cstmm,
)


def test_tmm_reference_column_matches_edger_raw_f75_rule():
    # edgeR's calcNormFactors selects the reference from the RAW 75th
    # percentile: f75 <- apply(counts, 2, quantile, probs = 0.75);
    # refColumn <- which.min(abs(f75 - mean(f75))). Library sizes must NOT
    # enter this selection. With heterogeneous library sizes the previous
    # implementation (f75 divided by lib.size) picked a different column
    # (391/500 randomized trials diverged from edgeR).
    counts = pandas.DataFrame(
        {
            'SP1': [100, 200, 50, 10, 5, 1],
            'SP2': [50, 400, 30, 8, 4, 2],
            'SP3': [200, 100, 60, 12, 6, 3],
            'SP4': [10, 20, 90, 5, 2, 1],
        },
        index=[f'G{i}' for i in range(1, 7)],
    )
    library_sizes = pandas.Series([1000.0, 2000.0, 500.0, 3000.0], index=list(counts.columns), dtype=float)

    matrix = counts.to_numpy(dtype=float)
    f75_raw = numpy.asarray(
        [
            float(numpy.quantile(matrix[:, idx], 0.75, method='linear'))
            for idx in range(matrix.shape[1])
        ],
        dtype=float,
    )
    expected_ref = int(numpy.argmin(numpy.abs(f75_raw - numpy.mean(f75_raw))))

    resolved = _resolve_tmm_reference_column(matrix, library_sizes.to_numpy(dtype=float))

    assert resolved == expected_ref
    # This fixture specifically exercises the divergence: the raw-f75 rule
    # (edgeR) picks column 1, while the old f75/lib.size rule picked column 0.
    assert expected_ref == 1
    assert f75_raw[0] != f75_raw[1]


def test_tmm_factors_golden_values_with_edger_reference_selection():
    # Golden-value regression: factor values produced with the edgeR reference
    # selection rule. Computed from a faithful reimplementation of edgeR's
    # calcFactorTMM (log-ratio, average-expression, variance weights,
    # quantile trimming floor(n*trim)+1 .. n+1-lo, 2^f) with the reference
    # column chosen by the raw-f75 rule.
    counts = pandas.DataFrame(
        {
            'SP1': [100, 200, 50, 10, 5, 1],
            'SP2': [50, 400, 30, 8, 4, 2],
            'SP3': [200, 100, 60, 12, 6, 3],
            'SP4': [10, 20, 90, 5, 2, 1],
        },
        index=[f'G{i}' for i in range(1, 7)],
    )
    library_sizes = pandas.Series([1000.0, 2000.0, 500.0, 3000.0], index=list(counts.columns), dtype=float)

    observed = calc_norm_factors_tmm(counts, lib_size=library_sizes, ref_column=None)

    expected = pandas.Series(
        [2.122576, 0.685809, 5.051604, 0.135989],
        index=list(counts.columns),
        dtype=float,
    )
    numpy.testing.assert_allclose(observed.to_numpy(dtype=float), expected.to_numpy(dtype=float), rtol=1e-5, atol=1e-5)


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
    assert len(observed.median_reference_columns) >= 1
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
