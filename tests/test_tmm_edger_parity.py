import numpy
import pandas

from amalgkit.normalization_tmm import (
    apply_tmm_factors,
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
