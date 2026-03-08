import numpy
import pandas
import amalgkit.batch_effect_sva as batch_effect_sva_module

from amalgkit.batch_effect_sva import (
    build_intercept_only_design_matrix,
    build_sample_group_design_matrix,
    clean_y_matrix,
    f_pvalue,
    irwsva_build,
    estimate_num_sv_be,
    estimate_num_sv_leek,
    estimate_num_sv_at_B,
    run_sva_backend,
)


def test_build_sample_group_design_matrix_matches_expected_treatment_coding():
    matrix, names = build_sample_group_design_matrix(['B', 'A', 'B'])
    assert names == ['(Intercept)', 'sample_groupB']
    assert matrix.tolist() == [
        [1.0, 1.0],
        [1.0, 0.0],
        [1.0, 1.0],
    ]


def test_build_intercept_only_design_matrix_returns_expected_shape():
    matrix, names = build_intercept_only_design_matrix(3)
    assert names == ['(Intercept)']
    assert matrix.shape == (3, 1)
    assert numpy.all(matrix == 1.0)


def test_clean_y_matrix_is_noop_when_sv_matrix_is_empty():
    y = numpy.array([[1.0, 2.0], [3.0, 4.0]])
    mod = numpy.array([[1.0], [1.0]])
    sv = numpy.zeros((2, 0))
    adjusted = clean_y_matrix(y, mod, sv)
    numpy.testing.assert_allclose(adjusted, y)


def test_run_sva_backend_returns_noop_for_explicit_nsv_zero():
    counts = pandas.DataFrame(
        {'RUN1': [1.0, 2.0], 'RUN2': [3.0, 4.0]},
        index=['G1', 'G2'],
    )
    metadata = pandas.DataFrame({
        'run': ['RUN1', 'RUN2'],
        'sample_group': ['A', 'B'],
    })
    corrected, sv_df, summary = run_sva_backend(
        counts_df=counts,
        metadata_df=metadata,
        nsv_setting='0',
        B_setting='auto',
        B_auto_max=100,
    )
    pandas.testing.assert_frame_equal(corrected, counts)
    assert list(sv_df.index) == ['RUN1', 'RUN2']
    assert sv_df.shape == (2, 0)
    assert summary['resolved_sva_nsv'] == 0
    assert summary['skip_reason'] == 'sva_nsv_zero'


def test_estimate_num_sv_be_returns_unresolved_for_constant_residual_matrix():
    data = numpy.array([
        [10.0, 10.0, 10.0, 10.0],
        [5.0, 5.0, 5.0, 5.0],
        [0.0, 0.0, 0.0, 0.0],
    ])
    mod, _ = build_sample_group_design_matrix(['A', 'A', 'B', 'B'])
    estimate = estimate_num_sv_be(
        data_matrix=data,
        mod_matrix=mod,
        B_value=5,
        max_nsv=1,
        random_seed=7,
    )
    assert estimate.method == 'be'
    assert estimate.nsv is None


def test_run_sva_backend_supports_auto_nsv_when_estimate_is_zero():
    counts = pandas.DataFrame(
        {
            'RUN1': [10.0, 5.0, 0.0],
            'RUN2': [10.0, 5.0, 0.0],
            'RUN3': [10.0, 5.0, 0.0],
            'RUN4': [10.0, 5.0, 0.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    metadata = pandas.DataFrame({
        'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
        'sample_group': ['A', 'A', 'B', 'B'],
    })
    corrected, sv_df, summary = run_sva_backend(
        counts_df=counts,
        metadata_df=metadata,
        nsv_setting='auto',
        B_setting='5',
        B_auto_max=100,
        random_seed=7,
    )
    pandas.testing.assert_frame_equal(corrected, counts)
    assert list(sv_df.index) == ['RUN1', 'RUN2', 'RUN3', 'RUN4']
    assert sv_df.shape == (4, 0)
    assert summary['resolved_sva_nsv'] == 0
    assert summary['resolved_sva_B'] == 5
    assert summary['skip_reason'] == 'sva_nsv_zero'


def test_estimate_num_sv_leek_returns_bounded_integer_estimate():
    data = numpy.array([
        [10.0, 12.0, 40.0, 42.0],
        [11.0, 13.0, 41.0, 43.0],
        [12.0, 14.0, 42.0, 44.0],
        [13.0, 15.0, 43.0, 45.0],
        [50.0, 48.0, 20.0, 18.0],
        [51.0, 49.0, 19.0, 17.0],
        [52.0, 50.0, 18.0, 16.0],
        [53.0, 51.0, 17.0, 15.0],
        [25.0, 26.0, 30.0, 31.0],
        [26.0, 27.0, 31.0, 32.0],
    ])
    mod, _ = build_sample_group_design_matrix(['A', 'A', 'B', 'B'])
    estimate = estimate_num_sv_leek(
        data_matrix=data,
        mod_matrix=mod,
        max_nsv=1,
    )
    assert estimate.method == 'leek'
    assert estimate.nsv in (0, 1)


def test_estimate_num_sv_at_B_falls_back_to_leek_when_be_fails(monkeypatch):
    data = numpy.array([
        [10.0, 12.0, 40.0, 42.0],
        [11.0, 13.0, 41.0, 43.0],
        [12.0, 14.0, 42.0, 44.0],
        [13.0, 15.0, 43.0, 45.0],
        [50.0, 48.0, 20.0, 18.0],
        [51.0, 49.0, 19.0, 17.0],
        [52.0, 50.0, 18.0, 16.0],
        [53.0, 51.0, 17.0, 15.0],
        [25.0, 26.0, 30.0, 31.0],
        [26.0, 27.0, 31.0, 32.0],
    ])
    mod, _ = build_sample_group_design_matrix(['A', 'A', 'B', 'B'])
    monkeypatch.setattr(
        batch_effect_sva_module,
        'estimate_num_sv_be',
        lambda *args, **kwargs: (_ for _ in ()).throw(ValueError('forced failure')),
    )
    estimate = estimate_num_sv_at_B(
        data_matrix=data,
        mod_matrix=mod,
        B_value=5,
        max_nsv=1,
        random_seed=7,
    )
    assert estimate.method == 'leek'


def test_f_pvalue_returns_valid_probabilities():
    data = numpy.array([
        [10.0, 11.0, 20.0, 21.0],
        [1.0, 2.0, 1.0, 2.0],
    ])
    mod, _ = build_sample_group_design_matrix(['A', 'A', 'B', 'B'])
    mod0, _ = build_intercept_only_design_matrix(4)
    p_values = f_pvalue(data, mod, mod0)
    assert p_values.shape == (2,)
    assert numpy.all((p_values >= 0.0) & (p_values <= 1.0))


def test_irwsva_build_returns_expected_shapes():
    data = numpy.array([
        [10.0, 12.0, 40.0, 42.0],
        [11.0, 13.0, 41.0, 43.0],
        [50.0, 48.0, 20.0, 18.0],
        [51.0, 49.0, 19.0, 17.0],
    ])
    mod, _ = build_sample_group_design_matrix(['A', 'A', 'B', 'B'])
    mod0, _ = build_intercept_only_design_matrix(4)
    out = irwsva_build(
        data_matrix=data,
        mod_matrix=mod,
        mod0_matrix=mod0,
        nsv=1,
        B_iterations=2,
    )
    assert out['sv'].shape == (4, 1)
    assert out['pprob_gam'].shape == (4,)
    assert out['pprob_b'].shape == (4,)
    assert out['n_svs'] == 1


def test_run_sva_backend_supports_positive_manual_nsv():
    counts = pandas.DataFrame(
        {
            'RUN1': [10.0, 11.0, 50.0, 51.0],
            'RUN2': [12.0, 13.0, 48.0, 49.0],
            'RUN3': [40.0, 41.0, 20.0, 19.0],
            'RUN4': [42.0, 43.0, 18.0, 17.0],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    metadata = pandas.DataFrame({
        'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
        'sample_group': ['A', 'A', 'B', 'B'],
    })
    corrected, sv_df, summary = run_sva_backend(
        counts_df=counts,
        metadata_df=metadata,
        nsv_setting='1',
        B_setting='5',
        B_auto_max=100,
        random_seed=7,
    )
    assert corrected.shape == counts.shape
    assert sv_df.shape == (4, 1)
    assert summary['resolved_sva_nsv'] == 1
    assert summary['resolved_sva_B'] == 5
    assert summary['skip_reason'] == ''
    assert summary['corrected_run_ids'] == ['RUN1', 'RUN2', 'RUN3', 'RUN4']


def test_run_sva_backend_supports_positive_manual_nsv_on_transformed_duplicate_groups():
    runs = ['RUN01', 'RUN02', 'RUN03', 'RUN04']
    genes = [f'G{i:03d}' for i in range(1, 51)]
    counts = pandas.DataFrame(index=genes)
    for run_id, sample_group in zip(runs, ['A', 'A', 'B', 'B']):
        if sample_group == 'A':
            values = [100 + i for i in range(25)] + [5 + i for i in range(25)]
        else:
            values = [5 + i for i in range(25)] + [100 + i for i in range(25)]
        counts[run_id] = values
    transformed = numpy.log2((counts.div(1000.0).div(counts.sum(axis=0), axis=1) * 1e9) + 1.0)
    metadata = pandas.DataFrame({
        'run': runs,
        'sample_group': ['A', 'A', 'B', 'B'],
    })
    corrected, sv_df, summary = run_sva_backend(
        counts_df=transformed,
        metadata_df=metadata,
        nsv_setting='1',
        B_setting='5',
        B_auto_max=80,
        random_seed=7,
    )
    assert corrected.shape == transformed.shape
    assert sv_df.shape == (4, 1)
    assert summary['resolved_sva_nsv'] == 1
    assert summary['resolved_sva_B'] == 5
    numpy.testing.assert_allclose(corrected.to_numpy(), transformed.to_numpy(), rtol=0.0, atol=1e-8)


def test_run_sva_backend_auto_transformed_duplicate_groups_falls_back_to_leek():
    runs = ['RUN01', 'RUN02', 'RUN03', 'RUN04']
    genes = [f'G{i:03d}' for i in range(1, 51)]
    counts = pandas.DataFrame(index=genes)
    for run_id, sample_group in zip(runs, ['A', 'A', 'B', 'B']):
        if sample_group == 'A':
            values = [100 + i for i in range(25)] + [5 + i for i in range(25)]
        else:
            values = [5 + i for i in range(25)] + [100 + i for i in range(25)]
        counts[run_id] = values
    transformed = numpy.log2((counts.div(1000.0).div(counts.sum(axis=0), axis=1) * 1e9) + 1.0)
    metadata = pandas.DataFrame({
        'run': runs,
        'sample_group': ['A', 'A', 'B', 'B'],
    })
    corrected, sv_df, summary = run_sva_backend(
        counts_df=transformed,
        metadata_df=metadata,
        nsv_setting='auto',
        B_setting='auto',
        B_auto_max=80,
        random_seed=7,
    )
    assert corrected.shape == transformed.shape
    assert sv_df.shape == (4, 0)
    assert summary['resolved_sva_nsv'] == 0
    assert summary['resolved_sva_B'] == 45
    assert summary['sva_estimation_method'] == 'auto_leek'
    assert summary['sva_stable'] is True
    assert summary['trace_B'] == [30, 45]
    assert summary['trace_nsv'] == [0, 0]
    assert summary['trace_method'] == ['leek', 'leek']
    numpy.testing.assert_allclose(corrected.to_numpy(), transformed.to_numpy(), rtol=0.0, atol=1e-8)
