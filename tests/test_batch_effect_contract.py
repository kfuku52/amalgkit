import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from amalgkit.batch_effect_combatseq import run_combatseq_backend
from amalgkit.batch_effect_contract import BatchModelError, build_biological_design
from amalgkit.batch_effect_latent_glm import run_latent_glm_backend
from amalgkit.batch_effect_ruvseq import run_ruvseq_backend, ruvr_correct_counts
from amalgkit.batch_effect_sva import _stabilize_surrogate_matrix, irwsva_build, run_sva_backend
from amalgkit.batch_effect_common import write_batch_effect_summary_tsv
from amalgkit.batch_effect_io import read_backend_summary_dcf, write_backend_summary_dcf


def fixture():
    rng = np.random.default_rng(94)
    counts = pd.DataFrame(rng.poisson(100, (200, 12)),
                          index=[f'g{i}' for i in range(200)], columns=[f'r{i}' for i in range(12)])
    metadata = pd.DataFrame({
        'run': counts.columns, 'sample_group': ['A'] * 6 + ['B'] * 6,
        'bioproject': ['P1', 'P2'] * 6, 'sex': ['M', 'M', 'F', 'F', 'M', 'F'] * 2,
        'age': np.arange(12),
    })
    return counts, metadata


@pytest.mark.parametrize('backend,options', [
    (run_combatseq_backend, {}), (run_sva_backend, {'nsv_setting': 0}),
    (run_latent_glm_backend, {'k_setting': 0}), (run_ruvseq_backend, {'k_setting': 0}),
])
def test_invalid_inputs_are_not_silently_skipped(backend, options):
    counts, metadata = fixture()
    missing = metadata.copy()
    missing.loc[0, 'sample_group'] = None
    with pytest.raises(ValueError, match='missing values'):
        backend(counts, missing, **options)
    with pytest.raises(ValueError, match='duplicate run'):
        backend(counts, pd.concat([metadata, metadata.iloc[:1]]), **options)
    damaged = counts.astype(float)
    damaged.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match='finite'):
        backend(damaged, metadata, **options)


def test_combat_failure_skips_whole_species_without_retry(monkeypatch):
    counts, metadata = fixture()
    calls = []

    def fail(**kwargs):
        calls.append(kwargs)
        raise ValueError('numerical fit failure')

    monkeypatch.setattr('amalgkit.batch_effect_combatseq._load_pycombat_seq', lambda: fail)
    corrected, summary = run_combatseq_backend(counts, metadata)
    pd.testing.assert_frame_equal(corrected, counts)
    assert len(calls) == 1 and 'covar_mod' in calls[0]
    assert summary['batch_failure_policy'] == 'skip'
    assert summary['status'] == 'skipped'
    assert summary['corrected_run_ids'] == []
    assert summary['group_fallback_used'] is False
    with pytest.raises(BatchModelError, match='numerical fit failure'):
        run_combatseq_backend(counts, metadata, failure_policy='error')


@pytest.mark.parametrize('samples', [1, 12])
def test_combat_fractional_input_is_error_even_when_model_cannot_fit(samples):
    counts, metadata = fixture()
    counts = counts.iloc[:, :samples].astype(float) + 0.5
    metadata = metadata.iloc[:samples].copy()
    metadata['bioproject'] = metadata['sample_group']
    with pytest.raises(ValueError, match='integer raw counts'):
        run_combatseq_backend(counts, metadata)


@pytest.mark.parametrize('problem', ['confounded', 'singleton'])
def test_combat_checks_actual_design_before_loading_backend(monkeypatch, problem):
    counts, metadata = fixture()
    if problem == 'confounded':
        metadata['bioproject'] = metadata['sample_group']
    else:
        metadata.loc[0, 'bioproject'] = 'singleton'

    def unexpected():
        pytest.fail('ineligible data reached ComBat backend')

    monkeypatch.setattr('amalgkit.batch_effect_combatseq._load_pycombat_seq', unexpected)
    corrected, summary = run_combatseq_backend(counts, metadata)
    pd.testing.assert_frame_equal(corrected, counts)
    assert summary['status'] == 'skipped'
    assert problem in summary['skip_reason']


def test_combat_batch_only_requires_explicit_choice(monkeypatch):
    counts, metadata = fixture()
    metadata['bioproject'] = metadata['sample_group']
    calls = []

    def backend(**kwargs):
        calls.append(kwargs)
        return kwargs['counts'] + 1

    monkeypatch.setattr('amalgkit.batch_effect_combatseq._load_pycombat_seq', lambda: backend)
    corrected, summary = run_combatseq_backend(counts, metadata, protect_group=False)
    assert 'covar_mod' not in calls[0]
    pd.testing.assert_frame_equal(corrected, counts + 1)
    assert summary['design']['group_protected'] is False


def test_unexpected_backend_error_is_not_recoverable(monkeypatch):
    counts, metadata = fixture()

    def backend(**kwargs):
        raise RuntimeError('implementation bug')

    monkeypatch.setattr('amalgkit.batch_effect_combatseq._load_pycombat_seq', lambda: backend)
    with pytest.raises(RuntimeError, match='implementation bug'):
        run_combatseq_backend(counts, metadata)


@pytest.mark.parametrize('backend', [run_latent_glm_backend, run_ruvseq_backend])
def test_uncalibrated_auto_is_explicitly_skipped(backend):
    counts, metadata = fixture()
    corrected, factors, summary = backend(counts, metadata)
    pd.testing.assert_frame_equal(corrected, counts)
    assert factors.empty
    assert summary['status'] == 'skipped'
    assert summary['skip_reason'].endswith('auto_not_calibrated')
    with pytest.raises(BatchModelError):
        backend(counts, metadata, failure_policy='error')


def test_sva_never_invents_directions_for_zero_or_collinear_candidates():
    design = np.column_stack([np.ones(6), [0, 0, 0, 1, 1, 1]])
    for candidates in (np.ones((6, 1)), design, np.zeros((6, 2))):
        assert _stabilize_surrogate_matrix(candidates, design).shape == (6, 0)
    fit = irwsva_build(np.zeros((20, 6)), design, nsv=2)
    assert fit['sv'].shape == (6, 0)
    assert fit['n_svs'] == 0
    assert fit['irw_iterations_completed'] == 0
    for scale in (1e-6, 1.0, 1e6):
        fitted_only = np.ones((20, 6)) * scale + np.arange(20)[:, None] * design[:, 1] * scale
        assert irwsva_build(fitted_only, design, nsv=2)['n_svs'] == 0


def test_latent_preserves_all_declared_biological_covariates_on_model_scale():
    counts, metadata = fixture()
    corrected, _, summary = run_latent_glm_backend(
        counts, metadata, k_setting=1, family='uniform',
        categorical_covariates=['sex'], continuous_covariates=['age'], failure_policy='error',
    )
    size_factors = counts.sum(axis=0) / counts.sum(axis=0).median()
    before = np.log(counts.div(size_factors, axis=1) + 0.5)
    after = np.log(corrected.div(size_factors, axis=1) + 0.5)
    design = np.array(summary['design']['design_matrix'])
    np.testing.assert_allclose(design.T @ (after - before).to_numpy().T, 0, atol=1e-10)
    assert summary['latent_objective_kind'] == 'weighted_log_residual_mse'
    assert summary['design']['design_encoding']['age']['type'] == 'continuous'
    shuffled, _, _ = run_latent_glm_backend(
        counts, metadata.sample(frac=1, random_state=4), k_setting=1, family='uniform',
        categorical_covariates=['sex'], continuous_covariates=['age'], failure_policy='error',
    )
    pd.testing.assert_frame_equal(corrected, shuffled)


def test_ruv_removes_only_design_orthogonal_component_but_returns_original_w():
    rng = np.random.default_rng(2)
    design = np.column_stack([np.ones(8), [0] * 4 + [1] * 4])
    nuisance = design[:, 1] + np.tile([-0.5, 0.5], 4)
    log_expression = pd.DataFrame(rng.normal(size=(30, 8)) + 3 * design[:, 1])
    residuals = pd.DataFrame(np.outer(rng.normal(size=30), nuisance))
    corrected, factors = ruvr_correct_counts(
        log_expression, np.ones(30, dtype=bool), 1, residuals,
        is_log=True, design_matrix=design,
    )
    assert np.linalg.norm(design.T @ factors.to_numpy()) > 0.1
    np.testing.assert_allclose(design.T @ (corrected - log_expression).to_numpy().T, 0, atol=1e-12)
    assert np.linalg.norm(corrected - log_expression) > 1


def test_ruv_controls_shortage_does_not_substitute_all_genes():
    counts, metadata = fixture()
    corrected, _, summary = run_ruvseq_backend(counts.iloc[:10], metadata, k_setting=1, min_controls=100)
    pd.testing.assert_frame_equal(corrected, counts.iloc[:10])
    assert summary['skip_reason'] == 'ruvseq_insufficient_controls'


def test_ruv_k_zero_bypasses_fitting(monkeypatch):
    counts, metadata = fixture()

    def unexpected(**kwargs):
        pytest.fail('k=0 must not run GLM or control selection')

    monkeypatch.setattr('amalgkit.batch_effect_ruvseq._compute_glm_pvalues_and_residuals', unexpected)
    corrected, factors, summary = run_ruvseq_backend(counts, metadata, k_setting=0)
    pd.testing.assert_frame_equal(corrected, counts)
    assert factors.empty
    assert summary['status'] == 'not_needed'


def test_ruv_nb_failure_does_not_revert_to_poisson(monkeypatch):
    import statsmodels.api as sm
    counts, metadata = fixture()
    original = sm.GLM
    families = []

    def glm(*args, **kwargs):
        families.append(type(kwargs['family']).__name__)
        if isinstance(kwargs['family'], sm.families.NegativeBinomial):
            raise ValueError('NB fit failed')
        return original(*args, **kwargs)

    monkeypatch.setattr(sm, 'GLM', glm)
    monkeypatch.setattr('amalgkit.batch_effect_ruvseq._estimate_nb_alpha_from_poisson_fit', lambda **kwargs: 0.1)
    corrected, factors, summary = run_ruvseq_backend(counts, metadata, k_setting=1, control_mode='all')
    pd.testing.assert_frame_equal(corrected, counts)
    assert factors.empty
    assert families == ['Poisson', 'NegativeBinomial']
    assert summary['skip_reason'] == 'ruvseq_nb_glm_failed'
    assert summary['ruv_fallback_used'] is not True


def test_explicit_control_ids_and_cli_defaults(tmp_path):
    from amalgkit.batch_effect_contract import read_control_gene_ids
    from amalgkit.batch_effect_runner import build_parser
    from amalgkit.main import build_main_parser
    public = build_main_parser().parse_args(['finalize'])
    internal = build_parser().parse_args(['--backend', 'sva', '--counts_tsv', 'c.tsv', '--metadata_tsv', 'm.tsv'])
    assert public.batch_failure_policy == internal.batch_failure_policy == 'skip'
    assert public.sva_irw_iterations == internal.sva_irw_iterations == 5
    path = tmp_path / 'controls.txt'
    path.write_text('g0\ng1\ng2\n')
    counts, metadata = fixture()
    _, _, summary = run_ruvseq_backend(counts, metadata, k_setting=1, control_mode='file',
                                       control_gene_ids=read_control_gene_ids(path), min_controls=2)
    assert summary['ruv_control_gene_ids'] == ['g0', 'g1', 'g2']
    assert summary['ruv_control_mode'] == 'file'


def test_diagnostics_survive_json_dcf_and_factor_exports(tmp_path):
    counts, metadata = fixture()
    _, factors, summary = run_latent_glm_backend(counts, metadata, k_setting=1)
    paths = write_batch_effect_summary_tsv(summary, 'Example species', 'Example_species', tmp_path)
    payload = json.loads((tmp_path / 'Example_species.latent_loglinear.batch_effect_diagnostics.json').read_text())
    assert payload['design']['design_run_ids'] == counts.columns.tolist()
    assert payload['factor_columns'] == factors.columns.tolist()
    assert paths['diagnostics_path'].endswith('.json')
    exported = pd.read_csv(tmp_path / 'Example_species.latent_loglinear.batch_effect_factors.tsv', sep='\t', index_col=0)
    np.testing.assert_allclose(exported, factors)
    dcf = tmp_path / 'summary.dcf'
    write_backend_summary_dcf(summary, dcf)
    decoded = read_backend_summary_dcf(dcf)
    assert decoded['design'] == summary['design']
    assert decoded['postprocessing'] == summary['postprocessing']


def test_design_rank_and_complete_batch_confounding_are_distinguished():
    _, metadata = fixture()
    metadata['bioproject'] = metadata['sample_group']
    design = build_biological_design(metadata)
    assert design.diagnostics['design_rank'] == 2
    assert design.diagnostics['batch_design_confounded'] is True
    assert design.diagnostics['batch_design_rank'] == 2


def test_finalizer_skip_preserves_requested_normalization():
    from amalgkit.per_species_finalize_python import _run_batch_effect_step
    counts, metadata = fixture()
    lengths = counts.astype(float) * 0 + 1000
    args = SimpleNamespace(norm='log2p1-none', batch_effect_alg='latent_loglinear')
    result = _run_batch_effect_step(counts, metadata, lengths, args)
    pd.testing.assert_frame_equal(result['tc'], np.log2(counts.loc[:, result['tc'].columns] + 1))
    assert result['batch_info']['status'] == 'skipped'
    assert result['batch_info']['gene_ids_fitted'] == []
    assert set(result['batch_info']['gene_ids_not_fitted']) == set(counts.index)


def test_manual_latent_dimension_failure_is_not_reported_as_requested_zero():
    counts, metadata = fixture()
    counts.iloc[:, :] = 100
    corrected, factors, summary = run_latent_glm_backend(counts, metadata, k_setting=1)
    pd.testing.assert_frame_equal(corrected, counts)
    assert factors.empty
    assert summary['status'] == 'skipped'
    assert summary['skip_reason'] == 'latent_degenerate_factors'
    with pytest.raises(BatchModelError, match='dimensions could not be estimated'):
        run_latent_glm_backend(counts, metadata, k_setting=1, failure_policy='error')


def test_public_copy_and_sanity_preserve_diagnostic_contract(tmp_path):
    from amalgkit.finalize import _copy_species_tables
    from amalgkit.sanity import _validate_finalize_batch_contract
    counts, metadata = fixture()
    _, _, summary = run_latent_glm_backend(counts, metadata, k_setting=1)
    source = tmp_path / 'per_species'
    tables = source / 'Example_species' / 'tables'
    write_batch_effect_summary_tsv(summary, 'Example species', 'Example_species', tables)
    destination = tmp_path / 'finalize'
    _copy_species_tables(source, destination, 'latent_loglinear')
    public = destination / 'Example_species'
    assert _validate_finalize_batch_contract(public, 'Example_species') == ''
    (public / 'Example_species_batch_effect_factors.tsv').unlink()
    assert 'missing' in _validate_finalize_batch_contract(public, 'Example_species')
