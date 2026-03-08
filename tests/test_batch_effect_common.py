import pandas

from amalgkit.batch_effect_common import (
    BatchEffectResult,
    annotate_metadata_with_batch_info,
    build_batch_effect_summary_dataframe,
    initialize_batch_info,
    normalize_run_ids,
    write_batch_effect_summary_tsv,
)


def test_normalize_run_ids_trims_deduplicates_and_skips_empty():
    assert normalize_run_ids([' RUN1 ', 'RUN1', '', None, 'RUN2']) == ['RUN1', 'RUN2']


def test_batch_effect_result_to_jsonable_normalizes_runs_and_merges_extra():
    result = BatchEffectResult(
        backend='sva',
        method='placeholder',
        corrected_run_ids=[' RUN1 ', 'RUN1'],
        uncorrected_run_ids=['RUN2'],
        resolved_sva_nsv=2,
        extra={'custom_flag': 'yes'},
    )
    payload = result.to_jsonable()
    assert payload['backend'] == 'sva'
    assert payload['corrected_run_ids'] == ['RUN1']
    assert payload['uncorrected_run_ids'] == ['RUN2']
    assert payload['resolved_sva_nsv'] == 2
    assert payload['custom_flag'] == 'yes'


def test_initialize_batch_info_matches_finalize_defaults():
    observed = initialize_batch_info(run_ids=[' RUN1 ', 'RUN1', 'RUN2'], batch_effect_alg='sva')
    assert observed['batch_effect_alg_requested'] == 'sva'
    assert observed['batch_effect_alg_applied'] == 'sva'
    assert observed['corrected_runs'] == []
    assert observed['uncorrected_runs'] == ['RUN1', 'RUN2']
    assert observed['skip_reason'] == 'not_run'
    assert observed['resolved_sva_nsv'] is None
    assert observed['resolved_ruv_k'] is None
    assert observed['resolved_latent_k'] is None


def test_annotate_metadata_with_batch_info_marks_corrected_runs():
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3'],
            'sample_group': ['A', 'A', 'B'],
        }
    )
    batch_info = initialize_batch_info(run_ids=metadata_df['run'], batch_effect_alg='sva')
    batch_info['batch_effect_alg_applied'] = 'sva'
    batch_info['corrected_runs'] = ['RUN1', 'RUN3']
    batch_info['uncorrected_runs'] = ['RUN2']

    observed = annotate_metadata_with_batch_info(metadata_df, batch_info)
    assert observed['batch_corrected'].tolist() == ['yes', 'no', 'yes']
    assert observed['batch_alg_used'].tolist() == ['sva', 'no', 'sva']


def test_annotate_metadata_with_batch_info_supports_batch_effect_result_shape():
    metadata_df = pandas.DataFrame({'run': ['RUN1', 'RUN2']})
    batch_result = BatchEffectResult(
        backend='combatseq',
        method='combatseq',
        corrected_run_ids=['RUN1'],
        uncorrected_run_ids=['RUN2'],
    )

    observed = annotate_metadata_with_batch_info(metadata_df, batch_result)
    assert observed['batch_corrected'].tolist() == ['yes', 'no']
    assert observed['batch_alg_used'].tolist() == ['combatseq', 'no']


def test_build_batch_effect_summary_dataframe_matches_expected_columns():
    batch_info = initialize_batch_info(run_ids=['RUN1', 'RUN2', 'RUN3'], batch_effect_alg='ruvseq')
    batch_info.update(
        {
            'batch_effect_alg_applied': 'ruvseq',
            'corrected_runs': ['RUN1', 'RUN2'],
            'uncorrected_runs': ['RUN3'],
            'resolved_ruv_k': 2,
            'resolved_ruv_controls': 150,
            'ruv_baseline_score': 1.25,
            'ruv_selected_score': 0.5,
            'ruv_selected_penalized_score': 0.6,
            'ruv_penalty': 0.1,
            'skip_reason': '',
        }
    )

    observed = build_batch_effect_summary_dataframe(
        batch_info=batch_info,
        scientific_name='Oryza sativa',
        random_seed_value=123,
    )
    assert observed.loc[0, 'scientific_name'] == 'Oryza sativa'
    assert observed.loc[0, 'batch_effect_alg_requested'] == 'ruvseq'
    assert observed.loc[0, 'batch_effect_alg_applied'] == 'ruvseq'
    assert observed.loc[0, 'random_seed'] == '123'
    assert int(observed.loc[0, 'corrected_run_count']) == 2
    assert observed.loc[0, 'corrected_runs'] == 'RUN1|RUN2'
    assert int(observed.loc[0, 'uncorrected_run_count']) == 1
    assert observed.loc[0, 'uncorrected_runs'] == 'RUN3'
    assert int(observed.loc[0, 'resolved_ruv_k']) == 2
    assert int(observed.loc[0, 'resolved_ruv_controls']) == 150


def test_build_batch_effect_summary_dataframe_accepts_batch_effect_result():
    batch_result = BatchEffectResult(
        backend='sva',
        method='sva',
        corrected_run_ids=['RUN1', 'RUN2'],
        uncorrected_run_ids=['RUN3'],
        resolved_sva_nsv=1,
        resolved_sva_B=10,
        stable=True,
        skip_reason='',
        extra={'sva_estimation_method': 'be'},
    )
    observed = build_batch_effect_summary_dataframe(
        batch_info=batch_result,
        scientific_name='Zea mays',
        random_seed_value=None,
    )
    assert observed.loc[0, 'batch_effect_alg_requested'] == 'sva'
    assert observed.loc[0, 'batch_effect_alg_applied'] == 'sva'
    assert observed.loc[0, 'random_seed'] == 'auto'
    assert int(observed.loc[0, 'resolved_sva_nsv']) == 1
    assert int(observed.loc[0, 'resolved_sva_B']) == 10
    assert observed.loc[0, 'sva_estimation_method'] == 'be'
    assert bool(observed.loc[0, 'sva_stable']) is True


def test_write_batch_effect_summary_tsv_writes_expected_file(tmp_path):
    batch_info = initialize_batch_info(run_ids=['RUN1', 'RUN2'], batch_effect_alg='combatseq')
    batch_info.update(
        {
            'batch_effect_alg_applied': 'combatseq',
            'corrected_runs': ['RUN1'],
            'uncorrected_runs': ['RUN2'],
            'skip_reason': 'combatseq_singleton_kept',
        }
    )
    out = write_batch_effect_summary_tsv(
        batch_info=batch_info,
        scientific_name='Glycine max',
        species_tag='Glycine_max',
        dir_tsv=str(tmp_path),
        random_seed_value=None,
    )
    out_path = tmp_path / 'Glycine_max.combatseq.batch_effect_summary.tsv'
    assert out['summary_path'] == str(out_path)
    assert out_path.exists()
    loaded = pandas.read_csv(out_path, sep='\t')
    pandas.testing.assert_frame_equal(loaded, out['summary_df'], check_dtype=False)


def test_build_batch_effect_summary_dataframe_supports_latent_glm_fields():
    batch_info = initialize_batch_info(run_ids=['RUN1', 'RUN2'], batch_effect_alg='latent_glm')
    batch_info.update(
        {
            'batch_effect_alg_applied': 'latent_glm',
            'corrected_runs': ['RUN1', 'RUN2'],
            'uncorrected_runs': [],
            'resolved_latent_k': 1,
            'latent_family': 'nb',
            'latent_iterations': 7,
            'latent_objective': 0.125,
            'latent_converged': True,
            'skip_reason': '',
        }
    )
    observed = build_batch_effect_summary_dataframe(
        batch_info=batch_info,
        scientific_name='Arabidopsis thaliana',
        random_seed_value=99,
    )
    assert observed.loc[0, 'batch_effect_alg_requested'] == 'latent_glm'
    assert observed.loc[0, 'batch_effect_alg_applied'] == 'latent_glm'
    assert int(observed.loc[0, 'resolved_latent_k']) == 1
    assert observed.loc[0, 'latent_family'] == 'nb'
    assert int(observed.loc[0, 'latent_iterations']) == 7
    assert float(observed.loc[0, 'latent_objective']) == 0.125
    assert bool(observed.loc[0, 'latent_converged']) is True
