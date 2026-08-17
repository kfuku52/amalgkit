import os
import sys
import types
from types import SimpleNamespace

import numpy
import pandas
import pytest

from amalgkit.command_context import PerSpeciesTableContext
from amalgkit.metadata_utils import Metadata
from amalgkit.per_species_finalize_python import (
    _apply_transformation_logic,
    _compute_corr_matrix,
    _compute_pca_coordinates,
    _compute_tsne_coordinates,
    _compute_distance_matrix,
    _resolve_scientific_name,
    _run_batch_effect_step,
    _transform_raw_to_fpkm,
    _transform_raw_to_tpm,
    load_quant_model_table,
    resolve_length_models,
)
from amalgkit import per_species_tables as per_species_tables_module
from tests.support.per_species import build_per_species_args


@pytest.mark.parametrize(('seed', 'expected'), [(37, 37), ('auto', None)])
def test_compute_tsne_coordinates_honors_finalize_seed(monkeypatch, seed, expected):
    observed = {}

    class FakeTSNE:
        def __init__(self, **kwargs):
            observed.update(kwargs)

        def fit_transform(self, values):
            return numpy.zeros((values.shape[0], 2), dtype=float)

    monkeypatch.setitem(sys.modules, 'sklearn.manifold', types.SimpleNamespace(TSNE=FakeTSNE))
    counts_df = pandas.DataFrame(
        numpy.arange(20, dtype=float).reshape(5, 4),
        columns=['RUN1', 'RUN2', 'RUN3', 'RUN4'],
    )

    coords = _compute_tsne_coordinates(counts_df, random_seed=seed)

    assert coords.shape == (4, 2)
    assert observed['random_state'] == expected


def test_resolve_scientific_name_honors_explicit_species_token():
    # Regression for #174: the tag->name reverse map must use the same token
    # resolution as merge/cstmm (explicit species_token honored), so a worker
    # dispatched with token "human" resolves back to "Homo sapiens".
    metadata = pandas.DataFrame({
        'run': ['R1', 'R2'],
        'scientific_name': ['Homo sapiens', 'Mus musculus'],
        'sample_group': ['g1', 'g2'],
        'exclusion': ['no', 'no'],
        'species_token': ['human', ''],
    })
    assert _resolve_scientific_name(metadata, 'human') == 'Homo sapiens'
    assert _resolve_scientific_name(metadata, 'Mus_musculus') == 'Mus musculus'
    # Fallback for a naive space->underscore tag remains available.
    assert _resolve_scientific_name(metadata, 'Homo_sapiens') == 'Homo sapiens'


def _write_species_input_fixture(tmp_path, species='Finalizus example', sample_groups=None, bioprojects=None):
    if sample_groups is None:
        sample_groups = ['A', 'A', 'B', 'B']
    if bioprojects is None:
        bioprojects = ['BP1', 'BP1', 'BP2', 'BP2']
    species_tag = species.replace(' ', '_')
    runs = ['RUN{:02d}'.format(i + 1) for i in range(len(sample_groups))]
    genes = ['G{:03d}'.format(i + 1) for i in range(20)]

    input_dir = tmp_path / 'input'
    species_dir = input_dir / species_tag
    species_dir.mkdir(parents=True, exist_ok=True)

    metadata_df = pandas.DataFrame(
        {
            'run': runs,
            'scientific_name': [species] * len(runs),
            'sample_group': sample_groups,
            'bioproject': bioprojects,
            'exclusion': ['no'] * len(runs),
            'mapping_rate': [100.0] * len(runs),
            'instrument': ['Illumina'] * len(runs),
            'lib_layout': ['PAIRED'] * len(runs),
            'lib_selection': ['cDNA'] * len(runs),
            'total_spots': [1_000_000] * len(runs),
            'total_bases': [150_000_000] * len(runs),
        }
    )
    metadata = Metadata.from_DataFrame(metadata_df)

    counts_df = pandas.DataFrame(index=genes)
    for run_id, sample_group in zip(runs, sample_groups):
        if sample_group == 'A':
            values = [100 + i for i in range(10)] + [2 + i for i in range(10)]
        else:
            values = [2 + i for i in range(10)] + [100 + i for i in range(10)]
        counts_df[run_id] = values
    eff_length_df = pandas.DataFrame(1000, index=genes, columns=runs)

    counts_df.reset_index(names='target_id').to_csv(species_dir / '{}_est_counts.tsv'.format(species_tag), sep='\t', index=False)
    eff_length_df.reset_index(names='target_id').to_csv(species_dir / '{}_eff_length.tsv'.format(species_tag), sep='\t', index=False)
    return {
        'metadata': metadata,
        'input_dir': str(input_dir),
        'species': species,
        'species_tag': species_tag,
    }


def _inject_latent_batch_signal(input_dir, species_tag):
    counts_path = os.path.join(input_dir, species_tag, '{}_est_counts.tsv'.format(species_tag))
    counts_df = pandas.read_csv(counts_path, sep='\t')
    run_columns = [column for column in counts_df.columns if column != 'target_id']
    values = {
        'RUN01': [120, 118, 122, 119, 12, 14, 15, 13, 80, 82, 78, 81, 10, 11, 12, 13, 14, 15, 16, 17],
        'RUN02': [121, 117, 123, 120, 13, 15, 14, 12, 32, 30, 34, 31, 10, 11, 12, 13, 14, 15, 16, 17],
        'RUN03': [14, 13, 12, 11, 131, 129, 133, 130, 79, 83, 77, 80, 18, 17, 16, 15, 14, 13, 12, 11],
        'RUN04': [13, 12, 11, 10, 132, 130, 134, 131, 33, 31, 35, 32, 18, 17, 16, 15, 14, 13, 12, 11],
    }
    for run_id in run_columns:
        counts_df.loc[:, run_id] = values[run_id]
    counts_df.to_csv(counts_path, sep='\t', index=False)


def _build_args(tmp_path):
    return build_per_species_args(
        tmp_path,
        internal_jobs='auto',
        mapping_rate=0.20,
        correlation_threshold=0.30,
        worker_mode='finalize',
    )


def test_compute_distance_matrix_clips_negative_roundoff():
    corr_df = pandas.DataFrame(
        [
            [1.0, 1.0000000001],
            [1.0000000001, 1.0],
        ],
        columns=['A', 'B'],
        index=['A', 'B'],
    )

    dist = _compute_distance_matrix(corr_df)

    assert numpy.all(dist >= 0.0)
    assert dist[0, 1] == 0.0
    assert dist[1, 0] == 0.0


def test_batch_step_aligns_effective_lengths_after_gene_partitioning():
    counts = pandas.DataFrame(
        {
            'RUN1': [100.0, 0.0, 100.0],
            'RUN2': [100.0, 0.0, 100.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    effective_lengths = pandas.DataFrame(
        {
            'RUN1': [100.0, 200.0, 1000.0],
            'RUN2': [100.0, 200.0, 1000.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2'],
            'sample_group': ['A', 'B'],
            'scientific_name': ['Example species', 'Example species'],
            'bioproject': ['BP1', 'BP2'],
            'exclusion': ['no', 'no'],
        }
    )
    args = SimpleNamespace(
        norm='none-fpkm',
        batch_effect_alg='combatseq',
        clip_negative=True,
    )

    result = _run_batch_effect_step(
        counts_df=counts,
        metadata_df=metadata,
        eff_length_df=effective_lengths,
        args=args,
    )

    assert list(result['tc'].index) == ['G1', 'G2', 'G3']
    assert result['tc'].loc['G3', 'RUN1'] == 500_000.0
    assert result['tc'].loc['G2', 'RUN1'] == 0.0


@pytest.mark.optional_dependency
def test_combatseq_batch_step_preserves_group_fallback_diagnostics():
    pytest.importorskip('inmoose.pycombat')
    counts = pandas.DataFrame(
        {
            'RUN1': [100, 101, 10, 11],
            'RUN2': [100, 101, 10, 11],
            'RUN3': [5, 6, 100, 101],
            'RUN4': [5, 6, 100, 101],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2'],
        }
    )
    effective_lengths = pandas.DataFrame(100.0, index=counts.index, columns=counts.columns)

    with pytest.warns(UserWarning, match='fell back to batch-only correction'):
        result = _run_batch_effect_step(
            counts_df=counts,
            metadata_df=metadata,
            eff_length_df=effective_lengths,
            args=SimpleNamespace(norm='none-none', batch_effect_alg='combatseq', clip_negative=True),
        )

    batch_info = result['batch_info']
    assert batch_info['group_model_used'] is False
    assert batch_info['group_fallback_used'] is True
    assert 'confounded with the batches' in batch_info['group_error_message']


def test_single_sample_batch_skip_still_applies_requested_normalization():
    counts = pandas.DataFrame({'RUN1': [100.0, 100.0]}, index=['G1', 'G2'])
    effective_lengths = pandas.DataFrame(
        {'RUN1': [100.0, 1000.0]},
        index=['G1', 'G2'],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1'],
            'sample_group': ['A'],
            'scientific_name': ['Example species'],
            'bioproject': ['BP1'],
            'exclusion': ['no'],
        }
    )

    for batch_effect_alg in ('combatseq', 'ruvseq', 'latent_glm'):
        result = _run_batch_effect_step(
            counts_df=counts,
            metadata_df=metadata,
            eff_length_df=effective_lengths,
            args=SimpleNamespace(
                norm='none-tpm',
                batch_effect_alg=batch_effect_alg,
                clip_negative=True,
            ),
        )
        numpy.testing.assert_allclose(
            result['tc'].loc[:, 'RUN1'].to_numpy(dtype=float),
            [909_090.9090909091, 90_909.09090909091],
        )
        assert result['batch_info']['skip_reason'] == 'single_sample'


def test_expression_transform_rejects_unknown_normalization_method():
    counts = pandas.DataFrame({'RUN1': [1.0]}, index=['G1'])
    effective_lengths = pandas.DataFrame({'RUN1': [100.0]}, index=['G1'])
    metadata = pandas.DataFrame({'run': ['RUN1']})

    with pytest.raises(ValueError, match='Unsupported expression normalization method'):
        _apply_transformation_logic(
            counts_df=counts,
            eff_length_df=effective_lengths,
            transform_method='banana',
            batch_effect_alg='no',
            step='before_batch',
            metadata_df=metadata,
        )


def test_ruvseq_all_nonexpressed_genes_returns_safe_noop():
    counts = pandas.DataFrame(
        {
            'RUN1': [0.0, 0.0],
            'RUN2': [0.0, 0.0],
        },
        index=['G1', 'G2'],
    )
    effective_lengths = pandas.DataFrame(
        {
            'RUN1': [100.0, 100.0],
            'RUN2': [100.0, 100.0],
        },
        index=['G1', 'G2'],
    )
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2'],
            'sample_group': ['A', 'B'],
            'scientific_name': ['Example species', 'Example species'],
            'bioproject': ['BP1', 'BP2'],
            'exclusion': ['no', 'no'],
        }
    )

    result = _run_batch_effect_step(
        counts_df=counts,
        metadata_df=metadata,
        eff_length_df=effective_lengths,
        args=SimpleNamespace(
            norm='none-tpm',
            batch_effect_alg='ruvseq',
            clip_negative=True,
            ruvseq_control_genes='auto',
            ruvseq_k='auto',
            ruvseq_k_max=2,
            ruvseq_control_top_n=1000,
            ruvseq_min_controls=2,
        ),
    )

    pandas.testing.assert_frame_equal(result['tc'], counts)
    assert result['batch_info']['skip_reason'] == 'no_expressed_genes'
    assert result['batch_info']['batch_effect_alg_applied'] == 'no'
    assert result['batch_info']['corrected_runs'] == []
    assert result['batch_info']['uncorrected_runs'] == ['RUN1', 'RUN2']


def test_expression_transform_allows_missing_or_zero_length_for_zero_counts():
    counts = pandas.DataFrame(
        {
            'RUN1': [0.0, 0.0, 5.0],
            'RUN2': [0.0, 2.0, 0.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    effective_lengths = pandas.DataFrame(
        {
            'RUN1': [0.0, numpy.nan, 10.0],
            'RUN2': [numpy.nan, 20.0, 0.0],
        },
        index=['G1', 'G2', 'G3'],
    )

    transformed = _transform_raw_to_tpm(
        counts_df=counts,
        eff_length_df=effective_lengths,
    )

    numpy.testing.assert_allclose(
        transformed.to_numpy(dtype=float),
        numpy.array(
            [
                [0.0, 0.0],
                [0.0, 1_000_000.0],
                [1_000_000.0, 0.0],
            ]
        ),
    )


def test_expression_transform_rejects_nonpositive_or_nonfinite_effective_lengths():
    counts = pandas.DataFrame({'RUN1': [1.0]}, index=['G1'])

    for invalid_value in (0.0, -1.0, numpy.nan, numpy.inf):
        effective_lengths = pandas.DataFrame(
            {'RUN1': [invalid_value]},
            index=['G1'],
        )
        with pytest.raises(ValueError, match='finite positive values'):
            _transform_raw_to_tpm(
                counts_df=counts,
                eff_length_df=effective_lengths,
            )


@pytest.mark.parametrize('invalid_value', [0.0, -1.0, numpy.nan, numpy.inf])
def test_fpkm_transform_rejects_invalid_tmm_library_size(invalid_value):
    counts = pandas.DataFrame({'RUN1': [1.0]}, index=['G1'])
    effective_lengths = pandas.DataFrame({'RUN1': [100.0]}, index=['G1'])
    metadata = pandas.DataFrame(
        {
            'run': ['RUN1'],
            'tmm_library_size': [invalid_value],
        }
    )

    with pytest.raises(ValueError, match='tmm_library_size'):
        _transform_raw_to_fpkm(
            counts_df=counts,
            eff_length_df=effective_lengths,
            metadata_df=metadata,
        )


@pytest.mark.integration
def test_generate_per_species_tables_uses_python_finalize_worker_for_skip_curation(tmp_path):
    fixture = _write_species_input_fixture(tmp_path)
    args = _build_args(tmp_path)
    args.skip_curation = True

    per_species_tables_module.generate_per_species_tables(
        args,
        context=PerSpeciesTableContext(metadata=fixture['metadata'], input_dir=fixture['input_dir']),
    )

    tables_dir = tmp_path / 'out' / 'per_species' / fixture['species_tag'] / 'tables'
    metadata_df = pandas.read_csv(tables_dir / '{}.metadata.tsv'.format(fixture['species_tag']), sep='\t')
    summary_df = pandas.read_csv(tables_dir / '{}.no.batch_effect_summary.tsv'.format(fixture['species_tag']), sep='\t')
    round_df = pandas.read_csv(tables_dir / '{}.no.curation_round_summary.tsv'.format(fixture['species_tag']), sep='\t')
    assert set(metadata_df['batch_corrected'].astype(str)) == {'no'}
    assert set(metadata_df['batch_alg_used'].astype(str)) == {'no'}
    assert summary_df.loc[0, 'skip_reason'] == 'skip_curation_requested'
    assert summary_df.loc[0, 'batch_effect_alg_applied'] == 'no'
    assert round_df.loc[0, 'step'] == 'skip_curation'
    assert (tmp_path / 'out' / 'per_species' / fixture['species_tag'] / 'per_species_completion_flag.txt').is_file()


@pytest.mark.integration
def test_generate_per_species_tables_uses_python_finalize_worker_for_disable_auto_outlier_filter(tmp_path):
    fixture = _write_species_input_fixture(tmp_path)
    args = _build_args(tmp_path)
    args.disable_auto_outlier_filter = True
    args.batch_effect_alg = 'sva'
    args.sva_nsv = '0'
    args.sva_B = '5'

    per_species_tables_module.generate_per_species_tables(
        args,
        context=PerSpeciesTableContext(metadata=fixture['metadata'], input_dir=fixture['input_dir']),
    )

    species_tag = fixture['species_tag']
    tables_dir = tmp_path / 'out' / 'per_species' / species_tag / 'tables'
    plots_dir = tmp_path / 'out' / 'per_species' / species_tag / 'plots'
    summary_df = pandas.read_csv(tables_dir / '{}.sva.batch_effect_summary.tsv'.format(species_tag), sep='\t')
    tau_df = pandas.read_csv(tables_dir / '{}.sva.tau.tsv'.format(species_tag), sep='\t')
    metadata_df = pandas.read_csv(tables_dir / '{}.metadata.tsv'.format(species_tag), sep='\t')
    assert int(summary_df.loc[0, 'resolved_sva_nsv']) == 0
    assert summary_df.loc[0, 'skip_reason'] == 'sva_nsv_zero'
    assert set(metadata_df['batch_corrected'].astype(str)) == {'no'}
    assert (tables_dir / '{}.sva.correlation_statistics.tsv'.format(species_tag)).exists()
    assert (plots_dir / '{}.before_after.sva.pdf'.format(species_tag)).is_file()
    assert (plots_dir / '{}.tau_hist.sva.pdf'.format(species_tag)).is_file()
    assert {'target_id', 'tau', 'highest', 'order'}.issubset(set(tau_df.columns))


@pytest.mark.integration
def test_generate_per_species_tables_supports_python_finalize_worker_for_latent_glm(tmp_path):
    fixture = _write_species_input_fixture(tmp_path)
    _inject_latent_batch_signal(fixture['input_dir'], fixture['species_tag'])
    args = _build_args(tmp_path)
    args.disable_auto_outlier_filter = True
    args.batch_effect_alg = 'latent_glm'
    args.latent_family = 'nb'
    args.latent_k = '1'
    args.latent_k_max = 3
    args.latent_max_iter = 50
    args.latent_tol = 1e-6

    per_species_tables_module.generate_per_species_tables(
        args,
        context=PerSpeciesTableContext(metadata=fixture['metadata'], input_dir=fixture['input_dir']),
    )

    species_tag = fixture['species_tag']
    tables_dir = tmp_path / 'out' / 'per_species' / species_tag / 'tables'
    plots_dir = tmp_path / 'out' / 'per_species' / species_tag / 'plots'
    summary_df = pandas.read_csv(tables_dir / '{}.latent_glm.batch_effect_summary.tsv'.format(species_tag), sep='\t')
    corrected_df = pandas.read_csv(tables_dir / '{}.latent_glm.tc.tsv'.format(species_tag), sep='\t', index_col=0)
    assert int(summary_df.loc[0, 'resolved_latent_k']) == 1
    assert summary_df.loc[0, 'latent_family'] == 'nb'
    assert (corrected_df.to_numpy(dtype=float) >= 0.0).all()
    assert (plots_dir / '{}.before_after.latent_glm.pdf'.format(species_tag)).is_file()
    assert (plots_dir / '{}.tau_hist.latent_glm.pdf'.format(species_tag)).is_file()


def test_transform_raw_to_tpm_is_cpm_when_effective_length_is_one():
    counts = pandas.DataFrame({'SRR001': [10.0, 10.0]}, index=['short', 'long'])
    eff = pandas.DataFrame({'SRR001': [1.0, 1.0]}, index=['short', 'long'])
    tpm = _transform_raw_to_tpm(counts, eff)
    assert numpy.allclose(tpm['SRR001'], [5.0e5, 5.0e5])


def test_apply_transformation_rejects_fpkm_without_effective_length():
    counts = pandas.DataFrame({'SRR001': [10.0, 10.0]}, index=['short', 'long'])
    eff = pandas.DataFrame({'SRR001': [1.0, 1.0]}, index=['short', 'long'])
    metadata = pandas.DataFrame({'run': ['SRR001']})
    with pytest.raises(ValueError, match='FPKM is undefined'):
        _apply_transformation_logic(
            counts,
            eff,
            'log2p1-fpkm',
            'no',
            'before_batch',
            metadata,
            {'SRR001': 'none'},
        )


def test_apply_transformation_keeps_tpm_for_unlength_normalized_runs():
    counts = pandas.DataFrame({'SRR001': [10.0, 10.0]}, index=['short', 'long'])
    eff = pandas.DataFrame({'SRR001': [1.0, 1.0]}, index=['short', 'long'])
    metadata = pandas.DataFrame({'run': ['SRR001']})
    transformed = _apply_transformation_logic(
        counts,
        eff,
        'none-tpm',
        'no',
        'before_batch',
        metadata,
        {'SRR001': 'none'},
    )
    assert numpy.allclose(transformed['SRR001'], [5.0e5, 5.0e5])


def test_quant_model_validation_rejects_duplicates_and_missing_runs(tmp_path):
    path = tmp_path / 'Species_A_quant_model.tsv'
    path.write_text(
        'run\tbackend\tlength_model\n'
        'SRR001\toarfish\tnone\n'
        'SRR001\toarfish\tnone\n'
    )
    with pytest.raises(ValueError, match='duplicate run values'):
        load_quant_model_table(path)

    model = pandas.DataFrame(
        {'run': ['SRR001'], 'backend': ['oarfish'], 'length_model': ['none']}
    )
    with pytest.raises(ValueError, match=r'missing run\(s\): SRR002'):
        resolve_length_models(['SRR001', 'SRR002'], model)


def test_compute_corr_matrix_preserves_undefined_constant_sample():
    counts = pandas.DataFrame(
        {
            'A': [1.0, 2.0, 3.0, 4.0],
            'B': [2.0, 4.0, 6.0, 8.0],
            'C': [5.0, 5.0, 5.0, 5.0],
        }
    )

    corr = _compute_corr_matrix(counts, 'pearson')
    pca = _compute_pca_coordinates(corr)

    assert corr.loc['C'].isna().all()
    assert corr.loc[:, 'C'].isna().all()
    assert numpy.isfinite(pca[0]).all()
    assert numpy.isfinite(pca[1]).all()
    assert numpy.isnan(pca[2]).all()
    assert numpy.isclose(corr.loc['A', 'B'], 1.0)
