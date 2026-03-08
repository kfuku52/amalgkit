import os
from types import SimpleNamespace

import pandas

from amalgkit.command_context import PerSpeciesTableContext
from amalgkit.metadata_utils import Metadata
from amalgkit import per_species_tables as per_species_tables_module


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


def _build_args(tmp_path):
    return SimpleNamespace(
        out_dir=str(tmp_path / 'out'),
        redo=False,
        metadata='inferred',
        input_dir='inferred',
        sample_group=None,
        sample_group_color='DEFAULT',
        batch=None,
        threads='auto',
        internal_jobs='auto',
        internal_cpu_budget='auto',
        dist_method='pearson',
        mapping_rate=0.20,
        correlation_threshold=0.30,
        plot_intermediate=False,
        one_outlier_per_iter=False,
        norm='log2p1-fpkm',
        clip_negative=True,
        maintain_zero=True,
        batch_effect_alg='no',
        skip_curation=False,
        disable_auto_outlier_filter=False,
        worker_mode='finalize',
        ruvseq_control_genes='auto',
        ruvseq_k='auto',
        ruvseq_k_max=5,
        ruvseq_control_top_n=1000,
        ruvseq_min_controls=100,
        seed='auto',
        sva_nsv='auto',
        sva_B='auto',
        sva_B_auto_max=100,
        sva_backend='python',
        combatseq_backend='python',
        ruvseq_backend='python',
        python_executable='python',
        latent_family='nb',
        latent_k='auto',
        latent_k_max=5,
        latent_max_iter=200,
        latent_tol=1e-5,
    )


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
    assert (plots_dir / '{}.batch_compare.sva.pdf'.format(species_tag)).is_file()
    assert (plots_dir / '{}.tau_hist.sva.pdf'.format(species_tag)).is_file()
    assert {'target_id', 'tau', 'highest', 'order'}.issubset(set(tau_df.columns))
