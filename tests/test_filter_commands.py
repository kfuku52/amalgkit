import os
from types import SimpleNamespace

import pandas

from amalgkit.command_context import CrossSpeciesFilterContext, PerSpeciesTableContext
from amalgkit.util import Metadata
from amalgkit import wsfilter as wsfilter_module
from amalgkit import csfilter as csfilter_module
from amalgkit import finalize as finalize_module
from amalgkit import per_species_tables as per_species_tables_module


def _base_metadata():
    return Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1', 'R2'],
        'scientific_name': ['Species A', 'Species A'],
        'sample_group': ['leaf', 'root'],
        'bioproject': ['BP1', 'BP1'],
        'exclusion': ['no', 'no'],
    }))


def _base_args(tmp_path):
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
        orthogroup_table=None,
        dir_busco=None,
        missing_strategy='em_pca',
        margin_threshold=0.0,
        robust_z_threshold=-2.5,
    )


def test_wsfilter_outputs_metadata_excluded_and_species_pdfs(tmp_path, monkeypatch):
    args = _base_args(tmp_path)
    metadata = _base_metadata()
    captured = {}

    def fake_resolve_per_species_input(_args):
        return metadata, str(tmp_path / 'input')

    def fake_generate_per_species_tables(per_species_args, context=None):
        captured['context'] = context
        tables_dir = os.path.join(per_species_args.out_dir, 'per_species', 'Species_A', 'tables')
        plots_dir = os.path.join(per_species_args.out_dir, 'per_species', 'Species_A', 'plots')
        os.makedirs(tables_dir, exist_ok=True)
        os.makedirs(plots_dir, exist_ok=True)
        pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['low_within_sample_group_correlation'],
            'ws_margin': [-0.3],
        }).to_csv(
            os.path.join(tables_dir, 'Species_A.metadata.tsv'),
            sep='\t',
            index=False,
        )
        with open(os.path.join(plots_dir, 'ws_plot.pdf'), 'wb') as handle:
            handle.write(b'%PDF-1.4\n')

    monkeypatch.setattr(wsfilter_module, 'resolve_per_species_input', fake_resolve_per_species_input)
    monkeypatch.setattr(wsfilter_module, 'generate_per_species_tables', fake_generate_per_species_tables)
    def fake_exclusion_plot(df_metadata, out_pdf_path, y_label='Sample count', font_size=8):
        _ = (df_metadata, y_label, font_size)
        os.makedirs(os.path.dirname(out_pdf_path), exist_ok=True)
        with open(out_pdf_path, 'wb') as handle:
            handle.write(b'%PDF-1.4\n')
    monkeypatch.setattr(wsfilter_module, 'save_exclusion_plot_pdf', fake_exclusion_plot)

    wsfilter_module.wsfilter_main(args)

    out_metadata = pandas.read_csv(tmp_path / 'out' / 'wsfilter' / 'metadata.tsv', sep='\t')
    out_excluded = pandas.read_csv(tmp_path / 'out' / 'wsfilter' / 'excluded.tsv', sep='\t')
    assert out_metadata.loc[out_metadata['run'] == 'R1', 'exclusion'].iloc[0] == 'low_within_sample_group_correlation'
    assert out_metadata.loc[out_metadata['run'] == 'R2', 'exclusion'].iloc[0] == 'no'
    assert out_excluded['run'].tolist() == ['R1']
    assert out_excluded['exclusion'].tolist() == ['low_within_sample_group_correlation']
    assert (tmp_path / 'out' / 'wsfilter' / 'wsfilter_exclusion.pdf').is_file()
    assert (tmp_path / 'out' / 'wsfilter' / 'Species_A' / 'Species_A_ws_plot.pdf').is_file()
    assert not (tmp_path / 'out' / 'wsfilter' / 'plots').exists()
    assert not (tmp_path / 'out' / 'wsfilter' / 'tables').exists()
    assert isinstance(captured['context'], PerSpeciesTableContext)
    assert captured['context'].metadata is metadata
    assert captured['context'].input_dir == str(tmp_path / 'input')


def test_wsfilter_uses_latest_filter_metadata_when_inferred(tmp_path, monkeypatch):
    args = _base_args(tmp_path)
    metadata = _base_metadata()
    latest_meta = tmp_path / 'out' / 'csfilter' / 'metadata.tsv'
    latest_meta.parent.mkdir(parents=True, exist_ok=True)
    latest_meta.write_text('run\texclusion\nR1\tno\n')
    captured = {}

    def fake_resolve_per_species_input(passed_args):
        captured['metadata'] = passed_args.metadata
        return metadata, str(tmp_path / 'input')

    def fake_generate_per_species_tables(per_species_args, context=None):
        captured['context'] = context
        tables_dir = os.path.join(per_species_args.out_dir, 'per_species', 'Species_A', 'tables')
        os.makedirs(tables_dir, exist_ok=True)
        pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
        }).to_csv(os.path.join(tables_dir, 'Species_A.metadata.tsv'), sep='\t', index=False)

    monkeypatch.setattr(wsfilter_module, 'resolve_per_species_input', fake_resolve_per_species_input)
    monkeypatch.setattr(wsfilter_module, 'generate_per_species_tables', fake_generate_per_species_tables)
    def fake_exclusion_plot(df_metadata, out_pdf_path, y_label='Sample count', font_size=8):
        _ = (df_metadata, y_label, font_size)
        os.makedirs(os.path.dirname(out_pdf_path), exist_ok=True)
        with open(out_pdf_path, 'wb') as handle:
            handle.write(b'%PDF-1.4\n')
    monkeypatch.setattr(wsfilter_module, 'save_exclusion_plot_pdf', fake_exclusion_plot)

    wsfilter_module.wsfilter_main(args)

    assert captured['metadata'] == os.path.realpath(str(latest_meta))
    assert isinstance(captured['context'], PerSpeciesTableContext)
    assert captured['context'].metadata is metadata
    assert captured['context'].input_dir == str(tmp_path / 'input')


def test_csfilter_outputs_metadata_excluded_and_root_pdfs(tmp_path, monkeypatch):
    args = _base_args(tmp_path)
    metadata = _base_metadata()
    captured = {}

    def fake_resolve_per_species_input(_args):
        return metadata, str(tmp_path / 'input')

    def fake_generate_per_species_tables(_per_species_args, context=None):
        captured['per_species_context'] = context
        return None

    def fake_run_cross_species_filter(cross_species_args, context=None):
        captured['sample_group'] = cross_species_args.sample_group
        captured['outlier_method'] = cross_species_args.outlier_method
        captured['batch_effect_alg'] = cross_species_args.batch_effect_alg
        captured['plot_mode'] = cross_species_args.plot_mode
        captured['cross_species_context'] = context
        cross_species_dir = os.path.join(cross_species_args.out_dir, 'cross_species')
        os.makedirs(cross_species_dir, exist_ok=True)
        pandas.DataFrame({
            'run': ['R2'],
            'exclusion': ['low_cross_species_group_correlation'],
            'cs_margin_corrected': [-0.2],
        }).to_csv(os.path.join(cross_species_dir, 'metadata.tsv'), sep='\t', index=False)
        with open(os.path.join(cross_species_dir, 'cross_species_overview.pdf'), 'wb') as handle:
            handle.write(b'%PDF-1.4\n')
        for name in [
            'cross_species_group_cor_scatter.pdf',
            'cross_species_sample_number_heatmap.pdf',
            'cross_species_SVA_heatmap.pdf',
            'cross_species_averaged_summary.pdf',
            'cross_species_unaveraged_pca_PC34_uncorrected.pdf',
        ]:
            with open(os.path.join(cross_species_dir, name), 'wb') as handle:
                handle.write(b'%PDF-1.4\n')

    monkeypatch.setattr(csfilter_module, 'resolve_per_species_input', fake_resolve_per_species_input)
    monkeypatch.setattr(csfilter_module, 'generate_per_species_tables', fake_generate_per_species_tables)
    monkeypatch.setattr(csfilter_module, 'run_cross_species_filter', fake_run_cross_species_filter)
    def fake_exclusion_plot(df_metadata, out_pdf_path, y_label='Sample count', font_size=8):
        _ = (df_metadata, y_label, font_size)
        os.makedirs(os.path.dirname(out_pdf_path), exist_ok=True)
        with open(out_pdf_path, 'wb') as handle:
            handle.write(b'%PDF-1.4\n')
    monkeypatch.setattr(csfilter_module, 'save_exclusion_plot_pdf', fake_exclusion_plot)

    csfilter_module.csfilter_main(args)

    out_metadata = pandas.read_csv(tmp_path / 'out' / 'csfilter' / 'metadata.tsv', sep='\t')
    out_excluded = pandas.read_csv(tmp_path / 'out' / 'csfilter' / 'excluded.tsv', sep='\t')
    assert out_metadata.loc[out_metadata['run'] == 'R1', 'exclusion'].iloc[0] == 'no'
    assert out_metadata.loc[out_metadata['run'] == 'R2', 'exclusion'].iloc[0] == 'low_cross_species_group_correlation'
    assert 'cs_margin' in out_metadata.columns
    assert 'cs_margin_corrected' not in out_metadata.columns
    assert out_metadata.loc[out_metadata['run'] == 'R2', 'cs_margin'].iloc[0] == -0.2
    assert out_excluded['run'].tolist() == ['R2']
    assert out_excluded['exclusion'].tolist() == ['low_cross_species_group_correlation']
    assert 'cs_margin' in out_excluded.columns
    assert 'cs_margin_corrected' not in out_excluded.columns
    assert (tmp_path / 'out' / 'csfilter' / 'csfilter_exclusion.pdf').is_file()
    assert (tmp_path / 'out' / 'csfilter' / 'csfilter_overview.pdf').is_file()
    assert (tmp_path / 'out' / 'csfilter' / 'csfilter_group_cor_scatter.pdf').is_file()
    assert (tmp_path / 'out' / 'csfilter' / 'csfilter_sample_number_heatmap.pdf').is_file()
    assert (tmp_path / 'out' / 'csfilter' / 'csfilter_SVA_heatmap.pdf').is_file()
    assert (tmp_path / 'out' / 'csfilter' / 'csfilter_averaged_summary.pdf').is_file()
    assert (tmp_path / 'out' / 'csfilter' / 'csfilter_unaveraged_pca_PC34_uncorrected.pdf').is_file()
    assert not (tmp_path / 'out' / 'csfilter' / 'Species_A').exists()
    assert not (tmp_path / 'out' / 'csfilter' / 'plots').exists()
    assert not (tmp_path / 'out' / 'csfilter' / 'tables').exists()
    assert captured['sample_group'] == 'leaf|root'
    assert captured['outlier_method'] == 'robust_margin'
    assert captured['batch_effect_alg'] == 'no'
    assert captured['plot_mode'] == 'single'
    assert isinstance(captured['per_species_context'], PerSpeciesTableContext)
    assert captured['per_species_context'].metadata is metadata
    assert captured['per_species_context'].input_dir == str(tmp_path / 'input')
    assert isinstance(captured['cross_species_context'], CrossSpeciesFilterContext)
    assert captured['cross_species_context'].metadata is metadata


def test_finalize_outputs_tables_and_merged_metadata(tmp_path, monkeypatch):
    args = _base_args(tmp_path)
    args.batch_effect_alg = 'sva'
    args.sva_backend = 'python'
    args.combatseq_backend = 'python'
    args.ruvseq_backend = 'python'
    args.latent_family = 'poisson'
    args.latent_k = 2
    args.latent_k_max = 6
    args.latent_max_iter = 111
    args.latent_tol = 1e-4
    metadata = _base_metadata()
    captured = {}

    def fake_resolve_per_species_input(_args):
        return metadata, str(tmp_path / 'input')

    def fake_generate_per_species_tables(per_species_args, context=None):
        captured['seed'] = getattr(per_species_args, 'seed')
        captured['sva_backend'] = getattr(per_species_args, 'sva_backend')
        captured['combatseq_backend'] = getattr(per_species_args, 'combatseq_backend')
        captured['ruvseq_backend'] = getattr(per_species_args, 'ruvseq_backend')
        captured['latent_family'] = getattr(per_species_args, 'latent_family')
        captured['latent_k'] = getattr(per_species_args, 'latent_k')
        captured['latent_k_max'] = getattr(per_species_args, 'latent_k_max')
        captured['latent_max_iter'] = getattr(per_species_args, 'latent_max_iter')
        captured['latent_tol'] = getattr(per_species_args, 'latent_tol')
        captured['context'] = context
        tables_dir = os.path.join(per_species_args.out_dir, 'per_species', 'Species_A', 'tables')
        plots_dir = os.path.join(per_species_args.out_dir, 'per_species', 'Species_A', 'plots')
        os.makedirs(tables_dir, exist_ok=True)
        os.makedirs(plots_dir, exist_ok=True)
        pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['low_cross_species_group_correlation'],
            'batch_corrected': ['yes'],
            'batch_alg_used': ['sva'],
        }).to_csv(os.path.join(tables_dir, 'Species_A.metadata.tsv'), sep='\t', index=False)
        pandas.DataFrame({'target_id': ['G1'], 'R1': [1.0]}).to_csv(
            os.path.join(tables_dir, 'Species_A.uncorrected.tc.tsv'),
            sep='\t',
            index=False,
        )
        pandas.DataFrame({'target_id': ['G1'], 'R1': [2.0]}).to_csv(
            os.path.join(tables_dir, 'Species_A.sva.tc.tsv'),
            sep='\t',
            index=False,
        )
        pandas.DataFrame({
            'scientific_name': ['Species A'],
            'batch_effect_alg_requested': ['sva'],
            'batch_effect_alg_applied': ['sva'],
            'corrected_run_count': [1],
            'skip_reason': [''],
        }).to_csv(
            os.path.join(tables_dir, 'Species_A.sva.batch_effect_summary.tsv'),
            sep='\t',
            index=False,
        )
        with open(os.path.join(plots_dir, 'Species_A.batch_compare.sva.pdf'), 'wb') as handle:
            handle.write(b'%PDF-1.4\n')

    monkeypatch.setattr(finalize_module, 'resolve_per_species_input', fake_resolve_per_species_input)
    monkeypatch.setattr(finalize_module, 'generate_per_species_tables', fake_generate_per_species_tables)
    def fake_exclusion_plot(df_metadata, out_pdf_path, y_label='Sample count', font_size=8):
        _ = (df_metadata, y_label, font_size)
        os.makedirs(os.path.dirname(out_pdf_path), exist_ok=True)
        with open(out_pdf_path, 'wb') as handle:
            handle.write(b'%PDF-1.4\n')
    monkeypatch.setattr(finalize_module, 'save_exclusion_plot_pdf', fake_exclusion_plot)

    finalize_module.finalize_main(args)

    out_metadata = pandas.read_csv(tmp_path / 'out' / 'finalize' / 'metadata.tsv', sep='\t')
    assert out_metadata.loc[out_metadata['run'] == 'R1', 'exclusion'].iloc[0] == 'low_cross_species_group_correlation'
    assert out_metadata.loc[out_metadata['run'] == 'R1', 'batch_corrected'].iloc[0] == 'yes'
    assert out_metadata.loc[out_metadata['run'] == 'R1', 'batch_alg_used'].iloc[0] == 'sva'
    assert captured['seed'] == 'auto'
    assert captured['sva_backend'] == 'python'
    assert captured['combatseq_backend'] == 'python'
    assert captured['ruvseq_backend'] == 'python'
    assert captured['latent_family'] == 'poisson'
    assert captured['latent_k'] == 2
    assert captured['latent_k_max'] == 6
    assert captured['latent_max_iter'] == 111
    assert captured['latent_tol'] == 1e-4
    assert (tmp_path / 'out' / 'finalize' / 'finalize_exclusion.pdf').is_file()
    assert (tmp_path / 'out' / 'finalize' / 'Species_A' / 'Species_A_expression_uncorrected.tsv').is_file()
    assert (tmp_path / 'out' / 'finalize' / 'Species_A' / 'Species_A_expression.tsv').is_file()
    assert (tmp_path / 'out' / 'finalize' / 'Species_A' / 'Species_A_batch_effect_summary.tsv').is_file()
    assert (tmp_path / 'out' / 'finalize' / 'Species_A' / 'Species_A_batch_compare_sva.pdf').is_file()
    assert not (tmp_path / 'out' / 'finalize' / 'Species_A' / 'tables').exists()
    assert isinstance(captured['context'], PerSpeciesTableContext)
    assert captured['context'].metadata is metadata
    assert captured['context'].input_dir == str(tmp_path / 'input')


def test_finalize_main_runs_without_rscript_for_supported_python_worker(tmp_path, monkeypatch):
    species = 'Species A'
    species_tag = 'Species_A'
    input_dir = tmp_path / 'input'
    species_dir = input_dir / species_tag
    species_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = tmp_path / 'metadata.tsv'
    metadata_df = pandas.DataFrame(
        {
            'run': ['R1', 'R2', 'R3', 'R4'],
            'scientific_name': [species] * 4,
            'sample_group': ['leaf', 'leaf', 'root', 'root'],
            'bioproject': ['BP1', 'BP1', 'BP2', 'BP2'],
            'exclusion': ['no'] * 4,
            'mapping_rate': [100.0] * 4,
            'instrument': ['Illumina'] * 4,
            'lib_layout': ['PAIRED'] * 4,
            'lib_selection': ['cDNA'] * 4,
            'total_spots': [1_000_000] * 4,
            'total_bases': [150_000_000] * 4,
        }
    )
    metadata_df.to_csv(metadata_path, sep='\t', index=False)

    counts_df = pandas.DataFrame(
        {
            'target_id': ['G1', 'G2', 'G3', 'G4'],
            'R1': [10, 20, 5, 0],
            'R2': [11, 19, 5, 0],
            'R3': [5, 3, 20, 0],
            'R4': [6, 4, 21, 0],
        }
    )
    eff_length_df = pandas.DataFrame(
        {
            'target_id': ['G1', 'G2', 'G3', 'G4'],
            'R1': [1000, 1000, 1000, 1000],
            'R2': [1000, 1000, 1000, 1000],
            'R3': [1000, 1000, 1000, 1000],
            'R4': [1000, 1000, 1000, 1000],
        }
    )
    counts_df.to_csv(species_dir / '{}_est_counts.tsv'.format(species_tag), sep='\t', index=False)
    eff_length_df.to_csv(species_dir / '{}_eff_length.tsv'.format(species_tag), sep='\t', index=False)

    args = _base_args(tmp_path)
    args.metadata = str(metadata_path)
    args.input_dir = str(input_dir)
    args.batch_effect_alg = 'sva'
    args.sva_backend = 'python'
    args.combatseq_backend = 'python'
    args.ruvseq_backend = 'python'
    args.sva_nsv = '0'
    args.sva_B = '5'
    args.sva_B_auto_max = 100
    args.seed = 'auto'
    args.ruvseq_control_genes = 'auto'
    args.ruvseq_k = 'auto'
    args.ruvseq_k_max = 5
    args.ruvseq_control_top_n = 1000
    args.ruvseq_min_controls = 100
    args.disable_auto_outlier_filter = True
    args.python_executable = 'python'
    args.latent_family = 'nb'
    args.latent_k = 'auto'
    args.latent_k_max = 5
    args.latent_max_iter = 200
    args.latent_tol = 1e-5

    finalize_module.finalize_main(args)

    out_root = tmp_path / 'out' / 'finalize' / species_tag
    out_metadata = pandas.read_csv(tmp_path / 'out' / 'finalize' / 'metadata.tsv', sep='\t')
    assert set(out_metadata['batch_corrected'].astype(str)) == {'no'}
    assert set(out_metadata['batch_alg_used'].astype(str)) == {'no'}
    assert (out_root / '{}_expression_uncorrected.tsv'.format(species_tag)).is_file()
    assert (out_root / '{}_expression.tsv'.format(species_tag)).is_file()
    assert (out_root / '{}_batch_effect_summary.tsv'.format(species_tag)).is_file()
    assert (out_root / '{}_batch_compare_sva.pdf'.format(species_tag)).is_file()
