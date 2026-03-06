import os
from types import SimpleNamespace

import pandas

from amalgkit.command_context import CscaContext, CurateContext
from amalgkit.util import Metadata
from amalgkit import wsfilter as wsfilter_module
from amalgkit import csfilter as csfilter_module
from amalgkit import finalize as finalize_module


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

    def fake_resolve_curate_input(_args):
        return metadata, str(tmp_path / 'input')

    def fake_curate_main(curate_args, context=None):
        captured['context'] = context
        tables_dir = os.path.join(curate_args.out_dir, 'curate', 'Species_A', 'tables')
        plots_dir = os.path.join(curate_args.out_dir, 'curate', 'Species_A', 'plots')
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

    monkeypatch.setattr(wsfilter_module, 'resolve_curate_input', fake_resolve_curate_input)
    monkeypatch.setattr(wsfilter_module, 'curate_main', fake_curate_main)
    def fake_exclusion_plot(df_metadata, out_pdf_path, r_util_path, y_label='Sample count', font_size=8):
        _ = (df_metadata, r_util_path, y_label, font_size)
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
    assert isinstance(captured['context'], CurateContext)
    assert captured['context'].metadata is metadata
    assert captured['context'].input_dir == str(tmp_path / 'input')


def test_wsfilter_uses_latest_filter_metadata_when_inferred(tmp_path, monkeypatch):
    args = _base_args(tmp_path)
    metadata = _base_metadata()
    latest_meta = tmp_path / 'out' / 'csfilter' / 'metadata.tsv'
    latest_meta.parent.mkdir(parents=True, exist_ok=True)
    latest_meta.write_text('run\texclusion\nR1\tno\n')
    captured = {}

    def fake_resolve_curate_input(passed_args):
        captured['metadata'] = passed_args.metadata
        return metadata, str(tmp_path / 'input')

    def fake_curate_main(curate_args, context=None):
        captured['context'] = context
        tables_dir = os.path.join(curate_args.out_dir, 'curate', 'Species_A', 'tables')
        os.makedirs(tables_dir, exist_ok=True)
        pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
        }).to_csv(os.path.join(tables_dir, 'Species_A.metadata.tsv'), sep='\t', index=False)

    monkeypatch.setattr(wsfilter_module, 'resolve_curate_input', fake_resolve_curate_input)
    monkeypatch.setattr(wsfilter_module, 'curate_main', fake_curate_main)
    def fake_exclusion_plot(df_metadata, out_pdf_path, r_util_path, y_label='Sample count', font_size=8):
        _ = (df_metadata, r_util_path, y_label, font_size)
        os.makedirs(os.path.dirname(out_pdf_path), exist_ok=True)
        with open(out_pdf_path, 'wb') as handle:
            handle.write(b'%PDF-1.4\n')
    monkeypatch.setattr(wsfilter_module, 'save_exclusion_plot_pdf', fake_exclusion_plot)

    wsfilter_module.wsfilter_main(args)

    assert captured['metadata'] == os.path.realpath(str(latest_meta))
    assert isinstance(captured['context'], CurateContext)
    assert captured['context'].metadata is metadata
    assert captured['context'].input_dir == str(tmp_path / 'input')


def test_csfilter_outputs_metadata_excluded_and_root_pdfs(tmp_path, monkeypatch):
    args = _base_args(tmp_path)
    metadata = _base_metadata()
    captured = {}

    def fake_resolve_curate_input(_args):
        return metadata, str(tmp_path / 'input')

    def fake_curate_main(_curate_args, context=None):
        captured['curate_context'] = context
        return None

    def fake_csca_main(csca_args, context=None):
        captured['sample_group'] = csca_args.sample_group
        captured['outlier_method'] = csca_args.outlier_method
        captured['batch_effect_alg'] = csca_args.batch_effect_alg
        captured['plot_mode'] = csca_args.plot_mode
        captured['csca_context'] = context
        csca_dir = os.path.join(csca_args.out_dir, 'csca')
        os.makedirs(csca_dir, exist_ok=True)
        pandas.DataFrame({
            'run': ['R2'],
            'exclusion': ['low_cross_species_group_correlation'],
            'cs_margin_corrected': [-0.2],
        }).to_csv(os.path.join(csca_dir, 'metadata.tsv'), sep='\t', index=False)
        with open(os.path.join(csca_dir, 'csca_overview.pdf'), 'wb') as handle:
            handle.write(b'%PDF-1.4\n')

    monkeypatch.setattr(csfilter_module, 'resolve_curate_input', fake_resolve_curate_input)
    monkeypatch.setattr(csfilter_module, 'curate_main', fake_curate_main)
    monkeypatch.setattr(csfilter_module, 'csca_main', fake_csca_main)
    def fake_exclusion_plot(df_metadata, out_pdf_path, r_util_path, y_label='Sample count', font_size=8):
        _ = (df_metadata, r_util_path, y_label, font_size)
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
    assert not (tmp_path / 'out' / 'csfilter' / 'Species_A').exists()
    assert not (tmp_path / 'out' / 'csfilter' / 'plots').exists()
    assert not (tmp_path / 'out' / 'csfilter' / 'tables').exists()
    assert captured['sample_group'] == 'leaf|root'
    assert captured['outlier_method'] == 'robust_margin'
    assert captured['batch_effect_alg'] == 'no'
    assert captured['plot_mode'] == 'single'
    assert isinstance(captured['curate_context'], CurateContext)
    assert captured['curate_context'].metadata is metadata
    assert captured['curate_context'].input_dir == str(tmp_path / 'input')
    assert isinstance(captured['csca_context'], CscaContext)
    assert captured['csca_context'].metadata is metadata


def test_finalize_outputs_tables_and_merged_metadata(tmp_path, monkeypatch):
    args = _base_args(tmp_path)
    args.batch_effect_alg = 'sva'
    metadata = _base_metadata()
    captured = {}

    def fake_resolve_curate_input(_args):
        return metadata, str(tmp_path / 'input')

    def fake_curate_main(curate_args, context=None):
        captured['seed'] = getattr(curate_args, 'seed')
        captured['context'] = context
        tables_dir = os.path.join(curate_args.out_dir, 'curate', 'Species_A', 'tables')
        plots_dir = os.path.join(curate_args.out_dir, 'curate', 'Species_A', 'plots')
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

    monkeypatch.setattr(finalize_module, 'resolve_curate_input', fake_resolve_curate_input)
    monkeypatch.setattr(finalize_module, 'curate_main', fake_curate_main)
    def fake_exclusion_plot(df_metadata, out_pdf_path, r_util_path, y_label='Sample count', font_size=8):
        _ = (df_metadata, r_util_path, y_label, font_size)
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
    assert (tmp_path / 'out' / 'finalize' / 'finalize_exclusion.pdf').is_file()
    assert (tmp_path / 'out' / 'finalize' / 'Species_A' / 'Species_A_expression_uncorrected.tsv').is_file()
    assert (tmp_path / 'out' / 'finalize' / 'Species_A' / 'Species_A_expression.tsv').is_file()
    assert (tmp_path / 'out' / 'finalize' / 'Species_A' / 'Species_A_batch_effect_summary.tsv').is_file()
    assert (tmp_path / 'out' / 'finalize' / 'Species_A' / 'Species_A_batch_compare_sva.pdf').is_file()
    assert not (tmp_path / 'out' / 'finalize' / 'Species_A' / 'tables').exists()
    assert isinstance(captured['context'], CurateContext)
    assert captured['context'].metadata is metadata
    assert captured['context'].input_dir == str(tmp_path / 'input')
