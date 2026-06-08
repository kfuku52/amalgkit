import json
import os
from types import SimpleNamespace

import pandas
import pytest

from amalgkit.exceptions import AmalgkitExit
from amalgkit.rerun import rerun_main


def _write_metadata(path_metadata):
    df = pandas.DataFrame(
        [
            {'scientific_name': 'Homo sapiens', 'run': 'SRR001', 'exclusion': 'no'},
            {'scientific_name': 'Mus musculus', 'run': 'SRR002', 'exclusion': 'no'},
            {'scientific_name': 'Homo sapiens', 'run': 'SRR003', 'exclusion': 'no'},
        ]
    )
    df.to_csv(path_metadata, sep='\t', index=False)


def _write_report(path_report, out_dir, metadata_path, summary, issues):
    payload = {
        'generated_at': '2026-03-25T00:00:00',
        'metadata_path': str(metadata_path),
        'out_dir': str(out_dir),
        'requested_checks': [row['check'] for row in summary],
        'summary': summary,
        'issues': issues,
    }
    with open(path_report, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _build_args(tmp_path, **overrides):
    data = {
        'out_dir': str(tmp_path),
        'report': 'inferred',
        'metadata': 'report',
        'manifest': 'inferred',
        'check': None,
        'run': None,
        'species': None,
        'redo': True,
        'dry_run': False,
        'include_warnings': False,
        'threads': 'auto',
        'internal_jobs': 'auto',
        'internal_cpu_budget': 'auto',
        'download_dir': 'inferred',
        'download_lock_dir': 'inferred',
    }
    data.update(overrides)
    return SimpleNamespace(**data)


def test_rerun_run_checks_use_report_targets(tmp_path, monkeypatch):
    metadata_dir = tmp_path / 'metadata'
    metadata_dir.mkdir()
    metadata_path = metadata_dir / 'metadata.tsv'
    _write_metadata(metadata_path)

    sanity_dir = tmp_path / 'sanity'
    sanity_dir.mkdir()
    report_path = sanity_dir / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[
            {'check': 'getfastq'},
            {'check': 'quant'},
        ],
        issues=[
            {'check': 'getfastq', 'severity': 'error', 'target_type': 'run', 'target_id': 'SRR002'},
            {'check': 'quant', 'severity': 'error', 'target_type': 'run', 'target_id': 'SRR003'},
        ],
    )

    observed = {}

    def fake_getfastq_main(args):
        df = pandas.read_csv(args.metadata, sep='\t')
        observed['getfastq_runs'] = df['run'].tolist()
        observed['getfastq_redo'] = args.redo

    def fake_quant_main(args):
        df = pandas.read_csv(args.metadata, sep='\t')
        observed['quant_runs'] = df['run'].tolist()
        observed['quant_redo'] = args.redo

    monkeypatch.setattr('amalgkit.rerun.getfastq_main', fake_getfastq_main)
    monkeypatch.setattr('amalgkit.rerun.quant_main', fake_quant_main)

    rerun_main(_build_args(tmp_path))

    assert observed['getfastq_runs'] == ['SRR002']
    assert observed['quant_runs'] == ['SRR003']
    assert observed['getfastq_redo'] is True
    assert observed['quant_redo'] is True


def test_rerun_run_checks_respect_redo_flag(tmp_path, monkeypatch):
    metadata_dir = tmp_path / 'metadata'
    metadata_dir.mkdir()
    metadata_path = metadata_dir / 'metadata.tsv'
    _write_metadata(metadata_path)

    sanity_dir = tmp_path / 'sanity'
    sanity_dir.mkdir()
    report_path = sanity_dir / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[
            {'check': 'getfastq'},
            {'check': 'quant'},
        ],
        issues=[
            {'check': 'getfastq', 'severity': 'error', 'target_type': 'run', 'target_id': 'SRR002'},
            {'check': 'quant', 'severity': 'error', 'target_type': 'run', 'target_id': 'SRR003'},
        ],
    )

    observed = {}

    def fake_getfastq_main(args):
        observed['getfastq_redo'] = args.redo

    def fake_quant_main(args):
        observed['quant_redo'] = args.redo

    monkeypatch.setattr('amalgkit.rerun.getfastq_main', fake_getfastq_main)
    monkeypatch.setattr('amalgkit.rerun.quant_main', fake_quant_main)

    rerun_main(_build_args(tmp_path, redo=False))

    assert observed['getfastq_redo'] is False
    assert observed['quant_redo'] is False


def test_rerun_busco_warning_repairs_summary_only(tmp_path, monkeypatch):
    metadata_dir = tmp_path / 'metadata'
    metadata_dir.mkdir()
    metadata_path = metadata_dir / 'metadata.tsv'
    _write_metadata(metadata_path)

    sanity_dir = tmp_path / 'sanity'
    sanity_dir.mkdir()
    report_path = sanity_dir / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[{'check': 'busco'}],
        issues=[
            {'check': 'busco', 'severity': 'warning', 'target_type': 'global', 'target_id': 'busco_completeness'},
        ],
    )

    observed = {'process_calls': 0, 'plot_calls': 0}

    def fake_process_species_busco(*args, **kwargs):
        observed['process_calls'] += 1

    def fake_generate_busco_species_plot(*args, **kwargs):
        observed['plot_calls'] += 1
        observed['plot_out_path'] = kwargs['out_path']

    monkeypatch.setattr('amalgkit.rerun.process_species_busco', fake_process_species_busco)
    monkeypatch.setattr('amalgkit.rerun.generate_busco_species_plot', fake_generate_busco_species_plot)

    rerun_main(_build_args(tmp_path, check='busco', include_warnings=True))

    assert observed['process_calls'] == 0
    assert observed['plot_calls'] == 1
    assert observed['plot_out_path'].endswith('busco/busco_completeness.pdf')


def test_rerun_merge_global_issue_respects_species_filter_and_rebuilds_summary(tmp_path, monkeypatch):
    metadata_dir = tmp_path / 'metadata'
    metadata_dir.mkdir()
    metadata_path = metadata_dir / 'metadata.tsv'
    _write_metadata(metadata_path)

    sanity_dir = tmp_path / 'sanity'
    sanity_dir.mkdir()
    report_path = sanity_dir / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[{'check': 'merge'}],
        issues=[
            {'check': 'merge', 'severity': 'error', 'target_type': 'global', 'target_id': 'metadata'},
        ],
    )

    observed = {}

    def fake_scan_quant_abundance_paths(quant_dir, target_runs=None):
        observed['quant_dir'] = quant_dir
        observed['target_runs'] = sorted(target_runs)
        return {'SRR001': 'abundance1', 'SRR003': 'abundance3'}

    def fake_merge_species_quant_tables(sp, metadata, quant_dir, merge_dir, run_abundance_paths):
        observed.setdefault('species_calls', []).append(sp)
        observed['merge_dir'] = merge_dir
        observed['run_abundance_paths'] = dict(run_abundance_paths)

    def fake_merge_fastp_stats_into_metadata(metadata, out_dir, max_workers=None):
        observed['fastp_out_dir'] = out_dir
        observed['fastp_runs'] = metadata.df['run'].tolist()
        return metadata

    def fake_write_updated_metadata(metadata, path, args, max_workers=None):
        observed['metadata_path'] = path
        observed['written_runs'] = metadata.df['run'].tolist()

    def fake_generate_merge_plot_pdfs(merge_dir, path_metadata_merge):
        observed['plot_merge_dir'] = merge_dir
        observed['plot_metadata_path'] = path_metadata_merge

    monkeypatch.setattr('amalgkit.rerun.scan_quant_abundance_paths', fake_scan_quant_abundance_paths)
    monkeypatch.setattr('amalgkit.rerun.merge_species_quant_tables', fake_merge_species_quant_tables)
    monkeypatch.setattr('amalgkit.rerun.merge_fastp_stats_into_metadata', fake_merge_fastp_stats_into_metadata)
    monkeypatch.setattr('amalgkit.rerun.write_updated_metadata', fake_write_updated_metadata)
    monkeypatch.setattr('amalgkit.rerun.generate_merge_plot_pdfs', fake_generate_merge_plot_pdfs)

    rerun_main(_build_args(tmp_path, check='merge', species='Homo sapiens'))

    assert observed['species_calls'] == ['Homo sapiens']
    assert observed['target_runs'] == ['SRR001', 'SRR002', 'SRR003']
    assert observed['run_abundance_paths'] == {'SRR001': 'abundance1', 'SRR003': 'abundance3'}
    assert os.path.basename(observed['metadata_path']) == 'metadata.tsv'
    assert observed['plot_metadata_path'] == observed['metadata_path']


def test_rerun_merge_subset_keeps_global_plot_context(tmp_path, monkeypatch):
    merge_dir = tmp_path / 'merge'
    for species_token in ['Homo_sapiens', 'Mus_musculus']:
        species_dir = merge_dir / species_token
        species_dir.mkdir(parents=True)
        (species_dir / f'{species_token}_est_counts.tsv').write_text('target_id\tRUN1\ng1\t1\n', encoding='utf-8')
    (merge_dir / 'merge_mapping_rate.pdf').write_text('Homo_sapiens,Mus_musculus', encoding='utf-8')

    metadata_dir = tmp_path / 'metadata'
    metadata_dir.mkdir()
    metadata_path = metadata_dir / 'metadata.tsv'
    _write_metadata(metadata_path)

    sanity_dir = tmp_path / 'sanity'
    sanity_dir.mkdir()
    report_path = sanity_dir / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[{'check': 'merge'}],
        issues=[{'check': 'merge', 'severity': 'error', 'target_type': 'species', 'target_id': 'Mus musculus'}],
    )

    monkeypatch.setattr('amalgkit.rerun.scan_quant_abundance_paths', lambda quant_dir, target_runs=None: {'SRR001': 'a', 'SRR002': 'b', 'SRR003': 'c'})

    def fake_merge_species_quant_tables(sp, metadata, quant_dir, merge_dir, run_abundance_paths):
        token = sp.replace(' ', '_')
        species_dir = os.path.join(merge_dir, token)
        os.makedirs(species_dir, exist_ok=True)
        with open(os.path.join(species_dir, '{}_est_counts.tsv'.format(token)), 'w', encoding='utf-8') as handle:
            handle.write('target_id\tRUN1\ng1\t2\n')

    monkeypatch.setattr('amalgkit.rerun.merge_species_quant_tables', fake_merge_species_quant_tables)
    monkeypatch.setattr('amalgkit.rerun.merge_fastp_stats_into_metadata', lambda metadata, out_dir, max_workers=None: metadata)
    monkeypatch.setattr(
        'amalgkit.rerun.write_updated_metadata',
        lambda metadata, path, args, max_workers=None: metadata.df.to_csv(path, sep='\t', index=False),
    )

    def fake_generate_merge_plot_pdfs(merge_dir, path_metadata_merge):
        species_dirs = sorted([
            name for name in os.listdir(merge_dir)
            if os.path.isdir(os.path.join(merge_dir, name))
        ])
        with open(os.path.join(merge_dir, 'merge_mapping_rate.pdf'), 'w', encoding='utf-8') as handle:
            handle.write(','.join(species_dirs))

    monkeypatch.setattr('amalgkit.rerun.generate_merge_plot_pdfs', fake_generate_merge_plot_pdfs)

    rerun_main(_build_args(tmp_path, check='merge', species='Mus musculus'))

    assert (merge_dir / 'merge_mapping_rate.pdf').read_text(encoding='utf-8') == 'Homo_sapiens,Mus_musculus'


def test_rerun_finalize_global_issue_respects_species_filter(tmp_path, monkeypatch):
    metadata_dir = tmp_path / 'metadata'
    metadata_dir.mkdir()
    metadata_path = metadata_dir / 'metadata.tsv'
    _write_metadata(metadata_path)

    sanity_dir = tmp_path / 'sanity'
    sanity_dir.mkdir()
    report_path = sanity_dir / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[{'check': 'finalize'}],
        issues=[
            {'check': 'finalize', 'severity': 'error', 'target_type': 'global', 'target_id': 'per_species_tables'},
        ],
    )

    observed = {}

    def fake_resolve_per_species_input(args):
        return None, str(tmp_path / 'merge')

    def fake_build_per_species_args(args, input_dir, tmp_out_dir):
        observed['input_dir'] = input_dir
        observed['tmp_out_dir'] = tmp_out_dir
        return SimpleNamespace(batch_effect_alg='none')

    def fake_generate_per_species_tables(args, context):
        observed['selected_species'] = context.metadata.df['scientific_name'].drop_duplicates().tolist()

    def fake_copy_species_tables(per_species_dir, finalize_dir, batch_effect_alg):
        observed['copied_per_species_dir'] = per_species_dir
        observed['finalize_dir'] = finalize_dir
        observed['batch_effect_alg'] = batch_effect_alg

    def fake_copy_per_species_pdfs(per_species_dir, dst_dir, species_subset=None):
        observed['pdf_src_dir'] = per_species_dir
        observed['pdf_dst_dir'] = dst_dir
        observed['pdf_species_subset'] = species_subset

    def fake_load_merged_per_species_metadata(per_species_dir):
        observed['loaded_per_species_dir'] = per_species_dir
        return pandas.DataFrame([{'run': 'SRR002', 'extra_metric': 1}])

    def fake_merge_metadata_by_run(source_df, updated_metadata):
        observed['source_runs'] = source_df['run'].tolist()
        observed['updated_runs'] = updated_metadata['run'].tolist()
        return source_df.loc[source_df['scientific_name'].eq('Mus musculus'), :].copy()

    def fake_save_exclusion_plot_pdf(df_metadata, out_pdf_path, y_label, font_size):
        observed['plot_runs'] = df_metadata['run'].tolist()
        observed['plot_path'] = out_pdf_path

    monkeypatch.setattr('amalgkit.rerun.resolve_per_species_input', fake_resolve_per_species_input)
    monkeypatch.setattr('amalgkit.rerun._build_per_species_args', fake_build_per_species_args)
    monkeypatch.setattr('amalgkit.rerun.generate_per_species_tables', fake_generate_per_species_tables)
    monkeypatch.setattr('amalgkit.rerun._copy_species_tables', fake_copy_species_tables)
    monkeypatch.setattr('amalgkit.rerun.copy_per_species_pdfs', fake_copy_per_species_pdfs)
    monkeypatch.setattr('amalgkit.rerun.load_merged_per_species_metadata', fake_load_merged_per_species_metadata)
    monkeypatch.setattr('amalgkit.rerun.merge_metadata_by_run', fake_merge_metadata_by_run)
    monkeypatch.setattr('amalgkit.rerun.save_exclusion_plot_pdf', fake_save_exclusion_plot_pdf)

    rerun_main(_build_args(tmp_path, check='finalize', species='Mus musculus'))

    assert observed['selected_species'] == ['Mus musculus']
    assert observed['updated_runs'] == ['SRR002']
    assert observed['plot_runs'] == ['SRR002']
    assert os.path.basename(observed['plot_path']) == 'finalize_exclusion.pdf'
    assert (tmp_path / 'finalize' / 'metadata.tsv').exists()


def test_rerun_finalize_global_issue_overrides_species_subset_selection(tmp_path, monkeypatch):
    metadata_dir = tmp_path / 'metadata'
    metadata_dir.mkdir()
    metadata_path = metadata_dir / 'metadata.tsv'
    _write_metadata(metadata_path)

    sanity_dir = tmp_path / 'sanity'
    sanity_dir.mkdir()
    report_path = sanity_dir / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[{'check': 'finalize'}],
        issues=[
            {'check': 'finalize', 'severity': 'error', 'target_type': 'global', 'target_id': 'finalize'},
            {'check': 'finalize', 'severity': 'error', 'target_type': 'species', 'target_id': 'Mus musculus'},
        ],
    )

    observed = {}

    def fake_resolve_per_species_input(args):
        return None, str(tmp_path / 'merge')

    def fake_build_per_species_args(args, input_dir, tmp_out_dir):
        return SimpleNamespace(batch_effect_alg='none')

    def fake_generate_per_species_tables(args, context):
        observed['selected_species'] = context.metadata.df['scientific_name'].drop_duplicates().tolist()

    def fake_copy_species_tables(per_species_dir, finalize_dir, batch_effect_alg):
        observed['finalize_dir'] = finalize_dir

    def fake_copy_per_species_pdfs(per_species_dir, dst_dir, species_subset=None):
        observed['pdf_dst_dir'] = dst_dir
        observed['pdf_species_subset'] = species_subset

    def fake_load_merged_per_species_metadata(per_species_dir):
        return pandas.DataFrame([
            {'run': 'SRR001', 'extra_metric': 1},
            {'run': 'SRR002', 'extra_metric': 2},
            {'run': 'SRR003', 'extra_metric': 3},
        ])

    def fake_merge_metadata_by_run(source_df, updated_metadata):
        return source_df.copy()

    def fake_save_exclusion_plot_pdf(df_metadata, out_pdf_path, y_label, font_size):
        observed['plot_path'] = out_pdf_path

    monkeypatch.setattr('amalgkit.rerun.resolve_per_species_input', fake_resolve_per_species_input)
    monkeypatch.setattr('amalgkit.rerun._build_per_species_args', fake_build_per_species_args)
    monkeypatch.setattr('amalgkit.rerun.generate_per_species_tables', fake_generate_per_species_tables)
    monkeypatch.setattr('amalgkit.rerun._copy_species_tables', fake_copy_species_tables)
    monkeypatch.setattr('amalgkit.rerun.copy_per_species_pdfs', fake_copy_per_species_pdfs)
    monkeypatch.setattr('amalgkit.rerun.load_merged_per_species_metadata', fake_load_merged_per_species_metadata)
    monkeypatch.setattr('amalgkit.rerun.merge_metadata_by_run', fake_merge_metadata_by_run)
    monkeypatch.setattr('amalgkit.rerun.save_exclusion_plot_pdf', fake_save_exclusion_plot_pdf)

    rerun_main(_build_args(tmp_path, check='finalize'))

    assert observed['selected_species'] == ['Homo sapiens', 'Mus musculus']
    assert os.path.basename(observed['plot_path']) == 'finalize_exclusion.pdf'


def test_rerun_dry_run_writes_manifest_without_creating_outputs(tmp_path):
    metadata_dir = tmp_path / 'metadata'
    metadata_dir.mkdir()
    metadata_path = metadata_dir / 'metadata.tsv'
    _write_metadata(metadata_path)

    sanity_dir = tmp_path / 'sanity'
    sanity_dir.mkdir()
    report_path = sanity_dir / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[
            {'check': 'getfastq'},
            {'check': 'merge'},
            {'check': 'busco'},
        ],
        issues=[
            {'check': 'getfastq', 'severity': 'error', 'target_type': 'run', 'target_id': 'SRR002'},
            {'check': 'merge', 'severity': 'error', 'target_type': 'global', 'target_id': 'metadata'},
            {'check': 'merge', 'severity': 'error', 'target_type': 'species', 'target_id': 'Mus musculus'},
            {'check': 'busco', 'severity': 'warning', 'target_type': 'global', 'target_id': 'busco_completeness'},
        ],
    )

    rerun_main(_build_args(tmp_path, dry_run=True, include_warnings=True))

    manifest_path = tmp_path / 'sanity' / 'rerun_manifest.json'
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text())
    assert manifest['dry_run'] is True
    checks = {row['check']: row for row in manifest['checks']}
    assert checks['getfastq']['target_runs'] == ['SRR002']
    assert checks['merge']['target_species'] == ['Homo sapiens', 'Mus musculus']
    assert checks['merge']['selection_basis'] == 'global_issue_all_species'
    assert checks['busco']['rebuild_summary_only'] is True
    assert checks['busco']['selection_basis'] == 'global_issue_summary_rebuild'
    assert not (tmp_path / 'merge').exists()
    assert not (tmp_path / 'busco').exists()


def test_rerun_merge_honors_redo_no_for_existing_outputs(tmp_path, monkeypatch):
    merge_dir = tmp_path / 'merge'
    species_dir = merge_dir / 'Homo_sapiens'
    species_dir.mkdir(parents=True)
    (species_dir / 'Homo_sapiens_est_counts.tsv').write_text('OLD\n', encoding='utf-8')
    (merge_dir / 'merge_mapping_rate.pdf').write_text('OLD\n', encoding='utf-8')

    metadata = pandas.DataFrame([
        {'scientific_name': 'Homo sapiens', 'run': 'SRR001', 'exclusion': 'no'},
    ])
    metadata_path = tmp_path / 'metadata.tsv'
    metadata.to_csv(metadata_path, sep='\t', index=False)
    report_path = tmp_path / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[{'check': 'merge'}],
        issues=[{'check': 'merge', 'severity': 'error', 'target_type': 'species', 'target_id': 'Homo sapiens'}],
    )

    monkeypatch.setattr('amalgkit.rerun.merge_species_quant_tables', lambda *args, **kwargs: pytest.fail('merge rerun should not execute'))

    with pytest.raises(AmalgkitExit, match='--redo is set to "no". Exiting.'):
        rerun_main(_build_args(tmp_path, report=str(report_path), metadata=str(metadata_path), check='merge', redo=False))


def test_rerun_merge_failure_keeps_existing_outputs_unchanged(tmp_path, monkeypatch):
    merge_dir = tmp_path / 'merge'
    species_dir = merge_dir / 'Homo_sapiens'
    species_dir.mkdir(parents=True)
    (species_dir / 'Homo_sapiens_est_counts.tsv').write_text('OLD\n')
    pandas.DataFrame([
        {'run': 'SRR001', 'scientific_name': 'Homo sapiens', 'state': 'old'},
    ]).to_csv(merge_dir / 'metadata.tsv', sep='\t', index=False)

    metadata = pandas.DataFrame([
        {'scientific_name': 'Homo sapiens', 'run': 'SRR001', 'exclusion': 'no'},
    ])
    metadata_path = tmp_path / 'metadata.tsv'
    metadata.to_csv(metadata_path, sep='\t', index=False)
    report_path = tmp_path / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[{'check': 'merge'}],
        issues=[{'check': 'merge', 'severity': 'error', 'target_type': 'species', 'target_id': 'Homo sapiens'}],
    )

    monkeypatch.setattr('amalgkit.rerun.scan_quant_abundance_paths', lambda quant_dir, target_runs=None: {'SRR001': 'abundance'})

    def fake_merge_species_quant_tables(sp, metadata, quant_dir, merge_dir, run_abundance_paths):
        species_dir = os.path.join(merge_dir, 'Homo_sapiens')
        os.makedirs(species_dir, exist_ok=True)
        with open(os.path.join(species_dir, 'Homo_sapiens_est_counts.tsv'), 'w', encoding='utf-8') as handle:
            handle.write('NEW\n')

    monkeypatch.setattr('amalgkit.rerun.merge_species_quant_tables', fake_merge_species_quant_tables)
    monkeypatch.setattr('amalgkit.rerun.merge_fastp_stats_into_metadata', lambda metadata, out_dir, max_workers=None: metadata)
    monkeypatch.setattr('amalgkit.rerun.write_updated_metadata', lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError('boom')))
    monkeypatch.setattr('amalgkit.rerun.generate_merge_plot_pdfs', lambda *args, **kwargs: None)

    with pytest.raises(RuntimeError, match='boom'):
        rerun_main(_build_args(tmp_path, report=str(report_path), metadata=str(metadata_path), check='merge'))

    assert (species_dir / 'Homo_sapiens_est_counts.tsv').read_text() == 'OLD\n'
    assert pandas.read_csv(merge_dir / 'metadata.tsv', sep='\t').loc[0, 'state'] == 'old'


def test_rerun_finalize_failure_keeps_existing_outputs_unchanged(tmp_path, monkeypatch):
    finalize_dir = tmp_path / 'finalize'
    species_dir = finalize_dir / 'Homo_sapiens'
    species_dir.mkdir(parents=True)
    (species_dir / 'Homo_sapiens_expression.tsv').write_text('OLD\n')
    pandas.DataFrame([
        {'run': 'SRR001', 'scientific_name': 'Homo sapiens', 'state': 'old'},
    ]).to_csv(finalize_dir / 'metadata.tsv', sep='\t', index=False)

    metadata = pandas.DataFrame([
        {'scientific_name': 'Homo sapiens', 'run': 'SRR001', 'exclusion': 'no'},
    ])
    metadata_path = tmp_path / 'metadata.tsv'
    metadata.to_csv(metadata_path, sep='\t', index=False)
    report_path = tmp_path / 'sanity_report.json'
    _write_report(
        path_report=report_path,
        out_dir=tmp_path,
        metadata_path=metadata_path,
        summary=[{'check': 'finalize'}],
        issues=[{'check': 'finalize', 'severity': 'error', 'target_type': 'species', 'target_id': 'Homo sapiens'}],
    )

    monkeypatch.setattr('amalgkit.rerun.resolve_per_species_input', lambda args: (None, str(tmp_path / 'merge')))
    monkeypatch.setattr(
        'amalgkit.rerun._build_per_species_args',
        lambda args, input_dir, tmp_out_dir: SimpleNamespace(batch_effect_alg='none', tmp_out_dir=tmp_out_dir),
    )

    def fake_generate_per_species_tables(args, context):
        tables_dir = os.path.join(args.tmp_out_dir, 'per_species', 'Homo_sapiens', 'tables')
        os.makedirs(tables_dir, exist_ok=True)
        with open(os.path.join(tables_dir, 'Homo_sapiens.metadata.tsv'), 'w', encoding='utf-8') as handle:
            handle.write('run\tstate\nSRR001\tnew\n')
        with open(os.path.join(tables_dir, 'Homo_sapiens.none.tc.tsv'), 'w', encoding='utf-8') as handle:
            handle.write('target_id\tSRR001\ng1\t1\n')

    monkeypatch.setattr('amalgkit.rerun.generate_per_species_tables', fake_generate_per_species_tables)
    monkeypatch.setattr(
        'amalgkit.rerun.load_merged_per_species_metadata',
        lambda per_species_dir: (_ for _ in ()).throw(RuntimeError('boom')),
    )

    with pytest.raises(RuntimeError, match='boom'):
        rerun_main(_build_args(tmp_path, report=str(report_path), metadata=str(metadata_path), check='finalize'))

    assert (species_dir / 'Homo_sapiens_expression.tsv').read_text() == 'OLD\n'
    assert pandas.read_csv(finalize_dir / 'metadata.tsv', sep='\t').loc[0, 'state'] == 'old'
