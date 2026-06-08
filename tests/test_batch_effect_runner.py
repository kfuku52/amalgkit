import subprocess
import sys

import pandas

from amalgkit import batch_effect_runner


def test_batch_effect_runner_writes_combatseq_outputs(tmp_path, capsys, monkeypatch):
    counts_path = tmp_path / 'counts.tsv'
    metadata_path = tmp_path / 'metadata.tsv'
    summary_path = tmp_path / 'summary.json'
    corrected_path = tmp_path / 'corrected.tsv'

    pandas.DataFrame({
        'target_id': ['G1', 'G2'],
        'RUN1': [1.0, 2.0],
        'RUN2': [3.0, 4.0],
    }).to_csv(counts_path, sep='\t', index=False)
    pandas.DataFrame({
        'run': ['RUN1', 'RUN2'],
        'sample_group': ['A', 'B'],
        'bioproject': ['BP1', 'BP1'],
    }).to_csv(metadata_path, sep='\t', index=False)

    def fake_run_combatseq_backend(counts_df, metadata_df, batch_column, sample_group_column):
        assert batch_column == 'bioproject'
        assert sample_group_column == 'sample_group'
        corrected = counts_df.copy()
        corrected.loc[:, 'RUN1'] = [10.0, 20.0]
        corrected.loc[:, 'RUN2'] = [30.0, 40.0]
        summary = {
            'backend': 'combatseq',
            'method': 'group',
            'skip_reason': '',
            'corrected_run_ids': ['RUN1', 'RUN2'],
            'uncorrected_run_ids': [],
        }
        return corrected, summary

    monkeypatch.setattr(batch_effect_runner, 'run_combatseq_backend', fake_run_combatseq_backend)

    code = batch_effect_runner.main([
        '--backend', 'combatseq',
        '--counts_tsv', str(counts_path),
        '--metadata_tsv', str(metadata_path),
        '--out_summary_json', str(summary_path),
        '--out_counts_tsv', str(corrected_path),
    ])
    captured = capsys.readouterr()

    assert code == 0
    assert captured.err == ''
    payload = batch_effect_runner.read_backend_summary_json(summary_path)
    assert payload['backend'] == 'combatseq'
    assert payload['method'] == 'group'
    assert payload['skip_reason'] == ''
    corrected = pandas.read_csv(corrected_path, sep='\t', index_col=0)
    assert corrected.loc['G1', 'RUN1'] == 10.0
    assert corrected.loc['G2', 'RUN2'] == 40.0


def test_batch_effect_runner_reports_combatseq_dependency_error(tmp_path, capsys, monkeypatch):
    counts_path = tmp_path / 'counts.tsv'
    metadata_path = tmp_path / 'metadata.tsv'
    summary_path = tmp_path / 'summary.json'

    pandas.DataFrame({
        'target_id': ['G1', 'G2'],
        'RUN1': [1.0, 2.0],
    }).to_csv(counts_path, sep='\t', index=False)
    pandas.DataFrame({
        'run': ['RUN1'],
        'sample_group': ['A'],
        'bioproject': ['BP1'],
    }).to_csv(metadata_path, sep='\t', index=False)

    def fake_run_combatseq_backend(*args, **kwargs):
        raise ImportError('inmoose is required')

    monkeypatch.setattr(batch_effect_runner, 'run_combatseq_backend', fake_run_combatseq_backend)

    code = batch_effect_runner.main([
        '--backend', 'combatseq',
        '--counts_tsv', str(counts_path),
        '--metadata_tsv', str(metadata_path),
        '--out_summary_json', str(summary_path),
    ])
    captured = capsys.readouterr()

    assert code == 1
    assert 'inmoose is required' in captured.err
    payload = batch_effect_runner.read_backend_summary_json(summary_path)
    assert payload['backend'] == 'combatseq'
    assert payload['method'] == 'error'
    assert payload['skip_reason'] == 'backend_dependency_missing'


def test_batch_effect_runner_writes_ruvseq_outputs(tmp_path, capsys, monkeypatch):
    counts_path = tmp_path / 'counts.tsv'
    metadata_path = tmp_path / 'metadata.tsv'
    summary_path = tmp_path / 'summary.json'
    corrected_path = tmp_path / 'corrected.tsv'
    w_path = tmp_path / 'w.tsv'

    pandas.DataFrame({
        'target_id': ['G1', 'G2'],
        'RUN1': [1.0, 2.0],
        'RUN2': [3.0, 4.0],
    }).to_csv(counts_path, sep='\t', index=False)
    pandas.DataFrame({
        'run': ['RUN1', 'RUN2'],
        'sample_group': ['A', 'B'],
        'bioproject': ['BP1', 'BP2'],
    }).to_csv(metadata_path, sep='\t', index=False)

    def fake_run_ruvseq_backend(counts_df, metadata_df, control_mode, k_setting, k_max, top_n, min_controls, batch_column, sample_group_column):
        assert control_mode == 'auto'
        assert k_setting == '1'
        assert str(k_max) == '3'
        assert str(top_n) == '10'
        assert str(min_controls) == '2'
        assert batch_column == 'bioproject'
        assert sample_group_column == 'sample_group'
        corrected = counts_df.copy()
        corrected.loc[:, 'RUN1'] = [10.0, 20.0]
        corrected.loc[:, 'RUN2'] = [30.0, 40.0]
        w_df = pandas.DataFrame({'W_1': [0.1, -0.1]}, index=['RUN1', 'RUN2'])
        summary = {
            'backend': 'ruvseq',
            'method': 'manual',
            'skip_reason': '',
            'corrected_run_ids': ['RUN1', 'RUN2'],
            'uncorrected_run_ids': [],
            'resolved_ruv_k': 1,
            'resolved_ruv_controls': 2,
        }
        return corrected, w_df, summary

    monkeypatch.setattr(batch_effect_runner, 'run_ruvseq_backend', fake_run_ruvseq_backend)

    code = batch_effect_runner.main([
        '--backend', 'ruvseq',
        '--counts_tsv', str(counts_path),
        '--metadata_tsv', str(metadata_path),
        '--out_summary_json', str(summary_path),
        '--out_counts_tsv', str(corrected_path),
        '--out_sv_tsv', str(w_path),
        '--ruvseq_k', '1',
        '--ruvseq_k_max', '3',
        '--ruvseq_control_top_n', '10',
        '--ruvseq_min_controls', '2',
    ])
    captured = capsys.readouterr()

    assert code == 0
    assert captured.err == ''
    payload = batch_effect_runner.read_backend_summary_json(summary_path)
    assert payload['backend'] == 'ruvseq'
    assert payload['resolved_ruv_k'] == 1
    corrected = pandas.read_csv(corrected_path, sep='\t', index_col=0)
    assert corrected.loc['G1', 'RUN1'] == 10.0
    w_df = pandas.read_csv(w_path, sep='\t', index_col=0)
    assert list(w_df.columns) == ['W_1']


def test_batch_effect_runner_supports_sva_zero_case_and_writes_outputs(tmp_path, capsys):
    counts_path = tmp_path / 'counts.tsv'
    metadata_path = tmp_path / 'metadata.tsv'
    summary_path = tmp_path / 'summary.dcf'
    corrected_path = tmp_path / 'corrected.tsv'
    sv_path = tmp_path / 'sv.tsv'

    pandas.DataFrame({
        'target_id': ['G1', 'G2'],
        'RUN1': [1.0, 2.0],
        'RUN2': [3.0, 4.0],
    }).to_csv(counts_path, sep='\t', index=False)
    pandas.DataFrame({
        'run': ['RUN1', 'RUN2'],
        'sample_group': ['A', 'B'],
    }).to_csv(metadata_path, sep='\t', index=False)

    code = batch_effect_runner.main([
        '--backend', 'sva',
        '--counts_tsv', str(counts_path),
        '--metadata_tsv', str(metadata_path),
        '--out_counts_tsv', str(corrected_path),
        '--out_summary_dcf', str(summary_path),
        '--out_sv_tsv', str(sv_path),
        '--sva_nsv', '0',
    ])
    captured = capsys.readouterr()

    assert code == 0
    assert captured.err == ''
    summary = batch_effect_runner.read_backend_summary_dcf(summary_path)
    assert summary['resolved_sva_nsv'] == '0'
    assert summary['skip_reason'] == 'sva_nsv_zero'
    corrected = pandas.read_csv(corrected_path, sep='\t', index_col=0)
    assert corrected.shape == (2, 2)
    sv = pandas.read_csv(sv_path, sep='\t', index_col=0)
    assert sv.shape == (2, 0)


def test_batch_effect_runner_module_entrypoint_executes_main(tmp_path):
    counts_path = tmp_path / 'counts.tsv'
    metadata_path = tmp_path / 'metadata.tsv'
    summary_path = tmp_path / 'summary.dcf'
    corrected_path = tmp_path / 'corrected.tsv'
    sv_path = tmp_path / 'sv.tsv'

    pandas.DataFrame({
        'target_id': ['G1', 'G2'],
        'RUN1': [1.0, 2.0],
        'RUN2': [3.0, 4.0],
    }).to_csv(counts_path, sep='\t', index=False)
    pandas.DataFrame({
        'run': ['RUN1', 'RUN2'],
        'sample_group': ['A', 'B'],
    }).to_csv(metadata_path, sep='\t', index=False)

    proc = subprocess.run(
        [
            sys.executable,
            '-m',
            'amalgkit.batch_effect_runner',
            '--backend',
            'sva',
            '--counts_tsv',
            str(counts_path),
            '--metadata_tsv',
            str(metadata_path),
            '--out_counts_tsv',
            str(corrected_path),
            '--out_summary_dcf',
            str(summary_path),
            '--out_sv_tsv',
            str(sv_path),
            '--sva_nsv',
            '0',
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert proc.returncode == 0, proc.stderr
    assert corrected_path.exists()
    assert summary_path.exists()
    assert sv_path.exists()


def test_batch_effect_runner_writes_latent_glm_outputs(tmp_path, capsys, monkeypatch):
    counts_path = tmp_path / 'counts.tsv'
    metadata_path = tmp_path / 'metadata.tsv'
    summary_path = tmp_path / 'summary.json'
    corrected_path = tmp_path / 'corrected.tsv'
    latent_path = tmp_path / 'latent.tsv'

    pandas.DataFrame({
        'target_id': ['G1', 'G2'],
        'RUN1': [10.0, 20.0],
        'RUN2': [30.0, 40.0],
    }).to_csv(counts_path, sep='\t', index=False)
    pandas.DataFrame({
        'run': ['RUN1', 'RUN2'],
        'sample_group': ['A', 'B'],
        'bioproject': ['BP1', 'BP2'],
    }).to_csv(metadata_path, sep='\t', index=False)

    def fake_run_latent_glm_backend(counts_df, metadata_df, family, k_setting, k_max, sample_group_column, max_iter, tol):
        assert family == 'poisson'
        assert k_setting == '1'
        assert str(k_max) == '4'
        assert sample_group_column == 'sample_group'
        assert str(max_iter) == '111'
        assert str(tol) == '0.0001'
        corrected = counts_df.copy()
        corrected.loc[:, 'RUN1'] = [11.0, 21.0]
        corrected.loc[:, 'RUN2'] = [31.0, 41.0]
        latent_df = pandas.DataFrame({'latent_1': [0.1, -0.1]}, index=['RUN1', 'RUN2'])
        summary = {
            'backend': 'latent_glm',
            'method': 'manual',
            'skip_reason': '',
            'resolved_latent_k': 1,
            'latent_family': 'poisson',
            'corrected_run_ids': ['RUN1', 'RUN2'],
            'uncorrected_run_ids': [],
        }
        return corrected, latent_df, summary

    monkeypatch.setattr(batch_effect_runner, 'run_latent_glm_backend', fake_run_latent_glm_backend)

    code = batch_effect_runner.main([
        '--backend', 'latent_glm',
        '--counts_tsv', str(counts_path),
        '--metadata_tsv', str(metadata_path),
        '--out_summary_json', str(summary_path),
        '--out_counts_tsv', str(corrected_path),
        '--out_sv_tsv', str(latent_path),
        '--latent_family', 'poisson',
        '--latent_k', '1',
        '--latent_k_max', '4',
        '--latent_max_iter', '111',
        '--latent_tol', '0.0001',
    ])
    captured = capsys.readouterr()

    assert code == 0
    assert captured.err == ''
    payload = batch_effect_runner.read_backend_summary_json(summary_path)
    assert payload['backend'] == 'latent_glm'
    assert payload['resolved_latent_k'] == 1
    corrected = pandas.read_csv(corrected_path, sep='\t', index_col=0)
    assert corrected.loc['G1', 'RUN1'] == 11.0
    latent_df = pandas.read_csv(latent_path, sep='\t', index_col=0)
    assert list(latent_df.columns) == ['latent_1']
