import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def run_cli(*args):
    return subprocess.run(
        [sys.executable, '-m', 'amalgkit'] + list(args),
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def test_help_command_exits_zero():
    out = run_cli('help')
    assert out.returncode == 0
    assert 'usage:' in out.stdout.lower()


def test_help_topic_metadata_exits_zero():
    out = run_cli('help', 'metadata')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--search_string' in merged


def test_help_topic_wsfilter_exits_zero():
    out = run_cli('help', 'wsfilter')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--input_dir' in merged
    assert '--margin_threshold' in merged


def test_help_topic_csfilter_exits_zero():
    out = run_cli('help', 'csfilter')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--orthogroup_table' in merged
    assert '--robust_z_threshold' in merged
    assert '--single_copy_threshold' in merged
    assert 'escape a literal comma or pipe' in merged


def test_help_topic_finalize_exits_zero():
    out = run_cli('help', 'finalize')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    compact = ' '.join(merged.split())
    assert '--batch_effect_alg' in merged
    assert '--seed' in merged
    assert '--sva_nsv' in merged
    assert '--sva_b' in merged
    assert '--sva_backend' in merged
    assert '--combatseq_backend' in merged
    assert '--ruvseq_backend' in merged
    assert '--latent_family' in merged
    assert 'latent_glm' in merged
    assert 'default=0: random seed for stochastic steps' in compact
    assert '"auto" uses os entropy and may produce non-reproducible results' in compact


def test_help_topic_cstmm_includes_redo():
    out = run_cli('help', 'cstmm')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--redo' in merged
    assert '--tmm_backend' in merged
    assert '--single_copy_threshold' in merged


def test_single_copy_threshold_rejects_values_outside_percentage_range():
    for command in ['cstmm', 'csfilter']:
        for value in ['0', '101', 'nan']:
            out = run_cli(command, '--single_copy_threshold', value)
            assert out.returncode != 0
            assert 'must be > 0 and <= 100' in out.stderr


def test_help_topic_quant_mentions_backend_specific_index_building():
    out = run_cli('help', 'quant')
    assert out.returncode == 0
    merged = ' '.join((out.stdout + '\n' + out.stderr).lower().split())
    assert '--build_index' in merged
    assert 'kallisto .idx or oarfish .mmi' in merged
    assert '.fa.gz' in merged
    assert '.fasta.gz' in merged


def test_help_rejects_legacy_csca_command():
    out = run_cli('help', 'csca')
    assert out.returncode != 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert 'invalid choice' in merged


def test_help_rejects_legacy_curate_command():
    out = run_cli('help', 'curate')
    assert out.returncode != 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert 'invalid choice' in merged


def test_dataset_list_skips_runtime_banner():
    out = run_cli('dataset', '--list')
    assert out.returncode == 0
    assert 'available datasets:' in out.stdout.lower()
    assert 'amalgkit dependency' not in out.stdout.lower()
    assert 'amalgkit tool' not in out.stdout.lower()


def test_help_topic_getfastq_mentions_filter_runtime_safeguards():
    out = run_cli('help', 'getfastq')
    assert out.returncode == 0
    merged = ' '.join((out.stdout + '\n' + out.stderr).lower().split())
    assert '--rrna_filter_sensitivity' in merged
    assert '--rrna_filter_max_seqs' in merged
    assert '--rrna_filter_jobs' in merged
    assert '--rrna_filter_chunk_spots' in merged
    assert '--rrna_filter_memory_limit' in merged
    assert '--sra_download_wait_timeout_seconds' in merged
    assert '--sra_download_transfer_timeout_seconds' in merged
    assert 'default=1.0' in merged
    assert 'default=20' in merged
    assert 'default=5000000' in merged
    assert 'default=32g' in merged
    assert 'default=86400' in merged
    assert 'default=21600' in merged
    assert '--contam_filter_sensitivity' in merged
    assert '--contam_filter_max_seqs' in merged
    assert 'default=5000000' in merged
    assert 'default=32g' in merged
    assert '32-128 gb ram' in merged
    assert 'default=superkingdom' in merged
    assert 'domain" is accepted as an alias for "superkingdom"' in merged
    assert 'builds both the silva sequence db and its reusable mmseqs search index' in merged
    assert 'first run also downloads/builds the db' in merged
    assert 'higher-priority sources are busy' in merged
    assert 'tries the next enabled source before waiting' in merged


def test_getfastq_rejects_non_positive_download_timeouts():
    for option in [
        '--sra_download_wait_timeout_seconds',
        '--sra_download_transfer_timeout_seconds',
    ]:
        out = run_cli('getfastq', option, '0')
        assert out.returncode != 0
        assert 'must be > 0' in out.stderr


def test_help_topic_integrate_mentions_download_dir():
    out = run_cli('help', 'integrate')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--download_dir' in merged
    assert '--output_metadata' in merged


def test_help_topic_sanity_mentions_new_selection_and_strictness_options():
    out = run_cli('help', 'sanity')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--threads' in merged
    assert '--strict' in merged
    assert '--strict_level' in merged
    assert '--check' in merged
    assert '--run' in merged
    assert '--species' in merged
    assert '--merge' in merged
    assert '--busco' in merged
    assert '--finalize' in merged
    assert '--all is assumed' in merged


def test_help_topic_rerun_mentions_report_and_warning_selection():
    out = run_cli('help', 'rerun')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--report' in merged
    assert '--metadata' in merged
    assert '--check' in merged
    assert '--manifest' in merged
    assert '--include_warnings' in merged
    assert 'sanity_report.json' in merged
    assert 'rerun_manifest.json' in merged


def test_help_topic_dataset_mentions_init_workspace_scaffold():
    out = run_cli('help', 'dataset')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--name' in merged
    assert 'init' in merged
    assert 'workspace scaffold' in merged
