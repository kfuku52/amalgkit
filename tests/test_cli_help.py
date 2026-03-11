import subprocess
import sys
from pathlib import Path


CLI_PATH = Path(__file__).resolve().parents[1] / 'amalgkit' / 'amalgkit'


def run_cli(*args):
    return subprocess.run(
        [sys.executable, str(CLI_PATH)] + list(args),
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


def test_help_topic_finalize_exits_zero():
    out = run_cli('help', 'finalize')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--batch_effect_alg' in merged
    assert '--seed' in merged
    assert '--sva_nsv' in merged
    assert '--sva_b' in merged
    assert '--sva_backend' in merged
    assert '--combatseq_backend' in merged
    assert '--ruvseq_backend' in merged
    assert '--latent_family' in merged
    assert 'latent_glm' in merged


def test_help_topic_cstmm_includes_redo():
    out = run_cli('help', 'cstmm')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--redo' in merged
    assert '--tmm_backend' in merged


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


def test_help_topic_getfastq_mentions_filter_runtime_costs():
    out = run_cli('help', 'getfastq')
    assert out.returncode == 0
    merged = ' '.join((out.stdout + '\n' + out.stderr).lower().split())
    assert '--rrna_filter_sensitivity' in merged
    assert '--rrna_filter_max_seqs' in merged
    assert 'default=1.0' in merged
    assert 'default=20' in merged
    assert '--contam_filter_sensitivity' in merged
    assert '--contam_filter_max_seqs' in merged
    assert '2-8 gb ram' in merged
    assert '32-128 gb ram' in merged
    assert 'default=superkingdom' in merged
    assert 'domain" is accepted as an alias for "superkingdom"' in merged
    assert 'first run also builds the silva db' in merged
    assert 'first run also downloads/builds the db' in merged


def test_help_topic_integrate_mentions_download_dir():
    out = run_cli('help', 'integrate')
    assert out.returncode == 0
    merged = (out.stdout + '\n' + out.stderr).lower()
    assert '--download_dir' in merged
