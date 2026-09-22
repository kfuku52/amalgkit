import subprocess
import sys
from pathlib import Path

import pytest


pytestmark = pytest.mark.slow


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


def test_root_help_documents_debug_tracebacks():
    out = run_cli('--help')
    assert out.returncode == 0
    assert '--debug' in out.stdout
    assert '--log_level' in out.stdout
    assert '--log_file' in out.stdout
    assert 'full traceback' in out.stdout.lower()


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


def test_dataset_list_skips_runtime_banner():
    out = run_cli('dataset', '--list')
    assert out.returncode == 0
    assert 'available datasets:' in out.stdout.lower()
    assert 'amalgkit dependency' not in out.stdout.lower()
    assert 'amalgkit tool' not in out.stdout.lower()
