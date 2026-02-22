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
