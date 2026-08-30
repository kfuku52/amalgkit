from types import SimpleNamespace
import json
import subprocess

import pytest

import amalgkit.main as main_module
from amalgkit.logging_utils import configure_logging
from amalgkit.subprocess_utils import run_logged_command


class _FailingParser:
    def __init__(self, debug):
        self.debug = debug

    def parse_args(self, _argv):
        def fail(_args):
            raise ValueError('diagnostic failure')

        return SimpleNamespace(debug=self.debug, handler=fail)

    def print_help(self):
        raise AssertionError('handler should have been called')


def test_main_debug_mode_prints_traceback(monkeypatch, capsys):
    monkeypatch.setattr(main_module, 'build_main_parser', lambda: _FailingParser(debug=True))

    assert main_module.main(['amalgkit', 'fake']) == 1

    stderr = capsys.readouterr().err
    assert 'Traceback (most recent call last)' in stderr
    assert 'ValueError: diagnostic failure' in stderr


def test_main_default_error_output_remains_concise(monkeypatch, capsys):
    monkeypatch.setattr(main_module, 'build_main_parser', lambda: _FailingParser(debug=False))

    assert main_module.main(['amalgkit', 'fake']) == 1

    stderr = capsys.readouterr().err
    assert stderr == 'ERROR: diagnostic failure\n'


@pytest.mark.parametrize('debug', [False, True])
def test_timeout_secrets_never_reach_cli_or_jsonl(tmp_path, monkeypatch, capsys, debug):
    log_path = tmp_path / 'log.jsonl'
    url = 'https://USER_SECRET:PASS_SECRET@example.org/data?token=QUERY_SECRET#FRAGMENT_SECRET'

    def runner(command, **_kwargs):
        # Execution still receives the complete, usable URL.
        assert command[-1] == url
        raise subprocess.TimeoutExpired(command, 2, stderr=('fetch failed: ' + url).encode())

    def handler(_args):
        run_logged_command(['curl', url], runner=runner, timeout_seconds=2)

    args = SimpleNamespace(debug=debug, log_level='DEBUG', log_file=str(log_path), handler=handler)
    monkeypatch.setattr(main_module, 'build_main_parser', lambda: SimpleNamespace(parse_args=lambda _argv: args))
    try:
        assert main_module.main(['amalgkit', 'fake']) == 1
        captured = capsys.readouterr()
        all_output = captured.out + captured.err + log_path.read_text(encoding='utf-8')
        for secret in ('USER_SECRET', 'PASS_SECRET', 'QUERY_SECRET', 'FRAGMENT_SECRET'):
            assert secret not in all_output
        assert 'https://example.org/data' in all_output
        records = [json.loads(line) for line in log_path.read_text(encoding='utf-8').splitlines()]
        assert any('TimeoutExpired' in record.get('exception', '') for record in records)
        if debug:
            assert 'Traceback (most recent call last)' in captured.err
    finally:
        configure_logging()
