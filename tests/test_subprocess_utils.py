import subprocess

import pytest

from amalgkit.subprocess_utils import (
    probe_dependency_command,
    run_checked_command,
    run_logged_check_call,
    run_logged_command,
)


def test_run_logged_command_decodes_non_utf8_output():
    result, stdout_txt, stderr_txt = run_logged_command(
        command=['dummy'],
        runner=lambda *_args, **_kwargs: subprocess.CompletedProcess(
            args=['dummy'],
            returncode=0,
            stdout=b'\xff',
            stderr=b'\xfe',
        ),
        print_command=False,
        print_output=False,
    )

    assert result.returncode == 0
    assert stdout_txt != ''
    assert stderr_txt != ''


def test_probe_dependency_command_reports_missing_executable():
    def fake_run(_cmd, stdout=None, stderr=None):
        raise FileNotFoundError('missing')

    with pytest.raises(FileNotFoundError, match='dummy executable not found'):
        probe_dependency_command(['dummy', '--help'], 'dummy', runner=fake_run)


def test_run_checked_command_uses_custom_failure_message():
    def fake_run(cmd, stdout=None, stderr=None):
        return subprocess.CompletedProcess(cmd, 2, stdout=b'', stderr=b'bad')

    with pytest.raises(RuntimeError, match='custom failure'):
        run_checked_command(
            command=['dummy'],
            runner=fake_run,
            print_command=False,
            failure_message=lambda result, _stdout, stderr, command_txt: (
                'custom failure {} {} {}'.format(result.returncode, command_txt, stderr)
            ),
        )


def test_run_logged_check_call_reports_missing_executable():
    def fake_check_call(_cmd):
        raise FileNotFoundError('missing')

    with pytest.raises(FileNotFoundError, match='dummy executable not found'):
        run_logged_check_call(
            command=['dummy', '--help'],
            runner=fake_check_call,
            print_command=False,
            not_found_label='dummy',
        )
