import subprocess
import sys
from types import SimpleNamespace

import pytest

from amalgkit.subprocess_utils import (
    DEPENDENCY_PROBE_TIMEOUT_SECONDS,
    clear_dependency_probe_cache,
    format_command,
    probe_dependency_command,
    resolve_timeout_seconds,
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


@pytest.mark.slow
def test_run_logged_command_applies_timeout_to_real_subprocess():
    # A genuinely slow command is killed after a short timeout instead of
    # blocking forever.
    with pytest.raises(TimeoutError, match='timed out'):
        run_logged_command(
            command=[sys.executable, '-c', 'import time; time.sleep(30)'],
            print_command=False,
            print_output=False,
            timeout_seconds=0.2,
        )


@pytest.mark.parametrize('timeout_stderr', [b'hung stderr', 'hung stderr'])
def test_run_logged_command_timeout_preserves_stderr(timeout_stderr):
    def fake_run(_cmd, stdout=None, stderr=None, timeout=5):
        raise subprocess.TimeoutExpired(cmd=['dummy'], timeout=float(timeout), stderr=timeout_stderr)

    with pytest.raises(TimeoutError, match='hung stderr'):
        run_logged_command(
            command=['dummy'],
            runner=fake_run,
            print_command=False,
            print_output=False,
            timeout_seconds=5,
        )


def test_run_logged_command_falls_back_for_runner_without_timeout_param():
    # Injected runners that accept only (command, stdout, stderr) must keep
    # working when a timeout is configured; the timeout is applied only when
    # the runner supports it.
    calls = {}

    def fake_run(cmd, stdout=None, stderr=None):
        calls['timeout_kwarg'] = None
        return subprocess.CompletedProcess(cmd, 0, stdout=b'ok', stderr=b'')

    result, stdout_txt, _stderr_txt = run_logged_command(
        command=['dummy'],
        runner=fake_run,
        print_command=False,
        print_output=False,
        timeout_seconds=30,
    )

    assert result.returncode == 0
    assert stdout_txt == 'ok'


def test_run_logged_check_call_falls_back_for_runner_without_timeout_param():
    calls = {}

    def fake_check_call(cmd):
        calls['timeout_kwarg'] = None
        return 0

    result = run_logged_check_call(
        command=['dummy'],
        runner=fake_check_call,
        print_command=False,
        timeout_seconds=30,
    )

    assert result == 0


def test_run_logged_command_does_not_retry_runner_internal_type_error():
    calls = []

    def fake_run(cmd, stdout=None, stderr=None, timeout=None):
        calls.append(timeout)
        raise TypeError('timeout calculation failed inside runner')

    with pytest.raises(TypeError, match='timeout calculation failed'):
        run_logged_command(
            command=['dummy'],
            runner=fake_run,
            print_command=False,
            print_output=False,
            timeout_seconds=30,
        )

    assert calls == [30.0]


def test_run_logged_check_call_does_not_retry_runner_internal_type_error():
    calls = []

    def fake_check_call(cmd, timeout=None):
        calls.append(timeout)
        raise TypeError('timeout calculation failed inside runner')

    with pytest.raises(TypeError, match='timeout calculation failed'):
        run_logged_check_call(
            command=['dummy'],
            runner=fake_check_call,
            print_command=False,
            timeout_seconds=30,
        )

    assert calls == [30.0]


# --- Timeout policy resolution -------------------------------------------------

def test_resolve_timeout_seconds_returns_default_without_override():
    args = SimpleNamespace()
    assert resolve_timeout_seconds(args, 'tool_timeout_seconds', 43200) == 43200.0


def test_resolve_timeout_seconds_honours_explicit_override():
    args = SimpleNamespace(tool_timeout_seconds=90)
    assert resolve_timeout_seconds(args, 'tool_timeout_seconds', 43200) == 90.0


def test_resolve_timeout_seconds_zero_disables_the_limit():
    # The escape hatch: a legitimately long job must be able to opt out rather
    # than trade an indefinite hang for an unavoidable kill.
    args = SimpleNamespace(tool_timeout_seconds=0)
    assert resolve_timeout_seconds(args, 'tool_timeout_seconds', 43200) is None


def test_resolve_timeout_seconds_falls_back_on_unusable_values():
    for bad_value in ('not-a-number', float('nan'), float('inf'), None):
        args = SimpleNamespace(tool_timeout_seconds=bad_value)
        assert resolve_timeout_seconds(args, 'tool_timeout_seconds', 120) == 120.0


def test_probe_dependency_command_has_a_short_default_timeout():
    # A hanging `--version` must not block startup.
    assert DEPENDENCY_PROBE_TIMEOUT_SECONDS == 30
    observed = {}

    def fake_run(cmd, stdout=None, stderr=None, timeout=None):
        observed['timeout'] = timeout
        return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

    probe_dependency_command(command=['dummy', '--version'], label='dummy', runner=fake_run)
    assert observed['timeout'] == float(DEPENDENCY_PROBE_TIMEOUT_SECONDS)


def test_successful_default_dependency_probe_is_cached(tmp_path):
    clear_dependency_probe_cache()
    counter_path = tmp_path / 'probe-count.txt'
    script = (
        'from pathlib import Path; '
        'p=Path({!r}); '
        'p.write_text(str(int(p.read_text()) + 1) if p.exists() else "1")'
    ).format(str(counter_path))
    command = [sys.executable, '-c', script]

    probe_dependency_command(command=command, label='python')
    probe_dependency_command(command=command, label='python')

    assert counter_path.read_text(encoding='utf-8') == '1'
    clear_dependency_probe_cache()


@pytest.mark.slow
def test_probe_dependency_command_kills_a_hanging_probe():
    with pytest.raises(TimeoutError, match='timed out'):
        probe_dependency_command(
            command=[sys.executable, '-c', 'import time; time.sleep(30)'],
            label='dummy',
            timeout_seconds=0.2,
        )

def test_format_command_redacts_signed_url_query_and_userinfo():
    url = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1/SRR1.sra?X-Amz-Signature=SECRET&X-Amz-Credential=AKIA'
    redacted = format_command(['curl', '-L', '-o', '/tmp/x.sra', url])
    assert redacted == (
        'curl -L -o /tmp/x.sra '
        'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1/SRR1.sra'
    )
    userinfo = format_command(['curl', 'https://user:hunter2@ftp.sra.ebi.ac.uk/x.sra'])
    assert userinfo == 'curl https://ftp.sra.ebi.ac.uk/x.sra'


def test_format_command_redacts_embedded_malformed_and_unknown_scheme_urls():
    assert format_command([
        'tool',
        '--url=https://user:SECRET@ftp.sra.ebi.ac.uk/x?token=SECRET',
    ]) == 'tool --url=https://ftp.sra.ebi.ac.uk/x'
    assert format_command(['tool', 'https://[broken?token=SECRET']) == (
        'tool https://<redacted>'
    )
    assert format_command(['tool', 'custom://user:SECRET@host/x?token=SECRET']) == (
        'tool custom://host/x'
    )
