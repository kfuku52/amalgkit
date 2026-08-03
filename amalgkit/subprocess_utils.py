import math
import subprocess


# Dependency probes run `--version`/`-h` and should answer within seconds. They
# get their own short default so a hanging executable cannot block startup for
# as long as the analysis-call limits allow.
DEPENDENCY_PROBE_TIMEOUT_SECONDS = 300


def resolve_timeout_seconds(args, attribute_name, default_seconds):
    """Resolve a wall-clock timeout in seconds from parsed CLI arguments.

    Returns ``default_seconds`` when the attribute is absent, ``None``, or not a
    finite number. An explicit value of ``0`` (or any negative value) disables
    the timeout and returns ``None``, so a legitimately long-running job can opt
    out of the limit rather than being killed by it.
    """
    default_value = None if default_seconds is None else float(default_seconds)
    raw_value = None if args is None else getattr(args, attribute_name, None)
    if raw_value is None:
        return default_value
    try:
        timeout_seconds = float(raw_value)
    except (TypeError, ValueError):
        return default_value
    if not math.isfinite(timeout_seconds):
        return default_value
    if timeout_seconds <= 0:
        return None
    return timeout_seconds


def format_command(command):
    return ' '.join([str(part) for part in command])


def decode_command_output(output):
    if output is None:
        return ''
    if isinstance(output, str):
        return output
    return output.decode('utf8', errors='replace')


def print_command_output(stdout_txt, stderr_txt, stdout_label=None, stderr_label=None):
    if (stdout_label is None) and (stderr_label is None):
        print(stdout_txt, flush=True)
        print(stderr_txt, flush=True)
        return
    if stdout_label is not None:
        print(stdout_label, flush=True)
        print(stdout_txt, flush=True)
    elif stdout_txt != '':
        print(stdout_txt, flush=True)
    if stderr_label is not None:
        print(stderr_label, flush=True)
        print(stderr_txt, flush=True)
    elif stderr_txt != '':
        print(stderr_txt, flush=True)


def _timeout_error_message(timeout_seconds, command_txt, exc):
    timeout_txt = '' if timeout_seconds is None else ' after {:,} sec'.format(int(float(timeout_seconds)))
    stderr_txt = getattr(exc, 'stderr', None)
    if stderr_txt is None:
        stderr_txt = b''
    return 'Command timed out{}: {}\n{}'.format(
        timeout_txt,
        command_txt,
        stderr_txt.decode('utf8', errors='replace').rstrip(),
    )


def _run_command_with_timeout(runner, command, timeout_seconds):
    if timeout_seconds is None:
        return runner(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        return runner(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=float(timeout_seconds),
        )
    except TypeError as exc:
        # Injected test runners (and some thin wrappers) accept only
        # (command, stdout, stderr); fall back without the timeout so the
        # timeout capability remains opt-in for real subprocess calls.
        if 'timeout' not in str(exc):
            raise
        return runner(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def run_logged_command(
    command,
    runner=subprocess.run,
    print_command=True,
    command_prefix='Command',
    print_output=False,
    stdout_label=None,
    stderr_label=None,
    not_found_label=None,
    timeout_seconds=None,
):
    command_txt = format_command(command)
    if print_command:
        print('{}: {}'.format(command_prefix, command_txt), flush=True)
    try:
        result = _run_command_with_timeout(
            runner=runner,
            command=command,
            timeout_seconds=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        raise TimeoutError(
            _timeout_error_message(timeout_seconds, command_txt, exc)
        ) from exc
    except FileNotFoundError as exc:
        if not_found_label is None:
            raise
        raise FileNotFoundError(
            '{} executable not found: {}'.format(not_found_label, command[0])
        ) from exc
    stdout_txt = decode_command_output(result.stdout)
    stderr_txt = decode_command_output(result.stderr)
    if print_output:
        print_command_output(
            stdout_txt=stdout_txt,
            stderr_txt=stderr_txt,
            stdout_label=stdout_label,
            stderr_label=stderr_label,
        )
    return result, stdout_txt, stderr_txt


def run_checked_command(
    command,
    runner=subprocess.run,
    print_command=True,
    command_prefix='Command',
    print_output=False,
    stdout_label=None,
    stderr_label=None,
    not_found_label=None,
    failure_message=None,
    timeout_seconds=None,
):
    result, stdout_txt, stderr_txt = run_logged_command(
        command=command,
        runner=runner,
        print_command=print_command,
        command_prefix=command_prefix,
        print_output=print_output,
        stdout_label=stdout_label,
        stderr_label=stderr_label,
        not_found_label=not_found_label,
        timeout_seconds=timeout_seconds,
    )
    if result.returncode != 0:
        command_txt = format_command(command)
        if callable(failure_message):
            message = failure_message(result, stdout_txt, stderr_txt, command_txt)
        elif failure_message is not None:
            message = failure_message
        else:
            message = 'Command failed with exit code {}: {}'.format(result.returncode, command_txt)
        raise RuntimeError(message)
    return result, stdout_txt, stderr_txt


def run_logged_check_call(
    command,
    runner=subprocess.check_call,
    print_command=True,
    command_prefix='Command',
    not_found_label=None,
    timeout_seconds=None,
):
    command_txt = format_command(command)
    if print_command:
        print('{}: {}'.format(command_prefix, command_txt), flush=True)
    try:
        if timeout_seconds is None:
            return runner(command)
        try:
            return runner(command, timeout=float(timeout_seconds))
        except TypeError as exc:
            if 'timeout' not in str(exc):
                raise
            return runner(command)
    except subprocess.TimeoutExpired as exc:
        raise TimeoutError(
            _timeout_error_message(timeout_seconds, command_txt, exc)
        ) from exc
    except FileNotFoundError as exc:
        if not_found_label is None:
            raise
        raise FileNotFoundError(
            '{} executable not found: {}'.format(not_found_label, command[0])
        ) from exc


def probe_dependency_command(
    command,
    label,
    runner=subprocess.run,
    timeout_seconds=DEPENDENCY_PROBE_TIMEOUT_SECONDS,
):
    return run_checked_command(
        command=command,
        runner=runner,
        print_command=False,
        print_output=False,
        not_found_label=label,
        timeout_seconds=timeout_seconds,
        failure_message=lambda result, _stdout, _stderr, command_txt: (
            '{} dependency probe failed with exit code {}: {}'.format(
                label,
                result.returncode,
                command_txt,
            )
        ),
    )
