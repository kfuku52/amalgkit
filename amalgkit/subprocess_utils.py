import math
import re
import shutil
import subprocess
import threading
import time
import urllib.parse
from typing import Any

from amalgkit.logging_utils import get_logger

# Dependency probes run `--version`/`-h` and should answer within seconds. They
# get their own short default so a hanging executable cannot block startup for
# as long as the analysis-call limits allow.
DEPENDENCY_PROBE_TIMEOUT_SECONDS = 30

_DEFAULT_SUBPROCESS_RUN = subprocess.run
_DEPENDENCY_PROBE_CACHE: dict[tuple[str, tuple[str, ...], float | None], Any] = {}
_DEPENDENCY_PROBE_CACHE_LOCK = threading.Lock()


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


_URL_TOKEN_PATTERN = re.compile(r"[A-Za-z][A-Za-z0-9+.-]*://[^\s]+")


def _redact_url_token(url):
    scheme = url.split("://", 1)[0].casefold()
    try:
        parts = urllib.parse.urlsplit(url)
        host = parts.hostname
        port = parts.port
    except ValueError:
        return "{}://<redacted>".format(scheme)
    if host is None:
        return "{}://<redacted>".format(scheme)
    if ":" in host:
        host = "[{}]".format(host)
    netloc = host if port is None else "{}:{}".format(host, port)
    return urllib.parse.urlunsplit((scheme, netloc, parts.path, "", ""))


def redact_url_for_logging(value):
    """Strip query, fragment, and userinfo from every URL-like command token.

    Signed/requester-pays URLs embed credentials in the query string; only
    scheme+host+path may be logged so tokens never reach stderr or the JSONL
    log. Malformed URLs fail closed instead of returning the secret-bearing
    input unchanged.
    """
    text = str(value)
    return _URL_TOKEN_PATTERN.sub(lambda match: _redact_url_token(match.group(0)), text)


def format_command(command):
    return " ".join([redact_url_for_logging(part) for part in command])


def decode_command_output(output):
    if output is None:
        return ""
    if isinstance(output, str):
        return output
    return output.decode("utf8", errors="replace")


def print_command_output(stdout_txt, stderr_txt, stdout_label=None, stderr_label=None):
    if (stdout_label is None) and (stderr_label is None):
        print(stdout_txt, flush=True)
        print(stderr_txt, flush=True)
        return
    if stdout_label is not None:
        print(stdout_label, flush=True)
        print(stdout_txt, flush=True)
    elif stdout_txt != "":
        print(stdout_txt, flush=True)
    if stderr_label is not None:
        print(stderr_label, flush=True)
        print(stderr_txt, flush=True)
    elif stderr_txt != "":
        print(stderr_txt, flush=True)


def _timeout_error_message(timeout_seconds, command_txt, exc):
    timeout_txt = "" if timeout_seconds is None else f" after {int(float(timeout_seconds)):,} sec"
    stderr_txt = decode_command_output(getattr(exc, "stderr", None))
    return f"Command timed out{timeout_txt}: {command_txt}\n{stderr_txt.rstrip()}"


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
        if "unexpected keyword argument 'timeout'" not in str(exc):
            raise
        return runner(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def run_logged_command(
    command,
    runner=subprocess.run,
    print_command=True,
    command_prefix="Command",
    print_output=False,
    stdout_label=None,
    stderr_label=None,
    not_found_label=None,
    timeout_seconds=None,
):
    command_txt = format_command(command)
    logger = get_logger("subprocess")
    started_at_ns = time.monotonic_ns()
    logger.debug(
        "external_command_start",
        extra={"event": "external_command_start", "command": command_txt},
    )
    if print_command:
        print(f"{command_prefix}: {command_txt}", flush=True)
    try:
        result = _run_command_with_timeout(
            runner=runner,
            command=command,
            timeout_seconds=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        logger.exception(
            "external_command_timeout",
            extra={
                "event": "external_command_timeout",
                "command": command_txt,
                "duration_seconds": round((time.monotonic_ns() - started_at_ns) / 1_000_000_000, 6),
            },
        )
        raise TimeoutError(_timeout_error_message(timeout_seconds, command_txt, exc)) from exc
    except FileNotFoundError as exc:
        logger.exception(
            "external_command_not_found",
            extra={"event": "external_command_not_found", "command": command_txt},
        )
        if not_found_label is None:
            raise
        raise FileNotFoundError(f"{not_found_label} executable not found: {command[0]}") from exc
    stdout_txt = decode_command_output(result.stdout)
    stderr_txt = decode_command_output(result.stderr)
    logger.debug(
        "external_command_end",
        extra={
            "event": "external_command_end",
            "command": command_txt,
            "duration_seconds": round((time.monotonic_ns() - started_at_ns) / 1_000_000_000, 6),
            "returncode": result.returncode,
        },
    )
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
    command_prefix="Command",
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
            message = f"Command failed with exit code {result.returncode}: {command_txt}"
        raise RuntimeError(message)
    return result, stdout_txt, stderr_txt


def run_logged_check_call(
    command,
    runner=subprocess.check_call,
    print_command=True,
    command_prefix="Command",
    not_found_label=None,
    timeout_seconds=None,
):
    command_txt = format_command(command)
    if print_command:
        print(f"{command_prefix}: {command_txt}", flush=True)
    try:
        if timeout_seconds is None:
            return runner(command)
        try:
            return runner(command, timeout=float(timeout_seconds))
        except TypeError as exc:
            if "unexpected keyword argument 'timeout'" not in str(exc):
                raise
            return runner(command)
    except subprocess.TimeoutExpired as exc:
        raise TimeoutError(_timeout_error_message(timeout_seconds, command_txt, exc)) from exc
    except FileNotFoundError as exc:
        if not_found_label is None:
            raise
        raise FileNotFoundError(f"{not_found_label} executable not found: {command[0]}") from exc


def probe_dependency_command(
    command,
    label,
    runner=subprocess.run,
    timeout_seconds=DEPENDENCY_PROBE_TIMEOUT_SECONDS,
):
    cache_key = None
    if runner is _DEFAULT_SUBPROCESS_RUN:
        executable = shutil.which(str(command[0])) or str(command[0])
        cache_key = (
            executable,
            tuple(str(part) for part in command[1:]),
            None if timeout_seconds is None else float(timeout_seconds),
        )
        with _DEPENDENCY_PROBE_CACHE_LOCK:
            cached = _DEPENDENCY_PROBE_CACHE.get(cache_key)
        if cached is not None:
            get_logger("subprocess").debug(
                "dependency_probe_cache_hit",
                extra={"event": "dependency_probe_cache_hit", "command": format_command(command)},
            )
            return cached

    result = run_checked_command(
        command=command,
        runner=runner,
        print_command=False,
        print_output=False,
        not_found_label=label,
        timeout_seconds=timeout_seconds,
        failure_message=lambda result, _stdout, _stderr, command_txt: (
            f"{label} dependency probe failed with exit code {result.returncode}: {command_txt}"
        ),
    )
    if cache_key is not None:
        with _DEPENDENCY_PROBE_CACHE_LOCK:
            _DEPENDENCY_PROBE_CACHE[cache_key] = result
    return result


def clear_dependency_probe_cache():
    """Clear successful dependency-probe results, primarily for long-lived callers."""
    with _DEPENDENCY_PROBE_CACHE_LOCK:
        _DEPENDENCY_PROBE_CACHE.clear()
