import importlib
import importlib.metadata
import math
import os
import re
import shutil
import subprocess
import sys
import time

from amalgkit.__init__ import __version__
from amalgkit.logging_utils import get_logger
from amalgkit.redaction import redact_url_for_logging

EXPRESSION_NORMALIZATION_METHODS = tuple(
    '{}-{}'.format(log_method, abundance_method)
    for log_method in ('logn', 'log2', 'lognp1', 'log2p1', 'none')
    for abundance_method in ('fpkm', 'tpm', 'none')
)

DEPENDENCY_SPECS = [
    ('numpy', 'numpy'),
    ('pandas', 'pandas'),
    ('scipy', 'scipy'),
    ('matplotlib', 'matplotlib'),
    ('statsmodels', 'statsmodels'),
    ('scikit-learn', 'sklearn'),
    ('biopython', 'Bio'),
    ('defusedxml', 'defusedxml'),
    ('inmoose', 'inmoose'),
]

EXTERNAL_TOOL_SPECS = [
    ('cat', 'cat', [['--version']]),
    ('seqkit', 'seqkit', [['version'], ['--help']]),
    ('kallisto', 'kallisto', [['version'], ['-h']]),
    ('oarfish', 'oarfish', [['--version'], ['-h']]),
    ('fasterq-dump', 'fasterq-dump', [['--version'], ['-h']]),
    ('fastp', 'fastp', [['--version'], ['-h']]),
    ('busco', 'busco', [['--version'], ['-h']]),
    ('compleasm', 'compleasm', [['--version'], ['-h']]),
    ('mmseqs', 'mmseqs', [['version'], ['--help']]),
]
COMMAND_EXTERNAL_TOOL_LABELS = {
    'getfastq': {'cat', 'seqkit', 'fasterq-dump', 'fastp', 'mmseqs'},
    'quant': {'kallisto', 'oarfish'},
    'busco': {'busco', 'compleasm'},
    'integrate': {'seqkit'},
    # rerun can dispatch any of these workflows depending on its manifest.
    'rerun': {
        'cat', 'seqkit', 'fasterq-dump', 'fastp', 'mmseqs',
        'kallisto', 'oarfish', 'busco', 'compleasm',
    },
}
SEMVER_PATTERN = re.compile(r'\b\d+\.\d+(?:\.\d+)*(?:[-+][0-9A-Za-z_.-]+)?\b')
RUNTIME_BANNER_COMMANDS = {
    'metadata',
    'select',
    'getfastq',
    'quant',
    'merge',
    'busco',
    'cstmm',
    'wsfilter',
    'csfilter',
    'finalize',
    'sanity',
    'rerun',
    'integrate',
}


def describe_os():
    if hasattr(os, 'uname'):
        uname = os.uname()
        return '{} {} {} {}'.format(uname.sysname, uname.release, uname.version, uname.machine)
    return sys.platform


def resolve_dependency_version(package_name, module_name):
    _ = module_name  # keep signature stable for call sites
    try:
        return importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        return 'MISSING'
    except Exception as exc:
        return 'UNKNOWN ({})'.format(exc.__class__.__name__)


def iter_nonempty_lines(text):
    if text is None:
        return []
    return [
        line.strip()
        for line in str(text).splitlines()
        if line.strip() != ''
    ]


def run_command_capture(command, timeout=5):
    try:
        out = subprocess.run(  # noqa: S603 - argv comes from fixed internal tool/version probes
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
        )
        return out.returncode, out.stdout, out.stderr, None
    except Exception as exc:
        return None, '', '', exc


def _line_mentions_executable(line, executable_name):
    executable_txt = os.path.basename(str(executable_name)).strip().lower()
    if executable_txt == '':
        return False
    line_lower = str(line).lower()
    boundary_pattern = r'(^|[^A-Za-z0-9_.+-]){}($|[^A-Za-z0-9_.+-])'.format(re.escape(executable_txt))
    return re.search(boundary_pattern, line_lower) is not None


def _score_version_line(line, executable_name):
    line = str(line).strip()
    if line == '':
        return 0
    line_lower = line.lower()
    if line_lower.startswith('usage:'):
        return 0
    has_semver = SEMVER_PATTERN.search(line) is not None
    has_version_word = 'version' in line_lower
    mentions_executable = _line_mentions_executable(line, executable_name)
    if has_version_word and has_semver:
        return 5
    if mentions_executable and has_semver:
        return 4
    if has_version_word and any(ch.isdigit() for ch in line):
        return 3
    if has_semver:
        return 2
    if mentions_executable and any(ch.isdigit() for ch in line):
        return 1
    return 0


def find_version_line(stdout_txt, stderr_txt, executable_name):
    best_line = ''
    best_score = 0
    for text in [stdout_txt, stderr_txt]:
        for line in iter_nonempty_lines(text):
            score = _score_version_line(line, executable_name)
            if score <= best_score:
                continue
            best_line = line
            best_score = score
            if best_score >= 5:
                return best_line
    return best_line


def resolve_external_tool_status(executable_name, version_commands):
    tool_path = shutil.which(executable_name)
    if tool_path is None:
        return 'MISSING'
    for version_args in version_commands:
        returncode, stdout_txt, stderr_txt, exc = run_command_capture([tool_path] + version_args)
        _ = returncode
        if exc is not None:
            continue
        candidate = find_version_line(
            stdout_txt=stdout_txt,
            stderr_txt=stderr_txt,
            executable_name=executable_name,
        )
        if candidate != '':
            return '{} ({})'.format(candidate, tool_path)
    return 'FOUND ({}; version unavailable)'.format(tool_path)


def resolve_external_tool_availability(executable_name):
    """Return a cheap PATH-only status without executing the external tool."""
    tool_path = shutil.which(executable_name)
    if tool_path is None:
        return 'MISSING'
    return 'FOUND ({})'.format(tool_path)


def print_runtime_banner(argv):
    print('AMALGKIT version: {}'.format(__version__))
    print('AMALGKIT command: {}'.format(redact_url_for_logging(' '.join(argv))))
    print('AMALGKIT bug report: https://github.com/kfuku52/amalgkit/issues')
    print('AMALGKIT os: {} (os.name={})'.format(describe_os(), os.name))
    print('AMALGKIT python: {}'.format(sys.version.split()[0]))
    for package_name, module_name in DEPENDENCY_SPECS:
        print(
            'AMALGKIT dependency {}: {}'.format(
                package_name,
                resolve_dependency_version(package_name, module_name),
            )
        )
    active_command = resolve_active_command(argv)
    relevant_labels = COMMAND_EXTERNAL_TOOL_LABELS.get(active_command, set())
    for label, executable_name, _version_commands in EXTERNAL_TOOL_SPECS:
        if label not in relevant_labels:
            continue
        print(
            'AMALGKIT tool {}: {}'.format(
                label,
                resolve_external_tool_availability(executable_name),
            )
        )


def resolve_active_command(argv):
    for arg in argv[1:]:
        if arg.startswith('-'):
            continue
        return arg
    return ''


def should_print_runtime_banner(active_command):
    return active_command in RUNTIME_BANNER_COMMANDS


def strtobool(val):
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")


def int_or_auto(val):
    if isinstance(val, str) and (val.strip().lower() == 'auto'):
        return 'auto'
    int_val = int(val)
    if int_val <= 0:
        raise ValueError('must be > 0 or "auto"')
    return int_val


def nonnegative_int_or_auto(val):
    if isinstance(val, str) and (val.strip().lower() == 'auto'):
        return 'auto'
    int_val = int(val)
    if int_val < 0:
        raise ValueError('must be >= 0 or "auto"')
    return int_val


def positive_float_or_auto(val):
    if isinstance(val, str) and (val.strip().lower() == 'auto'):
        return 'auto'
    float_val = float(val)
    if (not math.isfinite(float_val)) or float_val <= 0:
        raise ValueError('must be > 0 or "auto"')
    return float_val


def build_timed_command_handler(command_name, module_name, function_name):
    def command(args):
        sys.stdout.write('amalgkit {}: start\n'.format(command_name))
        start = time.time()
        logger = get_logger('command')
        logger.info(
            'command_start',
            extra={'event': 'command_start', 'command': command_name},
        )
        try:
            module = importlib.import_module(module_name)
            getattr(module, function_name)(args)
        except Exception:
            setattr(args, '_amalgkit_failure_logged', True)
            logger.exception(
                'command_failed',
                extra={
                    'event': 'command_failed',
                    'command': command_name,
                    'duration_seconds': round(time.time() - start, 6),
                },
            )
            raise
        elapsed = time.time() - start
        logger.info(
            'command_end',
            extra={
                'event': 'command_end',
                'command': command_name,
                'duration_seconds': round(elapsed, 6),
            },
        )
        print('Time elapsed: {:,} sec'.format(int(elapsed)))
        sys.stdout.write('amalgkit {}: end\n'.format(command_name))

    return command


def build_help_command_handler(parser):
    def command_help(args):
        topic = getattr(args, 'topic', None)
        if topic is None:
            parser.print_help()
            return
        parser.parse_args([topic, '--help'])

    return command_help
