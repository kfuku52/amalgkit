import pytest
from types import SimpleNamespace
from unittest.mock import Mock

import amalgkit.cli_utils as cli_utils


def test_resolve_external_tool_status_ignores_usage_and_extracts_busco_version(monkeypatch):
    monkeypatch.setattr(cli_utils.shutil, 'which', lambda _exe: '/opt/conda/bin/busco')

    def fake_run(command, timeout=5):
        _ = timeout
        args = tuple(command[1:])
        if args == ('--version',):
            return 0, 'usage: busco -i [SEQUENCE_FILE] -l [LINEAGE]\n', '', None
        if args == ('-h',):
            return (
                0,
                'usage: busco -i [SEQUENCE_FILE] -l [LINEAGE]\n'
                '\n'
                'Welcome to BUSCO 6.0.0: the Benchmarking Universal Single-Copy Ortholog assessment tool.\n',
                '',
                None,
            )
        raise AssertionError('Unexpected command: {}'.format(command))

    monkeypatch.setattr(cli_utils, 'run_command_capture', fake_run)

    status = cli_utils.resolve_external_tool_status('busco', [['--version'], ['-h']])

    assert status == (
        'Welcome to BUSCO 6.0.0: the Benchmarking Universal Single-Copy Ortholog assessment tool. '
        '(/opt/conda/bin/busco)'
    )


def test_resolve_external_tool_status_accepts_plain_semver_output(monkeypatch):
    monkeypatch.setattr(cli_utils.shutil, 'which', lambda _exe: '/opt/conda/bin/fasterq-dump')
    monkeypatch.setattr(
        cli_utils,
        'run_command_capture',
        lambda command, timeout=5: (0, '2.9.6\n', '', None),
    )

    status = cli_utils.resolve_external_tool_status('fasterq-dump', [['--version']])

    assert status == '2.9.6 (/opt/conda/bin/fasterq-dump)'


def test_runtime_banner_reports_all_scientific_python_dependencies(monkeypatch, capsys):
    monkeypatch.setattr(
        cli_utils,
        'resolve_dependency_version',
        lambda package_name, _module_name: 'version-for-{}'.format(package_name),
    )
    monkeypatch.setattr(cli_utils, 'resolve_external_tool_availability', lambda *_args: 'MISSING')

    cli_utils.print_runtime_banner(['amalgkit', 'finalize'])

    output = capsys.readouterr().out
    for package_name in (
        'numpy',
        'pandas',
        'scipy',
        'matplotlib',
        'statsmodels',
        'scikit-learn',
        'biopython',
        'defusedxml',
        'inmoose',
    ):
        assert 'AMALGKIT dependency {}: version-for-{}'.format(package_name, package_name) in output


def test_runtime_banner_skips_external_probes_for_metadata_only_command(monkeypatch, capsys):
    monkeypatch.setattr(cli_utils, 'resolve_dependency_version', lambda *_args: '1.0')

    def fail_external_probe(*_args):
        raise AssertionError('metadata must not execute unrelated external tools')

    monkeypatch.setattr(cli_utils, 'resolve_external_tool_availability', fail_external_probe)

    cli_utils.print_runtime_banner(['amalgkit', 'metadata'])

    assert 'AMALGKIT tool' not in capsys.readouterr().out


def test_runtime_banner_probes_only_tools_relevant_to_active_command(monkeypatch, capsys):
    monkeypatch.setattr(cli_utils, 'resolve_dependency_version', lambda *_args: '1.0')
    observed = []

    def record_probe(executable_name):
        observed.append(executable_name)
        return 'MISSING'

    monkeypatch.setattr(cli_utils, 'resolve_external_tool_availability', record_probe)

    cli_utils.print_runtime_banner(['amalgkit', 'quant'])

    assert observed == ['kallisto', 'oarfish']
    output = capsys.readouterr().out
    assert 'AMALGKIT tool kallisto:' in output
    assert 'AMALGKIT tool fastp:' not in output


def test_runtime_banner_does_not_execute_external_tools(monkeypatch, capsys):
    monkeypatch.setattr(cli_utils, 'resolve_dependency_version', lambda *_args: '1.0')
    monkeypatch.setattr(cli_utils.shutil, 'which', lambda executable: '/tools/' + executable)

    def fail_command(*_args, **_kwargs):
        raise AssertionError('runtime banner must not execute an external tool')

    monkeypatch.setattr(cli_utils.subprocess, 'run', fail_command)

    cli_utils.print_runtime_banner(['amalgkit', 'getfastq'])

    output = capsys.readouterr().out
    assert 'AMALGKIT tool fasterq-dump: FOUND (/tools/fasterq-dump)' in output


def test_runtime_banner_reports_explicit_executable_path(tmp_path, monkeypatch, capsys):
    from types import SimpleNamespace
    executable = tmp_path / 'seqkit'
    executable.write_text('#!/bin/sh\nexit 0\n')
    executable.chmod(0o700)
    monkeypatch.setenv('PATH', '')
    cli_utils.print_runtime_banner(['amalgkit', 'getfastq'], args=SimpleNamespace(seqkit_exe=str(executable)))
    assert 'AMALGKIT tool seqkit: FOUND ({})'.format(executable) in capsys.readouterr().out


@pytest.mark.parametrize('converter', [cli_utils.int_or_auto, cli_utils.nonnegative_int_or_auto, cli_utils.positive_float_or_auto])
def test_auto_converters_normalize_and_accept_positive_values(converter):
    assert converter(' AuTo ') == 'auto'
    assert converter('2') == 2
    with pytest.raises(ValueError):
        converter('-1')
    with pytest.raises(ValueError):
        converter('invalid')


@pytest.mark.parametrize('value', ['0', '-0.1', 'nan', 'inf', '-inf'])
def test_positive_float_or_auto_rejects_nonpositive_and_nonfinite(value):
    with pytest.raises(ValueError):
        cli_utils.positive_float_or_auto(value)


def test_integer_auto_zero_boundary():
    assert cli_utils.nonnegative_int_or_auto('0') == 0
    with pytest.raises(ValueError):
        cli_utils.int_or_auto('0')


@pytest.mark.parametrize('fails', [False, True])
def test_timed_handler_imports_lazily_and_logs_outcome(monkeypatch, capsys, fails):
    failure = RuntimeError('command failure')
    entry = Mock(side_effect=failure if fails else None)
    importer = Mock(return_value=SimpleNamespace(main=entry))
    logger = Mock()
    monkeypatch.setattr(cli_utils.importlib, 'import_module', importer)
    monkeypatch.setattr(cli_utils, 'get_logger', lambda _: logger)
    handler = cli_utils.build_timed_command_handler('quant', 'amalgkit.quant', 'main')
    importer.assert_not_called()
    args = SimpleNamespace()
    if fails:
        with pytest.raises(RuntimeError) as caught:
            handler(args)
        assert caught.value is failure
        assert args._amalgkit_failure_logged is True
        assert logger.exception.call_args.kwargs['extra']['event'] == 'command_failed'
        assert 'amalgkit quant: end' not in capsys.readouterr().out
    else:
        handler(args)
        logger.exception.assert_not_called()
        assert logger.info.call_args.kwargs['extra']['event'] == 'command_end'
        assert 'amalgkit quant: end' in capsys.readouterr().out
    importer.assert_called_once_with('amalgkit.quant')
    entry.assert_called_once_with(args)


@pytest.mark.parametrize('topic', [None, 'quant'])
def test_help_handler_routes_topic(topic):
    parser = Mock()
    cli_utils.build_help_command_handler(parser)(SimpleNamespace(topic=topic))
    if topic is None:
        parser.print_help.assert_called_once_with()
        parser.parse_args.assert_not_called()
    else:
        parser.parse_args.assert_called_once_with(['quant', '--help'])
        parser.print_help.assert_not_called()
