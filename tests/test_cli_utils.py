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
    monkeypatch.setattr(cli_utils, 'resolve_external_tool_status', lambda *_args: 'MISSING')

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

    monkeypatch.setattr(cli_utils, 'resolve_external_tool_status', fail_external_probe)

    cli_utils.print_runtime_banner(['amalgkit', 'metadata'])

    assert 'AMALGKIT tool' not in capsys.readouterr().out


def test_runtime_banner_probes_only_tools_relevant_to_active_command(monkeypatch, capsys):
    monkeypatch.setattr(cli_utils, 'resolve_dependency_version', lambda *_args: '1.0')
    observed = []

    def record_probe(executable_name, _version_commands):
        observed.append(executable_name)
        return 'MISSING'

    monkeypatch.setattr(cli_utils, 'resolve_external_tool_status', record_probe)

    cli_utils.print_runtime_banner(['amalgkit', 'quant'])

    assert observed == ['kallisto', 'oarfish']
    output = capsys.readouterr().out
    assert 'AMALGKIT tool kallisto:' in output
    assert 'AMALGKIT tool fastp:' not in output
