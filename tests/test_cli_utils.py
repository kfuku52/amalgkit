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
