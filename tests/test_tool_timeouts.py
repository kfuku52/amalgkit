import subprocess
from types import SimpleNamespace

import pytest

from amalgkit.busco import BUSCO_TOOL_TIMEOUT_SECONDS, resolve_busco_tool_timeout_seconds
from amalgkit.getfastq import (
    GETFASTQ_TOOL_TIMEOUT_SECONDS,
    concatenate_files_with_system_cat,
    resolve_getfastq_dependency_probe_timeout_seconds,
    resolve_getfastq_tool_timeout_seconds,
)
from amalgkit.integrate import INTEGRATE_SEQKIT_TIMEOUT_SECONDS, resolve_integrate_tool_timeout_seconds
from amalgkit.quant import (
    QUANT_TOOL_TIMEOUT_SECONDS,
    check_kallisto_dependency,
    resolve_dependency_probe_timeout_seconds,
    resolve_quant_tool_timeout_seconds,
)
from amalgkit.main import build_main_parser
from amalgkit.subprocess_utils import DEPENDENCY_PROBE_TIMEOUT_SECONDS




HELP_TIMEOUT_CONTRACTS = (
    ('metadata', ('--ncbi_metadata_timeout_seconds',), None),
    (
        'quant',
        ('--tool_timeout_seconds', '--dependency_probe_timeout_seconds'),
        QUANT_TOOL_TIMEOUT_SECONDS,
    ),
    ('busco', ('--tool_timeout_seconds',), BUSCO_TOOL_TIMEOUT_SECONDS),
    (
        'getfastq',
        (
            '--tool_timeout_seconds',
            '--dependency_probe_timeout_seconds',
            '--ncbi_metadata_timeout_seconds',
        ),
        GETFASTQ_TOOL_TIMEOUT_SECONDS,
    ),
    ('integrate', ('--tool_timeout_seconds',), INTEGRATE_SEQKIT_TIMEOUT_SECONDS),
)

RESOLVERS = (
    (resolve_quant_tool_timeout_seconds, QUANT_TOOL_TIMEOUT_SECONDS),
    (resolve_busco_tool_timeout_seconds, BUSCO_TOOL_TIMEOUT_SECONDS),
    (resolve_getfastq_tool_timeout_seconds, GETFASTQ_TOOL_TIMEOUT_SECONDS),
    (resolve_integrate_tool_timeout_seconds, INTEGRATE_SEQKIT_TIMEOUT_SECONDS),
)



@pytest.mark.parametrize('command, options, default_seconds', HELP_TIMEOUT_CONTRACTS)
def test_timeout_help_contract(command, options, default_seconds, capsys):
    with pytest.raises(SystemExit) as exc:
        build_main_parser().parse_args([command, '--help'])
    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    for option in options:
        assert option in help_text, (command, option)
    if default_seconds is not None:
        assert 'default={}'.format(default_seconds) in help_text, command


@pytest.mark.parametrize('resolver,default_seconds', RESOLVERS)
def test_tool_resolvers_preserve_default_and_honor_override(resolver, default_seconds):
    assert resolver(SimpleNamespace()) == float(default_seconds)
    assert resolver(SimpleNamespace(tool_timeout_seconds=123)) == 123.0


def test_probe_resolvers_default_and_override():
    for resolver in (resolve_dependency_probe_timeout_seconds,
                     resolve_getfastq_dependency_probe_timeout_seconds):
        assert resolver(SimpleNamespace()) == float(DEPENDENCY_PROBE_TIMEOUT_SECONDS)
        assert resolver(SimpleNamespace(dependency_probe_timeout_seconds=7)) == 7.0


def test_dependency_probe_override_reaches_the_subprocess_runner(monkeypatch):
    # End to end for the probe path: the CLI value must arrive at subprocess.run.
    observed = {}

    def fake_run(cmd, stdout=None, stderr=None, timeout=None):
        observed['timeout'] = timeout
        return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

    monkeypatch.setattr('amalgkit.quant.subprocess.run', fake_run)
    check_kallisto_dependency(SimpleNamespace(dependency_probe_timeout_seconds=11))
    assert observed['timeout'] == 11.0

    observed.clear()
    check_kallisto_dependency(SimpleNamespace())
    assert observed['timeout'] == float(DEPENDENCY_PROBE_TIMEOUT_SECONDS)


def test_tool_timeout_override_reaches_the_subprocess_runner(monkeypatch, tmp_path):
    # End to end for an analysis path: integrate's seqkit scan must receive the
    # resolved value, and the documented default when nothing is supplied.
    from amalgkit import integrate

    fastq_path = tmp_path / 'sample.fq'
    fastq_path.write_text('@r1\nACGT\n+\nIIII\n')
    header = 'file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n'
    row = '{}\tFASTQ\tDNA\t1\t4\t4\t4\t4\n'.format(fastq_path)
    observed = {}

    def fake_run(cmd, stdout=None, stderr=None, timeout=None):
        observed['timeout'] = timeout
        return subprocess.CompletedProcess(cmd, 0, stdout=(header + row).encode(), stderr=b'')

    monkeypatch.setattr('amalgkit.integrate.subprocess.run', fake_run)

    integrate.scan_fastq_stats_with_seqkit_batch(
        path_fastq_paths=[str(fastq_path)],
        timeout_seconds=integrate.resolve_integrate_tool_timeout_seconds(
            SimpleNamespace(tool_timeout_seconds=42)
        ),
    )
    assert observed['timeout'] == 42.0

    observed.clear()
    integrate.scan_fastq_stats_with_seqkit_batch(path_fastq_paths=[str(fastq_path)])
    assert observed['timeout'] == float(INTEGRATE_SEQKIT_TIMEOUT_SECONDS)


def test_getfastq_cat_receives_timeout_and_removes_partial_output_on_timeout(monkeypatch, tmp_path):
    infile_path = tmp_path / 'input.fastq.gz'
    outfile_path = tmp_path / 'output.fastq.gz'
    infile_path.write_bytes(b'content')
    observed = {}

    monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda _name: '/bin/cat')

    def fake_run(command, stdout=None, stderr=None, timeout=None):
        observed['timeout'] = timeout
        stdout.write(b'partial')
        raise subprocess.TimeoutExpired(command, timeout)

    monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

    with pytest.raises(TimeoutError, match='Command timed out after 17 sec'):
        concatenate_files_with_system_cat(
            [str(infile_path)],
            str(outfile_path),
            timeout_seconds=17,
        )

    assert observed['timeout'] == 17.0
    assert not outfile_path.exists()


def test_tool_timeout_zero_disables_the_limit_end_to_end(monkeypatch, tmp_path):
    from amalgkit import integrate

    fastq_path = tmp_path / 'sample.fq'
    fastq_path.write_text('@r1\nACGT\n+\nIIII\n')
    header = 'file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n'
    row = '{}\tFASTQ\tDNA\t1\t4\t4\t4\t4\n'.format(fastq_path)
    observed = {'called_with_timeout': False}

    def fake_run(cmd, stdout=None, stderr=None, **kwargs):
        observed['called_with_timeout'] = 'timeout' in kwargs
        return subprocess.CompletedProcess(cmd, 0, stdout=(header + row).encode(), stderr=b'')

    monkeypatch.setattr('amalgkit.integrate.subprocess.run', fake_run)
    integrate.scan_fastq_stats_with_seqkit_batch(
        path_fastq_paths=[str(fastq_path)],
        timeout_seconds=integrate.resolve_integrate_tool_timeout_seconds(
            SimpleNamespace(tool_timeout_seconds=0)
        ),
    )
    # 0 resolves to None, and no timeout kwarg is passed to the runner at all.
    assert observed['called_with_timeout'] is False
