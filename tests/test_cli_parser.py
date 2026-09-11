"""Validate numeric option contracts through the real subcommand parser."""

import pytest

from amalgkit.main import build_main_parser


# Explicit option lists catch a validator accidentally removed from one consumer.
CONTRACTS = [
    (
        'positive_int',
        [('getfastq', 'contam_filter_chunk_spots'),
         ('getfastq', 'sra_download_wait_timeout_seconds'),
         ('getfastq', 'sra_download_transfer_timeout_seconds'),
         ('cstmm', 'tmm_imputation_rank'), ('cstmm', 'tmm_imputation_max_iter'),
         ('wsfilter', 'max_filter_iterations')],
        int, ['1', '23'], ['0', '-1', '1.5', 'bad', 'nan', 'inf', '-inf'],
    ),
    (
        'single_copy_threshold',
        [('cstmm', 'single_copy_threshold'), ('csfilter', 'single_copy_threshold')],
        float, ['0.1', '100'], ['0', '-1', '100.1', 'nan', 'inf', '-inf', 'bad'],
    ),
    (
        'finite_float',
        [('wsfilter', 'robust_z_threshold'), ('csfilter', 'robust_z_threshold')],
        float, ['-2.5', '0', '2.5'], ['nan', 'inf', '-inf', 'bad'],
    ),
    (
        'mapping_rate_threshold', [('wsfilter', 'mapping_rate')],
        float, ['0', '100', '0.2'], ['-0.1', '100.1', 'nan', 'inf', '-inf', 'bad'],
    ),
    (
        'correlation_threshold', [('wsfilter', 'correlation_threshold')],
        float, ['-1', '0', '1'], ['-1.1', '1.1', 'nan', 'inf', '-inf', 'bad'],
    ),
    (
        'correlation_margin', [('wsfilter', 'margin_threshold'), ('csfilter', 'margin_threshold')],
        float, ['-2', '0', '2'], ['-2.1', '2.1', 'nan', 'inf', '-inf', 'bad'],
    ),
    (
        'positive_float', [('cstmm', 'tmm_imputation_tol')],
        float, ['1e-6', '1'], ['0', '-1', 'nan', 'inf', '-inf', 'bad'],
    ),
    (
        'common_gene_threshold', [('wsfilter', 'min_common_genes'), ('csfilter', 'min_common_genes')],
        int, ['0', '2', '50'], ['-1', '1', '1.5', 'bad', 'nan', 'inf', '-inf'],
    ),
]


@pytest.fixture(scope='module')
def parser():
    return build_main_parser()


@pytest.mark.parametrize(
    'command, option, value, expected_type',
    [pytest.param(command, option, value, expected_type, id=f'{command}-{option}-{value}')
     for _, consumers, expected_type, valid, _ in CONTRACTS
     for command, option in consumers for value in valid],
)
def test_numeric_options_accept_valid_values_and_boundaries(parser, command, option, value, expected_type, capsys):
    result = parser.parse_args([command, f'--{option}={value}'])
    parsed = getattr(result, option)
    assert type(parsed) is expected_type
    assert parsed == expected_type(value)
    assert capsys.readouterr().err == ''


@pytest.mark.parametrize(
    'command, option, value',
    [pytest.param(command, option, value, id=f'{command}-{option}-{value}')
     for _, consumers, _, _, invalid in CONTRACTS
     for command, option in consumers for value in invalid],
)
def test_numeric_options_reject_invalid_values_with_usage_error(parser, command, option, value, capsys):
    # The equals form ensures -inf reaches the type validator, rather than
    # argparse rejecting it as an apparent option with a missing value.
    with pytest.raises(SystemExit) as exc:
        parser.parse_args([command, f'--{option}={value}'])
    assert exc.value.code == 2
    stderr = capsys.readouterr().err
    assert 'usage:' in stderr.lower()
    assert 'error:' in stderr.lower()
    assert f'--{option}' in stderr
    assert 'Traceback' not in stderr
