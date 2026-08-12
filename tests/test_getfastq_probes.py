"""Direct tests for amalgkit.getfastq_probes parsing helpers.

Table-driven positive/negative cases for the fasterq-dump version parser,
the compatibility probe, and the written-spots / written-reads extractors.
"""

import pytest

from amalgkit.getfastq_probes import (
    detect_fasterq_spot_range_support,
    ensure_supported_fasterq_dump_version,
    parse_fasterq_dump_version,
    parse_fasterq_dump_written_reads,
    parse_fasterq_dump_written_spots,
)


@pytest.mark.parametrize(
    ('stdout_txt', 'stderr_txt', 'expected'),
    [
        ('fasterq-dump : 3.0.3', '', (3, 0, 3)),
        ('fasterq-dump : 3.0', '', (3, 0, 0)),
        ('version 2.10.9 of sra-tools', '', (2, 10, 9)),
        ('', 'fasterq-dump : 4.1.2', (4, 1, 2)),
        ('no version here', '', None),
        ('', '', None),
    ],
)
def test_parse_fasterq_dump_version(stdout_txt, stderr_txt, expected):
    assert parse_fasterq_dump_version(stdout_txt, stderr_txt) == expected


@pytest.mark.parametrize(
    ('stdout_txt', 'stderr_txt', 'expected'),
    [
        ('fasterq-dump : 3.0.3', '', (3, 0, 3)),
        ('fasterq-dump : 3.0', '', (3, 0, 0)),
        ('', 'fasterq-dump : 4.1.2', (4, 1, 2)),
    ],
)
def test_ensure_supported_fasterq_dump_version_accepts_v3(stdout_txt, stderr_txt, expected):
    assert ensure_supported_fasterq_dump_version(stdout_txt, stderr_txt) == expected


def test_ensure_supported_fasterq_dump_version_rejects_v2():
    with pytest.raises(RuntimeError, match='Unsupported fasterq-dump version'):
        ensure_supported_fasterq_dump_version('fasterq-dump : 2.10.9', '')


def test_ensure_supported_fasterq_dump_version_rejects_undetectable():
    with pytest.raises(RuntimeError, match='Could not determine fasterq-dump version'):
        ensure_supported_fasterq_dump_version('unhelpful output', '')


@pytest.mark.parametrize(
    ('stdout_txt', 'stderr_txt', 'expected'),
    [
        ('  spots written : 123,456', '', 123456),
        ('junk\nspots written: 5', '', 5),
        ('', 'SPOTS WRITTEN : 7,890,123', 7890123),
        ('spots written : 1\nspots written : 2', '', 2),
        ('nothing relevant', '', None),
        ('', '', None),
    ],
)
def test_parse_fasterq_dump_written_spots(stdout_txt, stderr_txt, expected):
    assert parse_fasterq_dump_written_spots(stdout_txt, stderr_txt) == expected


@pytest.mark.parametrize(
    ('stdout_txt', 'stderr_txt', 'expected'),
    [
        ('reads written: 10,000', '', 10000),
        ('', 'reads written : 42', 42),
        ('reads written : 7\nreads written : 8', '', 8),
        ('no reads line', '', None),
    ],
)
def test_parse_fasterq_dump_written_reads(stdout_txt, stderr_txt, expected):
    assert parse_fasterq_dump_written_reads(stdout_txt, stderr_txt) == expected


@pytest.mark.parametrize(
    ('help_txt', 'expected'),
    [
        # Empty legacy help is treated as supported (can't disprove support).
        ('', True),
        # Newer long flags.
        ('--min-spot-id X\n--max-spot-id Y', True),
        # Older sra-tools camelCase flags.
        ('--minSpotId\n--maxSpotId', True),
        # Short flags.
        ('  -N|  \n  -X|  ', True),
        # Only one bound present -> no spot-range support.
        ('--minSpotId', False),
        ('-N|', False),
        # No flag at all.
        ('generic help text', False),
    ],
)
def test_detect_fasterq_spot_range_support(help_txt, expected):
    assert detect_fasterq_spot_range_support(help_txt) is expected
