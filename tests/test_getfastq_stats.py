import pandas

from amalgkit.getfastq_stats import read_getfastq_stats_row


def test_read_getfastq_stats_row_returns_last_matching_run(tmp_path):
    pandas.DataFrame([
        {'run': 'SRR001', 'num_written': 10},
        {'run': 'SRR002', 'num_written': 20},
        {'run': 'SRR001', 'num_written': 30},
    ]).to_csv(tmp_path / 'getfastq_stats.tsv', sep='\t', index=False)

    row = read_getfastq_stats_row(str(tmp_path), 'SRR001')

    assert int(row['num_written']) == 30


def test_read_getfastq_stats_row_returns_none_for_missing_run(tmp_path):
    pandas.DataFrame([
        {'run': 'SRR002', 'num_written': 20},
    ]).to_csv(tmp_path / 'getfastq_stats.tsv', sep='\t', index=False)

    assert read_getfastq_stats_row(str(tmp_path), 'SRR001') is None


def test_read_getfastq_stats_row_reports_unreadable_file(tmp_path):
    (tmp_path / 'getfastq_stats.tsv').write_bytes(b'\xff')
    warnings = []

    row = read_getfastq_stats_row(str(tmp_path), 'SRR001', warning_writer=warnings.append)

    assert row is None
    assert len(warnings) == 1
    assert 'Failed to read getfastq stats file' in warnings[0]
