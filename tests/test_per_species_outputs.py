import numpy
import pandas

from amalgkit.per_species_outputs import (
    CORRELATION_STAT_COLUMNS,
    initialize_correlation_statistics,
    intersect_counts_and_metadata,
    save_correlation_statistics,
    save_state_overview_pdf,
    save_tau_histogram_pdf,
    write_table_with_index_name,
)


def test_intersect_counts_and_metadata_keeps_common_runs_only():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [1.0, 2.0],
            'RUN2': [3.0, 4.0],
            'RUN3': [5.0, 6.0],
        },
        index=['G1', 'G2'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN3', 'RUN1', 'RUN4'],
            'sample_group': ['A', 'A', 'B'],
        }
    )
    observed = intersect_counts_and_metadata(counts_df, metadata_df)
    assert list(observed['tc'].columns) == ['RUN1', 'RUN3']
    assert observed['sra']['run'].tolist() == ['RUN3', 'RUN1']


def test_write_table_with_index_name_sorts_and_writes_tsv(tmp_path):
    df = pandas.DataFrame({'RUN1': [2.0, 1.0]}, index=['G2', 'G1'])
    out_path = tmp_path / 'table.tsv'
    observed = write_table_with_index_name(df, str(out_path), index_name='target_id', sort=True)
    assert out_path.exists()
    assert observed['target_id'].tolist() == ['G1', 'G2']
    loaded = pandas.read_csv(out_path, sep='\t')
    pandas.testing.assert_frame_equal(loaded, observed)


def test_save_correlation_statistics_reports_expected_pair_counts():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [1.0, 2.0, 3.0, 4.0],
            'RUN2': [1.0, 2.1, 3.2, 4.1],
            'RUN3': [4.0, 3.0, 2.0, 1.0],
            'RUN4': [4.2, 3.1, 2.2, 1.1],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['root', 'root', 'leaf', 'leaf'],
            'bioproject': ['BP1', 'BP2', 'BP1', 'BP2'],
        }
    )
    observed = save_correlation_statistics(
        counts_df=counts_df,
        metadata_df=metadata_df,
        dist_method='pearson',
        round_value=3,
    )
    assert observed.index.tolist() == ['round_3']
    assert list(observed.columns) == CORRELATION_STAT_COLUMNS
    assert observed.loc['round_3', 'bwbw_n'] == 2.0
    assert observed.loc['round_3', 'wibw_n'] == 2.0
    assert observed.loc['round_3', 'bwwi_n'] == 2.0
    assert observed.loc['round_3', 'wiwi_n'] == 0.0
    assert observed.loc['round_3', 'bwbw_mean'] <= 1.0


def test_save_correlation_statistics_handles_single_sample_case():
    counts_df = pandas.DataFrame({'RUN1': [1.0, 2.0, 3.0]}, index=['G1', 'G2', 'G3'])
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1'],
            'sample_group': ['root'],
            'bioproject': ['BP1'],
        }
    )
    observed = save_correlation_statistics(
        counts_df=counts_df,
        metadata_df=metadata_df,
        dist_method='pearson',
        round_value=1,
        correlation_statistics=initialize_correlation_statistics(),
    )
    assert list(observed.columns) == CORRELATION_STAT_COLUMNS
    assert observed.index.tolist() == ['round_1']
    assert observed.loc['round_1', ['bwbw_n', 'wibw_n', 'bwwi_n', 'wiwi_n']].tolist() == [0.0, 0.0, 0.0, 0.0]
    assert numpy.isnan(observed.loc['round_1', 'bwbw_mean'])


def test_save_correlation_statistics_appends_rounds():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [1.0, 2.0, 3.0],
            'RUN2': [1.1, 2.1, 3.1],
        },
        index=['G1', 'G2', 'G3'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2'],
            'sample_group': ['root', 'root'],
            'bioproject': ['BP1', 'BP2'],
        }
    )
    stats_round_1 = save_correlation_statistics(
        counts_df=counts_df,
        metadata_df=metadata_df,
        dist_method='pearson',
        round_value=1,
    )
    stats_round_2 = save_correlation_statistics(
        counts_df=counts_df,
        metadata_df=metadata_df,
        dist_method='pearson',
        round_value=2,
        correlation_statistics=stats_round_1,
    )
    assert stats_round_2.index.tolist() == ['round_1', 'round_2']


def test_save_tau_histogram_pdf_writes_pdf(tmp_path):
    counts_df = pandas.DataFrame(
        {
            'RUN1': [0.0, 2.0, 4.0, 8.0],
            'RUN2': [0.0, 2.5, 3.5, 8.0],
            'RUN3': [0.0, 8.0, 1.0, 2.0],
            'RUN4': [0.0, 7.5, 1.5, 2.5],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['root', 'root', 'leaf', 'leaf'],
            'bioproject': ['BP1', 'BP2', 'BP1', 'BP2'],
            'exclusion': ['no', 'no', 'no', 'no'],
        }
    )
    out_path = tmp_path / 'tau_hist.pdf'
    result = save_tau_histogram_pdf(
        counts_df=counts_df,
        metadata_df=metadata_df,
        selected_sample_groups=['root', 'leaf'],
        out_pdf_path=str(out_path),
        font_size=8,
        transform_method='log2p1-fpkm',
    )
    assert out_path.exists()
    assert out_path.stat().st_size > 0
    assert result['num_total_genes'] == 4
    assert result['num_no_expression'] == 1


def test_save_state_overview_pdf_writes_pdf(tmp_path):
    counts_df = pandas.DataFrame(
        {
            'RUN1': [1.0, 2.0, 3.0, 4.0],
            'RUN2': [1.2, 2.2, 3.2, 4.2],
            'RUN3': [4.0, 3.0, 2.0, 1.0],
            'RUN4': [4.1, 3.1, 2.1, 1.1],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['root', 'root', 'leaf', 'leaf'],
            'bioproject': ['BP1', 'BP2', 'BP1', 'BP2'],
            'exclusion': ['no', 'no', 'no', 'no'],
        }
    )
    out_path = tmp_path / 'state_overview.pdf'
    observed = save_state_overview_pdf(
        counts_df=counts_df,
        metadata_df=metadata_df,
        selected_sample_groups=['root', 'leaf'],
        out_pdf_path=str(out_path),
        dist_method='pearson',
        transform_method='log2p1-fpkm',
        font_size=8,
    )
    assert observed == str(out_path)
    assert out_path.exists()
    assert out_path.stat().st_size > 0
