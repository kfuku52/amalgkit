import numpy
import pandas
import pytest

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


@pytest.mark.slow
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


def test_tau_histogram_uses_linear_run_means(tmp_path, stub_pdf_encoding):
    counts = pandas.DataFrame([[0., 100., 10., 10.]], columns=['a1', 'a2', 'b1', 'b2'])
    metadata = pandas.DataFrame({
        'run': counts.columns, 'sample_group': ['A', 'A', 'B', 'B'], 'exclusion': ['no'] * 4,
    })
    result = save_tau_histogram_pdf(
        numpy.log2(counts + 1), metadata, ['A', 'B'], str(tmp_path / 'tau.pdf'),
        transform_method='log2p1-none',
    )
    numpy.testing.assert_allclose(result['tau_df']['tau'], [0.8])


@pytest.mark.parametrize('comparison', [False, True])
def test_tau_plot_does_not_discard_missing_run_metadata(tmp_path, monkeypatch, comparison, stub_pdf_encoding):
    from amalgkit import per_species_outputs as overview
    from amalgkit import per_species_finalize_python as finalize
    from amalgkit.per_species_common import sample_group_to_tau

    counts = pandas.DataFrame({
        'a1': [0., 10., 2., 30.], 'a2': [100., 12., 3., 20.],
        'b1': [10., 20., 4., 10.], 'b2': [10., 21., 5., 11.],
    })
    metadata = pandas.DataFrame({
        'run': counts.columns, 'sample_group': ['A', 'A', 'B', 'B'],
        'bioproject': ['p1', 'p1', 'p2', 'p2'], 'exclusion': ['no'] * 4,
    })
    captured = []

    def capture(*args, **kwargs):
        result = sample_group_to_tau(*args, **kwargs)
        captured.append(result)
        return result

    incomplete = counts.drop(columns='a1')
    if comparison:
        monkeypatch.setattr(finalize, 'sample_group_to_tau', capture)
        finalize.save_quick_state_comparison_plot(
            counts, incomplete, metadata, 'pearson', str(tmp_path / 'comparison.pdf'),
            ['A', 'B'], 'none-none', 'no',
        )
        assert captured[0]['tau'].iloc[0] == pytest.approx(0.8)
    else:
        monkeypatch.setattr(overview, 'sample_group_to_tau', capture)
        overview.save_state_overview_pdf(
            incomplete, metadata, ['A', 'B'], str(tmp_path / 'overview.pdf'), transform_method='none-none',
        )
    assert captured[-1]['tau'].isna().all()


@pytest.mark.slow
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


@pytest.mark.parametrize('labels', [
    pandas.Series(['0001', '1', 'NA', '', '葉', 'leaf', '0001', 'NA']),
    pandas.Series([1, '1', True, 'True', 'é', 'e\u0301', 1.0, '1.0'], dtype=object),
    pandas.Series([None, numpy.nan, pandas.NA, 'nan', '<NA>', 'None', '', 'NA'], dtype=object),
    pandas.Series(['x', None, 'x', '', 'NA', '葉', None, 'nan'], dtype='str'),
    pandas.Series(['x', None, 'x', '', 'NA', '葉', None, '<NA>'], dtype='string'),
    pandas.Series(pandas.Categorical(['x', None, 'x', 'NA', '葉', 'y', None, 'NA'])),
    pandas.Series([1, None, 1, 2, 3, None, 5, 6], dtype='Int64'),
    pandas.Series([True, False, True, None, False, None, True, False], dtype='boolean'),
    pandas.Series(['same'] * 8),
    pandas.Series([f'unique_{i}' for i in range(8)]),
])
def test_correlation_label_equality_preserves_converted_labels(labels):
    from amalgkit.per_species_outputs import _label_equality_matrix

    converted = labels.astype(str).to_numpy()
    expected = numpy.equal.outer(converted, converted)
    numpy.testing.assert_array_equal(_label_equality_matrix(labels), expected)


@pytest.mark.parametrize('precomputed_frame', [False, True])
def test_correlation_statistics_preserve_pair_order_and_nonfinite_values(precomputed_frame):
    runs = ['0001', 'NA', 'RUN3', 'RUN4', 'RUN5', 'RUN6']
    counts = pandas.DataFrame(numpy.zeros((2, 6)), columns=runs)
    metadata = pandas.DataFrame({
        'run': runs,
        'bioproject': ['0001', '1', '0001', '1', '0001', '1'],
        'sample_group': ['NA', 'NA', '葉', '葉', 'NA', 'NA'],
    })
    matrix = numpy.array([
        [1., .75, .2, -.1, .9, .8],
        [-.9, 1., numpy.nan, numpy.inf, .8, -.3],
        [-.8, -.7, 1., -.5, .4, .6],
        [-.6, -.5, -.4, 1., -numpy.inf, .1],
        [-.3, -.2, -.1, 0., 1., -.2],
        [.1, .2, .3, .4, .5, 1.],
    ])
    precomputed = pandas.DataFrame(matrix, index=runs, columns=runs) if precomputed_frame else matrix
    # The public wrapper aligns DataFrame matrices by their run order and arrays
    # by metadata order. Both paths must preserve their existing upper triangle.
    ordered_metadata = metadata.iloc[::-1].reset_index(drop=True)
    observed = save_correlation_statistics(
        counts, ordered_metadata, 'pearson', 7, precomputed_tc_dist_matrix=precomputed,
    )
    aligned = metadata if precomputed_frame else ordered_metadata
    values_by_group = {label: [] for label in ['bwbw', 'wibw', 'bwwi', 'wiwi']}
    for i in range(6):
        for j in range(i + 1, 6):
            if not numpy.isfinite(matrix[i, j]):
                continue
            same_batch = aligned.iloc[i]['bioproject'] == aligned.iloc[j]['bioproject']
            same_group = aligned.iloc[i]['sample_group'] == aligned.iloc[j]['sample_group']
            label = ('wi' if same_batch else 'bw') + ('wi' if same_group else 'bw')
            values_by_group[label].append(matrix[i, j])
    expected = {}
    for label, values in values_by_group.items():
        expected[f'{label}_n'] = float(len(values))
        expected[f'{label}_mean'] = float(numpy.mean(values)) if values else numpy.nan
        expected[f'{label}_median'] = float(numpy.median(values)) if values else numpy.nan
        expected[f'{label}_variance'] = float(numpy.var(values, ddof=1)) if len(values) > 1 else numpy.nan
    expected = pandas.DataFrame([expected], index=['round_7'], columns=CORRELATION_STAT_COLUMNS)
    pandas.testing.assert_frame_equal(observed, expected, check_exact=True)
