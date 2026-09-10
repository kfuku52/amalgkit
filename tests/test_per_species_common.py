import numpy
import pandas
import pytest
import json

from amalgkit.per_species_common import (
    append_round_summary,
    build_curation_final_summary,
    initialize_round_summary,
    linear_sample_group_summary,
    write_tau_outputs,
    sample_group_mean,
    sample_group_to_tau,
    write_curation_summaries,
)


def _tau_fixture():
    counts = pandas.DataFrame([[0., 100., 10., 10.]], columns=['a1', 'a2', 'b1', 'b2'], index=['g'])
    metadata = pandas.DataFrame({
        'run': counts.columns, 'sample_group': ['A', 'A', 'B', 'B'],
        'exclusion': ['no'] * 4, 'bioproject': ['p1', 'p2', 'p1', 'p2'],
        'donor': ['d1', 'd2', 'd1', 'd2'], 'biosample': ['s1', 's2', 's3', 's4'],
    })
    return counts, metadata


@pytest.mark.parametrize('transform', ['none', 'log2', 'logn', 'log2p1', 'lognp1'])
def test_tau_inverse_transforms_runs_before_averaging(transform):
    counts, metadata = _tau_fixture()
    with numpy.errstate(divide='ignore'):
        encoded = {
            'none': counts, 'log2': numpy.log2(counts), 'logn': numpy.log(counts),
            'log2p1': numpy.log2(counts + 1), 'lognp1': numpy.log1p(counts),
        }[transform]
    result = linear_sample_group_summary(encoded, metadata, transform_method=transform + '-none')
    numpy.testing.assert_allclose(result['linear_mean'].loc['g'], [50, 10])
    tau = sample_group_to_tau(result['linear_mean'])
    assert tau.loc['g', 'tau'] == pytest.approx(0.8)
    assert tau.loc['g', 'highest'] == 'A'
    assert tau.loc['g', 'tau_status'] == 'ok'


def test_tau_keeps_tiny_log2p1_expression():
    counts, metadata = _tau_fixture()
    counts.loc['g'] = [1e-18, 1e-18, 0, 0]
    result = linear_sample_group_summary(counts, metadata, transform_method='log2p1-none')['linear_mean']
    assert result.loc['g', 'A'] == pytest.approx(numpy.log(2) * 1e-18, rel=1e-14, abs=0)
    assert sample_group_to_tau(result).loc['g', 'tau'] == 1


@pytest.mark.parametrize('balance', [False, True])
def test_tau_mean_does_not_overflow_finite_runs(balance):
    counts, metadata = _tau_fixture()
    counts.loc['g'] = [1e308, 1e308, 0, 0]
    result = linear_sample_group_summary(counts, metadata, transform_method='none-none', balance_projects=balance)
    assert result['linear_mean'].loc['g', 'A'] == 1e308
    assert sample_group_to_tau(result['linear_mean']).loc['g', 'tau'] == 1


@pytest.mark.parametrize('columns', [['A', 'B', 'A'], ['A', ' A '], ['A', ''], ['A', None]])
def test_tau_rejects_ambiguous_tissue_labels(columns):
    with pytest.raises(ValueError, match='sample_group labels'):
        sample_group_to_tau(pandas.DataFrame([[10.] * len(columns)], columns=columns))


@pytest.mark.parametrize('transform', ['typo', 'log2p1-typo', 'none-fpkm-extra'])
def test_tau_is_linear_only_and_rejects_unknown_input_transform(transform):
    counts, metadata = _tau_fixture()
    with pytest.raises(ValueError, match='Unsupported tau input transformation'):
        linear_sample_group_summary(counts, metadata, transform_method=transform)
    with pytest.raises(TypeError, match='transform_method'):
        sample_group_to_tau(pandas.DataFrame([[10., 0.]]), transform_method='log2p1-none')


def test_tau_annotation_lists_escape_group_delimiters():
    from amalgkit.text_utils import parse_sample_group_argument

    columns = ['brain|forebrain', 'liver,adult', r'root\tip']
    tau = sample_group_to_tau(pandas.DataFrame([[10., 10., 1.]], columns=columns)).iloc[0]
    assert parse_sample_group_argument(tau['order']) == columns
    assert parse_sample_group_argument(tau['highest_ties']) == columns[:2]
    assert tau['highest'] == columns[0]


@pytest.mark.parametrize('unit', ['run', 'biosample', 'donor'])
def test_tau_run_identity_and_unit_order_do_not_change_aggregation(unit):
    counts, metadata = _tau_fixture()
    # Run IDs are lexical keys even for nonstandard private data.
    counts.columns = [' a1 ', 'a2', 'b1', 'b2']
    metadata['run'] = counts.columns
    metadata['sample_group'] = [' A ', 'A', 'B', 'B']
    result = linear_sample_group_summary(counts, metadata, transform_method='none-none', unit=unit)
    permuted = linear_sample_group_summary(
        counts.iloc[:, ::-1], metadata.iloc[::-1], ['A', 'B'], 'none-none', unit=unit,
    )
    pandas.testing.assert_frame_equal(result['linear_mean'], permuted['linear_mean'])
    numpy.testing.assert_allclose(result['linear_mean'].loc['g'], [50, 10])


def test_tau_weights_match_explicit_unit_and_project_means():
    counts = pandas.DataFrame([[0., 0., 60., 100.]], columns=['r1', 'r2', 'r3', 'r4'])
    metadata = pandas.DataFrame({
        'run': counts.columns, 'sample_group': ['A'] * 4, 'exclusion': ['no'] * 4,
        'donor': ['d1', 'd1', 'd2', 'd3'], 'bioproject': ['p1', 'p1', 'p1', 'p2'],
    })
    for unit, balanced, expected in [('run', False, 40), ('donor', False, 160 / 3), ('donor', True, 65)]:
        result = linear_sample_group_summary(counts, metadata, transform_method='none-none', unit=unit, balance_projects=balanced)
        assert result['linear_mean'].iloc[0, 0] == pytest.approx(expected)
        weights = result['weights'].set_index('run')['weight']
        assert weights.sum() == pytest.approx(1)
        assert (counts.iloc[0] * weights).sum() == pytest.approx(expected)
    duplicated = pandas.concat([counts, counts[['r1']].rename(columns={'r1': 'r5'})], axis=1)
    duplicated_metadata = pandas.concat([metadata, metadata.iloc[[0]].assign(run='r5')], ignore_index=True)
    result = linear_sample_group_summary(duplicated, duplicated_metadata, transform_method='none-none', unit='donor')
    assert result['linear_mean'].iloc[0, 0] == pytest.approx(160 / 3)


def test_tau_never_drops_missing_groups_or_renormalizes_gene_weights():
    counts, metadata = _tau_fixture()
    counts.loc['g', 'a1'] = numpy.nan
    result = linear_sample_group_summary(counts, metadata, ['A', 'B', 'C'], 'none-none')
    assert result['linear_mean'].columns.tolist() == ['A', 'B', 'C']
    assert numpy.isnan(result['linear_mean'].loc['g', 'A'])
    assert result['coverage'].set_index('sample_group').loc['C', 'status'] == 'not_collected'
    tau = sample_group_to_tau(result['linear_mean'])
    assert numpy.isnan(tau.loc['g', 'tau'])
    assert numpy.isnan(tau.loc['g', 'highest'])
    assert tau.loc['g', 'num_groups_required'] == 3
    metadata.loc[metadata.sample_group.eq('A'), 'exclusion'] = 'manual'
    result = linear_sample_group_summary(counts, metadata, ['A', 'B'], 'none-none')
    assert result['coverage'].iloc[0]['status'] == 'all_excluded'
    assert result['linear_mean']['A'].isna().all()


@pytest.mark.parametrize('values,status', [
    ([10., 0., numpy.nan], 'missing_or_nonfinite_expression'),
    ([10., 0., numpy.inf], 'missing_or_nonfinite_expression'),
    ([10., -1., 0.], 'negative_expression'),
    ([0., 0., 0.], 'all_zero'),
])
def test_tau_invalid_rows_are_not_reported_as_specific(values, status):
    result = sample_group_to_tau(pandas.DataFrame([values]))
    assert result['tau'].isna().all()
    assert result['highest'].isna().all()
    assert result['tau_status'].iloc[0] == status


def test_tau_ties_are_explicit_and_panel_order_invariant():
    counts = pandas.DataFrame([[10., 10., 0.]], columns=['B', 'A', 'C'])
    for columns in [['B', 'A', 'C'], ['C', 'A', 'B']]:
        result = sample_group_to_tau(counts[columns])
        assert result['tau'].iloc[0] == pytest.approx(0.5)
        assert result['highest_ties'].iloc[0] == 'A|B'


@pytest.mark.parametrize('columns', [[], ['A', 'B']])
def test_tau_empty_gene_table_retains_schema(columns):
    result = sample_group_to_tau(pandas.DataFrame(columns=columns))
    assert result.empty
    assert {'tau', 'highest', 'order', 'tau_status', 'num_groups_required'}.issubset(result.columns)


def test_tau_missing_run_is_flagged_without_changing_weights():
    counts, metadata = _tau_fixture()
    summary = linear_sample_group_summary(counts.drop(columns='a1'), metadata, ['A', 'B'], 'none-none')
    assert summary['coverage'].iloc[0]['status'] == 'missing_runs'
    assert summary['linear_mean']['A'].isna().all()
    assert sample_group_to_tau(summary['linear_mean'])['tau'].isna().all()


def test_tau_preserves_literal_ids_and_does_not_count_unknown_projects():
    counts, metadata = _tau_fixture()
    metadata['donor'] = ['001', '1', 'NA', 'null']
    metadata['bioproject'] = 'not_provided'
    result = linear_sample_group_summary(counts, metadata, transform_method='none-none', unit='donor')
    assert set(result['weights']['unit_id']) == {'001', '1', 'NA', 'null'}
    assert result['coverage']['num_projects'].eq(0).all()
    with pytest.raises(ValueError, match='nonmissing bioproject'):
        linear_sample_group_summary(counts, metadata, unit='donor', balance_projects=True)


def test_tau_rejects_unknown_and_cross_project_units():
    counts, metadata = _tau_fixture()
    metadata.loc[0, 'donor'] = ''
    with pytest.raises(ValueError, match='nonmissing donor'):
        linear_sample_group_summary(counts, metadata, unit='donor')
    metadata['donor'] = ['shared', 'shared', 'd1', 'd2']
    with pytest.raises(ValueError, match='multiple projects'):
        linear_sample_group_summary(counts, metadata, unit='donor', balance_projects=True)


def test_tau_export_records_definition_and_reproducible_inputs(tmp_path):
    counts, metadata = _tau_fixture()
    linear = write_tau_outputs(numpy.log2(counts + 1), metadata, ['A', 'B'], 'log2p1-none', tmp_path, 'Species', 'no')
    definition = json.loads((tmp_path / 'Species.no.tau.definition.json').read_text())
    tau = pandas.read_csv(tmp_path / 'Species.no.tau.tsv', sep='\t', index_col='target_id')
    assert definition['representative'] == 'linear_arithmetic_mean'
    assert definition['unit'] == 'run'
    assert tau.loc['g', 'panel_id'] == definition['panel_id']
    assert tau.loc['g', 'tau'] == pytest.approx(sample_group_to_tau(linear).loc['g', 'tau'])
    assert (tmp_path / 'Species.no.tau.weights.tsv').is_file()


def test_sample_group_mean_and_tau_for_two_groups():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [1.0, 5.0, 0.0, 7.0],
            'RUN2': [3.0, 7.0, 2.0, 9.0],
            'RUN3': [10.0, 0.0, 4.0, 8.0],
            'RUN4': [14.0, 2.0, 6.0, 10.0],
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

    observed = sample_group_mean(
        counts_df=counts_df,
        metadata_df=metadata_df,
        selected_sample_groups=None,
        balance_bp=False,
    )
    expected_mean = pandas.DataFrame(
        {
            'root': [2.0, 6.0, 1.0, 8.0],
            'leaf': [12.0, 1.0, 5.0, 9.0],
        },
        index=['G1', 'G2', 'G3', 'G4'],
    )
    pandas.testing.assert_frame_equal(observed['tc_ave'], expected_mean)
    assert observed['selected_sample_groups'] == ['root', 'leaf']

    tau_df = sample_group_to_tau(
        tc_sample_group_df=observed['tc_ave'],
        rich_annotation=True,
    )
    numpy.testing.assert_allclose(
        tau_df['tau'].to_numpy(dtype=float),
        numpy.array([5.0 / 6.0, 5.0 / 6.0, 0.8, 1.0 / 9.0]),
        rtol=0.0,
        atol=1e-12,
    )
    assert tau_df['highest'].tolist() == ['leaf', 'root', 'leaf', 'leaf']
    assert tau_df['order'].tolist() == ['leaf|root', 'root|leaf', 'leaf|root', 'leaf|root']


def test_sample_group_mean_excludes_flagged_runs_from_group_average():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [10.0, 20.0],
            'RUN2': [100.0, 200.0],
        },
        index=['G1', 'G2'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2'],
            'sample_group': ['root', 'root'],
            'bioproject': ['BP1', 'BP1'],
            'exclusion': ['no', 'manual_removal'],
        }
    )

    expected_mean = pandas.DataFrame({'root': [10.0, 20.0]}, index=['G1', 'G2'])
    for balance_bp in (False, True):
        observed = sample_group_mean(
            counts_df=counts_df,
            metadata_df=metadata_df,
            balance_bp=balance_bp,
        )
        pandas.testing.assert_frame_equal(observed['tc_ave'], expected_mean)
        assert observed['selected_sample_groups'] == ['root']


def test_sample_group_mean_balance_bp_drops_fully_excluded_group():
    counts_df = pandas.DataFrame(
        {
            'RUN1': [2.0, 4.0, 8.0],
            'RUN2': [6.0, 8.0, 10.0],
            'RUN3': [3.0, 9.0, 12.0],
            'RUN4': [5.0, 7.0, 11.0],
        },
        index=['G1', 'G2', 'G3'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3', 'RUN4'],
            'sample_group': ['flower', 'flower', 'root', 'root'],
            'bioproject': ['BP1', 'BP2', 'BP1', 'BP2'],
            'exclusion': ['manual_removal', 'manual_removal', 'no', 'no'],
        }
    )

    observed = sample_group_mean(
        counts_df=counts_df,
        metadata_df=metadata_df,
        selected_sample_groups=['flower', 'root'],
        balance_bp=True,
    )
    expected_mean = pandas.DataFrame({'root': [4.0, 8.0, 11.5]}, index=['G1', 'G2', 'G3'])
    pandas.testing.assert_frame_equal(observed['tc_ave'], expected_mean)
    assert observed['selected_sample_groups'] == ['root']

    tau_df = sample_group_to_tau(
        tc_sample_group_df=observed['tc_ave'],
        rich_annotation=True,
    )
    assert tau_df['highest'].isna().all()
    assert tau_df['order'].isna().all()
    assert tau_df['tau_status'].eq('insufficient_groups').all()
    assert tau_df['tau'].isna().all()


def test_sample_group_to_tau_handles_empty_input():
    empty_df = pandas.DataFrame(index=['G1', 'G2'])
    observed = sample_group_to_tau(empty_df, rich_annotation=True)
    assert {'tau', 'highest', 'order', 'tau_status'}.issubset(observed.columns)
    assert list(observed.index) == ['G1', 'G2']
    assert observed[['tau', 'highest', 'order']].isna().all().all()
    assert observed['tau_status'].eq('insufficient_groups').all()


def test_round_summary_helpers_match_expected_shape_and_values():
    round_summary = initialize_round_summary()
    assert list(round_summary.columns) == [
        'step',
        'round',
        'reason',
        'num_runs_before',
        'num_runs_after',
        'num_runs_removed',
        'removed_runs',
    ]
    observed = append_round_summary(
        round_summary=round_summary,
        step='mapping_rate_cutoff',
        round_value=2,
        reason='low_mapping_rate',
        runs_before=['RUN1', 'RUN2', 'RUN3'],
        runs_after=['RUN1', 'RUN3'],
    )
    assert observed.shape == (1, 7)
    assert observed.loc[0, 'step'] == 'mapping_rate_cutoff'
    assert int(observed.loc[0, 'round']) == 2
    assert observed.loc[0, 'reason'] == 'low_mapping_rate'
    assert int(observed.loc[0, 'num_runs_before']) == 3
    assert int(observed.loc[0, 'num_runs_after']) == 2
    assert int(observed.loc[0, 'num_runs_removed']) == 1
    assert observed.loc[0, 'removed_runs'] == 'RUN2'


def test_build_and_write_curation_summaries(tmp_path):
    round_summary = append_round_summary(
        round_summary=initialize_round_summary(),
        step='correlation_iter',
        round_value=3,
        reason='low_within_sample_group_correlation',
        runs_before=['RUN1', 'RUN2', 'RUN3'],
        runs_after=['RUN1'],
    )
    metadata_df = pandas.DataFrame(
        {
            'run': ['RUN1', 'RUN2', 'RUN3'],
            'exclusion': ['no', 'manual_removal', 'low_mapping_rate'],
        }
    )
    final_summary = build_curation_final_summary(
        metadata_df=metadata_df,
        scientific_name='Arabidopsis thaliana',
        batch_effect_alg='sva',
        mapping_rate_cutoff=0.25,
        correlation_threshold=0.3,
        one_outlier_per_iteration=True,
        num_total_runs_species=5,
        num_runs_after_sample_group_filter=3,
        total_runtime_sec=12.3456789,
    )
    assert final_summary.loc[0, 'scientific_name'] == 'Arabidopsis thaliana'
    assert final_summary.loc[0, 'batch_effect_alg'] == 'sva'
    assert int(final_summary.loc[0, 'num_runs_final_kept']) == 1
    assert int(final_summary.loc[0, 'num_runs_final_excluded']) == 2
    assert final_summary.loc[0, 'final_kept_runs'] == 'RUN1'
    assert final_summary.loc[0, 'final_excluded_runs'] == 'RUN2 RUN3'
    assert float(final_summary.loc[0, 'total_runtime_sec']) == 12.345679

    out = write_curation_summaries(
        round_summary=round_summary,
        metadata_df=metadata_df,
        scientific_name='Arabidopsis thaliana',
        batch_effect_alg='sva',
        dir_tsv=str(tmp_path),
        mapping_rate_cutoff=0.25,
        correlation_threshold=0.3,
        one_outlier_per_iteration=True,
        num_total_runs_species=5,
        num_runs_after_sample_group_filter=3,
        total_runtime_sec=12.3456789,
    )
    round_path = tmp_path / 'Arabidopsis_thaliana.sva.curation_round_summary.tsv'
    final_path = tmp_path / 'Arabidopsis_thaliana.sva.curation_final_summary.tsv'
    assert out['round_path'] == str(round_path)
    assert out['final_path'] == str(final_path)
    assert round_path.exists()
    assert final_path.exists()
    loaded_round = pandas.read_csv(round_path, sep='\t')
    loaded_final = pandas.read_csv(final_path, sep='\t')
    pandas.testing.assert_frame_equal(loaded_round, round_summary, check_dtype=False)
    pandas.testing.assert_frame_equal(loaded_final, out['final_summary'], check_dtype=False)


def test_write_curation_summaries_honors_resolved_species_tag(tmp_path):
    metadata_df = pandas.DataFrame({'run': ['RUN1'], 'exclusion': ['no']})

    out = write_curation_summaries(
        round_summary=initialize_round_summary(),
        metadata_df=metadata_df,
        scientific_name='Homo sapiens',
        batch_effect_alg='no',
        dir_tsv=str(tmp_path),
        mapping_rate_cutoff=0.0,
        correlation_threshold=0.3,
        one_outlier_per_iteration=False,
        num_total_runs_species=1,
        num_runs_after_sample_group_filter=1,
        total_runtime_sec=0.0,
        species_tag='human',
    )

    assert out['round_path'] == str(tmp_path / 'human.no.curation_round_summary.tsv')
    assert out['final_path'] == str(tmp_path / 'human.no.curation_final_summary.tsv')
