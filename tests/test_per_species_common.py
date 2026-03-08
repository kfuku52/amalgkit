import numpy
import pandas

from amalgkit.per_species_common import (
    append_round_summary,
    build_curation_final_summary,
    initialize_round_summary,
    sample_group_mean,
    sample_group_to_tau,
    write_curation_summaries,
)


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
        transform_method='fpkm',
    )
    numpy.testing.assert_allclose(
        tau_df['tau'].to_numpy(dtype=float),
        numpy.array([5.0 / 6.0, 5.0 / 6.0, 0.8, 1.0 / 9.0]),
        rtol=0.0,
        atol=1e-12,
    )
    assert tau_df['highest'].tolist() == ['leaf', 'root', 'leaf', 'leaf']
    assert tau_df['order'].tolist() == ['leaf|root', 'root|leaf', 'leaf|root', 'leaf|root']


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
        transform_method='fpkm',
    )
    assert tau_df['highest'].tolist() == ['root', 'root', 'root']
    assert tau_df['order'].tolist() == ['root', 'root', 'root']
    assert tau_df['tau'].isna().all()


def test_sample_group_to_tau_handles_empty_input():
    empty_df = pandas.DataFrame(index=['G1', 'G2'])
    observed = sample_group_to_tau(empty_df, rich_annotation=True, transform_method='log2p1-fpkm')
    assert list(observed.columns) == ['tau', 'highest', 'order']
    assert list(observed.index) == ['G1', 'G2']
    assert observed.shape == (2, 3)
    assert observed.isna().all().all()


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
    pandas.testing.assert_frame_equal(loaded_round, round_summary)
    pandas.testing.assert_frame_equal(loaded_final, out['final_summary'])
