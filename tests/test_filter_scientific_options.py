"""Regression contracts for optional reference/support/iteration policies."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from amalgkit.cross_species_computation import calculate_correlation_within_group
from amalgkit.cross_species_filter import _apply_csfilter_outlier_flags
from amalgkit.csfilter import _normalize_csfilter_metadata_columns
from amalgkit.per_species_python import (
    _compute_sample_group_correlation_metrics,
    _apply_within_group_filter,
    _reduce_outlier_candidates,
    _should_stop_within_group_filter,
)


def fixture():
    metadata = pd.DataFrame(dict(
        run=['a1', 'a2', 'b1', 'b2', 'a3', 'b3'],
        species_tag=['S', 'S', 'S', 'S', 'T', 'T'],
        sample_group=['A', 'A', 'B', 'B', 'A', 'B'],
        bioproject=['p1', 'p2', 'p1', 'p2', 'p3', 'p3'], exclusion='no',
    ))
    counts = pd.DataFrame(dict(
        a1=[1., 2, 4, 3], a2=[1., 3, 2, 4], b1=[4., 3, 2, 1],
        b2=[3., 4, 1, 2], a3=[2., 1, 3, 4], b3=[4., 2, 3, 1],
    ))
    return counts, metadata


def cs_matrix(counts, metadata):
    result = counts.copy()
    result.columns = metadata.species_tag + '_' + metadata.run
    return result


def test_species_exclusion_removes_own_species_from_both_references():
    counts, metadata = fixture()
    result = calculate_correlation_within_group(metadata, cs_matrix(counts, metadata), 'corrected',
                                               reference_exclusion='species')
    assert result.loc[0, 'within_group_cor_corrected'] == pytest.approx(counts.a1.corr(counts.a3))
    assert result.loc[0, 'max_nongroup_cor_corrected'] == pytest.approx(counts.a1.corr(counts.b3))
    # Changing other runs from the target species cannot change the target's score.
    changed = counts.copy()
    changed.loc[:, ['a2', 'b1', 'b2']] *= -1
    updated = calculate_correlation_within_group(metadata, cs_matrix(changed, metadata), 'corrected',
                                                reference_exclusion='species')
    pd.testing.assert_series_equal(result.loc[0], updated.loc[0])


def test_project_exclusion_removes_project_from_both_references():
    counts, metadata = fixture()
    result = _compute_sample_group_correlation_metrics(counts, metadata, ['A', 'B'], 'pearson',
                                                       reference_exclusion='bioproject')
    assert result.loc[0, 'ws_within_group_cor'] == pytest.approx(counts.a1.corr(counts[['a2', 'a3']].mean(axis=1)))
    assert result.loc[0, 'ws_max_nongroup_cor'] == pytest.approx(counts.a1.corr(counts[['b2', 'b3']].mean(axis=1)))
    metadata.loc[0, 'bioproject'] = ''
    with pytest.raises(ValueError, match='nonempty'):
        _compute_sample_group_correlation_metrics(counts, metadata, ['A', 'B'], 'pearson',
                                                  reference_exclusion='bioproject')


@pytest.mark.parametrize('command', ['ws', 'cs'])
@pytest.mark.parametrize('threshold, scoreable', [(2, True), (3, True), (4, False)])
def test_pair_count_cutoff_is_inclusive_and_never_imputes(command, threshold, scoreable):
    counts, metadata = fixture()
    counts.loc[3, 'a1'] = np.nan
    if command == 'cs':
        result = calculate_correlation_within_group(metadata, cs_matrix(counts, metadata), 'corrected',
                                                   min_common_genes=threshold)
        result = _normalize_csfilter_metadata_columns(result)
        score, count = 'within_group_cor', 'cs_within_common_genes'
    else:
        result = _compute_sample_group_correlation_metrics(counts, metadata, ['A', 'B'], 'pearson',
                                                           min_common_genes=threshold)
        score, count = 'ws_within_group_cor', 'ws_within_common_genes'
    assert result.loc[0, count] == 3
    assert bool(np.isfinite(result.loc[0, score])) == scoreable


@pytest.mark.parametrize('command', ['ws', 'cs'])
def test_under_supported_competitor_cannot_be_silently_dropped(command):
    counts, metadata = fixture()
    metadata.loc[5, 'sample_group'] = 'C'
    counts.loc[2:, 'b3'] = np.nan
    if command == 'cs':
        result = calculate_correlation_within_group(metadata, cs_matrix(counts, metadata), 'corrected',
                                                   min_common_genes=3)
        assert np.isnan(result.loc[0, 'max_nongroup_cor_corrected'])
        assert result.loc[0, 'min_nongroup_common_genes_corrected'] == 2
    else:
        result = _compute_sample_group_correlation_metrics(counts, metadata, ['A', 'B', 'C'], 'pearson',
                                                           min_common_genes=3)
        assert np.isnan(result.loc[0, 'ws_margin'])
        assert result.loc[0, 'ws_min_nongroup_common_genes'] == 2


def test_species_z_preserves_consistent_species_shift_but_reports_scope():
    margins = [-.9, -1., -.95, .70, .75, .8, .85, .90, .95, 1.]
    metadata = pd.DataFrame(dict(
        species_tag=['A'] * 3 + ['B'] * 7, sample_group='leaf', exclusion='no',
        within_group_cor_corrected=margins, max_nongroup_cor_corrected=0.,
        within_group_cor_uncorrected=margins, max_nongroup_cor_uncorrected=0.,
    ))
    pooled = _apply_csfilter_outlier_flags(metadata, outlier_method='robust_margin')
    within = _apply_csfilter_outlier_flags(metadata, outlier_method='robust_margin', robust_z_scope='species_group')
    assert pooled.loc[:2, 'cs_outlier_candidate'].all()
    assert not within['cs_outlier_candidate'].any()
    assert within['cs_robust_z_scope'].eq('species_group').all()
    # Species labels can include delimiters without colliding with group labels.
    metadata['species_tag'] = ['A|leaf'] * 3 + ['A'] * 7
    metadata['sample_group'] = ['X'] * 3 + ['leaf|X'] * 7
    split = _apply_csfilter_outlier_flags(metadata, outlier_method='robust_margin', robust_z_scope='species_group')
    assert not split['cs_outlier_candidate'].any()


def test_species_exclusion_with_no_other_species_is_unscoreable_not_outlier():
    counts, metadata = fixture()
    metadata['species_tag'] = 'S'
    scores = calculate_correlation_within_group(metadata, cs_matrix(counts, metadata), 'corrected',
                                                reference_exclusion='species')
    scores['within_group_cor_uncorrected'] = scores['within_group_cor_corrected']
    scores['max_nongroup_cor_uncorrected'] = scores['max_nongroup_cor_corrected']
    result = _apply_csfilter_outlier_flags(scores, outlier_method='robust_margin')
    assert result['within_common_genes_corrected'].eq(0).all()
    assert result['cs_margin_corrected'].isna().all()
    assert result['exclusion'].eq('no').all()


def test_single_pass_stops_despite_remaining_candidates_and_validates_limit():
    before = pd.DataFrame(columns=['a', 'b', 'c'])
    after = pd.DataFrame(columns=['b', 'c'])
    assert not _should_stop_within_group_filter(before, after, ['a'], 1)
    assert _should_stop_within_group_filter(before, after, ['a'], 1, 1)
    with pytest.raises(ValueError, match='positive'):
        _should_stop_within_group_filter(before, after, ['a'], 1, 0)


def test_one_per_iteration_enforces_both_constraints_independent_of_order():
    candidates = pd.DataFrame(dict(run=['r1', 'r2', 'r3', 'r4'],
                                   sample_group=['A', 'A', 'B', 'B'],
                                   bioproject=['P', 'Q', 'P', 'Q'], ws_margin=[-1., -.8, -.7, -.6]))
    assert _reduce_outlier_candidates(candidates) == ['r1', 'r4']
    assert _reduce_outlier_candidates(candidates.iloc[::-1]) == ['r1', 'r4']


def test_support_threshold_retains_sample_and_preserves_manual_exclusion():
    counts, metadata = fixture()
    metadata.loc[5, 'exclusion'] = 'manual'
    counts = counts.drop(columns='b3')
    _, result, excluded = _apply_within_group_filter(counts, metadata, SimpleNamespace(min_common_genes=5), ['A', 'B'])
    assert excluded == []
    assert result.loc[5, 'exclusion'] == 'manual'
    assert result.loc[:4, 'ws_min_common_genes'].eq(5).all()
    assert pd.isna(result.loc[5, 'ws_min_common_genes'])


def test_new_options_reach_wrappers_without_species_policy_leaking_into_prepare(tmp_path):
    from amalgkit.cli_entry import build_main_parser
    from amalgkit.wsfilter import _build_per_species_args
    from amalgkit.csfilter import _build_prepare_per_species_args, _build_cross_species_args

    parser = build_main_parser()
    args = parser.parse_args(['wsfilter', '--min_common_genes', '50', '--reference_exclusion', 'bioproject',
                              '--small_group_policy', 'retain', '--max_filter_iterations', '1'])
    worker = _build_per_species_args(args, str(tmp_path), str(tmp_path))
    assert (worker.min_common_genes, worker.reference_exclusion, worker.max_filter_iterations) == (50, 'bioproject', 1)
    args = parser.parse_args(['csfilter', '--min_common_genes', '100', '--reference_exclusion', 'species',
                              '--robust_z_scope', 'species_group', '--small_group_policy', 'retain'])
    prepare = _build_prepare_per_species_args(args, str(tmp_path), str(tmp_path))
    cross = _build_cross_species_args(args, str(tmp_path), 'A,B', 50)
    assert prepare.reference_exclusion == 'run'
    assert (cross.min_common_genes, cross.reference_exclusion, cross.robust_z_scope) == (100, 'species', 'species_group')
    assert cross.small_group_policy == 'retain'


def test_rerun_metrics_override_stale_canonical_columns():
    metadata = pd.DataFrame(dict(cs_margin=[-1.], cs_margin_corrected=[.3],
                                 within_group_cor=[-.8], within_group_cor_corrected=[.7],
                                 cs_within_common_genes=[2], within_common_genes_corrected=[200]))
    output = _normalize_csfilter_metadata_columns(metadata)
    assert output.loc[0, 'cs_margin'] == .3
    assert output.loc[0, 'within_group_cor'] == .7
    assert output.loc[0, 'cs_within_common_genes'] == 200


def test_default_ws_margin_preserves_pandas_rounding_at_zero_boundary():
    counts = pd.DataFrame(dict(
        a1=[5.12588318783472, 15.364690481856902, np.nan],
        a2=[5.561957608630067, 14.981829500375296, 8.706798092670898],
        b1=[6.663754620053448, 13.07409070145881, 9.678133698831825],
        b2=[7.144597967048354, 12.607317411772504, 9.709964987151759],
    ))
    metadata = pd.DataFrame(dict(run=counts.columns, sample_group=['A', 'A', 'B', 'B'], exclusion='no'))
    expected = counts.a1.corr(counts.a2) - counts.a1.corr(counts[['b1', 'b2']].mean(axis=1))
    result = _compute_sample_group_correlation_metrics(counts, metadata, ['A', 'B'], 'pearson')
    assert expected < 0
    assert result.loc[0, 'ws_margin'] == expected


def test_ws_rescoring_preserves_excluded_run_evidence_and_its_settings():
    counts, metadata = fixture()
    metadata.loc[0, 'exclusion'] = 'manual'
    metadata['ws_margin'] = -.7
    metadata['ws_min_common_genes'] = 20
    metadata['ws_margin_threshold'] = -.2
    metadata['ws_reference_exclusion'] = 'run'
    metadata['ws_small_group_policy'] = 'margin_fallback'
    metadata.loc[0, 'bioproject'] = ''
    # Even when the excluded run is still present in the matrix, it must not
    # enter references, trigger project validation, or receive the new settings.
    kept, result, _ = _apply_within_group_filter(
        counts, metadata, SimpleNamespace(min_common_genes=50, reference_exclusion='bioproject',
                                          small_group_policy='retain'), ['A', 'B'])
    assert 'a1' not in kept
    assert result.loc[0, 'exclusion'] == 'manual'
    assert result.loc[0, 'ws_margin'] == -.7
    assert result.loc[0, 'ws_min_common_genes'] == 20
    assert result.loc[0, 'ws_margin_threshold'] == -.2
    assert result.loc[1, 'ws_margin_threshold'] == 0.
    assert result.loc[0, 'ws_reference_exclusion'] == 'run'
    assert result.loc[1, 'ws_min_common_genes'] == 50
    assert pd.isna(result.loc[1, 'ws_margin'])


def test_cs_rescoring_preserves_excluded_canonical_scores_and_settings():
    counts, metadata = fixture()
    metadata.loc[0, 'exclusion'] = 'low_cross_species_group_correlation'
    metadata['within_group_cor'] = -.5
    metadata['max_nongroup_cor'] = .2
    metadata['cs_margin'] = -.7
    metadata['cs_robust_z'] = -3.
    metadata['cs_min_common_genes'] = 20
    metadata['cs_margin_threshold'] = -.2
    metadata['cs_reference_exclusion'] = 'run'
    metadata['cs_robust_z_scope'] = 'sample_group'
    matrix = cs_matrix(counts, metadata)
    result = calculate_correlation_within_group(metadata, matrix, 'corrected', 'species', 50)
    result = _apply_csfilter_outlier_flags(result, outlier_method='robust_margin', robust_z_scope='species_group')
    result = _normalize_csfilter_metadata_columns(result)
    assert result.loc[0, 'cs_margin'] == -.7
    assert result.loc[0, 'cs_robust_z'] == -3.
    assert result.loc[0, 'cs_reference_exclusion'] == 'run'
    assert result.loc[0, 'cs_robust_z_scope'] == 'sample_group'
    assert result.loc[0, 'cs_min_common_genes'] == 20
    assert result.loc[0, 'cs_margin_threshold'] == -.2
    assert result.loc[1, 'cs_margin_threshold'] == 0.
    assert pd.isna(result.loc[1, 'cs_margin'])
    assert result.loc[1, 'cs_robust_z_scope'] == 'species_group'


def test_species_z_does_not_make_blank_groups_scoreable():
    frame = pd.DataFrame(dict(species_tag=['S', 'S'], sample_group=['', '  '], exclusion='no',
                              within_group_cor_corrected=[-.5, -.6], max_nongroup_cor_corrected=0.))
    result = _apply_csfilter_outlier_flags(frame, outlier_method='robust_margin', robust_z_scope='species_group')
    assert result['cs_robust_z'].isna().all()
    assert not result['cs_small_group'].any()
    assert result['exclusion'].eq('no').all()


@pytest.mark.parametrize('exclusion', ['no', ' No ', 'NO'])
def test_ws_active_flags_keep_established_case_and_whitespace_normalization(exclusion):
    counts, metadata = fixture()
    metadata['exclusion'] = exclusion
    kept, result, _ = _apply_within_group_filter(counts, metadata, SimpleNamespace(min_common_genes=50), ['A', 'B'])
    assert kept.shape[1] == 6
    assert result['ws_min_common_genes'].eq(50).all()


def test_standardized_missing_project_is_rejected_and_not_a_shared_project():
    from amalgkit.per_species_finalize_python import _standardize_metadata_all
    counts, metadata = fixture()
    metadata.loc[0, 'bioproject'] = None
    metadata = _standardize_metadata_all(metadata)
    with pytest.raises(ValueError, match='nonempty bioproject'):
        _compute_sample_group_correlation_metrics(counts, metadata, ['A', 'B'], 'pearson',
                                                  reference_exclusion='bioproject')
    candidates = pd.DataFrame(dict(run=['a', 'b'], sample_group=['A', 'B'],
                                   bioproject=['not_provided', 'not_provided'], ws_margin=[-.8, -.7]))
    assert _reduce_outlier_candidates(candidates) == ['a', 'b']


def test_mapping_filter_preserves_prior_exclusions_and_does_not_reintroduce_them():
    from amalgkit.per_species_python import _filter_low_mapping_rate
    counts, metadata = fixture()
    metadata['mapping_rate'] = [0., 1., 0., 1., 1., 1.]
    metadata.loc[0, 'exclusion'] = 'manual'
    metadata.loc[1, 'exclusion'] = 'previous_filter'
    kept, result, removed = _filter_low_mapping_rate(counts, metadata, .2)
    assert removed == ['b1']
    assert result.loc[0, 'exclusion'] == 'manual'
    assert result.loc[1, 'exclusion'] == 'previous_filter'
    assert set(kept.columns) == {'b2', 'a3', 'b3'}


@pytest.mark.parametrize('command', ['ws', 'cs'])
def test_pairwise_support_does_not_require_a_shared_gene_set(command):
    # Each comparison has two pairs, but their common intersection has only one.
    counts = pd.DataFrame(dict(a1=[1., 2., 3.], a2=[np.nan, 2., 3.],
                               b1=[1., np.nan, 3.], b2=[1., np.nan, 3.]))
    metadata = pd.DataFrame(dict(run=counts.columns, species_tag=['S', 'T', 'S', 'T'],
                                 sample_group=['A', 'A', 'B', 'B'], exclusion='no'))
    if command == 'cs':
        result = calculate_correlation_within_group(metadata, cs_matrix(counts, metadata), 'corrected',
                                                    min_common_genes=2)
        assert result.loc[0, 'within_common_genes_corrected'] == 2
        assert result.loc[0, 'min_nongroup_common_genes_corrected'] == 2
        assert np.isfinite(result.loc[0, 'within_group_cor_corrected'])
        assert np.isfinite(result.loc[0, 'max_nongroup_cor_corrected'])
    else:
        result = _compute_sample_group_correlation_metrics(counts, metadata, ['A', 'B'], 'pearson', min_common_genes=2)
        assert result.loc[0, 'ws_within_common_genes'] == 2
        assert result.loc[0, 'ws_min_nongroup_common_genes'] == 2
        assert np.isfinite(result.loc[0, 'ws_margin'])
