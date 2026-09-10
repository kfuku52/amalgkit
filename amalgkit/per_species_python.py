import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import pandas

from amalgkit.batch_effect_common import (
    annotate_metadata_with_batch_info,
    initialize_batch_info,
    write_batch_effect_summary_tsv,
)
from amalgkit.outlier_utils import flag_margin_outliers
from amalgkit.per_species_common import (
    _is_non_excluded_flag,
    append_round_summary,
    initialize_round_summary,
    sample_group_mean,
    tau_options_from_args,
    write_tau_outputs,
    write_curation_summaries,
)
from amalgkit.per_species_finalize_python import (
    _apply_transformation_logic,
    load_quant_model_table,
    resolve_length_models,
    _exclude_inappropriate_sample_from_eff_length,
    _exclude_inappropriate_sample_from_tc,
    _get_species_metadata,
    _normalize_dataframe_columns,
    _normalize_metadata_df,
    _read_expression_tsv,
    _resolve_scientific_name,
    _resolve_selected_sample_groups,
    _sort_tc_and_metadata,
    _standardize_metadata_all,
    run_finalize_python_worker,
    should_use_python_finalize_worker,
    record_expression_library_sizes,
)
from amalgkit.per_species_outputs import (
    initialize_correlation_statistics,
    save_correlation_statistics,
    save_state_overview_pdf,
    save_tau_histogram_pdf,
    write_table_with_index_name,
)


def should_use_python_per_species_worker(args):
    requested_mode = str(getattr(args, 'worker_mode', 'prepare_tables'))
    if requested_mode in {'prepare_tables', 'wsfilter'}:
        return True
    if requested_mode == 'finalize':
        return should_use_python_finalize_worker(args)
    return False


def _resolve_species_input_paths(input_dir_abs, species_tag):
    species_dir = os.path.join(input_dir_abs, species_tag)
    count_path_candidates = [
        os.path.join(species_dir, species_tag + '_cstmm_counts.tsv'),
        os.path.join(species_dir, species_tag + '_est_counts.tsv'),
    ]
    count_path = next((path for path in count_path_candidates if os.path.isfile(path)), count_path_candidates[0])
    eff_length_path = os.path.join(species_dir, species_tag + '_eff_length.tsv')
    return count_path, eff_length_path


def _filter_low_mapping_rate(tc, sra, mapping_rate_cutoff):
    if float(mapping_rate_cutoff) <= 0:
        return tc.copy(), sra.copy(), []
    if 'mapping_rate' not in sra.columns:
        warnings.warn('mapping_rate column is missing; low mapping-rate filtering will be skipped.')
        return tc.copy(), sra.copy(), []
    mapping_rate = pandas.to_numeric(sra.loc[:, 'mapping_rate'], errors='coerce')
    is_mapping_good = mapping_rate >= float(mapping_rate_cutoff)
    is_mapping_good = is_mapping_good.fillna(False)
    is_active = _is_non_excluded_flag(sra['exclusion']) & sra['run'].astype(str).isin(tc.columns)
    newly_excluded = is_active & ~is_mapping_good
    excluded_runs = sra.loc[newly_excluded, 'run'].astype(str).tolist()
    out_sra = sra.copy()
    if len(excluded_runs) > 0:
        out_sra.loc[newly_excluded, 'exclusion'] = 'low_mapping_rate'
    keep_runs = [run_id for run_id in out_sra.loc[is_active & is_mapping_good, 'run'].astype(str).tolist() if run_id in tc.columns]
    return tc.loc[:, keep_runs].copy(), out_sra, excluded_runs


def _finite_pair_correlation(left, right, method):
    valid = numpy.isfinite(left) & numpy.isfinite(right)
    left, right = left.loc[valid], right.loc[valid]
    if len(left) < 2 or left.nunique() <= 1 or right.nunique() <= 1:
        return numpy.nan
    # Preserve the established pandas arithmetic: even tiny rounding changes
    # near a zero margin can change the strict absolute-threshold decision.
    return left.corr(right, method=method)


def _compute_sample_group_correlation_metrics(tc, sra, selected_sample_groups, dist_method,
                                              min_common_genes=0, reference_exclusion='run'):
    dist_method = str(dist_method).strip().lower()
    if dist_method not in {'pearson', 'spearman', 'kendall'}:
        raise ValueError(
            'Unsupported correlation method: {}. Choose pearson, spearman, or kendall.'.format(
                dist_method
            )
        )
    is_active = _is_non_excluded_flag(sra['exclusion'])
    if min_common_genes < 0 or min_common_genes == 1:
        raise ValueError('min_common_genes must be 0 (disabled) or at least 2')
    if reference_exclusion not in {'run', 'bioproject'}:
        raise ValueError('reference_exclusion must be run or bioproject')
    if reference_exclusion == 'bioproject':
        if 'bioproject' not in sra.columns:
            raise ValueError('bioproject reference exclusion requires bioproject metadata')
        active = sra.loc[is_active & sra['run'].astype(str).isin(tc.columns.astype(str)), 'bioproject']
        if active.fillna('').astype(str).str.strip().isin(['', 'not_provided']).any():
            raise ValueError('bioproject reference exclusion requires nonempty bioproject metadata')
    active_ids = sra.loc[_is_non_excluded_flag(sra['exclusion']), 'run'].astype(str)
    tc = tc.loc[:, tc.columns.astype(str).isin(active_ids)]
    required_sample_groups = [group for group in selected_sample_groups if group in set(sra.loc[:, 'sample_group'].astype(str))]
    out = sra.copy()
    for column in [
        'ws_within_group_cor',
        'ws_max_nongroup_cor',
        'ws_margin',
        'ws_robust_z',
        'ws_outlier_candidate',
        'ws_small_group',
    ]:
        out.loc[:, column] = False if column in {'ws_outlier_candidate', 'ws_small_group'} else numpy.nan
    out['ws_within_common_genes'] = 0
    out['ws_min_nongroup_common_genes'] = 0
    if (len(required_sample_groups) <= 1) or (tc.shape[1] == 0) or (out.shape[0] == 0):
        return out
    tc_ave = sample_group_mean(tc, out, required_sample_groups)['tc_ave']
    if tc_ave.shape[1] <= 1:
        return out
    corr_by_group = pandas.DataFrame(
        numpy.nan,
        index=tc.columns.astype(str),
        columns=tc_ave.columns.astype(str),
        dtype=float,
    )
    tc_numeric = tc.apply(pandas.to_numeric, errors='coerce').to_numpy(dtype=float)
    tc_ave_numeric = tc_ave.apply(pandas.to_numeric, errors='coerce').to_numpy(dtype=float)
    common_by_group = pandas.DataFrame(
        numpy.isfinite(tc_numeric).astype(numpy.int64).T
        @ numpy.isfinite(tc_ave_numeric).astype(numpy.int64),
        index=corr_by_group.index, columns=corr_by_group.columns,
    )
    can_vectorize = (
        str(dist_method).lower() == 'pearson'
        and numpy.isfinite(tc_numeric).all()
        and numpy.isfinite(tc_ave_numeric).all()
    )
    if can_vectorize:
        tc_centered = tc_numeric - numpy.mean(tc_numeric, axis=0, keepdims=True)
        tc_ave_centered = tc_ave_numeric - numpy.mean(tc_ave_numeric, axis=0, keepdims=True)
        numerators = tc_centered.T @ tc_ave_centered
        denominators = numpy.sqrt(
            numpy.sum(tc_centered ** 2, axis=0).reshape(-1, 1)
            * numpy.sum(tc_ave_centered ** 2, axis=0).reshape(1, -1)
        )
        corr_by_group.loc[:, :] = numpy.divide(
            numerators,
            denominators,
            out=numpy.full_like(numerators, numpy.nan),
            where=denominators > 0,
        )
    else:
        for run_id in corr_by_group.index:
            sample_values = pandas.to_numeric(tc.loc[:, run_id], errors='coerce')
            for sample_group in corr_by_group.columns:
                corr_value = _finite_pair_correlation(
                    sample_values, pandas.to_numeric(tc_ave.loc[:, sample_group], errors='coerce'),
                    method=dist_method,
                )
                corr_by_group.loc[run_id, sample_group] = corr_value

    # A sample must not contribute to the reference profile used to judge that
    # same sample. Besides biasing correlations upward, self-inclusion can make
    # a two-sample group mean constant when the profiles oppose one another,
    # turning an actionable negative margin into NaN. Non-group references do
    # not contain the evaluated sample and can retain the precomputed means.
    group_members = {
        sample_group: [
            run_id
            for run_id in out.loc[
                out['sample_group'].astype(str).eq(sample_group),
                'run',
            ].astype(str).tolist()
            if run_id in tc.columns
        ]
        for sample_group in required_sample_groups
    }
    sample_group_by_run = dict(
        zip(out.loc[:, 'run'].astype(str), out.loc[:, 'sample_group'].astype(str))
    )
    project_by_run = (
        dict(zip(out['run'].astype(str), out['bioproject'].fillna('').astype(str).str.strip()))
        if reference_exclusion == 'bioproject' else {}
    )
    for run_id in corr_by_group.index:
        sample_group = sample_group_by_run.get(str(run_id))
        if sample_group not in corr_by_group.columns:
            continue
        reference_groups = list(corr_by_group.columns) if reference_exclusion == 'bioproject' else [sample_group]
        for reference_group in reference_groups:
            other_runs = [
                other_run for other_run in group_members.get(reference_group, [])
                if other_run != str(run_id)
                and (reference_exclusion == 'run' or project_by_run[other_run] != project_by_run[str(run_id)])
            ]
            if not other_runs:
                corr_by_group.loc[run_id, reference_group] = numpy.nan
                common_by_group.loc[run_id, reference_group] = 0
                continue
            reference = tc.loc[:, other_runs].apply(pandas.to_numeric, errors='coerce').mean(axis=1, skipna=True)
            sample = pandas.to_numeric(tc.loc[:, run_id], errors='coerce')
            valid = numpy.isfinite(sample) & numpy.isfinite(reference)
            common_by_group.loc[run_id, reference_group] = int(valid.sum())
            corr_by_group.loc[run_id, reference_group] = _finite_pair_correlation(sample, reference, method=dist_method)
    run_values = out.loc[:, 'run'].astype(str).tolist()
    sample_group_values = out.loc[:, 'sample_group'].astype(str).tolist()
    within_values = []
    nongroup_values = []
    margin_values = []
    within_counts = []
    nongroup_counts = []
    for run_id, sample_group in zip(run_values, sample_group_values):
        if (run_id not in corr_by_group.index) or (sample_group not in corr_by_group.columns):
            within_values.append(numpy.nan)
            nongroup_values.append(numpy.nan)
            margin_values.append(numpy.nan)
            within_counts.append(0)
            nongroup_counts.append(0)
            continue
        corr_row = pandas.to_numeric(corr_by_group.loc[run_id, :], errors='coerce')
        within_cor = corr_row.get(sample_group, numpy.nan)
        nongroup = pandas.to_numeric(corr_row.loc[corr_row.index != sample_group], errors='coerce').dropna()
        max_nongroup = float(nongroup.max()) if nongroup.shape[0] > 0 else numpy.nan
        within_count = int(common_by_group.loc[run_id, sample_group])
        other_counts = common_by_group.loc[run_id, common_by_group.columns != sample_group]
        other_count = int(other_counts.min()) if len(other_counts) else 0
        within_counts.append(within_count)
        nongroup_counts.append(other_count)
        if within_count < min_common_genes:
            within_cor = numpy.nan
        if other_count < min_common_genes:
            max_nongroup = numpy.nan
        margin_val = within_cor - max_nongroup if numpy.isfinite(within_cor) and numpy.isfinite(max_nongroup) else numpy.nan
        within_values.append(within_cor)
        nongroup_values.append(max_nongroup)
        margin_values.append(margin_val)
    out.loc[:, 'ws_within_group_cor'] = within_values
    out.loc[:, 'ws_max_nongroup_cor'] = nongroup_values
    out.loc[:, 'ws_margin'] = margin_values
    out['ws_within_common_genes'] = within_counts
    out['ws_min_nongroup_common_genes'] = nongroup_counts
    return out


def _reduce_outlier_candidates(candidate_df):
    if candidate_df.shape[0] == 0:
        return []
    candidates = candidate_df.copy()
    candidates['run'] = candidates['run'].fillna('').astype(str)
    sort_cols = ['ws_margin', 'run'] if 'ws_margin' in candidates else ['run']
    candidates = candidates.sort_values(sort_cols, kind='stable', na_position='last')
    selected = []
    used = {column: set() for column in ('bioproject', 'sample_group') if column in candidates}
    for _, row in candidates.iterrows():
        if not row['run']:
            continue
        keys = {column: str(row[column]).strip() if pandas.notna(row[column]) else '' for column in used}
        if keys.get('bioproject') == 'not_provided':
            keys['bioproject'] = ''
        if any(value and value in used[column] for column, value in keys.items()):
            continue
        selected.append(row['run'])
        for column, value in keys.items():
            if value:
                used[column].add(value)
    return selected


def _apply_within_group_filter(tc, sra, args, selected_sample_groups, min_dif=0.0):
    active_ids = sra.loc[_is_non_excluded_flag(sra['exclusion']), 'run'].astype(str)
    tc = tc.loc[:, tc.columns.astype(str).isin(active_ids)]
    out = _compute_sample_group_correlation_metrics(
        tc=tc,
        sra=sra,
        selected_sample_groups=selected_sample_groups,
        dist_method=str(getattr(args, 'dist_method', 'pearson')),
        min_common_genes=int(getattr(args, 'min_common_genes', 0)),
        reference_exclusion=str(getattr(args, 'reference_exclusion', 'run')),
    )
    filtered = flag_margin_outliers(
        df=out,
        margin_col='ws_margin',
        group_col='sample_group',
        margin_threshold=float(getattr(args, 'margin_threshold', 0.0)) + float(min_dif),
        robust_z_threshold=float(getattr(args, 'robust_z_threshold', -2.5)),
        robust_z_col='ws_robust_z',
        outlier_col='ws_outlier_candidate',
        small_group_policy=str(getattr(args, 'small_group_policy', 'margin_fallback')),
        small_group_col='ws_small_group',
    )
    candidate_runs = (
        filtered.loc[
            filtered['ws_outlier_candidate'].fillna(False).astype(bool)
            & filtered['run'].astype(str).isin(tc.columns.astype(str)),
            'run',
        ]
        .fillna('')
        .astype(str)
        .tolist()
    )
    candidate_runs = [run_id for run_id in candidate_runs if run_id != '']
    excluded_runs = list(dict.fromkeys(candidate_runs))
    if bool(getattr(args, 'one_outlier_per_iter', False)) and len(excluded_runs) > 0:
        candidate_df = filtered.loc[
            filtered['run'].astype(str).isin(excluded_runs),
            [col for col in ['run', 'sample_group', 'bioproject', 'ws_margin'] if col in filtered.columns],
        ].copy()
        excluded_runs = _reduce_outlier_candidates(candidate_df)
    out_sra = sra.copy()
    active_rows = out_sra['run'].astype(str).isin(tc.columns.astype(str))
    settings = {
        'ws_min_common_genes': int(getattr(args, 'min_common_genes', 0)),
        'ws_margin_threshold': float(getattr(args, 'margin_threshold', 0.0)) + float(min_dif),
        'ws_robust_z_threshold': float(getattr(args, 'robust_z_threshold', -2.5)),
        'ws_reference_exclusion': str(getattr(args, 'reference_exclusion', 'run')),
        'ws_small_group_policy': str(getattr(args, 'small_group_policy', 'margin_fallback')),
    }
    for column, value in settings.items():
        if column not in out_sra:
            out_sra[column] = pandas.Series(pandas.NA, index=out_sra.index, dtype='object')
        out_sra.loc[active_rows, column] = value
    boolean_metric_cols = {'ws_outlier_candidate', 'ws_small_group'}
    metric_cols = [
        'ws_within_common_genes',
        'ws_min_nongroup_common_genes',
        'ws_within_group_cor',
        'ws_max_nongroup_cor',
        'ws_margin',
        'ws_robust_z',
        'ws_outlier_candidate',
        'ws_small_group',
    ]
    for metric_col in metric_cols:
        if metric_col not in out_sra.columns:
            out_sra.loc[:, metric_col] = False if metric_col in boolean_metric_cols else numpy.nan
        run_map = filtered.set_index(filtered['run'].astype(str))[metric_col]
        # Preserve the scoring evidence from the removal round for inactive runs.
        out_sra.loc[active_rows, metric_col] = out_sra.loc[active_rows, 'run'].astype(str).map(run_map)
        if metric_col in boolean_metric_cols:
            out_sra.loc[:, metric_col] = out_sra.loc[:, metric_col].fillna(False).astype(bool)
    if len(excluded_runs) > 0:
        out_sra.loc[out_sra['run'].astype(str).isin(excluded_runs), 'exclusion'] = 'low_within_sample_group_correlation'
    out_tc = tc.loc[:, [run_id for run_id in tc.columns if run_id not in set(excluded_runs)]].copy()
    return out_tc, out_sra, excluded_runs


def _should_stop_within_group_filter(current_tc, next_tc, excluded_runs, completed_iterations=0, max_iterations=None):
    if max_iterations is not None:
        if max_iterations < 1:
            raise ValueError('max_filter_iterations must be positive')
        if completed_iterations >= max_iterations:
            return True
    if len(excluded_runs) == 0:
        return True
    if next_tc.shape[1] == 0:
        return True
    return list(next_tc.columns) == list(current_tc.columns)


def _save_ws_scatter_plot(metadata_df, out_pdf_path, font_size=8):
    required_cols = {'ws_within_group_cor', 'ws_max_nongroup_cor'}
    if not required_cols.issubset(metadata_df.columns):
        return None
    plot_df = metadata_df.copy()
    plot_df.loc[:, 'ws_within_group_cor'] = pandas.to_numeric(plot_df.loc[:, 'ws_within_group_cor'], errors='coerce')
    plot_df.loc[:, 'ws_max_nongroup_cor'] = pandas.to_numeric(plot_df.loc[:, 'ws_max_nongroup_cor'], errors='coerce')
    plot_df = plot_df.loc[
        plot_df['ws_within_group_cor'].notna() & plot_df['ws_max_nongroup_cor'].notna(),
        :,
    ].copy()
    if plot_df.shape[0] == 0:
        return None
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    groups = plot_df.loc[:, 'sample_group'].fillna('').astype(str).tolist()
    unique_groups = list(dict.fromkeys(groups))
    cmap = plt.get_cmap('tab20')
    color_map = {group: cmap(idx % max(1, cmap.N)) for idx, group in enumerate(unique_groups)}
    colors = [color_map[group] for group in groups]
    is_outlier = plot_df.get('ws_outlier_candidate', pandas.Series(False, index=plot_df.index)).fillna(False).astype(bool).to_numpy()
    ax.scatter(
        plot_df['ws_max_nongroup_cor'].to_numpy(dtype=float),
        plot_df['ws_within_group_cor'].to_numpy(dtype=float),
        c=colors,
        s=numpy.where(is_outlier, 60.0, 35.0),
        edgecolors=numpy.where(is_outlier, 'red', 'black'),
        linewidths=numpy.where(is_outlier, 1.2, 0.4),
        alpha=0.85,
    )
    ax.set_xlabel('ws_max_nongroup_cor', fontsize=font_size)
    ax.set_ylabel('ws_within_group_cor', fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size)
    ax.grid(color='#d0d0d0', linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_pdf_path)
    plt.close(fig)
    return out_pdf_path


def _write_prepare_outputs(
    args,
    species_tag,
    scientific_name,
    dir_tsv,
    dir_pdf,
    tc_before_filter,
    tc_final,
    sra_out,
    selected_sample_groups,
    round_summary,
    correlation_statistics,
    num_total_runs_species,
    num_runs_after_sample_group_filter,
):
    batch_effect_alg = 'no'
    tc_sample_group_uncorrected = sample_group_mean(tc_before_filter, sra_out, selected_sample_groups)['tc_ave']
    tc_sample_group_final = sample_group_mean(tc_final, sra_out, selected_sample_groups)['tc_ave']
    write_table_with_index_name(
        df=tc_before_filter,
        file_path=os.path.join(dir_tsv, '{}.uncorrected.tc.tsv'.format(species_tag)),
        index_name='target_id',
    )
    write_table_with_index_name(
        df=tc_sample_group_uncorrected,
        file_path=os.path.join(dir_tsv, '{}.uncorrected.sample_group.mean.tsv'.format(species_tag)),
        index_name='target_id',
    )
    write_table_with_index_name(
        df=tc_final,
        file_path=os.path.join(dir_tsv, '{}.{}.tc.tsv'.format(species_tag, batch_effect_alg)),
        index_name='target_id',
    )
    write_table_with_index_name(
        df=tc_sample_group_final,
        file_path=os.path.join(dir_tsv, '{}.{}.sample_group.mean.tsv'.format(species_tag, batch_effect_alg)),
        index_name='target_id',
    )
    tau_linear_mean = write_tau_outputs(
        tc_final, sra_out, selected_sample_groups,
        str(getattr(args, 'norm', 'log2p1-fpkm')), dir_tsv, species_tag, batch_effect_alg,
        **tau_options_from_args(args),
    )
    correlation_statistics.to_csv(
        os.path.join(dir_tsv, '{}.{}.correlation_statistics.tsv'.format(species_tag, batch_effect_alg)),
        sep='\t',
    )
    write_curation_summaries(
        round_summary=round_summary,
        metadata_df=sra_out,
        scientific_name=scientific_name,
        batch_effect_alg=batch_effect_alg,
        dir_tsv=dir_tsv,
        mapping_rate_cutoff=float(getattr(args, 'mapping_rate', 0.0)),
        correlation_threshold=float(getattr(args, 'correlation_threshold', 0.3)),
        one_outlier_per_iteration=bool(getattr(args, 'one_outlier_per_iter', False)),
        num_total_runs_species=num_total_runs_species,
        num_runs_after_sample_group_filter=num_runs_after_sample_group_filter,
        total_runtime_sec=0.0,
        species_tag=species_tag,
    )
    save_tau_histogram_pdf(
        counts_df=tc_final,
        metadata_df=sra_out,
        selected_sample_groups=selected_sample_groups,
        out_pdf_path=os.path.join(dir_pdf, '{}.tau_histogram.no.pdf'.format(species_tag)),
        transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
        linear_mean_df=tau_linear_mean,
    )
    _save_ws_scatter_plot(
        metadata_df=sra_out,
        out_pdf_path=os.path.join(dir_pdf, '{}.within_group_correlation.no.pdf'.format(species_tag)),
    )
    batch_info = initialize_batch_info(run_ids=sra_out.loc[:, 'run'].astype(str).tolist(), batch_effect_alg='no')
    batch_info['skip_reason'] = 'batch_effect_alg_no'
    batch_info['batch_effect_alg_applied'] = 'no'
    batch_info['corrected_runs'] = []
    batch_info['uncorrected_runs'] = list(tc_final.columns)
    write_batch_effect_summary_tsv(
        batch_info=batch_info,
        scientific_name=scientific_name,
        species_tag=species_tag,
        dir_tsv=dir_tsv,
        random_seed_value=getattr(args, 'seed', None),
    )


def _run_prepare_or_wsfilter_python_worker(args, metadata, species_tag, input_dir):
    input_dir_abs = os.path.abspath(input_dir)
    count_path, eff_length_path = _resolve_species_input_paths(input_dir_abs=input_dir_abs, species_tag=species_tag)
    needs_lengths = str(args.norm).split('-')[-1] in {'fpkm', 'tpm'}
    if not os.path.isfile(count_path) or (needs_lengths and not os.path.isfile(eff_length_path)):
        return 1

    counts_df = _normalize_dataframe_columns(_read_expression_tsv(count_path))
    eff_length_df = (_normalize_dataframe_columns(_read_expression_tsv(eff_length_path))
                     if needs_lengths else pandas.DataFrame())
    metadata_all = _standardize_metadata_all(_normalize_metadata_df(metadata.df))
    scientific_name = _resolve_scientific_name(metadata_all, species_tag)
    selected_sample_groups = _resolve_selected_sample_groups(args, metadata_all)
    num_total_runs_species = int(metadata_all.loc[:, 'scientific_name'].astype(str).eq(scientific_name).sum())
    sra = _get_species_metadata(metadata_all, scientific_name, selected_sample_groups, counts_df.columns)
    num_runs_after_sample_group_filter = int(sra.shape[0])
    sra = record_expression_library_sizes(count_path, counts_df, sra, args.norm)

    out_dir = os.path.realpath(args.out_dir)
    dir_per_species = os.path.join(out_dir, 'per_species')
    dir_pdf = os.path.join(dir_per_species, species_tag, 'plots')
    dir_tsv = os.path.join(dir_per_species, species_tag, 'tables')
    os.makedirs(dir_pdf, exist_ok=True)
    os.makedirs(dir_tsv, exist_ok=True)

    tc = _exclude_inappropriate_sample_from_tc(counts_df, sra)
    sorted_out = _sort_tc_and_metadata(tc, sra)
    tc = sorted_out['tc']
    sra = sorted_out['sra']
    eff_length_species = _exclude_inappropriate_sample_from_eff_length(eff_length_df, tc)
    length_models = resolve_length_models(
        list(tc.columns),
        load_quant_model_table(os.path.join(os.path.dirname(count_path), species_tag + '_quant_model.tsv')),
    )
    tc_original = _apply_transformation_logic(tc, eff_length_species, args.norm, 'no', 'before_batch', sra, length_models)
    correlation_statistics = save_correlation_statistics(
        counts_df=tc_original,
        metadata_df=sra,
        dist_method=str(getattr(args, 'dist_method', 'pearson')),
        round_value=0,
        correlation_statistics=initialize_correlation_statistics(),
    )
    save_state_overview_pdf(
        counts_df=tc_original,
        metadata_df=sra,
        selected_sample_groups=selected_sample_groups,
        out_pdf_path=os.path.join(dir_pdf, '{}.0.original.pdf'.format(species_tag)),
        dist_method=str(getattr(args, 'dist_method', 'pearson')),
        transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
        font_size=8,
        tau_options=tau_options_from_args(args),
    )

    tc, sra, mapping_excluded_runs = _filter_low_mapping_rate(
        tc=tc,
        sra=sra,
        mapping_rate_cutoff=float(getattr(args, 'mapping_rate', 0.0)),
    )
    tc = _apply_transformation_logic(tc, eff_length_species, args.norm, 'no', 'before_batch', sra, length_models)
    tc_before_filter = tc.copy()
    correlation_statistics = save_correlation_statistics(
        counts_df=tc_before_filter,
        metadata_df=sra,
        dist_method=str(getattr(args, 'dist_method', 'pearson')),
        round_value=1,
        correlation_statistics=correlation_statistics,
    )
    save_state_overview_pdf(
        counts_df=tc_before_filter,
        metadata_df=sra,
        selected_sample_groups=selected_sample_groups,
        out_pdf_path=os.path.join(dir_pdf, '{}.1.mapping_cutoff.pdf'.format(species_tag)),
        dist_method=str(getattr(args, 'dist_method', 'pearson')),
        transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
        font_size=8,
        tau_options=tau_options_from_args(args),
    )

    round_summary = initialize_round_summary()
    if len(mapping_excluded_runs) > 0:
        round_summary = append_round_summary(
            round_summary=round_summary,
            step='mapping_rate_filter',
            round_value=-1,
            reason='low_mapping_rate',
            runs_before=tc_before_filter.columns.tolist() + mapping_excluded_runs,
            runs_after=tc_before_filter.columns.tolist(),
        )
    if bool(getattr(args, 'skip_curation', False)):
        round_summary = append_round_summary(
            round_summary=round_summary,
            step='skip_curation',
            round_value=0,
            reason='skip_curation_requested',
            runs_before=tc_before_filter.columns,
            runs_after=tc_before_filter.columns,
        )
        batch_info = initialize_batch_info(run_ids=sra.loc[:, 'run'].astype(str).tolist(), batch_effect_alg='no')
        batch_info['skip_reason'] = 'skip_curation_requested'
        sra_out = annotate_metadata_with_batch_info(sra, batch_info)
        sra_out.to_csv(os.path.join(dir_tsv, '{}.metadata.tsv'.format(species_tag)), sep='\t', index=False)
        _write_prepare_outputs(
            args=args,
            species_tag=species_tag,
            scientific_name=scientific_name,
            dir_tsv=dir_tsv,
            dir_pdf=dir_pdf,
            tc_before_filter=tc_before_filter,
            tc_final=tc_before_filter,
            sra_out=sra_out,
            selected_sample_groups=selected_sample_groups,
            round_summary=round_summary,
            correlation_statistics=correlation_statistics,
            num_total_runs_species=num_total_runs_species,
            num_runs_after_sample_group_filter=num_runs_after_sample_group_filter,
        )
        return 0

    current_tc = tc.copy()
    current_sra = sra.copy()
    round_index = 0
    max_iterations = getattr(args, 'max_filter_iterations', None)
    if max_iterations is not None and max_iterations < 1:
        raise ValueError('max_filter_iterations must be positive')
    while True:
        next_tc, next_sra, excluded_runs = _apply_within_group_filter(
            tc=current_tc,
            sra=current_sra,
            args=args,
            selected_sample_groups=selected_sample_groups,
            min_dif=0.0,
        )
        reason = 'low_within_sample_group_correlation' if len(excluded_runs) > 0 else 'no_outlier_detected'
        round_summary = append_round_summary(
            round_summary=round_summary,
            step='within_group_filter',
            round_value=round_index,
            reason=reason,
            runs_before=current_tc.columns,
            runs_after=next_tc.columns,
        )
        limit_reached = max_iterations is not None and round_index + 1 >= max_iterations
        round_value = round_index + 2
        if limit_reached or (len(excluded_runs) == 0) or bool(getattr(args, 'plot_intermediate', False)):
            correlation_statistics = save_correlation_statistics(
                counts_df=next_tc,
                metadata_df=next_sra,
                dist_method=str(getattr(args, 'dist_method', 'pearson')),
                round_value=round_value,
                correlation_statistics=correlation_statistics,
            )
            save_state_overview_pdf(
                counts_df=next_tc,
                metadata_df=next_sra,
                selected_sample_groups=selected_sample_groups,
                out_pdf_path=os.path.join(dir_pdf, '{}.{}.correlation_cutoff.pdf'.format(species_tag, round_value)),
                dist_method=str(getattr(args, 'dist_method', 'pearson')),
                transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
                font_size=8,
                tau_options=tau_options_from_args(args),
            )
        should_stop = _should_stop_within_group_filter(
            current_tc=current_tc,
            next_tc=next_tc,
            excluded_runs=excluded_runs,
            completed_iterations=round_index + 1,
            max_iterations=max_iterations,
        )
        current_tc = next_tc
        current_sra = next_sra
        round_index += 1
        if should_stop:
            if limit_reached and excluded_runs and next_tc.shape[1] > 0:
                round_summary = append_round_summary(
                    round_summary=round_summary, step='within_group_filter_stop',
                    round_value=round_index, reason='iteration_limit_reached',
                    runs_before=current_tc.columns, runs_after=current_tc.columns,
                )
            break
    current_sra['ws_filter_stop_reason'] = (
        'no_outlier_detected' if not excluded_runs else
        'all_samples_excluded' if current_tc.shape[1] == 0 else
        'iteration_limit_reached' if limit_reached else 'no_progress'
    )
    current_sra['ws_filter_iterations'] = round_index
    current_sra['ws_max_filter_iterations'] = max_iterations if max_iterations is not None else 'until_stable'

    batch_info = initialize_batch_info(run_ids=current_sra.loc[:, 'run'].astype(str).tolist(), batch_effect_alg='no')
    batch_info['skip_reason'] = 'batch_effect_alg_no'
    batch_info['batch_effect_alg_applied'] = 'no'
    batch_info['uncorrected_runs'] = list(current_tc.columns)
    batch_info['corrected_runs'] = []
    sra_out = annotate_metadata_with_batch_info(current_sra, batch_info)
    sra_out.to_csv(os.path.join(dir_tsv, '{}.metadata.tsv'.format(species_tag)), sep='\t', index=False)
    _write_prepare_outputs(
        args=args,
        species_tag=species_tag,
        scientific_name=scientific_name,
        dir_tsv=dir_tsv,
        dir_pdf=dir_pdf,
        tc_before_filter=tc_before_filter,
        tc_final=current_tc,
        sra_out=sra_out,
        selected_sample_groups=selected_sample_groups,
        round_summary=round_summary,
        correlation_statistics=correlation_statistics,
        num_total_runs_species=num_total_runs_species,
        num_runs_after_sample_group_filter=num_runs_after_sample_group_filter,
    )
    return 0


def run_per_species_python_worker(args, metadata, species_tag, input_dir):
    requested_mode = str(getattr(args, 'worker_mode', 'prepare_tables'))
    if requested_mode == 'finalize':
        return run_finalize_python_worker(
            args=args,
            metadata=metadata,
            species_tag=species_tag,
            input_dir=input_dir,
        )
    if requested_mode in {'prepare_tables', 'wsfilter'}:
        return _run_prepare_or_wsfilter_python_worker(
            args=args,
            metadata=metadata,
            species_tag=species_tag,
            input_dir=input_dir,
        )
    raise ValueError('Unsupported per-species worker mode: {}'.format(requested_mode))


__all__ = [
    'run_per_species_python_worker',
    'should_use_python_per_species_worker',
]
