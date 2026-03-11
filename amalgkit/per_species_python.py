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
    append_round_summary,
    initialize_round_summary,
    sample_group_mean,
    sample_group_to_tau,
    write_curation_summaries,
)
from amalgkit.per_species_finalize_python import (
    _apply_transformation_logic,
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
    excluded_runs = sra.loc[~is_mapping_good, 'run'].astype(str).tolist()
    out_sra = sra.copy()
    if len(excluded_runs) > 0:
        out_sra.loc[~is_mapping_good, 'exclusion'] = 'low_mapping_rate'
    keep_runs = [run_id for run_id in out_sra.loc[is_mapping_good, 'run'].astype(str).tolist() if run_id in tc.columns]
    return tc.loc[:, keep_runs].copy(), out_sra, excluded_runs


def _compute_sample_group_correlation_metrics(tc, sra, selected_sample_groups, dist_method):
    required_sample_groups = [group for group in selected_sample_groups if group in set(sra.loc[:, 'sample_group'].astype(str))]
    out = sra.copy()
    for column in ['ws_within_group_cor', 'ws_max_nongroup_cor', 'ws_margin', 'ws_robust_z', 'ws_outlier_candidate']:
        if column not in out.columns:
            out.loc[:, column] = numpy.nan if column != 'ws_outlier_candidate' else False
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
    for run_id in corr_by_group.index:
        sample_values = pandas.to_numeric(tc.loc[:, run_id], errors='coerce')
        for sample_group in corr_by_group.columns:
            try:
                corr_value = sample_values.corr(
                    pandas.to_numeric(tc_ave.loc[:, sample_group], errors='coerce'),
                    method=str(dist_method),
                )
            except ValueError:
                corr_value = numpy.nan
            corr_by_group.loc[run_id, sample_group] = corr_value
    run_values = out.loc[:, 'run'].astype(str).tolist()
    sample_group_values = out.loc[:, 'sample_group'].astype(str).tolist()
    within_values = []
    nongroup_values = []
    margin_values = []
    for run_id, sample_group in zip(run_values, sample_group_values):
        if (run_id not in corr_by_group.index) or (sample_group not in corr_by_group.columns):
            within_values.append(numpy.nan)
            nongroup_values.append(numpy.nan)
            margin_values.append(numpy.nan)
            continue
        corr_row = pandas.to_numeric(corr_by_group.loc[run_id, :], errors='coerce')
        within_cor = corr_row.get(sample_group, numpy.nan)
        nongroup = pandas.to_numeric(corr_row.loc[corr_row.index != sample_group], errors='coerce').dropna()
        max_nongroup = float(nongroup.max()) if nongroup.shape[0] > 0 else numpy.nan
        margin_val = within_cor - max_nongroup if numpy.isfinite(within_cor) and numpy.isfinite(max_nongroup) else numpy.nan
        within_values.append(within_cor)
        nongroup_values.append(max_nongroup)
        margin_values.append(margin_val)
    out.loc[:, 'ws_within_group_cor'] = within_values
    out.loc[:, 'ws_max_nongroup_cor'] = nongroup_values
    out.loc[:, 'ws_margin'] = margin_values
    return out


def _reduce_outlier_candidates(candidate_df):
    if candidate_df.shape[0] == 0:
        return []
    keep_runs = []
    if 'bioproject' in candidate_df.columns:
        first_bp = candidate_df.drop_duplicates(subset=['bioproject'], keep='first')
        keep_runs.extend(first_bp.loc[:, 'run'].astype(str).tolist())
    if 'sample_group' in candidate_df.columns:
        first_group = candidate_df.drop_duplicates(subset=['sample_group'], keep='first')
        keep_runs.extend(first_group.loc[:, 'run'].astype(str).tolist())
    keep_runs = [run_id for run_id in keep_runs if run_id != '']
    return list(dict.fromkeys(keep_runs))


def _apply_within_group_filter(tc, sra, args, selected_sample_groups, min_dif=0.0):
    out = _compute_sample_group_correlation_metrics(
        tc=tc,
        sra=sra,
        selected_sample_groups=selected_sample_groups,
        dist_method=str(getattr(args, 'dist_method', 'pearson')),
    )
    filtered = flag_margin_outliers(
        df=out,
        margin_col='ws_margin',
        group_col='sample_group',
        margin_threshold=float(getattr(args, 'margin_threshold', 0.0)) + float(min_dif),
        robust_z_threshold=float(getattr(args, 'robust_z_threshold', -2.5)),
        robust_z_col='ws_robust_z',
        outlier_col='ws_outlier_candidate',
    )
    candidate_runs = (
        filtered.loc[filtered['ws_outlier_candidate'].fillna(False).astype(bool), 'run']
        .fillna('')
        .astype(str)
        .tolist()
    )
    candidate_runs = [run_id for run_id in candidate_runs if run_id != '']
    excluded_runs = list(dict.fromkeys(candidate_runs))
    if bool(getattr(args, 'one_outlier_per_iter', False)) and len(excluded_runs) > 0:
        candidate_df = filtered.loc[
            filtered['run'].astype(str).isin(excluded_runs),
            [col for col in ['run', 'sample_group', 'bioproject'] if col in filtered.columns],
        ].copy()
        excluded_runs = _reduce_outlier_candidates(candidate_df)
    out_sra = sra.copy()
    metric_cols = ['ws_within_group_cor', 'ws_max_nongroup_cor', 'ws_margin', 'ws_robust_z', 'ws_outlier_candidate']
    for metric_col in metric_cols:
        if metric_col not in out_sra.columns:
            out_sra.loc[:, metric_col] = numpy.nan if metric_col != 'ws_outlier_candidate' else False
        run_map = filtered.set_index(filtered['run'].astype(str))[metric_col]
        out_sra.loc[:, metric_col] = out_sra.loc[:, 'run'].astype(str).map(run_map)
    if len(excluded_runs) > 0:
        out_sra.loc[out_sra['run'].astype(str).isin(excluded_runs), 'exclusion'] = 'low_within_sample_group_correlation'
    out_tc = tc.loc[:, [run_id for run_id in tc.columns if run_id not in set(excluded_runs)]].copy()
    return out_tc, out_sra, excluded_runs


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
    tau_df = sample_group_to_tau(
        tc_sample_group_df=tc_sample_group_final,
        transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
    )
    write_table_with_index_name(
        df=tau_df,
        file_path=os.path.join(dir_tsv, '{}.{}.tau.tsv'.format(species_tag, batch_effect_alg)),
        index_name='target_id',
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
    )
    save_tau_histogram_pdf(
        counts_df=tc_final,
        metadata_df=sra_out,
        selected_sample_groups=selected_sample_groups,
        out_pdf_path=os.path.join(dir_pdf, '{}.tau_histogram.no.pdf'.format(species_tag)),
        transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
        tc_sample_group_df=tc_sample_group_final,
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
    if not os.path.isfile(count_path) or not os.path.isfile(eff_length_path):
        return 1

    counts_df = _normalize_dataframe_columns(_read_expression_tsv(count_path))
    eff_length_df = _normalize_dataframe_columns(_read_expression_tsv(eff_length_path))
    metadata_all = _standardize_metadata_all(_normalize_metadata_df(metadata.df))
    scientific_name = _resolve_scientific_name(metadata_all, species_tag)
    selected_sample_groups = _resolve_selected_sample_groups(args, metadata_all)
    num_total_runs_species = int(metadata_all.loc[:, 'scientific_name'].astype(str).eq(scientific_name).sum())
    sra = _get_species_metadata(metadata_all, scientific_name, selected_sample_groups, counts_df.columns)
    num_runs_after_sample_group_filter = int(sra.shape[0])

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
    tc_original = _apply_transformation_logic(tc, eff_length_species, args.norm, 'no', 'before_batch', sra)
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
    )

    tc, sra, mapping_excluded_runs = _filter_low_mapping_rate(
        tc=tc,
        sra=sra,
        mapping_rate_cutoff=float(getattr(args, 'mapping_rate', 0.0)),
    )
    tc = _apply_transformation_logic(tc, eff_length_species, args.norm, 'no', 'before_batch', sra)
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
        round_value = round_index + 2
        if (len(excluded_runs) == 0) or bool(getattr(args, 'plot_intermediate', False)):
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
            )
        current_tc = next_tc
        current_sra = next_sra
        round_index += 1
        if len(excluded_runs) == 0:
            break

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
