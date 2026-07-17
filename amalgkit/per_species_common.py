import os
import warnings

import numpy
import pandas


def _is_non_excluded_flag(exclusion_values):
    normalized = (
        pandas.Series(exclusion_values)
        .fillna('')
        .astype(str)
        .str.strip()
        .str.lower()
    )
    return normalized.eq('no').to_numpy(dtype=bool)


def sample_group_mean(
    counts_df,
    metadata_df,
    selected_sample_groups=None,
    balance_bp=False,
    sample_group_column='sample_group',
    run_column='run',
    exclusion_column='exclusion',
    batch_column='bioproject',
):
    if counts_df.shape[1] == 0:
        return {
            'tc_ave': pandas.DataFrame(index=counts_df.index),
            'selected_sample_groups': [] if selected_sample_groups is None else list(selected_sample_groups),
        }
    if sample_group_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(sample_group_column))
    if run_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(run_column))
    if exclusion_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(exclusion_column))
    if balance_bp and (batch_column not in metadata_df.columns):
        raise ValueError('Missing required metadata column for balance_bp: {}'.format(batch_column))

    sra = metadata_df.copy()
    sra.loc[:, run_column] = sra.loc[:, run_column].astype(str)
    sample_group_values = sra.loc[:, sample_group_column]
    if selected_sample_groups is None or all(pandas.isna(list(selected_sample_groups))):
        resolved_groups = (
            pandas.Series(sample_group_values)
            .drop_duplicates()
            .tolist()
        )
    else:
        resolved_groups = [
            group for group in list(selected_sample_groups)
            if group in set(sample_group_values.tolist())
        ]

    tc_ave = pandas.DataFrame(
        numpy.nan,
        index=counts_df.index,
        columns=list(resolved_groups),
        dtype=float,
    )
    sra_is_non_excluded = _is_non_excluded_flag(sra.loc[:, exclusion_column])
    is_run_in_tc = sra.loc[:, run_column].isin(list(counts_df.columns)).to_numpy(dtype=bool)
    if balance_bp:
        sra_bioproject = sra.loc[:, batch_column].astype(str).to_numpy(dtype=object)

    kept_groups = list(resolved_groups)
    for sample_group in list(resolved_groups):
        is_sample_group = sra.loc[:, sample_group_column].eq(sample_group).to_numpy(dtype=bool)
        exclusion_sample_group = sra_is_non_excluded[is_sample_group]
        run_sample_group = sra.loc[is_sample_group & is_run_in_tc, run_column].tolist()
        if (len(exclusion_sample_group) > 0) and bool(numpy.all(~exclusion_sample_group)):
            kept_groups = [group for group in kept_groups if group != sample_group]
            if sample_group in tc_ave.columns:
                tc_ave = tc_ave.drop(columns=[sample_group])
            warnings.warn(
                'All samples of sample_group {} are marked for exclusion. This sample_group will be omitted from further analysis.'.format(
                    sample_group
                )
            )
            continue
        if len(run_sample_group) == 0:
            continue
        if len(run_sample_group) == 1:
            exp_sample_group = counts_df.loc[:, run_sample_group[0]].to_numpy(dtype=float)
        else:
            if balance_bp:
                is_no_exclusion = sra_is_non_excluded
                bps = pandas.Series(
                    sra_bioproject[is_run_in_tc & is_sample_group & is_no_exclusion],
                    dtype=object,
                ).drop_duplicates().tolist()
                df_tmp = pandas.DataFrame(
                    numpy.nan,
                    index=counts_df.index,
                    columns=bps,
                    dtype=float,
                )
                for bp in bps:
                    sra_ids = sra.loc[
                        sra.loc[:, batch_column].eq(bp) &
                        sra.loc[:, sample_group_column].eq(sample_group) &
                        pandas.Series(sra_is_non_excluded, index=sra.index),
                        run_column,
                    ].tolist()
                    tc_bp = counts_df.loc[:, [run_id for run_id in sra_ids if run_id in counts_df.columns]]
                    if tc_bp.shape[1] == 0:
                        continue
                    if tc_bp.shape[1] == 1:
                        df_tmp.loc[:, bp] = tc_bp.iloc[:, 0].to_numpy(dtype=float)
                    else:
                        df_tmp.loc[:, bp] = tc_bp.mean(axis=1, skipna=True).to_numpy(dtype=float)
                exp_sample_group = df_tmp.mean(axis=1, skipna=True).to_numpy(dtype=float)
            else:
                exp_sample_group = counts_df.loc[:, run_sample_group].mean(axis=1, skipna=True).to_numpy(dtype=float)
        tc_ave.loc[:, sample_group] = exp_sample_group

    return {
        'tc_ave': tc_ave,
        'selected_sample_groups': kept_groups,
    }


def _inverse_transform_for_tau(tc_sample_group_df, transform_method):
    mat = tc_sample_group_df.astype(float).copy()
    method = str(transform_method)
    if 'logn-' in method:
        mat.loc[:, :] = numpy.exp(mat.to_numpy(dtype=float))
    elif 'log2-' in method:
        mat.loc[:, :] = numpy.power(2.0, mat.to_numpy(dtype=float))
    elif 'lognp1-' in method:
        mat.loc[:, :] = numpy.exp(mat.to_numpy(dtype=float)) - 1.0
    elif 'log2p1-' in method:
        mat.loc[:, :] = numpy.power(2.0, mat.to_numpy(dtype=float)) - 1.0
    mat.loc[:, :] = numpy.where(mat.to_numpy(dtype=float) < 0, 0.0, mat.to_numpy(dtype=float))
    return mat


def sample_group_to_tau(tc_sample_group_df, rich_annotation=True, transform_method='log2p1-fpkm'):
    columns = ['tau', 'highest', 'order'] if rich_annotation else ['tau']
    if (tc_sample_group_df.shape[0] == 0) or (tc_sample_group_df.shape[1] == 0):
        return pandas.DataFrame(index=tc_sample_group_df.index.copy(), columns=columns)

    transformed = _inverse_transform_for_tau(tc_sample_group_df=tc_sample_group_df, transform_method=transform_method)
    values = transformed.to_numpy(dtype=float)
    n_groups = transformed.shape[1]
    tau = numpy.full((transformed.shape[0],), numpy.nan, dtype=float)
    if n_groups > 1:
        xmax = numpy.nanmax(values, axis=1)
        with numpy.errstate(divide='ignore', invalid='ignore'):
            tau_terms = (1.0 - (values / xmax.reshape(-1, 1))) / float(n_groups - 1)
        invalid_rows = (~numpy.isfinite(tau_terms)).all(axis=1)
        tau = numpy.nansum(tau_terms, axis=1)
        tau[invalid_rows] = numpy.nan

    df_tau = pandas.DataFrame(index=transformed.index.copy(), data={'tau': tau})
    if not rich_annotation:
        return df_tau

    highest_values = []
    order_values = []
    filled = transformed.fillna(0.0)
    column_names = list(filled.columns)
    for row in filled.to_numpy(dtype=float):
        positive = row > 0
        if int(positive.sum()) == 0:
            highest_values.append(numpy.nan)
            order_values.append(numpy.nan)
            continue
        positive_names = [column_names[idx] for idx, is_positive in enumerate(positive) if bool(is_positive)]
        positive_values = row[positive]
        order_idx = numpy.argsort(-positive_values, kind='mergesort')
        ordered_groups = [positive_names[idx] for idx in order_idx]
        highest_values.append(ordered_groups[0])
        order_values.append('|'.join(ordered_groups))

    df_tau.loc[:, 'highest'] = highest_values
    df_tau.loc[:, 'order'] = order_values
    return df_tau


def initialize_round_summary():
    return pandas.DataFrame(
        {
            'step': pandas.Series(dtype=object),
            'round': pandas.Series(dtype='int64'),
            'reason': pandas.Series(dtype=object),
            'num_runs_before': pandas.Series(dtype='int64'),
            'num_runs_after': pandas.Series(dtype='int64'),
            'num_runs_removed': pandas.Series(dtype='int64'),
            'removed_runs': pandas.Series(dtype=object),
        }
    )


def append_round_summary(round_summary, step, round_value, reason, runs_before, runs_after):
    if round_summary is None:
        round_summary = initialize_round_summary()
    runs_before_list = [str(run_id) for run_id in list(runs_before) if str(run_id) != '']
    runs_after_set = {str(run_id) for run_id in list(runs_after) if str(run_id) != ''}
    removed_runs = []
    seen_removed = set()
    for run_id in runs_before_list:
        if (run_id in runs_after_set) or (run_id in seen_removed):
            continue
        seen_removed.add(run_id)
        removed_runs.append(run_id)
    new_row = pandas.DataFrame(
        [
            {
                'step': str(step),
                'round': int(round_value),
                'reason': str(reason),
                'num_runs_before': int(len(runs_before_list)),
                'num_runs_after': int(len([str(run_id) for run_id in list(runs_after) if str(run_id) != ''])),
                'num_runs_removed': int(len(removed_runs)),
                'removed_runs': '' if len(removed_runs) == 0 else ' '.join(removed_runs),
            }
        ]
    )
    return pandas.concat([round_summary, new_row], ignore_index=True)


def build_curation_final_summary(
    metadata_df,
    scientific_name,
    batch_effect_alg,
    mapping_rate_cutoff,
    correlation_threshold,
    one_outlier_per_iteration,
    num_total_runs_species,
    num_runs_after_sample_group_filter,
    total_runtime_sec,
    run_column='run',
    exclusion_column='exclusion',
):
    if run_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(run_column))
    if exclusion_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(exclusion_column))
    is_kept = _is_non_excluded_flag(metadata_df.loc[:, exclusion_column])
    run_values = metadata_df.loc[:, run_column].astype(str)
    kept_runs = run_values.loc[is_kept].tolist()
    excluded_runs = run_values.loc[~is_kept].tolist()
    return pandas.DataFrame(
        [
            {
                'scientific_name': str(scientific_name),
                'batch_effect_alg': str(batch_effect_alg),
                'mapping_rate_cutoff': float(mapping_rate_cutoff),
                'correlation_threshold': float(correlation_threshold),
                'one_outlier_per_iteration': bool(one_outlier_per_iteration),
                'total_runtime_sec': round(float(total_runtime_sec), 6),
                'num_total_runs_in_species': int(num_total_runs_species),
                'num_runs_after_sample_group_filter': int(num_runs_after_sample_group_filter),
                'num_runs_final_kept': int(is_kept.sum()),
                'num_runs_final_excluded': int((~is_kept).sum()),
                'final_kept_runs': '' if len(kept_runs) == 0 else ' '.join(kept_runs),
                'final_excluded_runs': '' if len(excluded_runs) == 0 else ' '.join(excluded_runs),
            }
        ]
    )


def write_curation_summaries(
    round_summary,
    metadata_df,
    scientific_name,
    batch_effect_alg,
    dir_tsv,
    mapping_rate_cutoff,
    correlation_threshold,
    one_outlier_per_iteration,
    num_total_runs_species,
    num_runs_after_sample_group_filter,
    total_runtime_sec,
):
    species_tag = str(scientific_name).replace(' ', '_')
    os.makedirs(dir_tsv, exist_ok=True)
    round_path = os.path.join(
        dir_tsv,
        '{}.{}.curation_round_summary.tsv'.format(species_tag, batch_effect_alg),
    )
    final_path = os.path.join(
        dir_tsv,
        '{}.{}.curation_final_summary.tsv'.format(species_tag, batch_effect_alg),
    )
    round_summary.to_csv(round_path, sep='\t', index=False)
    final_summary = build_curation_final_summary(
        metadata_df=metadata_df,
        scientific_name=scientific_name,
        batch_effect_alg=batch_effect_alg,
        mapping_rate_cutoff=mapping_rate_cutoff,
        correlation_threshold=correlation_threshold,
        one_outlier_per_iteration=one_outlier_per_iteration,
        num_total_runs_species=num_total_runs_species,
        num_runs_after_sample_group_filter=num_runs_after_sample_group_filter,
        total_runtime_sec=total_runtime_sec,
    )
    final_summary.to_csv(final_path, sep='\t', index=False)
    return {
        'round_path': round_path,
        'final_path': final_path,
        'final_summary': final_summary,
    }


__all__ = [
    'append_round_summary',
    'build_curation_final_summary',
    'initialize_round_summary',
    'sample_group_mean',
    'sample_group_to_tau',
    'write_curation_summaries',
]
