import os
import warnings
import json
import hashlib

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
        run_sample_group = sra.loc[
            is_sample_group & is_run_in_tc & sra_is_non_excluded,
            run_column,
        ].tolist()
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
    method = str(transform_method).split('-')[0]
    if method == 'logn':
        mat.loc[:, :] = numpy.exp(mat.to_numpy(dtype=float))
    elif method == 'log2':
        mat.loc[:, :] = numpy.power(2.0, mat.to_numpy(dtype=float))
    elif method == 'lognp1':
        mat.loc[:, :] = numpy.expm1(mat.to_numpy(dtype=float))
    elif method == 'log2p1':
        mat.loc[:, :] = numpy.power(2.0, mat.to_numpy(dtype=float)) - 1.0
    return mat


def sample_group_to_tau(tc_sample_group_df, rich_annotation=True, transform_method='none'):
    """Compute tau from complete, nonnegative tissue representatives.

    Pipeline callers supply linear arithmetic means, with no further transform.
    An explicit transform_method can decode an encoded representative table;
    it cannot recover linear means from means taken in log space.
    """
    transformed = _inverse_transform_for_tau(tc_sample_group_df=tc_sample_group_df, transform_method=transform_method)
    values = transformed.to_numpy(dtype=float)
    n_groups = transformed.shape[1]
    tau = numpy.full((transformed.shape[0],), numpy.nan, dtype=float)
    complete = numpy.isfinite(values).all(axis=1)
    nonnegative = (values >= 0).all(axis=1)
    xmax = values.max(axis=1) if n_groups else numpy.zeros(values.shape[0])
    valid = complete & nonnegative & (xmax > 0) & (n_groups > 1)
    if valid.any():
        tau[valid] = numpy.sum(1.0 - values[valid] / xmax[valid, None], axis=1) / (n_groups - 1)

    df_tau = pandas.DataFrame(index=transformed.index.copy(), data={'tau': tau})
    if not rich_annotation:
        return df_tau

    highest_values = []
    order_values = []
    highest_ties = []
    column_names = [str(column) for column in transformed.columns]
    for row, is_valid in zip(values, valid):
        positive = row > 0
        if not is_valid:
            highest_values.append(numpy.nan)
            order_values.append(numpy.nan)
            highest_ties.append(numpy.nan)
            continue
        positive_names = [column_names[idx] for idx, is_positive in enumerate(positive) if bool(is_positive)]
        positive_values = row[positive]
        order_idx = numpy.argsort(-positive_values, kind='mergesort')
        ordered_groups = [positive_names[idx] for idx in order_idx]
        highest_values.append(ordered_groups[0])
        order_values.append('|'.join(ordered_groups))
        highest_ties.append('|'.join(sorted(str(column_names[i]) for i in numpy.flatnonzero(row == row.max()))))

    df_tau.loc[:, 'highest'] = highest_values
    df_tau.loc[:, 'order'] = order_values
    df_tau.loc[:, 'highest_ties'] = highest_ties
    status = numpy.full(values.shape[0], 'ok', dtype=object)
    status[xmax == 0] = 'all_zero'
    status[~nonnegative & complete] = 'negative_expression'
    status[~complete] = 'missing_or_nonfinite_expression'
    if n_groups < 2:
        status[:] = 'insufficient_groups'
    df_tau.loc[:, 'tau_status'] = status
    df_tau['num_groups_required'] = n_groups
    df_tau.loc[:, 'num_groups_observed'] = numpy.isfinite(values).sum(axis=1)
    return df_tau


def tau_options_from_args(args):
    return {
        'unit': getattr(args, 'tau_unit', 'run'),
        'balance_projects': bool(getattr(args, 'tau_balance_projects', False)),
    }


def linear_sample_group_summary(
    counts_df, metadata_df, selected_sample_groups=None,
    transform_method='log2p1-fpkm', unit='run', balance_projects=False,
):
    """Inverse-transform each run before averaging; never drop requested groups.

    Units are averaged within each sample_group, optionally within projects
    first. Missing expression propagates rather than changing weights per gene.
    Donor IDs must already be curated, globally unambiguous IDs within a species.
    """
    if unit not in {'run', 'biosample', 'donor'}:
        raise ValueError('Unsupported tau unit: {}'.format(unit))
    if str(transform_method).split('-')[0] not in {'none', 'logn', 'log2', 'lognp1', 'log2p1'}:
        raise ValueError('Unsupported tau input transformation: {}'.format(transform_method))
    required = {'run', 'sample_group', 'exclusion', unit}
    if balance_projects:
        required.add('bioproject')
    missing = required.difference(metadata_df.columns)
    if missing:
        raise ValueError('Missing tau metadata columns: {}'.format(', '.join(sorted(missing))))
    metadata = metadata_df.copy()
    metadata['run'] = metadata['run'].astype(str)
    if metadata['run'].duplicated().any() or counts_df.columns.duplicated().any():
        raise ValueError('Tau aggregation requires unique run IDs.')
    if 'scientific_name' in metadata and metadata['scientific_name'].nunique() > 1:
        raise ValueError('Tau aggregation requires one species at a time.')
    groups = list(dict.fromkeys(
        metadata['sample_group'].dropna().tolist()
        if selected_sample_groups is None else selected_sample_groups
    ))
    retained = metadata.loc[_is_non_excluded_flag(metadata['exclusion'])].copy()
    retained = retained.loc[retained['sample_group'].isin(groups)]
    for column in {unit} | ({'bioproject'} if balance_projects else set()):
        retained[column] = retained[column].fillna('').astype(str).str.strip()
        if retained[column].isin({'', 'not_provided'}).any():
            raise ValueError('Tau {} aggregation requires nonmissing {} IDs.'.format(unit, column))
    if balance_projects and unit != 'run':
        projects_per_unit = retained.groupby(['sample_group', unit])['bioproject'].nunique()
        if projects_per_unit.gt(1).any():
            raise ValueError('Tau project balancing cannot count a shared {} in multiple projects.'.format(unit))
    with numpy.errstate(over='ignore', invalid='ignore'):
        linear = _inverse_transform_for_tau(counts_df, transform_method)
    representatives = pandas.DataFrame(numpy.nan, index=counts_df.index, columns=groups)
    coverage = []
    weights = []
    for group in groups:
        group_metadata = metadata.loc[metadata['sample_group'].eq(group)]
        eligible = retained.loc[retained['sample_group'].eq(group)]
        available = eligible.loc[eligible['run'].isin(linear.columns)]
        status = 'ok'
        if group_metadata.empty:
            status = 'not_collected'
        elif eligible.empty:
            status = 'all_excluded'
        elif len(available) != len(eligible):
            status = 'missing_runs'
        project_ids = available.get('bioproject', pandas.Series(dtype=object)).fillna('').astype(str).str.strip()
        project_ids = project_ids.loc[~project_ids.isin({'', 'not_provided'})]
        coverage.append({
            'sample_group': group, 'status': status,
            'num_runs_retained': len(eligible), 'num_runs_available': len(available),
            'num_units': available[unit].nunique(),
            'num_projects': project_ids.nunique(),
        })
        if status != 'ok':
            continue
        keys = ['bioproject', unit] if balance_projects else [unit]
        unit_means = []
        unit_projects = []
        for _, members in available.groupby(keys, sort=False):
            values = linear.loc[:, members['run'].tolist()]
            mean = values.mean(axis=1, skipna=False)
            valid = numpy.isfinite(values).all(axis=1) & values.ge(0).all(axis=1)
            unit_means.append(mean.where(valid))
            if balance_projects:
                unit_projects.append(members['bioproject'].iloc[0])
            denominator = available[unit].nunique()
            if balance_projects:
                project = members['bioproject'].iloc[0]
                denominator = available['bioproject'].nunique() * available.loc[
                    available['bioproject'].eq(project), unit,
                ].nunique()
            for run_id in members['run']:
                weights.append({
                    'sample_group': group, 'run': run_id,
                    'unit_id': members[unit].iloc[0],
                    'weight': 1.0 / (denominator * len(members)),
                })
        if not unit_means:
            continue
        unit_table = pandas.concat(unit_means, axis=1, ignore_index=True)
        if balance_projects:
            project_means = [
                unit_table.loc[:, [i for i, value in enumerate(unit_projects) if value == project]].mean(axis=1, skipna=False)
                for project in dict.fromkeys(unit_projects)
            ]
            unit_table = pandas.concat(project_means, axis=1, ignore_index=True)
        representatives[group] = unit_table.mean(axis=1, skipna=False)
    return {
        'linear_mean': representatives,
        'coverage': pandas.DataFrame(coverage, columns=[
            'sample_group', 'status', 'num_runs_retained', 'num_runs_available', 'num_units', 'num_projects',
        ]),
        'weights': pandas.DataFrame(weights, columns=['sample_group', 'run', 'unit_id', 'weight']),
    }


def write_tau_outputs(counts_df, metadata_df, selected_sample_groups, transform_method,
                      dir_tsv, species_tag, batch_effect_alg, **options):
    summary = linear_sample_group_summary(
        counts_df, metadata_df, selected_sample_groups, transform_method, **options,
    )
    prefix = os.path.join(dir_tsv, '{}.{}'.format(species_tag, batch_effect_alg))
    summary['linear_mean'].to_csv(prefix + '.tau.linear_mean.tsv', sep='\t', index_label='target_id')
    summary['coverage'].to_csv(prefix + '.tau.coverage.tsv', sep='\t', index=False)
    summary['weights'].to_csv(prefix + '.tau.weights.tsv', sep='\t', index=False)
    panel_id = hashlib.sha256(json.dumps(sorted(summary['linear_mean'].columns)).encode('utf-8')).hexdigest()
    tau = sample_group_to_tau(summary['linear_mean'])
    tau['panel_id'] = panel_id
    tau.to_csv(prefix + '.tau.tsv', sep='\t', index_label='target_id')
    with open(prefix + '.tau.definition.json', 'w', encoding='utf-8') as handle:
        json.dump({
            'schema_version': 1, 'representative': 'linear_arithmetic_mean',
            'input_transform': transform_method, 'tau_scale': 'linear',
            'unit': options.get('unit', 'run'),
            'balance_projects': options.get('balance_projects', False),
            'sample_groups': list(summary['linear_mean'].columns),
            'panel_id': panel_id,
            'missing_policy': 'require_complete',
            'source_stage': str(batch_effect_alg),
        }, handle, indent=2)
        handle.write('\n')
    return summary['linear_mean']


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
    species_tag=None,
):
    if species_tag is None:
        species_tag = str(scientific_name).replace(' ', '_')
    else:
        species_tag = str(species_tag)
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
    'linear_sample_group_summary',
    'tau_options_from_args',
    'write_tau_outputs',
    'sample_group_mean',
    'sample_group_to_tau',
    'write_curation_summaries',
]
