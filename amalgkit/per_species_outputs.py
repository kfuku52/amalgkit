import os
import warnings

import numpy
import pandas

from amalgkit.per_species_common import sample_group_mean, sample_group_to_tau


CORRELATION_STAT_COLUMNS = [
    'bwbw_n', 'bwbw_mean', 'bwbw_median', 'bwbw_variance',
    'wibw_n', 'wibw_mean', 'wibw_median', 'wibw_variance',
    'bwwi_n', 'bwwi_mean', 'bwwi_median', 'bwwi_variance',
    'wiwi_n', 'wiwi_mean', 'wiwi_median', 'wiwi_variance',
]


def write_table_with_index_name(df, file_path, index_name='target_id', sort=True):
    out_df = df.copy()
    index_values = pandas.Series(out_df.index, name=index_name)
    out_df = out_df.reset_index(drop=True)
    out_df.insert(0, index_name, index_values.astype(str).to_numpy())
    if sort and (index_name in out_df.columns):
        out_df = out_df.sort_values(by=index_name, kind='mergesort').reset_index(drop=True)
    os.makedirs(os.path.dirname(os.path.realpath(file_path)), exist_ok=True)
    out_df.to_csv(file_path, sep='\t', index=False)
    return out_df


def intersect_counts_and_metadata(counts_df, metadata_df, run_column='run'):
    if run_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(run_column))
    metadata_runs = metadata_df.loc[:, run_column].astype(str)
    counts = counts_df.loc[:, counts_df.columns.astype(str).isin(metadata_runs)].copy()
    metadata = metadata_df.loc[metadata_runs.isin(counts.columns.astype(str)), :].copy()
    return {
        'tc': counts,
        'sra': metadata,
    }


def initialize_correlation_statistics():
    return pandas.DataFrame(columns=CORRELATION_STAT_COLUMNS, dtype=float)


def _summarize_correlation_group(values):
    values = numpy.asarray(values, dtype=float)
    values = values[numpy.isfinite(values)]
    if values.size == 0:
        return {
            'n': 0.0,
            'mean': numpy.nan,
            'median': numpy.nan,
            'variance': numpy.nan,
        }
    variance = numpy.nan if values.size <= 1 else float(numpy.var(values, ddof=1))
    return {
        'n': float(values.size),
        'mean': float(numpy.mean(values)),
        'median': float(numpy.median(values)),
        'variance': variance,
    }


def _build_correlation_row(tc_dist_matrix, metadata_df, batch_column, sample_group_column):
    is_same_bp = numpy.equal.outer(
        metadata_df.loc[:, batch_column].astype(str).to_numpy(),
        metadata_df.loc[:, batch_column].astype(str).to_numpy(),
    )
    is_same_sample_group = numpy.equal.outer(
        metadata_df.loc[:, sample_group_column].astype(str).to_numpy(),
        metadata_df.loc[:, sample_group_column].astype(str).to_numpy(),
    )
    pair_mask = numpy.triu(numpy.ones(tc_dist_matrix.shape, dtype=bool), k=1)
    stats = {}
    group_masks = {
        'bwbw': (~is_same_bp) & (~is_same_sample_group),
        'wibw': is_same_bp & (~is_same_sample_group),
        'bwwi': (~is_same_bp) & is_same_sample_group,
        'wiwi': is_same_bp & is_same_sample_group,
    }
    for label, mask in group_masks.items():
        summary = _summarize_correlation_group(tc_dist_matrix[pair_mask & mask])
        stats['{}_n'.format(label)] = summary['n']
        stats['{}_mean'.format(label)] = summary['mean']
        stats['{}_median'.format(label)] = summary['median']
        stats['{}_variance'.format(label)] = summary['variance']
    return stats


def save_correlation_statistics(
    counts_df,
    metadata_df,
    dist_method,
    round_value,
    correlation_statistics=None,
    precomputed_tc_dist_matrix=None,
    run_column='run',
    batch_column='bioproject',
    sample_group_column='sample_group',
):
    if batch_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(batch_column))
    if sample_group_column not in metadata_df.columns:
        raise ValueError('Missing required metadata column: {}'.format(sample_group_column))
    out = intersect_counts_and_metadata(counts_df=counts_df, metadata_df=metadata_df, run_column=run_column)
    counts = out['tc']
    metadata = out['sra']
    if precomputed_tc_dist_matrix is not None:
        if isinstance(precomputed_tc_dist_matrix, pandas.DataFrame):
            run_order = [str(run_id) for run_id in precomputed_tc_dist_matrix.columns if str(run_id) in set(metadata.loc[:, run_column].astype(str))]
            tc_dist_matrix = precomputed_tc_dist_matrix.loc[run_order, run_order].to_numpy(dtype=float)
        else:
            tc_dist_matrix = numpy.asarray(precomputed_tc_dist_matrix, dtype=float)
            run_order = list(metadata.loc[:, run_column].astype(str))
        metadata = metadata.set_index(run_column, drop=False).loc[run_order, :].reset_index(drop=True)
    else:
        tc_dist_matrix = None

    if ((tc_dist_matrix is None) and (counts.shape[1] <= 1)) or ((tc_dist_matrix is not None) and (tc_dist_matrix.shape[1] <= 1)):
        row = {column: numpy.nan for column in CORRELATION_STAT_COLUMNS}
        for index in range(0, len(CORRELATION_STAT_COLUMNS), 4):
            row[CORRELATION_STAT_COLUMNS[index]] = 0.0
    else:
        if tc_dist_matrix is None:
            if str(dist_method) not in {'pearson', 'spearman', 'kendall'}:
                raise ValueError('Unsupported correlation method: {}'.format(dist_method))
            tc_dist_matrix = counts.corr(method=str(dist_method)).to_numpy(dtype=float)
        row = _build_correlation_row(
            tc_dist_matrix=tc_dist_matrix,
            metadata_df=metadata,
            batch_column=batch_column,
            sample_group_column=sample_group_column,
        )

    row_df = pandas.DataFrame([row], columns=CORRELATION_STAT_COLUMNS, index=['round_{}'.format(round_value)])
    if correlation_statistics is None:
        correlation_statistics = initialize_correlation_statistics()
    return pandas.concat([correlation_statistics, row_df], axis=0)


def save_tau_histogram_pdf(
    counts_df,
    metadata_df,
    selected_sample_groups,
    out_pdf_path,
    font_size=8,
    transform_method='log2p1-fpkm',
    tc_sample_group_df=None,
):
    if tc_sample_group_df is None:
        out = sample_group_mean(
            counts_df=counts_df,
            metadata_df=metadata_df,
            selected_sample_groups=selected_sample_groups,
            balance_bp=False,
        )
        tc_sample_group_df = out['tc_ave']
    df_tau = sample_group_to_tau(
        tc_sample_group_df=tc_sample_group_df,
        rich_annotation=False,
        transform_method=transform_method,
    )
    tau_values = pandas.to_numeric(df_tau.loc[:, 'tau'], errors='coerce')
    valid_values = tau_values.dropna().to_numpy(dtype=float)
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)

    try:
        import matplotlib
        matplotlib.use('Agg', force=True)
        from matplotlib import pyplot
    except ImportError as exc:
        warnings.warn('matplotlib is required to generate tau histogram {}: {}'.format(out_pdf_path, exc))
        return None

    breaks = numpy.arange(0.0, 1.000001, 0.05)
    fig, ax = pyplot.subplots(figsize=(6.0, 4.0))
    counts_hist, _bins, _patches = ax.hist(
        valid_values,
        bins=breaks,
        color='gray',
        edgecolor='black',
        linewidth=0.5,
    )
    ax.set_xlabel('Tau (expression specificity)', fontsize=font_size)
    ax.set_ylabel('Gene count', fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size)
    num_noexp = int(tau_values.isna().sum())
    num_all = int(df_tau.shape[0])
    max_count = float(numpy.nanmax(counts_hist)) if counts_hist.size > 0 else 0.0
    y_pos = max_count * 0.85 if max_count > 0 else 1.0
    ax.text(
        0.01,
        y_pos,
        'Excluded due to\nno expression:\n{}/{} genes'.format(num_noexp, num_all),
        ha='left',
        va='top',
        fontsize=font_size,
    )
    fig.tight_layout()
    fig.savefig(out_pdf_path)
    pyplot.close(fig)
    return {
        'tau_df': df_tau,
        'num_no_expression': num_noexp,
        'num_total_genes': num_all,
    }


__all__ = [
    'CORRELATION_STAT_COLUMNS',
    'initialize_correlation_statistics',
    'intersect_counts_and_metadata',
    'save_correlation_statistics',
    'save_tau_histogram_pdf',
    'write_table_with_index_name',
]
