import os
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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


def _plot_corr_method(dist_method):
    return str(dist_method) if str(dist_method) in {'pearson', 'spearman', 'kendall'} else 'pearson'


def _compute_corr_matrix(counts_df, dist_method):
    corr_df = counts_df.corr(method=_plot_corr_method(dist_method))
    return corr_df.fillna(0.0)


def _compute_distance_matrix(corr_df):
    dist = 1.0 - corr_df.to_numpy(dtype=float)
    dist = (dist + dist.T) / 2.0
    numpy.fill_diagonal(dist, 0.0)
    return dist


def _compute_pca_coords(corr_df):
    matrix = corr_df.to_numpy(dtype=float)
    n = matrix.shape[0]
    coords = numpy.zeros((n, 2), dtype=float)
    if n <= 1:
        return coords
    eigvals, eigvecs = numpy.linalg.eigh(matrix)
    order = numpy.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    for idx in range(min(2, eigvecs.shape[1])):
        value = max(float(eigvals[idx]), 0.0)
        coords[:, idx] = eigvecs[:, idx] * numpy.sqrt(value)
    return coords


def _compute_mds_coords(corr_df):
    dist = _compute_distance_matrix(corr_df)
    n = dist.shape[0]
    coords = numpy.zeros((n, 2), dtype=float)
    if n <= 1:
        return coords
    centering = numpy.eye(n) - (numpy.ones((n, n), dtype=float) / float(n))
    gram = -0.5 * centering.dot(dist ** 2).dot(centering)
    eigvals, eigvecs = numpy.linalg.eigh(gram)
    order = numpy.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    for idx in range(min(2, eigvecs.shape[1])):
        value = max(float(eigvals[idx]), 0.0)
        coords[:, idx] = eigvecs[:, idx] * numpy.sqrt(value)
    return coords


def _resolve_tsne_perplexity(num_samples):
    if int(num_samples) < 4:
        return None
    max_perplexity = int((int(num_samples) - 1) // 3)
    if max_perplexity < 1:
        return None
    return min(30, max_perplexity)


def _compute_tsne_coords(counts_df):
    num_samples = counts_df.shape[1]
    coords = numpy.full((num_samples, 2), numpy.nan, dtype=float)
    perplexity = _resolve_tsne_perplexity(num_samples)
    if perplexity is None:
        return coords
    try:
        from sklearn.manifold import TSNE
    except ImportError:
        return coords
    try:
        coords = TSNE(
            n_components=2,
            perplexity=float(perplexity),
            random_state=1,
            init='pca',
            learning_rate='auto',
            method='exact',
        ).fit_transform(counts_df.T.to_numpy(dtype=float))
    except ValueError:
        return numpy.full((num_samples, 2), numpy.nan, dtype=float)
    return coords


def _sample_group_color_map(sample_groups):
    groups = [str(value) for value in sample_groups]
    unique_groups = list(dict.fromkeys(groups))
    cmap = plt.get_cmap('tab20')
    return {group: cmap(idx % max(1, cmap.N)) for idx, group in enumerate(unique_groups)}


def _bioproject_color_map(bioprojects):
    values = [str(value) for value in bioprojects]
    unique_values = list(dict.fromkeys(values))
    cmap = plt.get_cmap('tab20b')
    return {value: cmap(idx % max(1, cmap.N)) for idx, value in enumerate(unique_values)}


def _draw_embedding_panel(ax, coords, metadata_df, title, x_label, y_label, font_size=8):
    if coords.shape[0] == 0 or numpy.isfinite(coords).sum() == 0:
        ax.text(0.5, 0.5, 'No finite coordinates', ha='center', va='center', fontsize=font_size)
        ax.set_axis_off()
        return
    group_colors = _sample_group_color_map(metadata_df.loc[:, 'sample_group'].fillna('').astype(str).tolist())
    bp_colors = _bioproject_color_map(metadata_df.loc[:, 'bioproject'].fillna('').astype(str).tolist())
    facecolors = [group_colors[str(value)] for value in metadata_df.loc[:, 'sample_group'].fillna('').astype(str)]
    edgecolors = [bp_colors[str(value)] for value in metadata_df.loc[:, 'bioproject'].fillna('').astype(str)]
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=facecolors,
        edgecolors=edgecolors,
        s=50.0,
        linewidths=0.8,
        alpha=0.9,
    )
    ax.set_title(title, fontsize=font_size)
    ax.set_xlabel(x_label, fontsize=font_size)
    ax.set_ylabel(y_label, fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size)


def _draw_dendrogram_panel(ax, corr_df, metadata_df, font_size=8):
    if corr_df.shape[0] <= 1:
        ax.text(0.5, 0.5, 'No dendrogram data', ha='center', va='center', fontsize=font_size)
        ax.set_axis_off()
        return
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import squareform
    except ImportError:
        ax.text(0.5, 0.5, 'SciPy not available', ha='center', va='center', fontsize=font_size)
        ax.set_axis_off()
        return
    dist = _compute_distance_matrix(corr_df)
    condensed = squareform(dist, checks=False)
    linkage_matrix = linkage(condensed, method='average')
    labels = metadata_df.loc[:, 'run'].astype(str).tolist()
    dendrogram(
        linkage_matrix,
        labels=labels,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=max(4, font_size - 2),
        color_threshold=None,
    )
    ax.set_ylabel('Distance', fontsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)


def _draw_boxplot_panel(ax, corr_df, metadata_df, font_size=8):
    is_same_bp = numpy.equal.outer(
        metadata_df.loc[:, 'bioproject'].astype(str).to_numpy(),
        metadata_df.loc[:, 'bioproject'].astype(str).to_numpy(),
    )
    is_same_group = numpy.equal.outer(
        metadata_df.loc[:, 'sample_group'].astype(str).to_numpy(),
        metadata_df.loc[:, 'sample_group'].astype(str).to_numpy(),
    )
    pair_mask = numpy.triu(numpy.ones(corr_df.shape, dtype=bool), k=1)
    corr = corr_df.to_numpy(dtype=float)
    values = [
        corr[pair_mask & (~is_same_bp) & (~is_same_group)],
        corr[pair_mask & is_same_bp & (~is_same_group)],
        corr[pair_mask & (~is_same_bp) & is_same_group],
        corr[pair_mask & is_same_bp & is_same_group],
    ]
    values = [numpy.asarray(group, dtype=float)[numpy.isfinite(group)] for group in values]
    ax.boxplot(values, patch_artist=True, widths=0.6)
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(['bw\nbw', 'bw\nwi', 'wi\nbw', 'wi\nwi'], fontsize=font_size)
    ax.set_ylabel("Pearson's correlation", fontsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)


def _draw_expression_histogram_panel(ax, counts_df, metadata_df, selected_sample_groups, transform_method, font_size=8):
    tc_sample_group = sample_group_mean(
        counts_df=counts_df,
        metadata_df=metadata_df,
        selected_sample_groups=selected_sample_groups,
        balance_bp=False,
    )['tc_ave']
    xmax = pandas.to_numeric(tc_sample_group.max(axis=1), errors='coerce').clip(lower=0, upper=15).dropna().to_numpy(dtype=float)
    ax.hist(xmax, bins=numpy.arange(0.0, 16.0, 1.0), color='gray', edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Max expression ({})'.format(transform_method), fontsize=font_size)
    ax.set_ylabel('Gene count', fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size)


def _draw_tau_histogram_panel(ax, counts_df, metadata_df, selected_sample_groups, transform_method, font_size=8):
    tc_sample_group = sample_group_mean(
        counts_df=counts_df,
        metadata_df=metadata_df,
        selected_sample_groups=selected_sample_groups,
        balance_bp=False,
    )['tc_ave']
    df_tau = sample_group_to_tau(
        tc_sample_group_df=tc_sample_group,
        rich_annotation=False,
        transform_method=transform_method,
    )
    tau_values = pandas.to_numeric(df_tau.loc[:, 'tau'], errors='coerce').dropna().to_numpy(dtype=float)
    counts_hist, _bins, _patches = ax.hist(
        tau_values,
        bins=numpy.arange(0.0, 1.000001, 0.05),
        color='gray',
        edgecolor='black',
        linewidth=0.5,
    )
    ax.set_xlabel('Tau', fontsize=font_size)
    ax.set_ylabel('Gene count', fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size)
    num_noexp = int(df_tau.loc[:, 'tau'].isna().sum())
    num_all = int(df_tau.shape[0])
    ymax = float(numpy.nanmax(counts_hist)) if counts_hist.size > 0 else 1.0
    ax.text(0.01, ymax * 0.85, 'Excluded due to\nno expression:\n{}/{} genes'.format(num_noexp, num_all), fontsize=font_size, va='top')


def _draw_legend_panel(ax, metadata_df, font_size=8):
    ax.set_axis_off()
    sample_group_colors = _sample_group_color_map(metadata_df.loc[:, 'sample_group'].fillna('').astype(str).tolist())
    bp_colors = _bioproject_color_map(metadata_df.loc[:, 'bioproject'].fillna('').astype(str).tolist())
    handles = [Line2D([], [], linestyle='none', label='Sample group')]
    handles.extend(
        Line2D([0], [0], marker='o', color='w', label=str(group), markerfacecolor=color, markeredgecolor='black', markersize=7)
        for group, color in sample_group_colors.items()
    )
    handles.append(Line2D([], [], linestyle='none', label='BioProject'))
    handles.extend(
        Line2D([0], [0], marker='o', color='white', label=str(bp), markerfacecolor='white', markeredgecolor=color, markersize=7)
        for bp, color in bp_colors.items()
    )
    ax.legend(handles=handles, frameon=False, loc='center', fontsize=font_size, ncol=2)


def save_state_overview_pdf(
    counts_df,
    metadata_df,
    selected_sample_groups,
    out_pdf_path,
    dist_method='pearson',
    transform_method='log2p1-fpkm',
    font_size=8,
):
    out = intersect_counts_and_metadata(counts_df=counts_df, metadata_df=metadata_df)
    counts = out['tc']
    metadata = out['sra'].copy()
    if counts.shape[1] <= 1:
        return None
    run_order = [run_id for run_id in metadata.loc[:, 'run'].astype(str).tolist() if run_id in counts.columns]
    counts = counts.loc[:, run_order].copy()
    metadata = metadata.set_index('run', drop=False).loc[run_order, :].reset_index(drop=True)
    corr_df = _compute_corr_matrix(counts, dist_method=dist_method)
    pca_coords = _compute_pca_coords(corr_df)
    mds_coords = _compute_mds_coords(corr_df)
    tsne_coords = _compute_tsne_coords(counts)

    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    fig, axes = plt.subplots(5, 2, figsize=(12.0, 18.0))
    _draw_dendrogram_panel(axes[0, 0], corr_df, metadata, font_size=font_size)
    axes[0, 0].set_title('Dendrogram', fontsize=font_size)
    image = axes[0, 1].imshow(corr_df.to_numpy(dtype=float), vmin=-1.0, vmax=1.0, cmap='coolwarm')
    axes[0, 1].set_title('Sample correlation', fontsize=font_size)
    axes[0, 1].set_xticks(range(len(run_order)))
    axes[0, 1].set_xticklabels(run_order, rotation=90, fontsize=max(4, font_size - 2))
    axes[0, 1].set_yticks(range(len(run_order)))
    axes[0, 1].set_yticklabels(run_order, fontsize=max(4, font_size - 2))
    fig.colorbar(image, ax=axes[0, 1], fraction=0.046, pad=0.04)
    _draw_embedding_panel(axes[1, 0], pca_coords, metadata, 'PCA', 'PC1', 'PC2', font_size=font_size)
    _draw_embedding_panel(axes[1, 1], mds_coords, metadata, 'MDS', 'Axis 1', 'Axis 2', font_size=font_size)
    _draw_embedding_panel(axes[2, 0], tsne_coords, metadata, 't-SNE', 't-SNE 1', 't-SNE 2', font_size=font_size)
    _draw_boxplot_panel(axes[2, 1], corr_df, metadata, font_size=font_size)
    axes[2, 1].set_title('Correlation boxplot', fontsize=font_size)
    _draw_expression_histogram_panel(
        axes[3, 0],
        counts,
        metadata,
        selected_sample_groups=selected_sample_groups,
        transform_method=transform_method,
        font_size=font_size,
    )
    axes[3, 0].set_title('Expression histogram', fontsize=font_size)
    _draw_tau_histogram_panel(
        axes[3, 1],
        counts,
        metadata,
        selected_sample_groups=selected_sample_groups,
        transform_method=transform_method,
        font_size=font_size,
    )
    axes[3, 1].set_title('Tau histogram', fontsize=font_size)
    _draw_legend_panel(axes[4, 0], metadata, font_size=font_size)
    axes[4, 1].set_axis_off()
    fig.tight_layout()
    fig.savefig(out_pdf_path)
    plt.close(fig)
    return out_pdf_path


__all__ = [
    'CORRELATION_STAT_COLUMNS',
    'initialize_correlation_statistics',
    'intersect_counts_and_metadata',
    'save_correlation_statistics',
    'save_state_overview_pdf',
    'save_tau_histogram_pdf',
    'write_table_with_index_name',
]
