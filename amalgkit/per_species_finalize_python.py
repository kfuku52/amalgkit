import os
import warnings

import numpy
import pandas

from amalgkit.batch_effect_common import (
    annotate_metadata_with_batch_info,
    initialize_batch_info,
    write_batch_effect_summary_tsv,
)
from amalgkit.batch_effect_combatseq import run_combatseq_backend
from amalgkit.batch_effect_ruvseq import compute_factor_r2
from amalgkit.batch_effect_latent_glm import run_latent_glm_backend
from amalgkit.batch_effect_ruvseq import run_ruvseq_backend
from amalgkit.batch_effect_sva import run_sva_backend
from amalgkit.filter_utils import _format_genus_species_label
from amalgkit.per_species_common import (
    append_round_summary,
    initialize_round_summary,
    sample_group_mean,
    sample_group_to_tau,
    write_curation_summaries,
)
from amalgkit.per_species_outputs import (
    initialize_correlation_statistics,
    save_tau_histogram_pdf,
    write_table_with_index_name,
)


def should_use_python_finalize_worker(args):
    requested_mode = str(getattr(args, 'worker_mode', 'prepare_tables'))
    if requested_mode != 'finalize':
        return False
    if bool(getattr(args, 'skip_curation', False)):
        return True
    if not bool(getattr(args, 'disable_auto_outlier_filter', False)):
        return False
    batch_effect_alg = str(getattr(args, 'batch_effect_alg', 'no')).lower()
    if batch_effect_alg == 'no':
        return True
    if batch_effect_alg == 'sva':
        return str(getattr(args, 'sva_backend', 'python')).lower() == 'python'
    if batch_effect_alg == 'combatseq':
        return str(getattr(args, 'combatseq_backend', 'python')).lower() == 'python'
    if batch_effect_alg == 'ruvseq':
        return str(getattr(args, 'ruvseq_backend', 'python')).lower() == 'python'
    if batch_effect_alg == 'latent_glm':
        return True
    return False


def _read_expression_tsv(path):
    return pandas.read_csv(path, sep='\t', index_col=0, low_memory=False)


def _normalize_dataframe_columns(df):
    out = df.copy()
    out.index = out.index.map(str)
    out.columns = out.columns.map(str)
    return out


def _normalize_metadata_df(metadata_df):
    out = metadata_df.copy()
    for column in out.columns:
        if out[column].dtype == object:
            out.loc[:, column] = out.loc[:, column].where(~out.loc[:, column].isna(), '')
    if 'run' in out.columns:
        out.loc[:, 'run'] = out.loc[:, 'run'].fillna('').astype(str)
    return out


def _resolve_scientific_name(metadata_df, species_tag):
    scientific_name_series = metadata_df.get('scientific_name', pandas.Series(dtype=object)).fillna('').astype(str)
    normalized = scientific_name_series.str.replace(' ', '_', regex=False)
    matched = scientific_name_series.loc[normalized == species_tag]
    if matched.shape[0] > 0:
        return str(matched.iloc[0])
    return str(species_tag).replace('_', ' ')


def _resolve_selected_sample_groups(args, metadata_df):
    sample_group_arg = getattr(args, 'sample_group', None)
    if sample_group_arg is None:
        if 'sample_group' not in metadata_df.columns:
            raise ValueError('The "sample_group" column was not found in metadata.')
        groups = metadata_df.loc[:, 'sample_group'].fillna('').astype(str).str.strip().tolist()
    else:
        groups = str(sample_group_arg).replace(',', '|').split('|')
    resolved = []
    seen = set()
    for value in groups:
        normalized = str(value).strip()
        if normalized == '':
            continue
        if normalized in seen:
            continue
        seen.add(normalized)
        resolved.append(normalized)
    if len(resolved) == 0:
        raise ValueError('No sample_group values were resolved for per-species finalize worker.')
    return resolved


def _standardize_metadata_all(metadata_df):
    out = metadata_df.copy()
    for column in ('instrument', 'bioproject'):
        if column not in out.columns:
            continue
        values = out.loc[:, column]
        is_missing = values.isna() | values.astype(str).eq('')
        out.loc[is_missing, column] = 'not_provided'
    return out


def _get_species_metadata(metadata_df, scientific_name, selected_sample_groups, count_columns):
    out = metadata_df.copy()
    is_sp = out.loc[:, 'scientific_name'].astype(str).eq(str(scientific_name))
    is_sample_group = out.loc[:, 'sample_group'].astype(str).isin(list(selected_sample_groups))
    out = out.loc[is_sp & is_sample_group, :].copy()
    if 'exclusion' not in out.columns:
        out.loc[:, 'exclusion'] = 'no'
    conditions = (
        out.loc[:, 'exclusion'].fillna('').astype(str).str.strip().str.lower().eq('no')
        & ~out.loc[:, 'run'].astype(str).isin(list(count_columns))
    )
    out.loc[conditions, 'exclusion'] = 'failed_quantification'
    return out


def _exclude_inappropriate_sample_from_tc(counts_df, metadata_df):
    is_not_excluded = metadata_df.loc[:, 'exclusion'].fillna('').astype(str).str.strip().str.lower().eq('no')
    run_ids = metadata_df.loc[is_not_excluded, 'run'].astype(str).tolist()
    keep_runs = [run_id for run_id in run_ids if run_id in counts_df.columns]
    return counts_df.loc[:, keep_runs].copy()


def _exclude_inappropriate_sample_from_eff_length(eff_length_df, counts_df):
    return eff_length_df.loc[:, [run_id for run_id in counts_df.columns if run_id in eff_length_df.columns]].copy()


def _sort_tc_and_metadata(counts_df, metadata_df, sort_columns=('sample_group', 'scientific_name', 'bioproject')):
    sorted_metadata = metadata_df.copy()
    present_columns = [column for column in sort_columns if column in sorted_metadata.columns]
    if len(present_columns) > 0:
        sorted_metadata = sorted_metadata.sort_values(by=present_columns, kind='mergesort').reset_index(drop=True)
    run_ids = [run_id for run_id in sorted_metadata.loc[:, 'run'].astype(str).tolist() if run_id in counts_df.columns]
    sorted_counts = counts_df.loc[:, run_ids].copy()
    return {
        'tc': sorted_counts,
        'sra': sorted_metadata,
    }


def _transform_raw_to_fpkm(counts_df, eff_length_df, metadata_df):
    if 'tmm_library_size' in metadata_df.columns:
        metadata_indexed = metadata_df.copy().set_index('run')
        library_sizes = metadata_indexed.loc[list(counts_df.columns), 'tmm_library_size'].astype(float).to_numpy()
    else:
        library_sizes = counts_df.sum(axis=0).to_numpy(dtype=float)
    effective_lengths = eff_length_df.loc[:, counts_df.columns].to_numpy(dtype=float)
    counts = counts_df.to_numpy(dtype=float)
    with numpy.errstate(divide='ignore', invalid='ignore'):
        values = counts / effective_lengths / library_sizes.reshape(1, -1) * 1e9
    values[~numpy.isfinite(values)] = 0.0
    return pandas.DataFrame(values, index=counts_df.index, columns=counts_df.columns)


def _transform_raw_to_tpm(counts_df, eff_length_df):
    counts = counts_df.to_numpy(dtype=float)
    effective_lengths = eff_length_df.loc[:, counts_df.columns].to_numpy(dtype=float)
    with numpy.errstate(divide='ignore', invalid='ignore'):
        x = counts / effective_lengths
        values = x * 1e6 / numpy.nansum(x, axis=0, keepdims=True)
    values[~numpy.isfinite(values)] = 0.0
    return pandas.DataFrame(values, index=counts_df.index, columns=counts_df.columns)


def _apply_transformation_logic(counts_df, eff_length_df, transform_method, batch_effect_alg, step, metadata_df):
    transform_method = str(transform_method)
    batch_effect_alg = str(batch_effect_alg)
    if batch_effect_alg in {'no', 'sva'}:
        bool_fpkm_tpm = step == 'before_batch'
        bool_log = step == 'before_batch'
    elif batch_effect_alg in {'ruvseq', 'combatseq', 'latent_glm'}:
        bool_fpkm_tpm = step in {'before_batch_plot', 'after_batch'}
        bool_log = step in {'before_batch_plot', 'after_batch'}
    else:
        raise ValueError('Unsupported batch effect algorithm for Python finalize worker: {}'.format(batch_effect_alg))

    transformed = counts_df.copy()
    if bool_fpkm_tpm:
        if 'fpkm' in transform_method:
            transformed = _transform_raw_to_fpkm(transformed, eff_length_df, metadata_df)
        elif 'tpm' in transform_method:
            transformed = _transform_raw_to_tpm(transformed, eff_length_df)
    if bool_log:
        values = transformed.to_numpy(dtype=float)
        with numpy.errstate(divide='ignore', invalid='ignore'):
            if 'logn-' in transform_method:
                values = numpy.log(values)
            elif 'log2-' in transform_method:
                values = numpy.log2(values)
            elif 'lognp1-' in transform_method:
                values = numpy.log(values + 1.0)
            elif 'log2p1-' in transform_method:
                values = numpy.log2(values + 1.0)
        transformed = pandas.DataFrame(values, index=transformed.index, columns=transformed.columns)
    return transformed


def _remove_nonexpressed_gene(counts_df):
    gene_sum = counts_df.sum(axis=1)
    return {
        'tc_ex': counts_df.loc[gene_sum > 0, :].copy(),
        'tc_ne': counts_df.loc[gene_sum == 0, :].copy(),
    }


def _run_batch_effect_step(counts_df, metadata_df, eff_length_df, args):
    transform_method = str(getattr(args, 'norm', 'log2p1-fpkm'))
    batch_effect_alg = str(getattr(args, 'batch_effect_alg', 'no'))
    clip_negative = bool(getattr(args, 'clip_negative', True))
    if batch_effect_alg == 'no':
        out = _sort_tc_and_metadata(counts_df, metadata_df)
        batch_info = initialize_batch_info(run_ids=out['tc'].columns, batch_effect_alg=batch_effect_alg)
        batch_info['skip_reason'] = 'batch_effect_alg_no'
        tc_batch_corrected = _apply_transformation_logic(
            counts_df=out['tc'],
            eff_length_df=eff_length_df,
            transform_method=transform_method,
            batch_effect_alg=batch_effect_alg,
            step='after_batch',
            metadata_df=out['sra'],
        )
        return {
            'tc': tc_batch_corrected,
            'sva': None,
            'batch_info': batch_info,
        }

    if counts_df.shape[1] == 1:
        batch_info = initialize_batch_info(run_ids=counts_df.columns, batch_effect_alg=batch_effect_alg)
        batch_info['skip_reason'] = 'single_sample'
        batch_info['batch_effect_alg_applied'] = 'no'
        return {
            'tc': counts_df.copy(),
            'sva': None,
            'batch_info': batch_info,
        }

    out = _sort_tc_and_metadata(counts_df, metadata_df)
    counts_sorted = out['tc']
    metadata_sorted = out['sra']
    nonexpressed = _remove_nonexpressed_gene(counts_sorted)
    counts_expressed = nonexpressed['tc_ex']
    counts_nonexpressed = nonexpressed['tc_ne']
    run_all = list(counts_expressed.columns)
    batch_info = initialize_batch_info(run_ids=run_all, batch_effect_alg=batch_effect_alg)

    if batch_effect_alg == 'sva':
        if metadata_sorted.loc[:, 'sample_group'].astype(str).nunique() <= 1:
            batch_info['skip_reason'] = 'sva_design_failed'
            corrected = counts_expressed.copy()
            sv_info = None
        else:
            corrected, sv_df, summary = run_sva_backend(
                counts_df=counts_expressed,
                metadata_df=metadata_sorted,
                nsv_setting=str(getattr(args, 'sva_nsv', 'auto')),
                B_setting=str(getattr(args, 'sva_B', 'auto')),
                B_auto_max=int(getattr(args, 'sva_B_auto_max', 100)),
                sample_group_column='sample_group',
                random_seed=getattr(args, 'seed', 'auto'),
            )
            batch_info['resolved_sva_nsv'] = summary.get('resolved_sva_nsv')
            batch_info['resolved_sva_B'] = summary.get('resolved_sva_B')
            batch_info['sva_estimation_method'] = summary.get('sva_estimation_method')
            batch_info['sva_stable'] = summary.get('sva_stable')
            batch_info['skip_reason'] = summary.get('skip_reason', '')
            batch_info['corrected_runs'] = summary.get('corrected_run_ids', [])
            sv_info = sv_df
        corrected_full = pandas.concat([corrected.loc[:, run_all], counts_nonexpressed.loc[:, run_all]], axis=0)
    elif batch_effect_alg == 'combatseq':
        corrected, summary = run_combatseq_backend(
            counts_df=counts_expressed,
            metadata_df=metadata_sorted,
            batch_column='bioproject',
            sample_group_column='sample_group',
        )
        batch_info['skip_reason'] = summary.get('skip_reason', '')
        batch_info['corrected_runs'] = summary.get('corrected_run_ids', [])
        sv_info = None
        corrected_full = pandas.concat([corrected.loc[:, run_all], counts_nonexpressed.loc[:, run_all]], axis=0)
    elif batch_effect_alg == 'ruvseq':
        corrected, w_df, summary = run_ruvseq_backend(
            counts_df=counts_expressed,
            metadata_df=metadata_sorted,
            control_mode=str(getattr(args, 'ruvseq_control_genes', 'auto')),
            k_setting=str(getattr(args, 'ruvseq_k', 'auto')),
            k_max=int(getattr(args, 'ruvseq_k_max', 5)),
            top_n=int(getattr(args, 'ruvseq_control_top_n', 1000)),
            min_controls=int(getattr(args, 'ruvseq_min_controls', 100)),
            batch_column='bioproject',
            sample_group_column='sample_group',
        )
        batch_info['resolved_ruv_k'] = summary.get('resolved_ruv_k')
        batch_info['resolved_ruv_controls'] = summary.get('resolved_ruv_controls')
        batch_info['ruv_baseline_score'] = summary.get('ruv_baseline_score')
        batch_info['ruv_selected_score'] = summary.get('ruv_selected_score')
        batch_info['ruv_selected_penalized_score'] = summary.get('ruv_selected_penalized_score')
        batch_info['ruv_penalty'] = summary.get('ruv_penalty')
        batch_info['skip_reason'] = summary.get('skip_reason', '')
        batch_info['corrected_runs'] = summary.get('corrected_run_ids', [])
        sv_info = w_df
        corrected_full = pandas.concat([corrected.loc[:, run_all], counts_nonexpressed.loc[:, run_all]], axis=0)
    elif batch_effect_alg == 'latent_glm':
        corrected, latent_df, summary = run_latent_glm_backend(
            counts_df=counts_expressed,
            metadata_df=metadata_sorted,
            family=str(getattr(args, 'latent_family', 'nb')),
            k_setting=str(getattr(args, 'latent_k', 'auto')),
            k_max=int(getattr(args, 'latent_k_max', 5)),
            sample_group_column='sample_group',
            max_iter=int(getattr(args, 'latent_max_iter', 200)),
            tol=float(getattr(args, 'latent_tol', 1e-5)),
        )
        batch_info['resolved_latent_k'] = summary.get('resolved_latent_k')
        batch_info['latent_family'] = summary.get('latent_family')
        batch_info['latent_iterations'] = summary.get('latent_iterations')
        batch_info['latent_objective'] = summary.get('latent_objective')
        batch_info['latent_converged'] = summary.get('latent_converged', summary.get('stable'))
        batch_info['skip_reason'] = summary.get('skip_reason', '')
        batch_info['corrected_runs'] = summary.get('corrected_run_ids', [])
        sv_info = latent_df
        corrected_full = pandas.concat([corrected.loc[:, run_all], counts_nonexpressed.loc[:, run_all]], axis=0)
    else:
        raise ValueError('Unsupported batch effect algorithm for Python finalize worker: {}'.format(batch_effect_alg))

    corrected_runs = [run_id for run_id in batch_info['corrected_runs'] if run_id in run_all]
    batch_info['corrected_runs'] = corrected_runs
    batch_info['uncorrected_runs'] = [run_id for run_id in run_all if run_id not in corrected_runs]
    batch_info['batch_effect_alg_applied'] = batch_effect_alg if len(corrected_runs) > 0 else 'no'
    if clip_negative and str(transform_method).startswith(('lognp1-', 'log2p1-')):
        negative_mask = corrected_full.to_numpy(dtype=float) < 0
        if negative_mask.any():
            corrected_full = corrected_full.copy()
            corrected_full.loc[:, :] = numpy.where(negative_mask, 0.0, corrected_full.to_numpy(dtype=float))

    corrected_after = _apply_transformation_logic(
        counts_df=corrected_full,
        eff_length_df=eff_length_df,
        transform_method=transform_method,
        batch_effect_alg=batch_effect_alg,
        step='after_batch',
        metadata_df=metadata_sorted,
    )
    return {
        'tc': corrected_after,
        'sva': sv_info,
        'batch_info': batch_info,
    }


def _compute_corr_matrix(counts_df, dist_method):
    matrix = counts_df.corr(method=str(dist_method))
    matrix = matrix.fillna(0.0)
    return matrix


def _compute_distance_matrix(corr_df):
    dist = 1.0 - corr_df.to_numpy(dtype=float)
    dist = (dist + dist.T) / 2.0
    numpy.fill_diagonal(dist, 0.0)
    return dist


def _compute_pca_coordinates(corr_df):
    matrix = corr_df.to_numpy(dtype=float)
    n = matrix.shape[0]
    if n <= 1:
        return numpy.zeros((n, 2), dtype=float)
    eigvals, eigvecs = numpy.linalg.eigh(matrix)
    order = numpy.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    coords = numpy.zeros((n, 2), dtype=float)
    for idx in range(min(2, n)):
        value = max(0.0, float(eigvals[idx]))
        coords[:, idx] = eigvecs[:, idx] * numpy.sqrt(value)
    return coords


def _compute_mds_coordinates(corr_df):
    dist = _compute_distance_matrix(corr_df)
    n = dist.shape[0]
    if n <= 1:
        return numpy.zeros((n, 2), dtype=float)
    centering = numpy.eye(n) - (numpy.ones((n, n), dtype=float) / float(n))
    gram = -0.5 * centering.dot(dist ** 2).dot(centering)
    eigvals, eigvecs = numpy.linalg.eigh(gram)
    order = numpy.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    coords = numpy.zeros((n, 2), dtype=float)
    for idx in range(min(2, n)):
        value = max(0.0, float(eigvals[idx]))
        coords[:, idx] = eigvecs[:, idx] * numpy.sqrt(value)
    return coords


def _compute_tsne_coordinates(counts_df):
    num_samples = counts_df.shape[1]
    coords = numpy.full((num_samples, 2), numpy.nan, dtype=float)
    if num_samples < 4:
        return coords
    max_perplexity = int((int(num_samples) - 1) // 3)
    if max_perplexity < 1:
        return coords
    try:
        from sklearn.manifold import TSNE
    except ImportError:
        return coords
    try:
        coords = TSNE(
            n_components=2,
            perplexity=float(min(30, max_perplexity)),
            random_state=1,
            init='pca',
            learning_rate='auto',
            method='exact',
        ).fit_transform(counts_df.transpose().to_numpy(dtype=float))
    except ValueError:
        return numpy.full((num_samples, 2), numpy.nan, dtype=float)
    return coords


def _sample_group_color_map(sample_groups):
    try:
        import matplotlib
        matplotlib.use('Agg', force=True)
        from matplotlib import pyplot
    except ImportError:
        return {group: '#1f77b4' for group in sample_groups}
    palette = pyplot.get_cmap('tab20')
    unique_groups = list(dict.fromkeys(sample_groups))
    colors = {}
    for idx, group in enumerate(unique_groups):
        colors[group] = palette(idx % max(1, palette.N))
    return colors


def _draw_embedding_panel(ax, coords, colors, labels, title, font_size, x_label, y_label):
    if coords.shape[0] == 0:
        ax.text(0.5, 0.5, 'No samples', ha='center', va='center', fontsize=font_size)
        ax.set_axis_off()
        return
    ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=50)
    for idx, label in enumerate(labels):
        ax.text(coords[idx, 0], coords[idx, 1], str(label), fontsize=max(4, font_size - 2))
    ax.set_title(title, fontsize=font_size)
    ax.set_xlabel(x_label, fontsize=font_size)
    ax.set_ylabel(y_label, fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size)


def _draw_dendrogram_panel(ax, corr_df, labels, font_size, title):
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
    dendrogram(
        linkage_matrix,
        labels=labels,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=max(4, font_size - 2),
        color_threshold=None,
    )
    ax.set_title(title, fontsize=font_size)
    ax.set_ylabel('Distance', fontsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)


def _compute_numeric_r2(values, covariate):
    y = pandas.to_numeric(pandas.Series(values), errors='coerce').to_numpy(dtype=float)
    x = pandas.to_numeric(pandas.Series(covariate), errors='coerce').to_numpy(dtype=float)
    valid = numpy.isfinite(y) & numpy.isfinite(x)
    if valid.sum() < 2:
        return numpy.nan
    y = y[valid]
    x = x[valid]
    if (numpy.nanmax(x) - numpy.nanmin(x)) <= 0:
        return numpy.nan
    design = numpy.column_stack([numpy.ones((x.shape[0],), dtype=float), x])
    beta, _, _, _ = numpy.linalg.lstsq(design, y, rcond=None)
    fitted = design @ beta
    rss = float(numpy.sum((y - fitted) ** 2))
    sst = float(numpy.sum((y - numpy.mean(y)) ** 2))
    if (not numpy.isfinite(sst)) or (sst <= 0):
        return numpy.nan
    return float(1.0 - (rss / sst))


def _compute_sva_summary_table(sv_df, metadata_df):
    if (sv_df is None) or (sv_df.shape[1] == 0):
        return pandas.DataFrame()
    aligned = metadata_df.copy().set_index('run', drop=False)
    aligned = aligned.loc[[run_id for run_id in sv_df.index if run_id in aligned.index], :].copy()
    if aligned.shape[0] == 0:
        return pandas.DataFrame()
    sv_aligned = sv_df.loc[aligned.index, :].copy()
    covariates = []
    if 'sample_group' in aligned.columns:
        covariates.append(('Sample group', 'factor', aligned['sample_group']))
    if 'bioproject' in aligned.columns:
        covariates.append(('BioProject', 'factor', aligned['bioproject']))
    if 'mapping_rate' in aligned.columns:
        covariates.append(('Mapping rate', 'numeric', aligned['mapping_rate']))
    if 'total_spots' in aligned.columns:
        covariates.append(('Log10 total reads', 'numeric', numpy.log10(pandas.to_numeric(aligned['total_spots'], errors='coerce'))))
    if 'total_bases' in aligned.columns:
        covariates.append(('Log10 total bases', 'numeric', numpy.log10(pandas.to_numeric(aligned['total_bases'], errors='coerce'))))
    if 'tmm_normalization_factor' in aligned.columns:
        covariates.append(('TMM normalization factor', 'numeric', aligned['tmm_normalization_factor']))
    rows = []
    index = []
    for label, mode, covariate in covariates:
        stats = []
        for sv_name in sv_aligned.columns:
            sv_values = pandas.to_numeric(sv_aligned.loc[:, sv_name], errors='coerce')
            if mode == 'factor':
                stats.append(compute_factor_r2(sv_values, covariate))
            else:
                stats.append(_compute_numeric_r2(sv_values, covariate))
        rows.append(stats)
        index.append(label)
    return pandas.DataFrame(rows, index=index, columns=sv_aligned.columns, dtype=float)


def _draw_sva_summary_panel(ax, sv_df, metadata_df, font_size, title):
    summary_df = _compute_sva_summary_table(sv_df=sv_df, metadata_df=metadata_df)
    if summary_df.shape[0] == 0:
        ax.text(0.5, 0.5, 'No SV summary data', ha='center', va='center', fontsize=font_size)
        ax.set_axis_off()
        return
    image = ax.imshow(summary_df.to_numpy(dtype=float), vmin=0.0, vmax=1.0, cmap='viridis', aspect='auto')
    ax.set_title(title, fontsize=font_size)
    ax.set_xticks(range(summary_df.shape[1]))
    ax.set_xticklabels(summary_df.columns.tolist(), rotation=90, fontsize=max(4, font_size - 2))
    ax.set_yticks(range(summary_df.shape[0]))
    ax.set_yticklabels(summary_df.index.tolist(), fontsize=max(4, font_size - 2))
    ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)


def save_quick_state_comparison_plot(
    tc_before,
    tc_after,
    metadata_df,
    dist_method,
    out_pdf_path,
    selected_sample_groups,
    transform_method,
    batch_effect_alg,
    sv_info=None,
    font_size=8,
):
    if (tc_before.shape[1] <= 1) or (tc_after.shape[1] <= 1):
        return None
    common_runs = [run_id for run_id in tc_before.columns if run_id in set(tc_after.columns)]
    common_runs = [run_id for run_id in common_runs if run_id in set(metadata_df.loc[:, 'run'].astype(str))]
    if len(common_runs) <= 1:
        return None
    before = tc_before.loc[:, common_runs].copy()
    after = tc_after.loc[:, common_runs].copy()
    metadata = metadata_df.set_index('run', drop=False).loc[common_runs, :].reset_index(drop=True)
    out = _sort_tc_and_metadata(before, metadata)
    before = out['tc']
    metadata = out['sra']
    after = after.loc[:, before.columns].copy()

    try:
        import matplotlib
        matplotlib.use('Agg', force=True)
        from matplotlib import pyplot
    except ImportError as exc:
        warnings.warn('matplotlib is required to generate batch comparison plot {}: {}'.format(out_pdf_path, exc))
        return None

    corr_before = _compute_corr_matrix(before, dist_method)
    corr_after = _compute_corr_matrix(after, dist_method)
    pca_before = _compute_pca_coordinates(corr_before)
    pca_after = _compute_pca_coordinates(corr_after)
    coords_before = _compute_mds_coordinates(corr_before)
    coords_after = _compute_mds_coordinates(corr_after)
    tsne_before = _compute_tsne_coordinates(before)
    tsne_after = _compute_tsne_coordinates(after)
    tau_before = sample_group_to_tau(
        tc_sample_group_df=sample_group_mean(before, metadata, selected_sample_groups)['tc_ave'],
        rich_annotation=False,
        transform_method=transform_method,
    )
    tau_after = sample_group_to_tau(
        tc_sample_group_df=sample_group_mean(after, metadata, selected_sample_groups)['tc_ave'],
        rich_annotation=False,
        transform_method=transform_method,
    )
    color_map = _sample_group_color_map(metadata.loc[:, 'sample_group'].astype(str).tolist())
    colors = [color_map[str(group)] for group in metadata.loc[:, 'sample_group'].astype(str)]
    labels = [_format_genus_species_label(str(run_id)).replace('\n', ' ') for run_id in before.columns]

    fig, axes = pyplot.subplots(7, 2, figsize=(12.0, 26.0))
    heatmaps = [
        (corr_before, axes[0, 0], 'Before {}'.format(batch_effect_alg)),
        (corr_after, axes[0, 1], 'After {}'.format(batch_effect_alg)),
    ]
    for corr_df, ax, title in heatmaps:
        im = ax.imshow(corr_df.to_numpy(dtype=float), vmin=-1.0, vmax=1.0, cmap='coolwarm')
        ax.set_title(title, fontsize=font_size)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=max(4, font_size - 2))
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=max(4, font_size - 2))
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    pca_panels = [
        (pca_before, axes[1, 0], 'PCA before'),
        (pca_after, axes[1, 1], 'PCA after'),
    ]
    for coords, ax, title in pca_panels:
        _draw_embedding_panel(ax, coords, colors, before.columns.tolist(), title, font_size, 'PC1', 'PC2')

    scatter_panels = [
        (coords_before, axes[2, 0], 'MDS before'),
        (coords_after, axes[2, 1], 'MDS after'),
    ]
    for coords, ax, title in scatter_panels:
        _draw_embedding_panel(ax, coords, colors, before.columns.tolist(), title, font_size, 'Axis 1', 'Axis 2')

    tsne_panels = [
        (tsne_before, axes[3, 0], 't-SNE before'),
        (tsne_after, axes[3, 1], 't-SNE after'),
    ]
    for coords, ax, title in tsne_panels:
        _draw_embedding_panel(ax, coords, colors, before.columns.tolist(), title, font_size, 't-SNE 1', 't-SNE 2')

    dendrogram_panels = [
        (corr_before, axes[4, 0], 'Dendrogram before'),
        (corr_after, axes[4, 1], 'Dendrogram after'),
    ]
    for corr_df, ax, title in dendrogram_panels:
        _draw_dendrogram_panel(ax, corr_df, labels, font_size, title)

    hist_panels = [
        (tau_before, axes[5, 0], 'Tau before'),
        (tau_after, axes[5, 1], 'Tau after'),
    ]
    breaks = numpy.arange(0.0, 1.000001, 0.05)
    for tau_df, ax, title in hist_panels:
        tau_values = pandas.to_numeric(tau_df.loc[:, 'tau'], errors='coerce').dropna().to_numpy(dtype=float)
        ax.hist(tau_values, bins=breaks, color='gray', edgecolor='black', linewidth=0.5)
        ax.set_title(title, fontsize=font_size)
        ax.set_xlabel('Tau', fontsize=font_size)
        ax.set_ylabel('Gene count', fontsize=font_size)
        ax.tick_params(axis='both', labelsize=font_size)

    if str(batch_effect_alg).lower() == 'sva':
        _draw_sva_summary_panel(axes[6, 0], sv_info, metadata, font_size, 'SV summary')
    else:
        axes[6, 0].text(0.5, 0.5, 'No SV summary for {}'.format(batch_effect_alg), ha='center', va='center', fontsize=font_size)
        axes[6, 0].set_axis_off()
    axes[6, 1].set_axis_off()
    axes[6, 1].text(
        0.0,
        1.0,
        'Sample groups are encoded by fill color.\nBioProjects are encoded by edge color.',
        ha='left',
        va='top',
        fontsize=font_size,
    )

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    fig.savefig(out_pdf_path)
    pyplot.close(fig)
    return out_pdf_path


def _write_empty_correlation_statistics(path):
    os.makedirs(os.path.dirname(os.path.realpath(path)), exist_ok=True)
    initialize_correlation_statistics().to_csv(path, sep='\t')


def run_finalize_python_worker(args, metadata, species_tag, input_dir):
    input_dir_abs = os.path.abspath(input_dir)
    species_dir = os.path.join(input_dir_abs, species_tag)
    count_path_candidates = [
        os.path.join(species_dir, species_tag + '_cstmm_counts.tsv'),
        os.path.join(species_dir, species_tag + '_est_counts.tsv'),
    ]
    count_path = next((path for path in count_path_candidates if os.path.isfile(path)), count_path_candidates[0])
    eff_length_path = os.path.join(species_dir, species_tag + '_eff_length.tsv')
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

    round_summary = initialize_round_summary()
    batch_info_current = initialize_batch_info(run_ids=sra.loc[:, 'run'].astype(str).tolist(), batch_effect_alg=str(getattr(args, 'batch_effect_alg', 'no')))

    tc = _exclude_inappropriate_sample_from_tc(counts_df, sra)
    sorted_out = _sort_tc_and_metadata(tc, sra)
    tc = sorted_out['tc']
    sra = sorted_out['sra']
    eff_length_species = _exclude_inappropriate_sample_from_eff_length(eff_length_df, tc)
    tc = _apply_transformation_logic(tc, eff_length_species, args.norm, args.batch_effect_alg, 'before_batch', sra)
    tc_tmp = _apply_transformation_logic(tc, eff_length_species, args.norm, args.batch_effect_alg, 'before_batch_plot', sra)
    is_input_zero = tc_tmp.eq(0)

    write_table_with_index_name(
        df=tc_tmp,
        file_path=os.path.join(dir_tsv, '{}.uncorrected.tc.tsv'.format(species_tag)),
        index_name='target_id',
    )
    original_sample_groups = list(selected_sample_groups)
    sample_group_out = sample_group_mean(tc_tmp, sra, selected_sample_groups)
    tc_sample_group_uncorrected = sample_group_out['tc_ave']
    selected_sample_groups = sample_group_out['selected_sample_groups']
    write_table_with_index_name(
        df=tc_sample_group_uncorrected,
        file_path=os.path.join(dir_tsv, '{}.uncorrected.sample_group.mean.tsv'.format(species_tag)),
        index_name='target_id',
    )

    if bool(getattr(args, 'skip_curation', False)):
        batch_info_current['skip_reason'] = 'skip_curation_requested'
        batch_info_current['batch_effect_alg_applied'] = 'no'
        batch_info_current['corrected_runs'] = []
        batch_info_current['uncorrected_runs'] = list(tc_tmp.columns)
        sra_out = annotate_metadata_with_batch_info(sra, batch_info_current)
        sra_out.to_csv(os.path.join(dir_tsv, '{}.metadata.tsv'.format(species_tag)), sep='\t', index=False)
        write_table_with_index_name(
            df=tc_tmp,
            file_path=os.path.join(dir_tsv, '{}.{}.tc.tsv'.format(species_tag, args.batch_effect_alg)),
            index_name='target_id',
        )
        write_table_with_index_name(
            df=tc_sample_group_uncorrected,
            file_path=os.path.join(dir_tsv, '{}.{}.sample_group.mean.tsv'.format(species_tag, args.batch_effect_alg)),
            index_name='target_id',
        )
        round_summary = append_round_summary(
            round_summary=round_summary,
            step='skip_curation',
            round_value=-1,
            reason='skip_curation_requested',
            runs_before=tc.columns,
            runs_after=tc.columns,
        )
        write_curation_summaries(
            round_summary=round_summary,
            metadata_df=sra_out,
            scientific_name=scientific_name,
            batch_effect_alg=str(args.batch_effect_alg),
            dir_tsv=dir_tsv,
            mapping_rate_cutoff=float(getattr(args, 'mapping_rate', 0.0)),
            correlation_threshold=float(getattr(args, 'correlation_threshold', 0.3)),
            one_outlier_per_iteration=bool(getattr(args, 'one_outlier_per_iter', False)),
            num_total_runs_species=num_total_runs_species,
            num_runs_after_sample_group_filter=num_runs_after_sample_group_filter,
            total_runtime_sec=0.0,
        )
        write_batch_effect_summary_tsv(
            batch_info=batch_info_current,
            scientific_name=scientific_name,
            species_tag=species_tag,
            dir_tsv=dir_tsv,
            random_seed_value=getattr(args, 'seed', None),
        )
        return 0

    if not bool(getattr(args, 'disable_auto_outlier_filter', False)):
        return None

    out = _run_batch_effect_step(tc, sra, eff_length_species, args)
    tc_batch_corrected = out['tc']
    batch_info_current = out['batch_info']
    if str(getattr(args, 'batch_effect_alg', 'no')).lower() in {'sva', 'combatseq', 'ruvseq', 'latent_glm'}:
        save_quick_state_comparison_plot(
            tc_before=tc_tmp,
            tc_after=tc_batch_corrected,
            metadata_df=sra,
            dist_method=str(getattr(args, 'dist_method', 'pearson')),
            out_pdf_path=os.path.join(dir_pdf, '{}.before_after.{}.pdf'.format(species_tag, args.batch_effect_alg)),
            selected_sample_groups=selected_sample_groups,
            transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
            batch_effect_alg=str(getattr(args, 'batch_effect_alg', 'no')),
            sv_info=out.get('sva'),
            font_size=8,
        )
    if bool(getattr(args, 'maintain_zero', True)):
        tc_batch_corrected = tc_batch_corrected.copy()
        aligned_zero = is_input_zero.reindex(index=tc_batch_corrected.index, columns=tc_batch_corrected.columns, fill_value=False)
        tc_batch_corrected = tc_batch_corrected.mask(aligned_zero, 0.0)

    round_summary = append_round_summary(
        round_summary=round_summary,
        step='auto_outlier_filter',
        round_value=0,
        reason='disabled',
        runs_before=tc.columns,
        runs_after=tc.columns,
    )
    sra_out = annotate_metadata_with_batch_info(sra, batch_info_current)
    sra_out.to_csv(os.path.join(dir_tsv, '{}.metadata.tsv'.format(species_tag)), sep='\t', index=False)
    write_table_with_index_name(
        df=tc_batch_corrected,
        file_path=os.path.join(dir_tsv, '{}.{}.tc.tsv'.format(species_tag, args.batch_effect_alg)),
        index_name='target_id',
    )
    corrected_sample_group = sample_group_mean(tc_batch_corrected, sra, selected_sample_groups)['tc_ave']
    write_table_with_index_name(
        df=corrected_sample_group,
        file_path=os.path.join(dir_tsv, '{}.{}.sample_group.mean.tsv'.format(species_tag, args.batch_effect_alg)),
        index_name='target_id',
    )
    _write_empty_correlation_statistics(os.path.join(dir_tsv, '{}.{}.correlation_statistics.tsv'.format(species_tag, args.batch_effect_alg)))
    save_tau_histogram_pdf(
        counts_df=tc_batch_corrected,
        metadata_df=sra,
        selected_sample_groups=selected_sample_groups,
        out_pdf_path=os.path.join(dir_pdf, '{}.tau_hist.{}.pdf'.format(species_tag, args.batch_effect_alg)),
        font_size=8,
        transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
        tc_sample_group_df=corrected_sample_group,
    )
    tc_tau = sample_group_to_tau(
        tc_sample_group_df=corrected_sample_group,
        rich_annotation=True,
        transform_method=str(getattr(args, 'norm', 'log2p1-fpkm')),
    )
    write_table_with_index_name(
        df=tc_tau,
        file_path=os.path.join(dir_tsv, '{}.{}.tau.tsv'.format(species_tag, args.batch_effect_alg)),
        index_name='target_id',
    )
    write_curation_summaries(
        round_summary=round_summary,
        metadata_df=sra_out,
        scientific_name=scientific_name,
        batch_effect_alg=str(args.batch_effect_alg),
        dir_tsv=dir_tsv,
        mapping_rate_cutoff=float(getattr(args, 'mapping_rate', 0.0)),
        correlation_threshold=float(getattr(args, 'correlation_threshold', 0.3)),
        one_outlier_per_iteration=bool(getattr(args, 'one_outlier_per_iter', False)),
        num_total_runs_species=num_total_runs_species,
        num_runs_after_sample_group_filter=num_runs_after_sample_group_filter,
        total_runtime_sec=0.0,
    )
    write_batch_effect_summary_tsv(
        batch_info=batch_info_current,
        scientific_name=scientific_name,
        species_tag=species_tag,
        dir_tsv=dir_tsv,
        random_seed_value=getattr(args, 'seed', None),
    )
    return 0


__all__ = [
    'run_finalize_python_worker',
    'should_use_python_finalize_worker',
    'save_quick_state_comparison_plot',
]
