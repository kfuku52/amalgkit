import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import pandas

from amalgkit.command_context import CrossSpeciesFilterContext
from amalgkit.filter_utils import save_exclusion_plot_pdf, staged_output_dir
from amalgkit.metadata_utils import load_metadata
from amalgkit.orthology_utils import (
    check_ortholog_parameter_compatibility,
    generate_multisp_busco_table,
    orthogroup2genecount,
)
from amalgkit.outlier_utils import flag_margin_outliers
from amalgkit.runtime_utils import cleanup_tmp_amalgkit_files


def _normalize_sample_groups(values):
    groups = []
    for value in values:
        if value is None:
            continue
        if pandas.isna(value):
            continue
        normalized = str(value).strip()
        if normalized == '':
            continue
        groups.append(normalized)
    return list(dict.fromkeys(groups))


def _parse_sample_group_argument(sample_group_arg):
    tokens = re.split(r'[,\|]+', str(sample_group_arg))
    return _normalize_sample_groups(tokens)


def _iter_visible_subdirs(path_dir):
    subdirs = []
    with os.scandir(path_dir) as entries:
        for entry in entries:
            if (not entry.is_dir()) or entry.name.startswith('.') or entry.name.startswith('tmp.'):
                continue
            subdirs.append(entry.name)
    return sorted(subdirs)


def _collect_tsv_files(path_tables_dir):
    files = []
    with os.scandir(path_tables_dir) as entries:
        for entry in entries:
            if (not entry.is_file()) or (not entry.name.endswith('.tsv')):
                continue
            files.append((entry.name, entry.path))
    return sorted(files)


def get_sample_group_string(args, context=None):
    metadata = None
    if isinstance(context, CrossSpeciesFilterContext):
        metadata = context.metadata
    if args.sample_group is None:
        if metadata is None:
            metadata = load_metadata(args)
        if 'sample_group' not in metadata.df.columns:
            raise ValueError(
                'The "sample_group" column was not found in metadata. '
                'Please add this column or provide --sample_group.'
            )
        sample_group = _normalize_sample_groups(metadata.df.loc[:, 'sample_group'].tolist())
    else:
        sample_group = _parse_sample_group_argument(args.sample_group)
    if len(sample_group) == 0:
        raise ValueError('No sample_group was selected. Provide --sample_group or fill metadata sample_group column.')
    print('sample_groups to be included: {}'.format(', '.join(sample_group)))
    return '|'.join(sample_group)


def get_species_from_dir(per_species_dir):
    return _iter_visible_subdirs(per_species_dir)


def generate_input_symlinks(input_table_dir, per_species_dir, spp):
    if os.path.lexists(input_table_dir):
        if os.path.islink(input_table_dir):
            os.remove(input_table_dir)
        elif os.path.isdir(input_table_dir):
            import shutil
            shutil.rmtree(input_table_dir)
        else:
            raise NotADirectoryError(
                'Cross-species input path exists but is not a directory: {}'.format(input_table_dir)
            )
    os.makedirs(input_table_dir, exist_ok=True)
    src_by_filename = {}
    for sp in spp:
        path_tables_dir = os.path.join(per_species_dir, sp, 'tables')
        if not os.path.isdir(path_tables_dir):
            raise FileNotFoundError(
                'Per-species tables directory not found for species {}: {}'.format(sp, path_tables_dir)
            )
        tsv_files = _collect_tsv_files(path_tables_dir)
        if len(tsv_files) == 0:
            raise FileNotFoundError(
                'No TSV table file was found in per-species tables directory for species {}: {}'.format(
                    sp,
                    path_tables_dir,
                )
            )
        for file, path_src in tsv_files:
            existing_src = src_by_filename.get(file)
            if (existing_src is not None) and (os.path.realpath(existing_src) != os.path.realpath(path_src)):
                raise ValueError(
                    'Duplicate table filename across species in per-species output: {} ({} vs {})'.format(
                        file, existing_src, path_src
                    )
                )
            src_by_filename[file] = path_src
            path_dst = os.path.join(input_table_dir, file)
            if os.path.lexists(path_dst):
                if os.path.isdir(path_dst) and (not os.path.islink(path_dst)):
                    raise IsADirectoryError(
                        'Cross-species input destination exists but is a directory: {}'.format(path_dst)
                    )
                os.remove(path_dst)
            os.symlink(path_src, path_dst)


def _species_tag(scientific_name):
    return str(scientific_name).strip().replace(' ', '_')


def _normalize_exclusion(series):
    normalized = pandas.Series(series).fillna('').astype(str).str.strip().str.lower()
    normalized = normalized.replace('', 'no')
    return normalized


def _normalize_cross_species_metadata_table(df_metadata):
    required_cols = ['run', 'scientific_name', 'sample_group', 'exclusion']
    missing = [col for col in required_cols if col not in df_metadata.columns]
    if len(missing) > 0:
        raise ValueError('Required metadata columns are missing for cross-species filtering: {}'.format(', '.join(missing)))
    out = df_metadata.copy()
    for column in out.columns:
        if out[column].dtype == object:
            out.loc[:, column] = out.loc[:, column].fillna('').astype(str).str.strip()
    out.loc[:, 'run'] = out.loc[:, 'run'].astype(str).str.strip()
    out.loc[:, 'scientific_name'] = out.loc[:, 'scientific_name'].astype(str).str.strip()
    out.loc[:, 'sample_group'] = out.loc[:, 'sample_group'].astype(str).str.strip()
    out.loc[:, 'exclusion'] = _normalize_exclusion(out.loc[:, 'exclusion'])
    out.loc[:, 'species_tag'] = out.loc[:, 'scientific_name'].map(_species_tag)
    return out


def _prepare_metadata_table(dir_cross_species_input_table, selected_sample_groups, spp):
    metadata_paths = []
    for name in sorted(os.listdir(dir_cross_species_input_table)):
        if not name.endswith('.metadata.tsv'):
            continue
        path = os.path.join(dir_cross_species_input_table, name)
        if os.path.isfile(path):
            metadata_paths.append(path)
    if len(metadata_paths) == 0:
        raise FileNotFoundError('No metadata files found in the cross-species input table directory.')
    frames = [pandas.read_csv(path, sep='\t', low_memory=False) for path in metadata_paths]
    df_metadata = pandas.concat(frames, axis=0, ignore_index=True, sort=False)
    df_metadata = _normalize_cross_species_metadata_table(df_metadata)
    return df_metadata.loc[
        df_metadata['sample_group'].isin(selected_sample_groups)
        & df_metadata['species_tag'].isin(spp),
        :,
    ].copy()


def _load_expression_tables(dir_cross_species_input_table, spp_filled, batch_effect_alg):
    all_files = [
        name
        for name in sorted(os.listdir(dir_cross_species_input_table))
        if name.endswith('.tc.tsv') and os.path.isfile(os.path.join(dir_cross_species_input_table, name))
    ]
    uncorrected_suffix = '.uncorrected.tc.tsv'
    corrected_suffix = '.{}.tc.tsv'.format(batch_effect_alg)
    out = {'uncorrected': {}, 'corrected': {}}
    for sp in spp_filled:
        species_prefix = '{}.'.format(sp)
        uncorrected_matches = [name for name in all_files if name.startswith(species_prefix) and name.endswith(uncorrected_suffix)]
        corrected_matches = [name for name in all_files if name.startswith(species_prefix) and name.endswith(corrected_suffix) and not name.endswith(uncorrected_suffix)]
        if len(uncorrected_matches) > 1:
            raise ValueError('Multiple uncorrected tc tables matched species {}: {}'.format(sp, ', '.join(uncorrected_matches)))
        if len(corrected_matches) > 1:
            raise ValueError('Multiple corrected tc tables matched species {}: {}'.format(sp, ', '.join(corrected_matches)))
        if (len(uncorrected_matches) == 0) or (len(corrected_matches) == 0):
            continue
        uncorrected_path = os.path.join(dir_cross_species_input_table, uncorrected_matches[0])
        corrected_path = os.path.join(dir_cross_species_input_table, corrected_matches[0])
        out['uncorrected'][sp] = pandas.read_csv(uncorrected_path, sep='\t', index_col=0, low_memory=False)
        out['corrected'][sp] = pandas.read_csv(corrected_path, sep='\t', index_col=0, low_memory=False)
        out['uncorrected'][sp].columns = out['uncorrected'][sp].columns.map(str)
        out['corrected'][sp].columns = out['corrected'][sp].columns.map(str)
    return out


def _select_single_copy_orthogroups(file_orthogroup_table, file_genecount, spp_filled):
    df_gc = pandas.read_csv(file_genecount, sep='\t', low_memory=False).set_index('orthogroup_id')
    df_og = pandas.read_csv(file_orthogroup_table, sep='\t', low_memory=False).set_index('busco_id')
    present_species = [sp for sp in spp_filled if (sp in df_gc.columns) and (sp in df_og.columns)]
    if len(present_species) == 0:
        raise ValueError('No species columns overlapped between orthogroup inputs and per-species tables.')
    is_singlecopy = df_gc.loc[:, present_species].eq(1).all(axis=1)
    df_singleog = df_og.loc[is_singlecopy, present_species].copy()
    for sp in present_species:
        df_singleog = df_singleog.loc[df_singleog.loc[:, sp].fillna('').astype(str).str.strip() != '', :]
    return df_singleog


def _extract_ortholog_unaveraged_expression_table(df_singleog, unaveraged_tcs):
    orthologs = {'uncorrected': pandas.DataFrame(index=df_singleog.index), 'corrected': pandas.DataFrame(index=df_singleog.index)}
    for correction in ['uncorrected', 'corrected']:
        slices = []
        for sp in df_singleog.columns:
            tc = unaveraged_tcs[correction].get(sp)
            if tc is None:
                continue
            tc_prefixed = tc.copy()
            tc_prefixed.columns = ['{}_{}'.format(sp, col) for col in tc_prefixed.columns]
            row_idx = df_singleog.loc[:, sp].astype(str).tolist()
            tc_slice = tc_prefixed.reindex(row_idx)
            tc_slice.index = df_singleog.index
            slices.append(tc_slice)
        if len(slices) > 0:
            orthologs[correction] = pandas.concat(slices, axis=1)
        else:
            orthologs[correction] = pandas.DataFrame(index=df_singleog.index)
    return orthologs


def _safe_corr(left_values, right_values, method):
    left = pandas.to_numeric(pandas.Series(left_values), errors='coerce')
    right = pandas.to_numeric(pandas.Series(right_values), errors='coerce')
    valid = left.notna() & right.notna()
    if int(valid.sum()) <= 1:
        return numpy.nan
    try:
        return float(left.loc[valid].corr(right.loc[valid], method=method))
    except ValueError:
        return numpy.nan


def _calculate_correlation_within_group(df_metadata, ortholog_matrix, correction_label):
    target_col = 'within_group_cor_{}'.format(correction_label)
    nongroup_col = 'max_nongroup_cor_{}'.format(correction_label)
    out = df_metadata.copy()
    out.loc[:, target_col] = numpy.nan
    out.loc[:, nongroup_col] = numpy.nan
    if ortholog_matrix.shape[1] == 0:
        return out
    kept = out.loc[out['exclusion'].eq('no'), :].copy()
    kept.loc[:, 'sample_id'] = kept['species_tag'].astype(str) + '_' + kept['run'].astype(str)
    kept = kept.loc[kept['sample_id'].isin(ortholog_matrix.columns), :].copy()
    if kept.shape[0] == 0:
        return out
    group_keys = (kept['species_tag'].astype(str) + '_' + kept['sample_group'].astype(str)).tolist()
    kept.loc[:, 'species_sample_group'] = group_keys
    group_order = list(dict.fromkeys(group_keys))
    ortholog_med = pandas.DataFrame(index=ortholog_matrix.index, columns=group_order, dtype=float)
    for group_key in group_order:
        sample_ids = kept.loc[kept['species_sample_group'].eq(group_key), 'sample_id'].astype(str).tolist()
        tc_group = ortholog_matrix.loc[:, [sample_id for sample_id in sample_ids if sample_id in ortholog_matrix.columns]].copy()
        if tc_group.shape[1] == 0:
            continue
        ortholog_med.loc[:, group_key] = tc_group.median(axis=1, skipna=True).to_numpy(dtype=float)
    for _, row in kept.iterrows():
        sample_id = str(row['sample_id'])
        group_key = str(row['species_sample_group'])
        if (sample_id not in ortholog_matrix.columns) or (group_key not in ortholog_med.columns):
            continue
        sample_values = ortholog_matrix.loc[:, sample_id]
        within_cor = _safe_corr(sample_values, ortholog_med.loc[:, group_key], method='pearson')
        other_keys = [key for key in group_order if key != group_key]
        nongroup_values = [
            _safe_corr(sample_values, ortholog_med.loc[:, other_key], method='pearson')
            for other_key in other_keys
        ]
        nongroup_values = [value for value in nongroup_values if numpy.isfinite(value)]
        max_nongroup = max(nongroup_values) if len(nongroup_values) > 0 else numpy.nan
        row_idx = out.index[out['species_tag'].eq(row['species_tag']) & out['run'].astype(str).eq(str(row['run']))]
        out.loc[row_idx, target_col] = within_cor
        out.loc[row_idx, nongroup_col] = max_nongroup
    return out


def _fill_missing_by_row_mean(df):
    if df.shape[0] == 0 or df.shape[1] == 0:
        return df.copy()
    filled = df.copy()
    row_means = filled.mean(axis=1, skipna=True).fillna(0.0)
    return filled.T.fillna(row_means).T


def _resolve_matrix_for_embedding(matrix_df, missing_strategy):
    if str(missing_strategy).lower() == 'strict':
        out = matrix_df.dropna(axis=0, how='any').copy()
    else:
        out = _fill_missing_by_row_mean(matrix_df)
    return out


def _compute_pca_coordinates(matrix_df, missing_strategy):
    if matrix_df.shape[1] <= 1:
        return pandas.DataFrame(index=matrix_df.columns, columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
    filled = _resolve_matrix_for_embedding(matrix_df, missing_strategy=missing_strategy)
    if filled.shape[0] == 0:
        return pandas.DataFrame(index=matrix_df.columns, columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
    tc_corr = filled.corr(method='pearson').fillna(0.0)
    eigvals, eigvecs = numpy.linalg.eigh(tc_corr.to_numpy(dtype=float))
    order = numpy.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    out = pandas.DataFrame(index=tc_corr.index, columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'], dtype=float)
    for idx in range(min(5, eigvecs.shape[1])):
        scale = numpy.sqrt(max(float(eigvals[idx]), 0.0))
        out.iloc[:, idx] = eigvecs[:, idx] * scale
    return out


def _assign_pca_to_metadata(df_metadata, pca_df, suffix):
    out = df_metadata.copy()
    sample_ids = out['species_tag'].astype(str) + '_' + out['run'].astype(str)
    for pc in ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']:
        column_name = '{}_{}'.format(pc, suffix)
        out.loc[:, column_name] = sample_ids.map(pca_df.loc[:, pc] if pc in pca_df.columns else pandas.Series(dtype=float))
        if suffix == 'corrected':
            out.loc[:, pc] = out.loc[:, column_name]
    return out


def _apply_csfilter_outlier_flags(df_metadata, outlier_method='none', margin_threshold=0.0, robust_z_threshold=-2.5):
    out = df_metadata.copy()
    out.loc[:, 'cs_margin_uncorrected'] = pandas.to_numeric(out.get('within_group_cor_uncorrected'), errors='coerce') - pandas.to_numeric(out.get('max_nongroup_cor_uncorrected'), errors='coerce')
    out.loc[:, 'cs_margin_corrected'] = pandas.to_numeric(out.get('within_group_cor_corrected'), errors='coerce') - pandas.to_numeric(out.get('max_nongroup_cor_corrected'), errors='coerce')
    out.loc[:, 'cs_robust_z'] = numpy.nan
    out.loc[:, 'cs_outlier_candidate'] = False
    out.loc[:, 'cs_outlier_reason'] = ''
    if str(outlier_method) != 'robust_margin':
        return out
    idx = out.index[out['exclusion'].eq('no') & pandas.to_numeric(out['cs_margin_corrected'], errors='coerce').notna()]
    if len(idx) == 0:
        return out
    flagged = flag_margin_outliers(
        df=out.loc[idx, ['sample_group', 'cs_margin_corrected']].copy(),
        margin_col='cs_margin_corrected',
        group_col='sample_group',
        margin_threshold=float(margin_threshold),
        robust_z_threshold=float(robust_z_threshold),
        robust_z_col='cs_robust_z',
        outlier_col='cs_outlier_candidate',
    )
    out.loc[idx, 'cs_robust_z'] = flagged['cs_robust_z'].to_numpy()
    out.loc[idx, 'cs_outlier_candidate'] = flagged['cs_outlier_candidate'].fillna(False).astype(bool).to_numpy()
    outlier_idx = idx[out.loc[idx, 'cs_outlier_candidate'].fillna(False).astype(bool).to_numpy()]
    if len(outlier_idx) > 0:
        out.loc[outlier_idx, 'cs_outlier_reason'] = 'low_cross_species_group_correlation'
        out.loc[outlier_idx, 'exclusion'] = 'low_cross_species_group_correlation'
    return out


def _save_csfilter_scatter_plot(df_metadata, out_pdf_path):
    required_cols = {'within_group_cor_corrected', 'max_nongroup_cor_corrected'}
    if not required_cols.issubset(df_metadata.columns):
        return None
    plot_df = df_metadata.copy()
    plot_df.loc[:, 'within_group_cor_corrected'] = pandas.to_numeric(plot_df['within_group_cor_corrected'], errors='coerce')
    plot_df.loc[:, 'max_nongroup_cor_corrected'] = pandas.to_numeric(plot_df['max_nongroup_cor_corrected'], errors='coerce')
    plot_df = plot_df.loc[
        plot_df['within_group_cor_corrected'].notna() & plot_df['max_nongroup_cor_corrected'].notna(),
        :,
    ].copy()
    if plot_df.shape[0] == 0:
        return None
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    groups = plot_df['sample_group'].fillna('').astype(str).tolist()
    unique_groups = list(dict.fromkeys(groups))
    cmap = plt.get_cmap('tab20')
    color_map = {group: cmap(idx % max(1, cmap.N)) for idx, group in enumerate(unique_groups)}
    colors = [color_map[group] for group in groups]
    is_outlier = plot_df.get('cs_outlier_candidate', pandas.Series(False, index=plot_df.index)).fillna(False).astype(bool).to_numpy()
    ax.scatter(
        plot_df['max_nongroup_cor_corrected'].to_numpy(dtype=float),
        plot_df['within_group_cor_corrected'].to_numpy(dtype=float),
        c=colors,
        s=numpy.where(is_outlier, 60.0, 35.0),
        edgecolors=numpy.where(is_outlier, 'red', 'black'),
        linewidths=numpy.where(is_outlier, 1.2, 0.4),
        alpha=0.85,
    )
    ax.set_xlabel('max_nongroup_cor', fontsize=8)
    ax.set_ylabel('within_group_cor', fontsize=8)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(color='#d0d0d0', linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_pdf_path)
    plt.close(fig)
    return out_pdf_path


def _save_within_group_histogram(df_metadata, out_pdf_path):
    if 'within_group_cor_corrected' not in df_metadata.columns:
        return None
    plot_df = df_metadata.copy()
    plot_df.loc[:, 'within_group_cor_corrected'] = pandas.to_numeric(plot_df['within_group_cor_corrected'], errors='coerce')
    plot_df = plot_df.loc[plot_df['within_group_cor_corrected'].notna(), :].copy()
    if plot_df.shape[0] == 0:
        return None
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    bins = numpy.linspace(-1.0, 1.0, 41)
    groups = list(dict.fromkeys(plot_df['sample_group'].fillna('').astype(str).tolist()))
    cmap = plt.get_cmap('tab20')
    for idx, group in enumerate(groups):
        values = plot_df.loc[plot_df['sample_group'].astype(str).eq(group), 'within_group_cor_corrected'].to_numpy(dtype=float)
        ax.hist(values, bins=bins, alpha=0.45, label=group, color=cmap(idx % max(1, cmap.N)))
    ax.set_xlabel('within_group_cor', fontsize=8)
    ax.set_ylabel('Sample count', fontsize=8)
    ax.tick_params(axis='both', labelsize=8)
    if len(groups) > 0:
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(out_pdf_path)
    plt.close(fig)
    return out_pdf_path


def _save_heatmap_pdf(matrix_df, out_pdf_path):
    if matrix_df.shape[1] == 0:
        return None
    corr = _resolve_matrix_for_embedding(matrix_df, missing_strategy='row_mean').corr(method='pearson').fillna(0.0)
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(4.0, 0.25 * corr.shape[0]), max(4.0, 0.25 * corr.shape[0])))
    im = ax.imshow(corr.to_numpy(dtype=float), vmin=-1.0, vmax=1.0, cmap='coolwarm')
    ax.set_xticks(range(corr.shape[0]))
    ax.set_xticklabels(corr.columns.tolist(), rotation=90, fontsize=6)
    ax.set_yticks(range(corr.shape[0]))
    ax.set_yticklabels(corr.index.tolist(), fontsize=6)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_pdf_path)
    plt.close(fig)
    return out_pdf_path


def _save_pca_pdf(df_metadata, out_pdf_path):
    required_cols = {'PC1_corrected', 'PC2_corrected'}
    if not required_cols.issubset(df_metadata.columns):
        return None
    plot_df = df_metadata.copy()
    plot_df.loc[:, 'PC1_corrected'] = pandas.to_numeric(plot_df['PC1_corrected'], errors='coerce')
    plot_df.loc[:, 'PC2_corrected'] = pandas.to_numeric(plot_df['PC2_corrected'], errors='coerce')
    plot_df = plot_df.loc[plot_df['PC1_corrected'].notna() & plot_df['PC2_corrected'].notna(), :].copy()
    if plot_df.shape[0] == 0:
        return None
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    groups = plot_df['sample_group'].fillna('').astype(str).tolist()
    unique_groups = list(dict.fromkeys(groups))
    cmap = plt.get_cmap('tab20')
    color_map = {group: cmap(idx % max(1, cmap.N)) for idx, group in enumerate(unique_groups)}
    colors = [color_map[group] for group in groups]
    ax.scatter(
        plot_df['PC1_corrected'].to_numpy(dtype=float),
        plot_df['PC2_corrected'].to_numpy(dtype=float),
        c=colors,
        s=35.0,
        edgecolors='black',
        linewidths=0.4,
        alpha=0.85,
    )
    ax.set_xlabel('PC1', fontsize=8)
    ax.set_ylabel('PC2', fontsize=8)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(color='#d0d0d0', linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_pdf_path)
    plt.close(fig)
    return out_pdf_path


def _save_overview_pdf(matrix_df, df_metadata, out_pdf_path):
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.2))
    heatmap_path = out_pdf_path + '.heatmap.tmp.pdf'
    _ = heatmap_path
    corr = _resolve_matrix_for_embedding(matrix_df, missing_strategy='row_mean').corr(method='pearson').fillna(0.0) if matrix_df.shape[1] > 0 else pandas.DataFrame()
    if corr.shape[0] > 0:
        im = axes[0, 0].imshow(corr.to_numpy(dtype=float), vmin=-1.0, vmax=1.0, cmap='coolwarm')
        axes[0, 0].set_title('Sample correlation', fontsize=8)
        axes[0, 0].set_xticks([])
        axes[0, 0].set_yticks([])
        fig.colorbar(im, ax=axes[0, 0], fraction=0.046, pad=0.04)
    else:
        axes[0, 0].text(0.5, 0.5, 'No heatmap data', ha='center', va='center', fontsize=8)
        axes[0, 0].set_axis_off()

    scatter_df = df_metadata.copy()
    scatter_df.loc[:, 'within_group_cor_corrected'] = pandas.to_numeric(scatter_df.get('within_group_cor_corrected'), errors='coerce')
    scatter_df.loc[:, 'max_nongroup_cor_corrected'] = pandas.to_numeric(scatter_df.get('max_nongroup_cor_corrected'), errors='coerce')
    scatter_df = scatter_df.loc[
        scatter_df['within_group_cor_corrected'].notna() & scatter_df['max_nongroup_cor_corrected'].notna(),
        :,
    ].copy()
    if scatter_df.shape[0] > 0:
        axes[0, 1].scatter(
            scatter_df['max_nongroup_cor_corrected'].to_numpy(dtype=float),
            scatter_df['within_group_cor_corrected'].to_numpy(dtype=float),
            s=30.0,
            color='#1f77b4',
            edgecolors='black',
            linewidths=0.3,
            alpha=0.85,
        )
        axes[0, 1].set_xlabel('max_nongroup_cor', fontsize=8)
        axes[0, 1].set_ylabel('within_group_cor', fontsize=8)
        axes[0, 1].tick_params(axis='both', labelsize=8)
        axes[0, 1].set_title('Cross-species margin', fontsize=8)
    else:
        axes[0, 1].text(0.5, 0.5, 'No scatter data', ha='center', va='center', fontsize=8)
        axes[0, 1].set_axis_off()

    hist_values = pandas.to_numeric(df_metadata.get('within_group_cor_corrected'), errors='coerce').dropna().to_numpy(dtype=float)
    if hist_values.size > 0:
        axes[1, 0].hist(hist_values, bins=numpy.linspace(-1.0, 1.0, 41), color='gray', edgecolor='black', linewidth=0.4)
        axes[1, 0].set_xlabel('within_group_cor', fontsize=8)
        axes[1, 0].set_ylabel('Sample count', fontsize=8)
        axes[1, 0].tick_params(axis='both', labelsize=8)
        axes[1, 0].set_title('Within-group correlation', fontsize=8)
    else:
        axes[1, 0].text(0.5, 0.5, 'No histogram data', ha='center', va='center', fontsize=8)
        axes[1, 0].set_axis_off()

    pca_df = df_metadata.copy()
    pca_df.loc[:, 'PC1_corrected'] = pandas.to_numeric(pca_df.get('PC1_corrected'), errors='coerce')
    pca_df.loc[:, 'PC2_corrected'] = pandas.to_numeric(pca_df.get('PC2_corrected'), errors='coerce')
    pca_df = pca_df.loc[pca_df['PC1_corrected'].notna() & pca_df['PC2_corrected'].notna(), :].copy()
    if pca_df.shape[0] > 0:
        axes[1, 1].scatter(
            pca_df['PC1_corrected'].to_numpy(dtype=float),
            pca_df['PC2_corrected'].to_numpy(dtype=float),
            s=30.0,
            color='#ff7f0e',
            edgecolors='black',
            linewidths=0.3,
            alpha=0.85,
        )
        axes[1, 1].set_xlabel('PC1', fontsize=8)
        axes[1, 1].set_ylabel('PC2', fontsize=8)
        axes[1, 1].tick_params(axis='both', labelsize=8)
        axes[1, 1].set_title('Corrected PCA', fontsize=8)
    else:
        axes[1, 1].text(0.5, 0.5, 'No PCA data', ha='center', va='center', fontsize=8)
        axes[1, 1].set_axis_off()
    fig.tight_layout()
    fig.savefig(out_pdf_path)
    plt.close(fig)
    return out_pdf_path


def run_cross_species_filter(args, context=None):
    orthology_params = check_ortholog_parameter_compatibility(args)
    if orthology_params is None:
        orthogroup_table = getattr(args, 'orthogroup_table', None)
        dir_busco = getattr(args, 'dir_busco', None)
    else:
        orthogroup_table, dir_busco = orthology_params
    dir_out = os.path.realpath(args.out_dir)
    if os.path.exists(dir_out) and (not os.path.isdir(dir_out)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(dir_out))
    per_species_dir = os.path.join(dir_out, 'per_species')
    if os.path.exists(per_species_dir) and (not os.path.isdir(per_species_dir)):
        raise NotADirectoryError('Per-species output path exists but is not a directory: {}'.format(per_species_dir))
    if not os.path.isdir(per_species_dir):
        raise FileNotFoundError(
            'Per-species table output directory not found: {}. Per-species table generation must run before cross-species analysis.'
            .format(per_species_dir)
        )
    cross_species_dir = os.path.join(dir_out, 'cross_species')
    if os.path.exists(cross_species_dir) and (not os.path.isdir(cross_species_dir)):
        raise NotADirectoryError('Cross-species path exists but is not a directory: {}'.format(cross_species_dir))
    spp = get_species_from_dir(per_species_dir)
    if len(spp) == 0:
        raise ValueError('No per-species directories were found in: {}'.format(per_species_dir))
    selected_sample_groups = get_sample_group_string(args, context=context).split('|')
    with staged_output_dir(
        cross_species_dir,
        redo=bool(getattr(args, 'redo', False)),
        prefix='amalgkit_cross_species_stage_',
    ) as stage_dir:
        input_table_dir = os.path.join(stage_dir, 'cross_species_input_symlinks')
        generate_input_symlinks(input_table_dir, per_species_dir, spp)
        if dir_busco is not None:
            file_orthogroup_table = os.path.join(stage_dir, 'multispecies_busco_table.tsv')
            generate_multisp_busco_table(dir_busco=dir_busco, outfile=file_orthogroup_table)
        else:
            file_orthogroup_table = os.path.realpath(orthogroup_table)
            if not os.path.exists(file_orthogroup_table):
                raise FileNotFoundError('Orthogroup table not found: {}'.format(file_orthogroup_table))
            if not os.path.isfile(file_orthogroup_table):
                raise IsADirectoryError('Orthogroup table path exists but is not a file: {}'.format(file_orthogroup_table))
        file_genecount = os.path.join(stage_dir, 'multispecies_genecount.tsv')
        orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
        df_metadata = _prepare_metadata_table(input_table_dir, selected_sample_groups, spp)
        batch_effect_alg = str(getattr(args, 'batch_effect_alg', 'no'))
        unaveraged_tcs = _load_expression_tables(input_table_dir, spp, batch_effect_alg=batch_effect_alg)
        spp_filled = sorted(set(unaveraged_tcs['uncorrected']).intersection(unaveraged_tcs['corrected']))
        if len(spp_filled) == 0:
            raise FileNotFoundError('No matching cross-species tc tables were found for batch_effect_alg={}'.format(batch_effect_alg))
        df_singleog = _select_single_copy_orthogroups(file_orthogroup_table, file_genecount, spp_filled)
        orthologs = _extract_ortholog_unaveraged_expression_table(df_singleog, unaveraged_tcs)
        df_metadata = _calculate_correlation_within_group(df_metadata, orthologs['uncorrected'], 'uncorrected')
        df_metadata = _calculate_correlation_within_group(df_metadata, orthologs['corrected'], 'corrected')
        pca_corrected = _compute_pca_coordinates(
            orthologs['corrected'],
            missing_strategy=str(getattr(args, 'missing_strategy', 'em_pca')),
        )
        pca_uncorrected = _compute_pca_coordinates(
            orthologs['uncorrected'],
            missing_strategy=str(getattr(args, 'missing_strategy', 'em_pca')),
        )
        df_metadata = _assign_pca_to_metadata(df_metadata, pca_corrected, 'corrected')
        df_metadata = _assign_pca_to_metadata(df_metadata, pca_uncorrected, 'uncorrected')
        df_metadata = _apply_csfilter_outlier_flags(
            df_metadata,
            outlier_method=str(getattr(args, 'outlier_method', 'none')),
            margin_threshold=float(getattr(args, 'margin_threshold', 0.0)),
            robust_z_threshold=float(getattr(args, 'robust_z_threshold', -2.5)),
        )
        df_metadata.to_csv(os.path.join(stage_dir, 'metadata.tsv'), sep='\t', index=False)
        _save_overview_pdf(orthologs['corrected'], df_metadata, os.path.join(stage_dir, 'cross_species_overview.pdf'))
        _save_csfilter_scatter_plot(df_metadata, os.path.join(stage_dir, 'cross_species_csfilter_scatter.pdf'))
        _save_heatmap_pdf(orthologs['corrected'], os.path.join(stage_dir, 'cross_species_heatmap.pdf'))
        _save_within_group_histogram(df_metadata, os.path.join(stage_dir, 'cross_species_within_group_cor.pdf'))
        _save_pca_pdf(df_metadata, os.path.join(stage_dir, 'cross_species_unaveraged_pca_PC12.pdf'))
        save_exclusion_plot_pdf(
            df_metadata=df_metadata,
            out_pdf_path=os.path.join(stage_dir, 'cross_species_exclusion.pdf'),
            y_label='Sample count',
            font_size=8,
        )
    cleanup_tmp_amalgkit_files(work_dir='.')
