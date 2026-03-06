import os
import shutil
import tempfile
from types import SimpleNamespace

import pandas

from amalgkit.command_context import CrossSpeciesFilterContext, PerSpeciesTableContext
from amalgkit.cross_species_filter import run_cross_species_filter
from amalgkit.filter_utils import (
    infer_latest_filter_metadata,
    merge_metadata_by_run,
    save_exclusion_plot_pdf,
    staged_output_dir,
)
from amalgkit.per_species_tables import generate_per_species_tables, resolve_per_species_input


def _copy_csfilter_pdf_outputs(src_dir, dst_dir):
    if not os.path.isdir(src_dir):
        return
    os.makedirs(dst_dir, exist_ok=True)
    for name in sorted(os.listdir(src_dir)):
        if not name.lower().endswith('.pdf'):
            continue
        src_path = os.path.join(src_dir, name)
        if not os.path.isfile(src_path):
            continue
        # cross_species_exclusion.pdf duplicates csfilter_exclusion.pdf generated in Python.
        if name == 'cross_species_exclusion.pdf':
            continue
        # group_cor_scatter is intentionally suppressed in csfilter output.
        if name == 'cross_species_group_cor_scatter.pdf':
            continue
        if name == 'cross_species_csfilter_scatter.pdf':
            dst_name = 'csfilter_outlier_scatter.pdf'
        elif name.startswith('cross_species_'):
            dst_name = 'csfilter_' + name[len('cross_species_'):]
        else:
            dst_name = name
        shutil.copy2(src_path, os.path.join(dst_dir, dst_name))


def _resolve_sample_group_arg(args, metadata_df):
    if args.sample_group is not None:
        return args.sample_group
    if 'sample_group' not in metadata_df.columns:
        raise ValueError('Column "sample_group" is required in metadata when --sample_group is not specified.')
    values = (
        metadata_df['sample_group']
        .fillna('')
        .astype(str)
        .str.strip()
    )
    values = values.loc[values != ''].drop_duplicates().tolist()
    if len(values) == 0:
        raise ValueError('No sample_group values were found in metadata.')
    return '|'.join(values)


def _build_prepare_per_species_args(args, input_dir, tmp_out_dir):
    data = vars(args).copy()
    data['out_dir'] = tmp_out_dir
    data['input_dir'] = input_dir
    data['batch_effect_alg'] = 'no'
    data['r_script_name'] = 'prepare_tables.r'
    data['skip_curation'] = True
    data['disable_auto_outlier_filter'] = True
    data['outlier_method'] = 'legacy'
    data.setdefault('dist_method', 'pearson')
    data.setdefault('mapping_rate', 0.0)
    data.setdefault('correlation_threshold', 0.30)
    data.setdefault('plot_intermediate', False)
    data.setdefault('one_outlier_per_iter', False)
    data.setdefault('norm', 'log2p1-fpkm')
    data.setdefault('clip_negative', True)
    data.setdefault('maintain_zero', True)
    data.setdefault('batch', None)
    data.setdefault('threads', 'auto')
    data.setdefault('internal_jobs', 'auto')
    data.setdefault('internal_cpu_budget', 'auto')
    data.setdefault('sample_group_color', 'DEFAULT')
    return SimpleNamespace(**data)


def _build_cross_species_args(args, tmp_out_dir, sample_group_arg):
    data = vars(args).copy()
    data['out_dir'] = tmp_out_dir
    data['sample_group'] = sample_group_arg
    data['batch_effect_alg'] = 'no'
    data['plot_mode'] = 'single'
    data.setdefault('sample_group_color', 'DEFAULT')
    data.setdefault('missing_strategy', 'em_pca')
    data['outlier_method'] = 'robust_margin'
    data.setdefault('margin_threshold', 0.0)
    data.setdefault('robust_z_threshold', -2.5)
    cross_species_args = SimpleNamespace(**data)
    return cross_species_args


def _write_excluded_table(df_metadata, out_path):
    reason = 'low_cross_species_group_correlation'
    exclusion_values = df_metadata['exclusion'].fillna('').astype(str).str.strip()
    excluded = df_metadata.loc[exclusion_values == reason, :].copy()
    preferred_cols = [
        'run',
        'scientific_name',
        'sample_group',
        'bioproject',
        'exclusion',
        'within_group_cor',
        'max_nongroup_cor',
        'cs_margin',
        'cs_robust_z',
    ]
    cols = [col for col in preferred_cols if col in excluded.columns]
    if len(cols) == 0:
        cols = ['run', 'exclusion']
    excluded = excluded.loc[:, cols]
    excluded = excluded.sort_values(by=[col for col in ['scientific_name', 'sample_group', 'run'] if col in excluded.columns])
    excluded.to_csv(out_path, sep='\t', index=False)


def _normalize_csfilter_metadata_columns(df_metadata):
    df = df_metadata.copy()

    def _coalesce_from_candidates(target_col, candidate_cols):
        if target_col in df.columns:
            return
        for candidate in candidate_cols:
            if candidate in df.columns:
                df[target_col] = df[candidate]
                return

    _coalesce_from_candidates(
        target_col='within_group_cor',
        candidate_cols=['within_group_cor_corrected', 'within_group_cor_uncorrected'],
    )
    _coalesce_from_candidates(
        target_col='max_nongroup_cor',
        candidate_cols=['max_nongroup_cor_corrected', 'max_nongroup_cor_uncorrected'],
    )
    _coalesce_from_candidates(
        target_col='cs_margin',
        candidate_cols=['cs_margin_corrected', 'cs_margin_uncorrected'],
    )
    if ('cs_margin' not in df.columns) and {'within_group_cor', 'max_nongroup_cor'}.issubset(df.columns):
        within = pandas.to_numeric(df['within_group_cor'], errors='coerce')
        nongroup = pandas.to_numeric(df['max_nongroup_cor'], errors='coerce')
        df['cs_margin'] = within - nongroup

    for pc in ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']:
        corrected_col = f'{pc}_corrected'
        if (pc not in df.columns) and (corrected_col in df.columns):
            df[pc] = df[corrected_col]

    drop_cols = [
        'within_group_cor_uncorrected',
        'within_group_cor_corrected',
        'max_nongroup_cor_uncorrected',
        'max_nongroup_cor_corrected',
        'cs_margin_uncorrected',
        'cs_margin_corrected',
        'corrected_per_uncorrected_group_cor',
        'corrected_per_uncorrected_max_nongroup_cor',
        'PC1_uncorrected',
        'PC2_uncorrected',
        'PC3_uncorrected',
        'PC4_uncorrected',
        'PC5_uncorrected',
        'PC1_corrected',
        'PC2_corrected',
        'PC3_corrected',
        'PC4_corrected',
        'PC5_corrected',
    ]
    drop_cols = [col for col in drop_cols if col in df.columns]
    if len(drop_cols) > 0:
        df = df.drop(columns=drop_cols)
    return df


def csfilter_main(args):
    resolve_args = args
    if getattr(args, 'metadata', 'inferred') == 'inferred':
        latest_metadata = infer_latest_filter_metadata(args.out_dir)
        if latest_metadata is not None:
            data = vars(args).copy()
            data['metadata'] = latest_metadata
            resolve_args = SimpleNamespace(**data)
            print('Using latest filter metadata: {}'.format(latest_metadata))
    metadata, input_dir = resolve_per_species_input(resolve_args)
    out_root = os.path.realpath(args.out_dir)
    dir_cs = os.path.join(out_root, 'csfilter')
    tmp_out_dir = tempfile.mkdtemp(prefix='amalgkit_csfilter_')
    try:
        per_species_args = _build_prepare_per_species_args(args=resolve_args, input_dir=input_dir, tmp_out_dir=tmp_out_dir)
        generate_per_species_tables(
            per_species_args,
            context=PerSpeciesTableContext(metadata=metadata, input_dir=input_dir),
        )
        sample_group_arg = _resolve_sample_group_arg(args=resolve_args, metadata_df=metadata.df)
        cross_species_args = _build_cross_species_args(args=resolve_args, tmp_out_dir=tmp_out_dir, sample_group_arg=sample_group_arg)
        run_cross_species_filter(cross_species_args, context=CrossSpeciesFilterContext(metadata=metadata))
        cross_species_metadata_path = os.path.join(tmp_out_dir, 'cross_species', 'metadata.tsv')
        if not os.path.isfile(cross_species_metadata_path):
            raise FileNotFoundError('csfilter metadata.tsv was not generated: {}'.format(cross_species_metadata_path))
        cross_species_metadata = pandas.read_csv(cross_species_metadata_path, sep='\t', low_memory=False)
        merged_metadata = merge_metadata_by_run(metadata.df, cross_species_metadata)
        merged_metadata = _normalize_csfilter_metadata_columns(merged_metadata)
        with staged_output_dir(dir_cs, redo=args.redo, prefix='amalgkit_csfilter_stage_') as stage_dir:
            merged_metadata.to_csv(os.path.join(stage_dir, 'metadata.tsv'), sep='\t', index=False)
            _write_excluded_table(
                df_metadata=merged_metadata,
                out_path=os.path.join(stage_dir, 'excluded.tsv'),
            )
            r_util_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'util.r')
            save_exclusion_plot_pdf(
                df_metadata=merged_metadata,
                out_pdf_path=os.path.join(stage_dir, 'csfilter_exclusion.pdf'),
                r_util_path=r_util_path,
                y_label='Sample count',
                font_size=8,
            )
            _copy_csfilter_pdf_outputs(
                src_dir=os.path.join(tmp_out_dir, 'cross_species'),
                dst_dir=stage_dir,
            )
    finally:
        if os.path.isdir(tmp_out_dir):
            shutil.rmtree(tmp_out_dir, ignore_errors=True)
