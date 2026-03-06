import os
import shutil
import tempfile
from types import SimpleNamespace

from amalgkit.curate import curate_main, resolve_curate_input
from amalgkit.filter_utils import (
    copy_curate_species_pdfs,
    infer_latest_filter_metadata,
    load_merged_species_metadata,
    merge_metadata_by_run,
    save_exclusion_plot_pdf,
    staged_output_dir,
)


def _build_curate_args(args, input_dir, tmp_out_dir):
    data = vars(args).copy()
    data['out_dir'] = tmp_out_dir
    data['input_dir'] = input_dir
    data['r_script_name'] = 'wsfilter.r'
    data.setdefault('margin_threshold', 0.0)
    data.setdefault('robust_z_threshold', -2.5)
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
    curate_args = SimpleNamespace(**data)
    curate_args._resolved_input_dir = input_dir
    return curate_args


def _write_excluded_table(df_metadata, out_path):
    reason = 'low_within_sample_group_correlation'
    exclusion_values = df_metadata['exclusion'].fillna('').astype(str).str.strip()
    excluded = df_metadata.loc[exclusion_values == reason, :].copy()
    preferred_cols = [
        'run',
        'scientific_name',
        'sample_group',
        'bioproject',
        'exclusion',
        'ws_within_group_cor',
        'ws_max_nongroup_cor',
        'ws_margin',
        'ws_robust_z',
    ]
    cols = [col for col in preferred_cols if col in excluded.columns]
    if len(cols) == 0:
        cols = ['run', 'exclusion']
    excluded = excluded.loc[:, cols]
    excluded = excluded.sort_values(by=[col for col in ['scientific_name', 'sample_group', 'run'] if col in excluded.columns])
    excluded.to_csv(out_path, sep='\t', index=False)


def wsfilter_main(args):
    resolve_args = args
    if getattr(args, 'metadata', 'inferred') == 'inferred':
        latest_metadata = infer_latest_filter_metadata(args.out_dir)
        if latest_metadata is not None:
            data = vars(args).copy()
            data['metadata'] = latest_metadata
            resolve_args = SimpleNamespace(**data)
            print('Using latest filter metadata: {}'.format(latest_metadata))
    metadata, input_dir = resolve_curate_input(resolve_args)
    out_root = os.path.realpath(args.out_dir)
    dir_ws = os.path.join(out_root, 'wsfilter')
    tmp_out_dir = tempfile.mkdtemp(prefix='amalgkit_wsfilter_')
    try:
        curate_args = _build_curate_args(args=resolve_args, input_dir=input_dir, tmp_out_dir=tmp_out_dir)
        curate_args._resolved_metadata = metadata
        curate_main(curate_args)
        merged_species_metadata = load_merged_species_metadata(curate_dir=os.path.join(tmp_out_dir, 'curate'))
        merged_metadata = merge_metadata_by_run(metadata.df, merged_species_metadata)
        with staged_output_dir(dir_ws, redo=args.redo, prefix='amalgkit_wsfilter_stage_') as stage_dir:
            out_metadata_path = os.path.join(stage_dir, 'metadata.tsv')
            merged_metadata.to_csv(out_metadata_path, sep='\t', index=False)
            _write_excluded_table(
                df_metadata=merged_metadata,
                out_path=os.path.join(stage_dir, 'excluded.tsv'),
            )
            r_util_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'util.r')
            save_exclusion_plot_pdf(
                df_metadata=merged_metadata,
                out_pdf_path=os.path.join(stage_dir, 'wsfilter_exclusion.pdf'),
                r_util_path=r_util_path,
                y_label='Sample count',
                font_size=8,
            )
            copy_curate_species_pdfs(
                curate_dir=os.path.join(tmp_out_dir, 'curate'),
                dst_dir=stage_dir,
            )
    finally:
        if os.path.isdir(tmp_out_dir):
            shutil.rmtree(tmp_out_dir, ignore_errors=True)
