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
    prepare_output_dir,
    save_exclusion_plot_pdf,
)


def _build_curate_args(args, input_dir, tmp_out_dir):
    data = vars(args).copy()
    data['out_dir'] = tmp_out_dir
    data['input_dir'] = input_dir
    data['r_script_name'] = 'finalize.r'
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
    data.setdefault('ruvseq_control_genes', 'auto')
    data.setdefault('ruvseq_k', 'auto')
    data.setdefault('ruvseq_k_max', 5)
    data.setdefault('ruvseq_control_top_n', 1000)
    data.setdefault('ruvseq_min_controls', 100)
    data.setdefault('seed', 'auto')
    data.setdefault('sva_nsv', 'auto')
    data.setdefault('sva_B', 'auto')
    data.setdefault('sva_B_auto_max', 100)
    return SimpleNamespace(**data)


def _simplify_table_filename(filename, species, batch_effect_alg):
    mapping = {
        '{}.metadata.tsv'.format(species): '{}_metadata.tsv'.format(species),
        '{}.uncorrected.tc.tsv'.format(species): '{}_expression_uncorrected.tsv'.format(species),
        '{}.uncorrected.sample_group.mean.tsv'.format(species): '{}_sample_group_mean_uncorrected.tsv'.format(species),
        '{}.{}.tc.tsv'.format(species, batch_effect_alg): '{}_expression.tsv'.format(species),
        '{}.{}.sample_group.mean.tsv'.format(species, batch_effect_alg): '{}_sample_group_mean.tsv'.format(species),
        '{}.{}.tau.tsv'.format(species, batch_effect_alg): '{}_tau.tsv'.format(species),
        '{}.{}.correlation_statistics.tsv'.format(species, batch_effect_alg): '{}_correlation_statistics.tsv'.format(species),
        '{}.{}.curation_round_summary.tsv'.format(species, batch_effect_alg): '{}_curation_round_summary.tsv'.format(species),
        '{}.{}.curation_final_summary.tsv'.format(species, batch_effect_alg): '{}_curation_final_summary.tsv'.format(species),
        '{}.{}.batch_effect_summary.tsv'.format(species, batch_effect_alg): '{}_batch_effect_summary.tsv'.format(species),
    }
    if filename in mapping:
        return mapping[filename]
    prefix = '{}.'.format(species)
    if filename.startswith(prefix):
        remaining = filename[len(prefix):]
    else:
        remaining = filename
    stem, ext = os.path.splitext(remaining)
    stem = stem.replace('.', '_')
    return '{}_{}{}'.format(species, stem, ext)


def _copy_species_tables(curate_dir, finalize_dir, batch_effect_alg):
    for species in sorted(os.listdir(curate_dir)):
        src_tables = os.path.join(curate_dir, species, 'tables')
        if not os.path.isdir(src_tables):
            continue
        dst_species_dir = os.path.join(finalize_dir, species)
        if os.path.lexists(dst_species_dir):
            if os.path.islink(dst_species_dir) or os.path.isfile(dst_species_dir):
                os.remove(dst_species_dir)
            else:
                shutil.rmtree(dst_species_dir)
        os.makedirs(dst_species_dir, exist_ok=True)
        for name in sorted(os.listdir(src_tables)):
            src_path = os.path.join(src_tables, name)
            if not os.path.isfile(src_path):
                continue
            dst_name = _simplify_table_filename(
                filename=name,
                species=species,
                batch_effect_alg=batch_effect_alg,
            )
            shutil.copy2(src_path, os.path.join(dst_species_dir, dst_name))


def finalize_main(args):
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
    dir_finalize = os.path.join(out_root, 'finalize')
    prepare_output_dir(dir_finalize, redo=args.redo)
    tmp_out_dir = tempfile.mkdtemp(prefix='amalgkit_finalize_')
    try:
        curate_args = _build_curate_args(args=resolve_args, input_dir=input_dir, tmp_out_dir=tmp_out_dir)
        curate_main(curate_args)
        curate_dir = os.path.join(tmp_out_dir, 'curate')
        _copy_species_tables(
            curate_dir=curate_dir,
            finalize_dir=dir_finalize,
            batch_effect_alg=curate_args.batch_effect_alg,
        )
        copy_curate_species_pdfs(curate_dir=curate_dir, dst_dir=dir_finalize)
        merged_species_metadata = load_merged_species_metadata(curate_dir=curate_dir)
        merged_metadata = merge_metadata_by_run(metadata.df, merged_species_metadata)
        merged_metadata.to_csv(os.path.join(dir_finalize, 'metadata.tsv'), sep='\t', index=False)
        r_util_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'util.r')
        save_exclusion_plot_pdf(
            df_metadata=merged_metadata,
            out_pdf_path=os.path.join(dir_finalize, 'finalize_exclusion.pdf'),
            r_util_path=r_util_path,
            y_label='Sample count',
            font_size=8,
        )
    finally:
        if os.path.isdir(tmp_out_dir):
            shutil.rmtree(tmp_out_dir, ignore_errors=True)
