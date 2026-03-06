import os
import re
import shutil
import subprocess

import pandas

from amalgkit.command_context import CrossSpeciesFilterContext
from amalgkit.filter_utils import staged_output_dir
from amalgkit.metadata_utils import load_metadata
from amalgkit.orthology_utils import (
    check_ortholog_parameter_compatibility,
    generate_multisp_busco_table,
    orthogroup2genecount,
)
from amalgkit.runtime_utils import check_rscript, cleanup_tmp_amalgkit_files
from amalgkit.subprocess_utils import run_logged_check_call

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
    # Preserve order while de-duplicating.
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
    sample_group_string = '|'.join(sample_group)
    return sample_group_string

def get_species_from_dir(per_species_dir):
    return _iter_visible_subdirs(per_species_dir)

def generate_input_symlinks(input_table_dir, per_species_dir, spp):
    if os.path.lexists(input_table_dir):
        if os.path.islink(input_table_dir):
            os.remove(input_table_dir)
        elif os.path.isdir(input_table_dir):
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

def run_cross_species_filter(args, context=None):
    check_rscript()
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
    sample_group_string = get_sample_group_string(args, context=context)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    cross_species_r_script_path = os.path.join(dir_amalgkit_script, 'cross_species_filter.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    from amalgkit.r_config import temporary_r_config
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
        elif orthogroup_table is not None:
            file_orthogroup_table = os.path.realpath(orthogroup_table)
            if not os.path.exists(file_orthogroup_table):
                raise FileNotFoundError('Orthogroup table not found: {}'.format(file_orthogroup_table))
            if not os.path.isfile(file_orthogroup_table):
                raise IsADirectoryError('Orthogroup table path exists but is not a file: {}'.format(file_orthogroup_table))
        file_genecount = os.path.join(stage_dir, 'multispecies_genecount.tsv')
        orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
        config_map = {
            'selected_sample_groups': sample_group_string,
            'sample_group_colors': args.sample_group_color,
            'dir_work': dir_out,
            'dir_cross_species_input_table': input_table_dir,
            'file_orthogroup': file_orthogroup_table,
            'file_genecount': file_genecount,
            'r_util_path': r_util_path,
            'dir_cross_species': stage_dir,
            'batch_effect_alg': args.batch_effect_alg,
            'missing_strategy': args.missing_strategy,
            'cross_species_outlier_method': str(getattr(args, 'outlier_method', 'none')),
            'cross_species_margin_threshold': str(getattr(args, 'margin_threshold', 0.0)),
            'cross_species_robust_z_threshold': str(getattr(args, 'robust_z_threshold', -2.5)),
            'cross_species_plot_mode': str(getattr(args, 'plot_mode', 'dual')),
        }
        with temporary_r_config(config_map, prefix='amalgkit_cross_species_r_') as config_path:
            call_list = ['Rscript', cross_species_r_script_path, config_path]
            run_logged_check_call(
                command=call_list,
                runner=subprocess.check_call,
                command_prefix='Rscript command',
                not_found_label='Rscript',
            )
    cleanup_tmp_amalgkit_files(work_dir='.')
