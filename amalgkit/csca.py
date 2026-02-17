from amalgkit.util import *

import re
import subprocess
import os
import shutil
import sys

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

def get_sample_group_string(args):
    if args.sample_group is None:
        metadata = load_metadata(args)
        sample_group = _normalize_sample_groups(metadata.df.loc[:, 'sample_group'].tolist())
    else:
        sample_group = _parse_sample_group_argument(args.sample_group)
    if len(sample_group) == 0:
        raise ValueError('No sample_group was selected. Provide --sample_group or fill metadata sample_group column.')
    print('sample_groups to be included: {}'.format(', '.join(sample_group)))
    sample_group_string = '|'.join(sample_group)
    return sample_group_string

def get_spp_from_dir(dir_curate):
    return _iter_visible_subdirs(dir_curate)

def generate_csca_input_symlinks(dir_csca_input_table, dir_curate, spp):
    if os.path.exists(dir_csca_input_table):
        shutil.rmtree(dir_csca_input_table)
    os.makedirs(dir_csca_input_table, exist_ok=True)
    for sp in spp:
        path_tables_dir = os.path.join(dir_curate, sp, 'tables')
        for file, path_src in _collect_tsv_files(path_tables_dir):
            path_dst = os.path.join(dir_csca_input_table, file)
            if os.path.lexists(path_dst):
                os.remove(path_dst)
            os.symlink(path_src, path_dst)

def csca_main(args):
    check_rscript()
    check_ortholog_parameter_compatibility(args)
    dir_out = os.path.realpath(args.out_dir)
    dir_curate = os.path.join(dir_out, 'curate')
    dir_csca = os.path.join(dir_out, 'csca')
    dir_csca_input_table = os.path.join(dir_csca, 'csca_input_symlinks')
    spp = get_spp_from_dir(dir_curate)
    generate_csca_input_symlinks(dir_csca_input_table, dir_curate, spp)
    sample_group_string = get_sample_group_string(args)
    os.makedirs(dir_csca, exist_ok=True)
    if args.dir_busco is not None:
        file_orthogroup_table = os.path.join(dir_csca, 'multispecies_busco_table.tsv')
        generate_multisp_busco_table(dir_busco=args.dir_busco, outfile=file_orthogroup_table)
    elif args.orthogroup_table is not None:
        file_orthogroup_table = os.path.realpath(args.orthogroup_table)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    csca_r_script_path = os.path.join(dir_amalgkit_script, 'csca.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    file_genecount = os.path.join(dir_csca, 'multispecies_genecount.tsv')
    orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
    call_list = ['Rscript',
                 csca_r_script_path,
                 sample_group_string,
                 args.sample_group_color,
                 dir_out,
                 dir_csca_input_table,
                 file_orthogroup_table,
                 file_genecount,
                 r_util_path,
                 dir_csca,
                 args.batch_effect_alg,
                 args.missing_strategy,
                 ]
    print(f"Rscript command: {' '.join(call_list)}")
    subprocess.check_call(call_list)
    cleanup_tmp_amalgkit_files(work_dir='.')
