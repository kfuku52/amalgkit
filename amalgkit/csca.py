from amalgkit.util import *

import re
import subprocess
import os
import shutil
import sys

def get_curate_group(args):
    if args.curate_group is None:
        metadata = load_metadata(args)
        curate_group = metadata.df.loc[:, 'curate_group'].dropna().unique()
    else:
        curate_group = re.findall(r"[\w]+", args.curate_group)
    print('curate_groups to be included: {}'.format(', '.join(curate_group)))
    curate_group = '|'.join(curate_group)
    return curate_group

def get_spp_from_dir(dir_curate):
    files = os.listdir(dir_curate)
    spp = [ f for f in files if (not f.startswith('.')) & (not f.startswith('tmp.')) ]
    return spp

def generate_csca_input_symlinks(dir_csca_input_table, dir_curate, spp):
    if os.path.exists(dir_csca_input_table):
        shutil.rmtree(dir_csca_input_table)
    os.makedirs(dir_csca_input_table)
    for sp in spp:
        files = os.listdir(os.path.join(dir_curate, sp, 'tables'))
        files = [ f for f in files if f.endswith('.tsv') ]
        for file in files:
            path_src = os.path.join(dir_curate, sp, 'tables', file)
            path_dst = os.path.join(dir_csca_input_table, file)
            if os.path.exists(path_dst):
                os.remove(path_dst)
            os.symlink(path_src, path_dst)
    return None

def csca_main(args):
    check_rscript()
    check_ortholog_parameter_compatibility(args)
    dir_out = os.path.realpath(args.out_dir)
    dir_curate = os.path.join(dir_out, 'curate')
    dir_csca_input_table = os.path.join(dir_out, 'tmp.amalgkit.csca_input_tables')
    spp = get_spp_from_dir(dir_curate)
    generate_csca_input_symlinks(dir_csca_input_table, dir_curate, spp)
    curate_group = get_curate_group(args)
    dir_csca = os.path.join(dir_out, 'csca')
    if not os.path.exists(dir_csca):
        os.makedirs(dir_csca)
    if args.dir_busco is not None:
        file_orthogroup_table = os.path.join(dir_csca, 'multispecies_busco_table.tsv')
        generate_multisp_busco_table(dir_busco=args.dir_busco, outfile=file_orthogroup_table)
    elif args.orthogroup_table is not None:
        file_orthogroup_table = os.path.realpath(args.orthogroup_table)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    csca_r_script_path = os.path.join(dir_amalgkit_script, 'csca.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    file_genecount = os.path.join(dir_out, 'tmp.amalgkit.orthogroup.genecount.tsv')
    orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
    subprocess.call(['Rscript',
                     csca_r_script_path,
                     curate_group,
                     dir_out,
                     dir_csca_input_table,
                     file_orthogroup_table,
                     file_genecount,
                     r_util_path,
                     dir_csca,
                     args.batch_effect_alg,
                     ])
    shutil.rmtree(dir_csca_input_table)
    for f in glob.glob("tmp.amalgkit.*"):
        os.remove(f)
