from amalgkit.util import *

import re
import subprocess
import os
import shutil
import sys

def get_sample_group(args):
    if args.sample_group is None:
        metadata = load_metadata(args)
        sample_group = metadata.df.loc[:, 'sample_group'].dropna().unique()
    else:
        sample_group = re.findall(r"[\w]+", args.sample_group)
    print('sample_groups to be included: {}'.format(', '.join(sample_group)))
    sample_group = '|'.join(sample_group)
    return sample_group

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
    dir_csca = os.path.join(dir_out, 'csca')
    dir_csca_input_table = os.path.join(dir_csca, 'csca_input_symlinks')
    spp = get_spp_from_dir(dir_curate)
    generate_csca_input_symlinks(dir_csca_input_table, dir_curate, spp)
    sample_group = get_sample_group(args)
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
    file_genecount = os.path.join(dir_csca, 'multispecies_genecount.tsv')
    orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
    subprocess.call(['Rscript',
                     csca_r_script_path,
                     sample_group,
                     dir_out,
                     dir_csca_input_table,
                     file_orthogroup_table,
                     file_genecount,
                     r_util_path,
                     dir_csca,
                     args.batch_effect_alg,
                     ])
    for f in glob.glob("tmp.amalgkit.*"):
        os.remove(f)
