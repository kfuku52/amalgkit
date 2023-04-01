from amalgkit.util import *

import re
import subprocess
import os
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

def dir_table2spp(dir_curate_group_mean):
    files = os.listdir(dir_curate_group_mean)
    files = [ f for f in files if f.endswith('.uncorrected.curate_group.mean.tsv') ]
    spp = [ f.replace('.uncorrected.curate_group.mean.tsv', '') for f in files ]
    return spp

def csca_main(args):
    check_rscript()
    check_ortholog_parameter_compatibility(args)
    dir_out = os.path.realpath(args.out_dir)
    if args.dir_tc == "inferred":
        dir_tc = os.path.join(dir_out, 'curate', 'tables')
    else:
        dir_tc = os.path.realpath(args.dir_tc)

    if args.dir_uncorrected_curate_group_mean == "inferred":
        dir_uncorrected_curate_group_mean = os.path.join(dir_out, 'curate', 'tables')
    else:
        dir_uncorrected_curate_group_mean = os.path.realpath(args.dir_uncorrected_curate_group_mean)

    if args.dir_curate_group_mean == "inferred":
        dir_curate_group_mean = os.path.join(dir_out, 'curate', 'tables')
    else:
        dir_curate_group_mean = os.path.realpath(args.dir_curate_group_mean)

    if args.dir_sra == "inferred":
        dir_sra = os.path.join(dir_out, 'curate', 'tables')
    else:
        dir_sra = os.path.realpath(args.dir_sra)
    # norm method may not be needed, curate outputs are all in the same directory => can be inferred
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
    spp = dir_table2spp(dir_curate_group_mean)
    orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)

    subprocess.call(['Rscript',
                     csca_r_script_path,
                     curate_group,
                     dir_out,
                     dir_tc,
                     dir_uncorrected_curate_group_mean,
                     dir_curate_group_mean,
                     dir_sra,
                     file_orthogroup_table,
                     file_genecount,
                     r_util_path,
                     dir_csca,
                     ])
    for f in glob.glob("tmp.amalgkit.*"):
        os.remove(f)
