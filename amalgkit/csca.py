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


def csca_main(args):
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print(",<ERROR> Rscript is not installed.")
        sys.exit(1)

    if args.dir_tc == "inferred":
        dir_tc = os.path.join(args.out_dir, 'curate/tables/')
    else:
        dir_tc = os.path.realpath(args.dir_tc)

    if args.dir_uncorrected_curate_group_mean == "inferred":
        dir_uncorrected_curate_group_mean = os.path.join(args.out_dir, 'curate/tables/')
    else:
        dir_uncorrected_curate_group_mean = os.path.realpath(args.dir_uncorrected_curate_group_mean)

    if args.dir_curate_group_mean == "inferred":
        dir_curate_group_mean = os.path.join(args.out_dir, 'curate/tables/')
    else:
        dir_curate_group_mean = os.path.realpath(args.dir_curate_group_mean)

    if args.dir_sra == "inferred":
        dir_sra = os.path.join(args.out_dir, 'curate/tables/')
    else:
        dir_sra = os.path.realpath(args.dir_sra)

    # norm method may not be needed, curate outputs are all in the same directory => can be inferred
    curate_group = get_curate_group(args)
    dir_out = args.out_dir
    file_species_tree = args.file_species_tree
    file_singlecopy = args.file_singlecopy
    file_orthogroup = args.file_orthogroup
    csca_path = os.path.dirname(os.path.realpath(__file__))
    r_script_path = csca_path + '/csca.r'

    subprocess.call(['Rscript',
                     r_script_path,
                     curate_group,
                     dir_out,
                     dir_tc,
                     dir_uncorrected_curate_group_mean,
                     dir_curate_group_mean,
                     dir_sra,
                     file_species_tree,
                     file_singlecopy,
                     file_orthogroup,
                     ])
