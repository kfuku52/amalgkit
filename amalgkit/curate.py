import datetime
import os
import re
import shutil
import subprocess
import sys
import warnings

from amalgkit.util import *


def get_sample_group(args, metadata):
    if args.sample_group is None:
        sample_group = metadata.df.loc[:, 'sample_group'].dropna().unique()
    else:
        sample_group = re.findall(r"[\w]+", args.sample_group)
    if (len(sample_group)==0):
        txt = 'The "sample_group" column in --metadata ({}) is not filled. Please manually edit the file.\n'
        txt += '`amalgkit curate` recognizes samples with the same string in this columns to belong the same group.\n'
        txt += 'Exiting.\n'
        sys.stderr.write(txt.format(args.metadata))
        sys.exit(1)
    print('Tissues to be included: {}'.format(', '.join(sample_group)))
    sample_group = '|'.join(sample_group)
    return sample_group

def run_curate_r_script(args, metadata, sp, input_dir):
    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    correlation_threshold = args.correlation_threshold
    intermediate = args.plot_intermediate
    sample_group = get_sample_group(args, metadata)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    r_script_path = os.path.join(dir_amalgkit_script, 'curate.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    path_curate_input_metadata = os.path.join(input_dir, 'metadata.tsv')
    len_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_eff_length.tsv')
    if 'cstmm' in input_dir:
        count_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_cstmm_counts.tsv')
    else:
        count_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_est_counts.tsv')
    if os.path.exists(count_file) and os.path.exists(len_file):
        print("Both counts and effective length files found: {}".format(sp), flush=True)
    else:
        if not os.path.exists(count_file):
            sys.stderr.write('Expected but undetected PATH of the count file: {}\n'.format(count_file))
        if not os.path.exists(len_file):
            sys.stderr.write('Expected but undetected PATH of the effective length file: {}\n'.format(len_file))
        sys.stderr.write('Skipping {}\n'.format(sp))
        print('Skipping {}'.format(sp), flush=True)
        return 1
    print("Starting Rscript to obtain curated {} values.".format(args.norm), flush=True)
    curate_r_exit_code = subprocess.call([
            'Rscript',
            r_script_path,
            count_file,
            path_curate_input_metadata,
            os.path.realpath(args.out_dir),
            len_file,
            dist_method,
            str(mr_cut),
            '0',
            str(intermediate),
            sample_group,
            args.sample_group_color,
            str(args.norm),
            str(args.one_outlier_per_iter),
            str(correlation_threshold),
            str(args.batch_effect_alg),
            str(args.clip_negative),
            str(args.maintain_zero),
            os.path.realpath(r_util_path),
            str(args.skip_curation)
         ])
    return curate_r_exit_code

def curate_main(args):
    check_rscript()
    if args.input_dir=='inferred':
        dir_merge = os.path.realpath(os.path.join(args.out_dir, 'merge'))
        dir_cstmm = os.path.realpath(os.path.join(args.out_dir, 'cstmm'))
        if os.path.exists(dir_cstmm):
            print('Subdirectory for amalgkit cstmm will be used as input: {}'.format(dir_cstmm))
            metadata = load_metadata(args, dir_subcommand='cstmm')
            input_dir = dir_cstmm
        else:
            print('Subdirectory for amalgkit merge will be used as input: {}'.format(dir_merge))
            metadata = load_metadata(args, dir_subcommand='merge')
            input_dir = dir_merge
    else:
        print('Input_directory: {}'.format(args.input_dir))
        metadata = load_metadata(args, dir_subcommand=os.path.basename(args.input_dir))
        input_dir = args.input_dir
    if ('tpm' in args.norm) & ('cstmm' in input_dir):
            txt = ("TPM and TMM are incompatible. "
                   "If input data are CSTMM-normalized, "
                   "please switch --norm to any of the 'fpkm' normalization methods instead.")
            sys.stderr.write(txt)
    is_selected = (metadata.df['exclusion']=='no')
    spp = metadata.df.loc[is_selected, 'scientific_name'].drop_duplicates().values
    curate_dir = os.path.join(args.out_dir, 'curate')
    if not os.path.exists(curate_dir):
        os.mkdir(curate_dir)
    print('Number of species in the selected metadata table ("exclusion"=="no"): {}'.format(len(spp)), flush=True)
    for sp in spp:
        sp = sp.replace(" ", "_")
        file_curate_completion_flag = os.path.join(curate_dir, sp, 'curate_completion_flag.txt')
        if os.path.exists(file_curate_completion_flag):
            if args.redo:
                print('Output file detected. Will be overwritten: {}'.format(sp), flush=True)
                shutil.rmtree(os.path.join(curate_dir, sp))
            else:
                print('Skipping. Output file detected: {}'.format(sp), flush=True)
                continue
        print('Starting: {}'.format(sp), flush=True)
        exit_status = run_curate_r_script(args, metadata, sp, input_dir)
        if exit_status == 0:
            with open(file_curate_completion_flag, 'w') as f:
                f.write('amalgkit curate completed at {}\n'.format(datetime.datetime.now()))
