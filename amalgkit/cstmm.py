import pandas

import subprocess
import os
import re
import sys

def check_cstmm_dependency():
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        raise Exception("R (Rscript) is not installed.")

def get_count_files(dir_count):
    sciname_dirs = os.listdir(dir_count)
    count_files = list()
    for sciname_dir in sciname_dirs:
        sciname_path = os.path.join(dir_count, sciname_dir)
        if not os.path.isdir(sciname_path):
            continue
        files = os.listdir(sciname_path)
        for file in files:
            sciname_file = os.path.join(sciname_path, file)
            if not os.path.isfile(sciname_file):
                continue
            if not file.endswith('est_counts.tsv'):
                continue
            count_files.append(sciname_file)
        sciname_count_files = [ f for f in count_files if '/'+sciname_dir+'/' in f ]
        num_sciname_count_file = len(sciname_count_files)
        if (num_sciname_count_file==0):
            sys.stderr.write('No est_counts.tsv file found in: {}\n'.format(sciname_path))
        elif (num_sciname_count_file==1):
            print('One est_counts.tsv file found in {}'.format(sciname_path), flush=True)
        elif (num_sciname_count_file>=2):
            raise Exception('Multiple est_counts.tsv files found in: {}\n'.format(sciname_path))
    if (len(count_files)==0):
        raise Exception('No est_counts.tsv file was detected.')
    return count_files

def filepath2spp(file_paths):
    spp = [ os.path.basename(cf) for cf in file_paths]
    spp = [ re.sub('_', 'PLACEHOLDER', sp, count=1) for sp in spp ]
    spp = [ re.sub('_.*', '', sp) for sp in spp ]
    spp = [ re.sub('PLACEHOLDER', '_', sp) for sp in spp ]
    return spp

def orthogroup2genecount(file_orthogroup, file_genecount, spp):
    df = pandas.read_csv(file_orthogroup, sep='\t', header=0, low_memory=False)
    is_spp = df.columns.isin(spp)
    df = df.loc[:,is_spp]
    df[df.isnull()] = ''
    gc = df.copy()
    for col in gc.columns:
        is_comma = (df[col].str.contains(','))
        gc[col] = 0
        gc.loc[(df[col] != '') & ~is_comma, col] = 1
        gc.loc[is_comma, col] = df.loc[is_comma, col].str.count(',') + 1
    gc.to_csv(file_genecount, index=False, sep='\t')

def cstmm_main(args):
    check_cstmm_dependency()
    dir_work = os.path.realpath(args.out_dir)
    dir_cstmm = os.path.join(dir_work, 'cstmm')
    if not os.path.exists(dir_cstmm):
        os.makedirs(dir_cstmm)
    if args.dir_count=='inferred':
        dir_count = os.path.join(dir_work, 'merge')
    else:
        dir_count = os.path.realpath(args.dir_count)
    file_orthogroup_table = os.path.realpath(args.orthogroup_table)
    count_files = get_count_files(dir_count)
    if (len(count_files)==1):
        txt = 'Only one species was detected. Standard TMM normalization will be applied.'
        print(txt, flush=True)
        mode_tmm = 'single_species'
    else:
        txt = 'Multiple species were detected. ' \
              'Cross-species TMM normalization will be applied with single-copy orthologs.'
        print(txt, flush=True)
        mode_tmm = 'multi_species'
    file_genecount = 'tmp.amalgkit.orthogroup.genecount.tsv'
    spp = filepath2spp(count_files)
    orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    r_cstmm_path = os.path.join(dir_amalgkit_script, 'cstmm.r')
    r_command = ['Rscript', r_cstmm_path, dir_count, file_orthogroup_table, file_genecount, dir_cstmm, mode_tmm]
    print('')
    print('Starting R script: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)
