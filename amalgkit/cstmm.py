import pandas

import glob
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

def check_cstmm_parameter_compatibility(args):
    if (args.orthogroup_table is None)&(args.dir_busco is None):
        raise Exception('One of --orthogroup_table and --dir_busco should be specified.')
    if (args.orthogroup_table is not None)&(args.dir_busco is not None):
        raise Exception('Only one of --orthogroup_table and --dir_busco should be specified.')

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

def generate_multisp_busco_table(dir_busco, outfile):
    print('Generating multi-species BUSCO table.', flush=True)
    col_names = ['busco_id', 'status', 'sequence', 'score', 'length', 'orthodb_url', 'description']
    species_infiles = [f for f in os.listdir(path=dir_busco) if f.endswith('.tsv')]
    print('BUSCO full tables for {} species were detected at: {}'.format(len(species_infiles), dir_busco), flush=True)
    for species_infile in species_infiles:
        print('Working on {}'.format(species_infile), flush=True)
        path_to_table = os.path.join(dir_busco, species_infile)
        if not os.path.exists(path_to_table):
            warnings.warn('full_table.tsv does not exist. Skipping: '.format(species_infile))
            continue
        tmp_table = pandas.read_table(path_to_table, sep='\t', header=None, comment='#', names=col_names)
        tmp_table.loc[:, 'sequence'] = tmp_table.loc[:, 'sequence'].str.replace(':[-\.0-9]*$', '', regex=True)
        tmp_table.loc[(tmp_table.loc[:, 'sequence'].isnull()), 'sequence'] = '-'
        tmp_table.loc[(tmp_table.loc[:, 'orthodb_url'].isnull()), 'orthodb_url'] = '-'
        tmp_table.loc[(tmp_table.loc[:, 'description'].isnull()), 'description'] = '-'
        if species_infile == species_infiles[0]:
            merged_table = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']]
            merged_table = merged_table.drop_duplicates(keep='first', inplace=False, ignore_index=True)
        else:
            is_mt_missing = (merged_table.loc[:, 'orthodb_url'] == '-')
            if is_mt_missing.sum() > 0:
                tmp_table2 = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']]
                tmp_table2 = tmp_table2.drop_duplicates(keep='first', inplace=False, ignore_index=True)
                merged_table.loc[is_mt_missing, 'orthodb_url'] = tmp_table2.loc[is_mt_missing, 'orthodb_url']
                merged_table.loc[is_mt_missing, 'description'] = tmp_table2.loc[is_mt_missing, 'description']
        tmp_table = tmp_table.loc[:, ['busco_id', 'sequence']].groupby(['busco_id'])['sequence'].apply(
            lambda x: ','.join(x))
        tmp_table = tmp_table.reset_index()
        species_colname = species_infile
        species_colname = re.sub('_', 'PLACEHOLDER', species_colname)
        species_colname = re.sub('[-\._].*', '',  species_colname)
        species_colname = re.sub('PLACEHOLDER', '_', species_colname)
        tmp_table = tmp_table.rename(columns={'sequence': species_colname})
        merged_table = merged_table.merge(tmp_table, on='busco_id', how='outer')
    merged_table.to_csv(outfile, sep='\t', index=None)

def cstmm_main(args):
    check_cstmm_dependency()
    check_cstmm_parameter_compatibility(args)
    dir_work = os.path.realpath(args.out_dir)
    dir_cstmm = os.path.join(dir_work, 'cstmm')
    if not os.path.exists(dir_cstmm):
        os.makedirs(dir_cstmm)
    if args.dir_count=='inferred':
        dir_count = os.path.join(dir_work, 'merge')
    else:
        dir_count = os.path.realpath(args.dir_count)
    if args.dir_busco is not None:
        file_orthogroup_table = os.path.join(dir_cstmm, 'multispecies_busco_table.tsv')
        generate_multisp_busco_table(dir_busco=args.dir_busco, outfile=file_orthogroup_table)
    elif args.orthogroup_table is not None:
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
    for f in glob.glob("tmp.amalgkit.*"):
        os.remove(f)