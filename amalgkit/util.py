import json
import numpy
import pandas

import glob
import inspect
import os
import re
import subprocess
import sys

from distutils.util import strtobool
from amalgkit.metadata import Metadata

def load_metadata(args):
    if args.metadata=='inferred':
        relative_path = os.path.join(args.out_dir, 'metadata', 'metadata', 'metadata.tsv')
        real_path = os.path.realpath(relative_path)
    else:
        real_path = os.path.realpath(args.metadata)
    print('Loading metadata from: {}'.format(real_path), flush=True)
    df = pandas.read_csv(real_path, sep='\t', header=0)
    metadata = Metadata.from_DataFrame(df)
    if 'batch' in dir(args):
        if args.batch is None:
            return metadata
        # --batch must be handled species-wise in curate.py
        # so we need to find out where the call came from
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        if mod.__name__ == 'amalgkit.curate':
            print('Entering --batch mode for amalgkit curate. processing 1 species')
            txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table.'
            spp = metadata.df.loc[:, 'scientific_name'].drop_duplicates().values
            print(txt.format(args.batch, len(spp)))
            sp = spp[args.batch - 1]
            print('processing species number ', args.batch, ' : ', sp)
            is_sampled = numpy.array([strtobool(yn) for yn in df.loc[:, 'is_sampled']], dtype=bool)
            metadata.df = metadata.df.loc[is_sampled, :]
            metadata.df = metadata.df.reset_index()
            metadata.df = metadata.df.loc[metadata.df['scientific_name'] == sp]
            return metadata

        print('--batch is specified. Processing one SRA per job.')
        is_sampled = numpy.array([strtobool(yn) for yn in df.loc[:, 'is_sampled']], dtype=bool)
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table. {:,} '
        txt += 'SRAs were excluded from the table (is_sampled==no).'
        print(txt.format(args.batch, sum(is_sampled), len(numpy.where(is_sampled == False)[0])))
        if args.batch>sum(is_sampled):
            sys.stderr.write('--batch {} is too large. Exiting.\n'.format(args.batch))
            sys.exit(0)
        if is_sampled.sum()==0:
            print('No sample is "sampled". Please check the "is_sampled" column in the metadata. Exiting.')
            sys.exit(1)
        metadata.df = metadata.df.loc[is_sampled,:]
        metadata.df = metadata.df.reset_index()
        metadata.df = metadata.df.loc[[args.batch-1,],:]
    return metadata

def get_sra_stat(sra_id, metadata, num_bp_per_sra=None):
    sra_stat = dict()
    sra_stat['sra_id'] = sra_id
    is_sra = (metadata.df.loc[:,'run']==sra_id)
    assert is_sra.sum()==1, 'There are multiple metadata rows with the same SRA ID: '+sra_id
    sra_stat['layout'] = metadata.df.loc[is_sra,'lib_layout'].values[0]
    sra_stat['total_spot'] = int(metadata.df.loc[is_sra,'total_spots'].values[0])
    original_spot_len = metadata.df.loc[is_sra,'spot_length'].values[0]
    if (numpy.isnan(original_spot_len) | (original_spot_len==0)):
        inferred_spot_len = int(metadata.df.loc[is_sra,'total_bases'].values[0]) / int(sra_stat['total_spot'])
        sra_stat['spot_length'] = int(inferred_spot_len)
        print('spot_length cannot be obtained directly from the metadata.')
        print('Using total_bases/total_spots instead: {:,}'.format(sra_stat['spot_length']))
    else:
        sra_stat['spot_length'] = int(original_spot_len)
    if num_bp_per_sra is not None:
        sra_stat['num_read_per_sra'] = int(num_bp_per_sra/sra_stat['spot_length'])
    return sra_stat

def get_newest_intermediate_file_extension(sra_stat, work_dir):
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz']
    if sra_stat['layout']=='single':
        subext = ''
    elif sra_stat['layout']=='paired':
        subext = '_1'
    files = os.listdir(work_dir)
    for ext in extensions:
        if any([ f==sra_stat['sra_id']+subext+ext for f in files ]):
            ext_out = ext
            break
    if not 'ext_out' in locals():
        safe_delete_files = glob.glob(os.path.join(work_dir, sra_stat['sra_id']+"*.safely_removed"))
        if len(safe_delete_files):
            txt = 'getfastq safely_removed flag was detected. `amalgkit quant` has been completed in this sample: {}\n'
            sys.stdout.write(txt.format(work_dir))
            for safe_delete_file in safe_delete_files:
                sys.stdout.write('{}\n'.format(safe_delete_file))
            return '.safely_removed'
        sys.stderr.write('getfastq output not found in: {}, layout = {}\n'.format(work_dir, sra_stat['layout']))
        txt = 'Skipping. If you wish to obtain the .fastq file(s), run: getfastq --id {}\n'
        sys.stderr.write(txt.format(sra_stat['sra_id']))
        raise FileNotFoundError
    return ext_out

def write_updated_metadata(metadata, outpath, args):
    try:
        overwrite_intermediate_metadata = args.overwrite_intermediate_metadata
    except AttributeError:
        overwrite_intermediate_metadata = 'yes'
    if os.path.exists(outpath):
        if not overwrite_intermediate_metadata:
            print('Intermediate metadata from previous run was detected and will not be overwritten.')
            return None
        else:
            print('Intermediate metadata from previous run was detected.')
            print('Intermediate metadata will be overwritten.')
            print('Preparing...')
    else:
        print('Intermediate metadata file was not detected. Preparing...')
    quant_dir = os.path.join(args.out_dir, 'quant')
    metadata = get_mapping_rate(metadata, quant_dir)
    print('Writing curate metadata containing mapping rate: {}'.format(outpath))
    metadata.df.to_csv(outpath, sep='\t', index=False)

def get_mapping_rate(metadata, quant_dir):
    if os.path.exists(quant_dir):
        print('quant directory found: {}'.format(quant_dir))
        metadata.df.loc[:, 'mapping_rate'] = numpy.nan
        sra_ids = metadata.df.loc[:, 'run'].values
        sra_dirs = [d for d in os.listdir(quant_dir) if d in sra_ids]
        print('Number of quant sub-directories that matched to metadata: {:,}'.format(len(sra_dirs)))
        for sra_id in sra_dirs:
            run_info_path = os.path.join(quant_dir, sra_id, sra_id + '_run_info.json')
            if not os.path.exists(run_info_path):
                sys.stderr.write('run_info.json not found. Skipping {}.\n'.format(sra_id))
                continue
            is_sra = (metadata.df.loc[:, 'run'] == sra_id)
            with open(run_info_path) as f:
                run_info = json.load(f)
            metadata.df.loc[is_sra, 'mapping_rate'] = run_info['p_pseudoaligned']
    else:
        txt = 'quant directory not found. Mapping rate cutoff will not be applied: {}\n'
        sys.stderr.write(txt.format(quant_dir))
    return metadata

def check_rscript():
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print("R (Rscript) is not installed. Exiting.")
        sys.exit(1)

def orthogroup2genecount(file_orthogroup, file_genecount, spp):
    df = pandas.read_csv(file_orthogroup, sep='\t', header=0, low_memory=False)
    is_spp = df.columns.isin(spp)
    df = df.loc[:,is_spp]
    df[df.isnull()] = ''
    df[df=='-'] = ''
    gc = df.copy()
    for col in gc.columns:
        is_comma = (df[col].str.contains(','))
        gc[col] = 0
        gc.loc[(df[col] != '') & ~is_comma, col] = 1
        gc.loc[is_comma, col] = df.loc[is_comma, col].str.count(',') + 1
    gc.to_csv(file_genecount, index=False, sep='\t')

def check_ortholog_parameter_compatibility(args):
    if (args.orthogroup_table is None)&(args.dir_busco is None):
        raise Exception('One of --orthogroup_table and --dir_busco should be specified.')
    if (args.orthogroup_table is not None)&(args.dir_busco is not None):
        raise Exception('Only one of --orthogroup_table and --dir_busco should be specified.')

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