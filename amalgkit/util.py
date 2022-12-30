import json
import numpy
import pandas

import glob
import os
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
    if os.path.exists(outpath):
        if not args.overwrite_metadata:
            print('Updated metadata was detected and will not be overwritten.')
            return None
        else:
            print('Updated metadata was detected.')
            print('--overwrite_metadata option was set to yes. Metadata will be overwritten.')
            print('Preparing...')
    else:
        print('Updated metadata file was not detected. Preparing...')
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