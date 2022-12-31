import pandas as pd

import glob
import os
import platform
import re
import subprocess
import time

from amalgkit.sanity import check_getfastq_outputs
from amalgkit.util import *

def check_seqkit_dependency():
    print("checking SeqKit dependency")
    try:
        subprocess.run(['seqkit','-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("SeqKit dependency satisfied. Moving on.")
    except FileNotFoundError:
        raise FileNotFoundError("SeqKit not found. Please make sure SeqKit is installed properly.")


def get_fastq_stats(args):
    print("Starting integration of fastq-file metadata...")
    if not os.path.exists(args.fastq_dir):
        raise ValueError("PATH to fastq directory does not exist: {}".format(args.fastq_dir))
    all_files = os.listdir(args.fastq_dir)
    fastq_extension_regex = r'(\.fq$|\.fastq$|\.fq.gz$|\.fastq.gz$)'
    all_fastq_files = [ f for f in all_files if re.search(fastq_extension_regex, f) ]
    if len(all_fastq_files)==0:
        txt = 'No detected fastq files (with regex "{}") in: {}'.format(fastq_extension_regex, args.fastq_dir)
        raise ValueError(txt)
    id_list = list(map(os.path.basename, all_fastq_files))
    id_list = [basename.split(os.extsep)[0] for basename in id_list] # This is necessary to correctly parse single-end fastq files
    # split off last (!) underscore in case sample name includes multiple underscores.
    # Assuming naming format: sample1_1.fastq, sample1_2.fastq
    # sample_1_1.fastq, sample_1_2.fastq is also possible
    id_list = [basename.rsplit('_', 1)[0] for basename in id_list]
    # duplicates (i.e. doublets) should indicate paired-end library
    num_fastq_files = {id: id_list.count(id) for id in id_list}
    column_names = ['scientific_name', 'curate_group', 'run', 'read1_path','read2_path', 'is_sampled',
                    'is_qualified','exclusion', 'lib_layout', 'spot_length', 'total_spots', 'total_bases', 'size', 'private_file']
    tmp_metadata = pd.DataFrame(columns = column_names)
    row = 0
    for id in num_fastq_files:
        if num_fastq_files[id] == 1:
            lib_layout = 'single'
        if num_fastq_files[id] == 2:
            lib_layout = 'paired'
        if num_fastq_files[id] > 2:
            raise ValueError("found more than 2 files (", num_fastq_files[id], ") for id", id,". Please check filenames.")
        fastq_files = [ os.path.join(args.fastq_dir, f) for f in all_fastq_files if re.search(id+'[\._]', f) ]
        print("Found {} file(s) for ID {}. Lib-layout: {}".format(num_fastq_files[id], id,lib_layout), flush=True)
        print("Getting sequence statistics.", flush=True)
        tmp_file = os.path.join(args.out_dir, id+'_seqkit_stats.tmp')
        if args.accurate_size:
            print('--accurate_size set to yes. Running accurate sequence scan.')
            seqkit_command = ['seqkit', 'stats', '-T', '-j', str(args.threads), fastq_files[0]]
            seqkit_stdout = open(tmp_file, 'w')
            subprocess.run(seqkit_command, stdout=seqkit_stdout)
            seqkit_stdout.close()
            tmp_stat_df = pandas.read_csv(tmp_file, sep='\t', header=0)
            total_spots = tmp_stat_df.loc[0, 'num_seqs']
        else:
            OS = platform.system()
            if OS == 'Darwin':
                zcat_command = 'zcat < '
            elif OS == 'Linux':
                zcat_command = 'zcat '
            else:
                zcat_command = 'zcat '
                sys.stderr.write('zcat may not be supported by this OS: {}\n'.format(OS))
            seqkit_command = zcat_command + fastq_files[0] + ' | head -n 4000 | seqkit stats -T -j ' + str(
                args.threads)
            total_lines_command = 'echo $(' + zcat_command + str(fastq_files[0]) + ' | wc -l)'
            seqkit_stdout = open(tmp_file, 'w')
            subprocess.run(seqkit_command, shell=True, stdout=seqkit_stdout)
            seqkit_stdout.close()
            tmp_stat_df = pandas.read_csv(tmp_file, sep='\t', header=0)
            total_lines_bytes = subprocess.check_output(total_lines_command, shell=True)
            total_lines = int(total_lines_bytes.decode().replace('\n', ''))
            total_spots = int(total_lines / 4)
        if args.remove_tmp:
            os.remove(tmp_file)
        tmp_stat_df.loc[0,'id'] = id
        try:
            tmp_stat_df.loc[0, 'file2'] = fastq_files[1]
        except IndexError:
            tmp_stat_df.loc[0, 'file2'] = 'unavailable'
        tmp_metadata.loc[row, 'scientific_name'] = 'Please add in format: Genus species'
        tmp_metadata.loc[row,'curate_group'] = 'Please add'
        tmp_metadata.loc[row,'run'] = tmp_stat_df.loc[0,'id']
        tmp_metadata.loc[row,'read1_path'] = os.path.abspath(fastq_files[0])
        if tmp_stat_df.loc[0, 'file2'] != 'unavailable':
            tmp_metadata.loc[row,'read2_path'] = os.path.abspath(tmp_stat_df.loc[0,'file2'])
        else:
            tmp_metadata.loc[row, 'read2_path'] = 'unavailable'
        tmp_metadata.loc[row,'is_sampled'] = 'yes'
        tmp_metadata.loc[row,'is_qualified'] = 'yes'
        tmp_metadata.loc[row, 'exclusion'] = 'no'
        tmp_metadata.loc[row,'lib_layout'] = lib_layout
        tmp_metadata.loc[row,'total_spots'] = total_spots
        tmp_metadata.loc[row,'size'] = os.path.getsize(fastq_files[0])
        tmp_metadata.loc[row, 'private_file'] = 'yes'
        tmp_metadata.loc[row,'spot_length'] = int(tmp_stat_df.loc[0,'avg_len'])
        total_bases = total_spots * int(tmp_stat_df.loc[0,'avg_len'])
        if lib_layout == 'paired':
            total_bases = total_bases * 2
        tmp_metadata.loc[row,'total_bases'] = total_bases
        row += 1
    if not os.path.exists(os.path.join(args.out_dir, 'metadata')):
        os.makedirs(os.path.join(args.out_dir, 'metadata'))
    tmp_metadata.to_csv(os.path.join(args.out_dir,'metadata','metadata_private_fastq.tsv'), sep='\t', index=False)
    return tmp_metadata

def integrate_main(args):
    check_seqkit_dependency()
    if args.metadata=='inferred':
        relative_path = os.path.join(args.out_dir, 'metadata', 'metadata', 'metadata.tsv')
        metadata_path = os.path.realpath(relative_path)
    else:
        metadata_path = os.path.realpath(args.metadata)
    if os.path.exists(metadata_path):
        print('Merging existing metadata and private fastq info: {}'.format(metadata_path))
        metadata = load_metadata(args)
        metadata.df.loc[:,'private_file'] = 'no'
        print("scanning for getfastq output")
        sra_ids = metadata.df.loc[:,'run']
        data_available, data_unavailable = check_getfastq_outputs(args, sra_ids, metadata, args.out_dir)
        metadata.df.loc[metadata.df['run'].isin(data_available), 'data_available'] = 'yes'
        metadata.df.loc[metadata.df['run'].isin(data_unavailable), 'data_available'] = 'no'
        tmp_metadata = get_fastq_stats(args)
        df = pandas.concat([metadata.df, tmp_metadata])
        df.to_csv(os.path.join(args.out_dir,'metadata','metadata_updated_for_private_fastq.tsv'), sep='\t', index=False)
    else:
        print('Generating a new metadata table.')
        get_fastq_stats(args)

