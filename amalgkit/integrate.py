import pandas as pd
from amalgkit.sanity import check_getfastq_outputs
from amalgkit.util import *
import glob
import subprocess
import time

def check_seqkit_dependency():
    print("checking SeqKit dependency")
    try:
        subprocess.run(['seqkit','-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("SeqKit dependency satisfied. Moving on.")
    except FileNotFoundError:
        raise FileNotFoundError("SeqKit not found. Please make sure SeqKit is installed properly.")


def get_fastq_stats(args):
    print("Starting integration of fastq-file metadata...")
    if os.path.exists(args.fastq_dir):
        fastq_files = glob.glob(os.path.join(args.fastq_dir, '*.fastq*'))
        fastq_files.extend(glob.glob(os.path.join(args.fastq_dir, '*.fq*')))
        if fastq_files:
                # Infer unique ID and lib-layout from filenames
                # Get basename
                id_list = list(map(os.path.basename, fastq_files))
                # split off extensions
                id_list = [basename.split(os.extsep)[0] for basename in id_list]
                # split off last (!) underscore in case sample name includes multiple underscores.
                # Assuming naming format: sample1_1.fastq, sample1_2.fastq
                # sample_1_1.fastq, sample_1_2.fastq is also possible
                id_list = [basename.rsplit('_', 1)[0] for basename in id_list]
                # duplicates (i.e. doublets) should indicate paired-end library
                id_dict = {id: id_list.count(id) for id in id_list}
                column_names = ['scientific_name', 'curate_group', 'run', 'read1_path','read2_path', 'is_sampled',
                                'is_qualified','exclusion', 'lib_layout', 'spot_length', 'total_spots', 'total_bases', 'size', 'private_file']
                tmp_metadata = pd.DataFrame(columns = column_names)
                row = 0
                for id in id_dict:
                    if id_dict[id] == 1:
                        lib_layout = 'single'
                    if id_dict[id] == 2:
                        lib_layout = 'paired'
                    if id_dict[id] > 2:
                        raise ValueError("found more than 2 files (", id_dict[id], ") for id", id,". Please check filenames.")

                    fastq_files = glob.glob(os.path.join(args.fastq_dir, id + "*"))

                    print("Found ",id_dict[id], " file[s] for ID ", id,". Lib-layout: ", lib_layout)
                    print("Getting sequence statistics.")
                    tmp_file = os.path.join(args.out_dir, 'seqkit_stats.tmp')
                    seqkit_stdout = open(tmp_file, 'w')
                    subprocess.run(['seqkit','stats', '-T', '-j', str(args.threads), fastq_files[0]], stdout=seqkit_stdout)

                    tmp_stat_df = pandas.read_csv(tmp_file, sep='\t', header=0)
                    tmp_stat_df.loc[0,'id'] = id
                    try:
                        tmp_stat_df.loc[0, 'file2'] = fastq_files[1]
                    except IndexError:
                        tmp_stat_df.loc[0, 'file2'] = 'no path'

                    subprocess.run(['rm', tmp_file])

                    tmp_metadata.loc[row, 'scientific_name'] = 'Please add in format: Genus species'
                    tmp_metadata.loc[row,'curate_group'] = 'Please add'
                    tmp_metadata.loc[row,'run'] = tmp_stat_df.loc[0,'id']
                    tmp_metadata.loc[row,'read1_path']= tmp_stat_df.loc[0,'file']
                    tmp_metadata.loc[row,'read2_path']= tmp_stat_df.loc[0,'file2']
                    tmp_metadata.loc[row,'is_sampled']= 'yes'
                    tmp_metadata.loc[row,'is_qualified']= 'yes'
                    tmp_metadata.loc[row, 'exclusion'] = 'no'
                    tmp_metadata.loc[row,'lib_layout']= lib_layout
                    tmp_metadata.loc[row,'spot_length']= tmp_stat_df.loc[0,'avg_len']
                    tmp_metadata.loc[row,'total_spots']= tmp_stat_df.loc[0,'num_seqs']
                    tmp_metadata.loc[row,'total_bases']= tmp_stat_df.loc[0,'sum_len']
                    tmp_metadata.loc[row,'size']= os.path.getsize(fastq_files[0])
                    tmp_metadata.loc[row, 'private_file'] = 'yes'

                    row += 1

        else:
            raise ValueError("ERROR: No files in the fastq directory. Please check PATH and presence of fastq files.")
    else:
        raise ValueError("ERROR: Path to fastq directory does not exist.")
    if not os.path.exists(os.path.join(args.out_dir, 'metadata')):
        os.makedirs(os.path.join(args.out_dir, 'metadata'))
    tmp_metadata.to_csv(os.path.join(args.out_dir,'metadata','metadata_private_fastq_' + time.strftime("%Y-%m-%d") + '.tsv'), sep='\t', index=False)
    return tmp_metadata

def integrate_main(args):
    check_seqkit_dependency()
    if args.metadata:
        if not os.path.exists(args.metadata):
            raise ValueError("Path to metadata table does not exist.")
        print("found metadata \n")
        print("reading metadata from:", args.metadata)
        metadata = load_metadata(args)
        metadata.df.loc[:,'private_file'] = 'no'
        print("scanning for getfastq output")
        sra_ids = metadata.df.loc[:,'run']
        data_available, data_unavailable = check_getfastq_outputs(args, sra_ids, metadata, args.out_dir)
        metadata.df.loc[metadata.df['run'].isin(data_available), 'data_available'] = 'yes'
        metadata.df.loc[metadata.df['run'].isin(data_unavailable), 'data_available'] = 'no'
        tmp_metadata = get_fastq_stats(args)
        df = pandas.concat([metadata.df, tmp_metadata])

        df.to_csv(os.path.join(args.out_dir,'metadata','metadata_updated_for_private_fastq_' + time.strftime("%Y-%m-%d") + '.tsv'), sep='\t', index=False)
    else:
        get_fastq_stats(args)

