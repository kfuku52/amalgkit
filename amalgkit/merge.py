import pandas
import numpy

import os
import warnings
from amalgkit.util import *

def merge_fastp_stats_into_metadata(metadata, out_dir):
    cols = ['fastp_duplication_rate', 'fastp_insert_size_peak']
    for col in cols:
        if col not in metadata.df.columns:
            metadata.df.loc[:, col] = numpy.nan
    getfastq_dir = os.path.realpath(os.path.join(out_dir, 'getfastq'))
    if not os.path.exists(getfastq_dir):
        print('getfastq directory not found. Skipping fastp stats import: {}'.format(getfastq_dir), flush=True)
        return metadata
    num_detected = 0
    for i in metadata.df.index:
        sra_id = metadata.df.at[i, 'run']
        fastp_stats_path = os.path.join(getfastq_dir, sra_id, 'fastp_stats.tsv')
        if not os.path.exists(fastp_stats_path):
            continue
        try:
            df_fastp = pandas.read_csv(fastp_stats_path, sep='\t', header=0)
        except Exception as e:
            warnings.warn('Failed to read fastp stats. Skipping {}: {}'.format(fastp_stats_path, e))
            continue
        if df_fastp.shape[0] == 0:
            continue
        for col in cols:
            if col in df_fastp.columns:
                metadata.df.at[i, col] = df_fastp.loc[0, col]
        num_detected += 1
    print('{:,} fastp stats files were detected and merged into metadata.'.format(num_detected), flush=True)
    return metadata

def merge_main(args):
    quant_dir = os.path.realpath(os.path.join(args.out_dir, 'quant'))
    merge_dir = os.path.realpath(os.path.join(args.out_dir, 'merge'))
    if not os.path.exists(merge_dir):
        os.makedirs(os.path.join(merge_dir))
    metadata = load_metadata(args)
    spp = metadata.df.loc[:,'scientific_name'].dropna().unique()
    for sp in spp:
        print('processing: {}'.format(sp), flush=True)
        sp_filled = sp.replace(' ', '_')
        merge_species_dir = os.path.join(os.path.join(merge_dir, sp_filled))
        is_sp = (metadata.df.loc[:,'scientific_name']==sp)
        sra_ids = metadata.df.loc[is_sp,'run'].values
        is_sampled = (metadata.df.loc[:,'exclusion']=='no')
        sampled_sra_ids = metadata.df.loc[is_sampled,'run'].values
        if len(sra_ids)==0:
            warnings.warn('No SRA Run ID found. Skipping: {}'.format(sp))
            continue
        quant_out_paths = list()
        detected_sra_ids = list()
        for sra_id in sra_ids:
            quant_out_path = os.path.join(quant_dir, sra_id, sra_id+'_abundance.tsv')
            if os.path.exists(quant_out_path):
                if not os.path.exists(merge_species_dir):
                    os.makedirs(merge_species_dir)
                quant_out_paths.append(quant_out_path)
                detected_sra_ids.append(sra_id)
            else:
                if sra_id in sampled_sra_ids:
                    print('quant outfile not found: {}'.format(quant_out_path))
        print('{:,} quant outfiles were detected.'.format(len(quant_out_paths)))
        if len(quant_out_paths) == 0:
            continue
        for col in ['eff_length','est_counts','tpm']:
            out = pandas.read_csv(quant_out_paths[0], header=0, sep='\t').loc[:,['target_id',]]
            values = [ pandas.read_csv(p, header=0, sep='\t').loc[:,[col,]] for p in quant_out_paths ]
            for sra_id,value in zip(detected_sra_ids,values):
                value.columns = [sra_id,]
            out = pandas.concat([out,]+values, axis=1, ignore_index=False, sort=False)
            outfile_name = sp_filled+'_'+col+'.tsv'
            outfile = os.path.join(merge_species_dir, outfile_name)
            print('Writing output file:', outfile)
            out.to_csv(outfile, sep='\t', index=False)
    print('Getting mapping rate from quant output and write new metadata file into merge directory.', flush=True)
    metadata = merge_fastp_stats_into_metadata(metadata, args.out_dir)
    path_metadata_merge = os.path.realpath(os.path.join(args.out_dir, 'merge', 'metadata.tsv'))
    write_updated_metadata(metadata, path_metadata_merge, args)
    r_merge_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'merge.r')
    r_util_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'util.r')
    r_command = ['Rscript', r_merge_path, merge_dir, path_metadata_merge, r_util_path]
    print('Starting R script for plot generation: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)
