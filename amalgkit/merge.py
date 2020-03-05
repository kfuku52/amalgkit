import pandas
from glob import glob
import os


def read_4th(fn, column, mode):

    return pandas.read_csv(fn, delim_whitespace=1, usecols=[0,column], index_col=0, header=0,names=[mode,'length_'+fn,'eff-length_'+fn, 'est-counts_'+fn, 'tpm_'+fn])


def exec_merge(args,mode):

    if mode == 'est_counts':
        column = 3
    if mode == 'eff_length':
        column = 2
    files = glob('*.tsv')

    sra_file = pandas.concat([read_4th(fn, column, mode) for fn in files], axis=1)
    sra_file = sra_file.rename(columns={col: col.split('_')[1] for col in sra_file.columns})
    sra_file.to_csv('./merged_file/'+mode+'_'+args.out_name+'.tsv', sep='\t')


def merge_main(args):

    if not os.path.exists(os.path.join(args.work_dir+'merged_file')):
        os.makedirs(os.path.join(args.work_dir+'merged_file'))
    exec_merge(args, "est_counts")
    exec_merge(args, "eff_length")
