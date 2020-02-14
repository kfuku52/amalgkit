import pandas
import glob
import os


def make_calculations(args, mode):
    if mode == 'est_counts':
        column = 3
    if mode == 'eff_length':
        column = 2
    in_files = glob.glob(args.work_dir + '/*.tsv')

    base = os.path.basename(in_files[0])
    sra_id = os.path.splitext(base)[0]
    sra_file = pandas.read_csv(in_files.pop(0), sep='\t')

    sra_file = sra_file[['target_id', mode]]
    sra_file = sra_file.rename(index=str, columns={sra_file.columns[1]: sra_id})

    for file_name in in_files:
        base = os.path.basename(file_name)
        merge_id = os.path.splitext(base)[0]
        file_to_merge = pandas.read_csv(file_name, sep='\t')

        # merge RAW-counts based on transcript-ID
        file_to_merge = file_to_merge.rename(index=str, columns={file_to_merge.columns[column]: merge_id})

        sra_file = sra_file.merge(file_to_merge[['target_id', merge_id]], left_on='target_id', right_on='target_id')

    if not os.path.exists(os.path.join(args.work_dir, 'merged_file')):
        os.makedirs(os.path.join(args.work_dir, 'merged_file'))

    # clip off "_abundance" from col names - might become obsolete
    sra_file = sra_file.rename(columns={col: col.split('_')[0] for col in sra_file.columns})

    sra_file.to_csv(os.path.join(args.work_dir, 'merged_file', args.out_name+mode+'.tsv'), sep='\t', index=False)


def merge_main(args):
    make_calculations(args, 'est_counts')
    make_calculations(args, 'eff_length')
