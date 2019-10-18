import pandas
import glob
import os


def merge_main(args):

    in_files = glob.glob(args.work_dir + '/*.tsv')

    base = os.path.basename(in_files[0])
    sra_id = os.path.splitext(base)[0]
    sra_file = pandas.read_csv(in_files.pop(0), sep='\t')
    sra_file = sra_file[['target_id', sra_file.columns[4]]]

    sra_file = sra_file.rename(index=str, columns={sra_file.columns[1]: sra_id})

    for file_name in in_files:
        base = os.path.basename(file_name)
        merge_id = os.path.splitext(base)[0]
        file_to_merge = pandas.read_csv(file_name, sep='\t')
        file_to_merge = file_to_merge.rename(index=str, columns={file_to_merge.columns[4]: merge_id})

        sra_file = sra_file.merge(file_to_merge[['target_id', merge_id]], left_on='target_id', right_on='target_id')

    if not os.path.exists(os.path.join(args.work_dir, 'merged_file')):
        os.makedirs(os.path.join(args.work_dir, 'merged_file'))

    sra_file.to_csv(os.path.join(args.work_dir, 'merged_file', args.out_name), sep='\t')
