import pandas
from glob import glob
import os

def merge_main(args):
    quant_dir = os.path.join(args.work_dir, 'quant')
    merge_dir = os.path.join(args.work_dir, 'merge')
    if not os.path.exists(merge_dir):
        os.makedirs(os.path.join(merge_dir))
    sra_ids = os.listdir(quant_dir)

    quant_out_paths = list()
    for sra_id in sra_ids:
        quant_out_path = os.path.join(quant_dir, sra_id, sra_id+'_abundance.tsv')
        if os.path.exists(quant_out_path):
            quant_out_paths.append(quant_out_path)
    print('{:,} quant outfiles were detected.'.format(len(quant_out_paths)))
    for col in ['eff_length','est_counts','tpm']:
        out = pandas.read_csv(quant_out_paths[0], header=0, sep='\t').loc[:,['target_id',]]
        values = [ pandas.read_csv(p, header=0, sep='\t').loc[:,[col,]] for p in quant_out_paths ]
        for sra_id,value in zip(sra_ids,values):
            value.columns = [sra_id,]
        out = pandas.concat([out,]+values, axis=1, ignore_index=False, sort=False)
        outfile = os.path.join(merge_dir, col+'.tsv')
        print('Writing output file:', outfile)
        out.to_csv(outfile, sep='\t', index=False)