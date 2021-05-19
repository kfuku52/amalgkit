import pandas
import os
import warnings

from amalgkit.util import *


def merge_main(args):
    quant_dir = os.path.join(args.work_dir, 'quant')
    merge_dir = os.path.join(args.work_dir, 'merge')
    if not os.path.exists(merge_dir):
        os.makedirs(os.path.join(merge_dir))

    if args.metadata is not None:
        print('Loading metadata from: {}'.format(args.metadata), flush=True)
        metadata = load_metadata(args)
        spp = metadata.df.loc[:,'scientific_name'].unique()
    else:
        raise Exception("If getfastq outputs are restructured like /getfastq/species_name/SRR00000,"
              "species names can be obtained from the directory names. Such change can drop "
              "the --metadata option, but not done yet. Use --metadata for now.")

    for sp in spp:
        print('processing: {}'.format(sp), flush=True)
        sp_filled = sp.replace(' ', '_')
        is_sp = (metadata.df.loc[:,'scientific_name']==sp)
        sra_ids = metadata.df.loc[is_sp,'run'].values
        if len(sra_ids)==0:
            warnings.warn('No SRA Run ID was found in {}.'.format(sp))
        quant_out_paths = list()
        detected_sra_ids = list()
        for sra_id in sra_ids:
            quant_out_path = os.path.join(quant_dir, sra_id, sra_id+'_abundance.tsv')
            if os.path.exists(quant_out_path):
                quant_out_paths.append(quant_out_path)
                detected_sra_ids.append(sra_id)
            else:
                print('quant outfile not found:', quant_out_path)
        print('{:,} quant outfiles were detected.'.format(len(quant_out_paths)))
        for col in ['eff_length','est_counts','tpm']:
            out = pandas.read_csv(quant_out_paths[0], header=0, sep='\t').loc[:,['target_id',]]
            values = [ pandas.read_csv(p, header=0, sep='\t').loc[:,[col,]] for p in quant_out_paths ]
            for sra_id,value in zip(detected_sra_ids,values):
                value.columns = [sra_id,]
            out = pandas.concat([out,]+values, axis=1, ignore_index=False, sort=False)
            outfile_name = sp_filled+'_'+col+'.tsv'
            outfile = os.path.join(merge_dir, outfile_name)
            print('Writing output file:', outfile)
            out.to_csv(outfile, sep='\t', index=False)