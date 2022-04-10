import pandas

import json
import re
import os
import subprocess
import sys
import warnings

from amalgkit.util import *

def get_curate_group(args):
    if args.curate_group is None:
        metadata = load_metadata(args)
        curate_group = metadata.df.loc[:, 'curate_group'].dropna().unique()
    else:
        curate_group = re.findall(r"[\w]+", args.curate_group)
    print('Tissues to be included: {}'.format(', '.join(curate_group)))
    curate_group = '|'.join(curate_group)
    return curate_group

def write_updated_metadata(metadata, outpath, args):
    if os.path.exists(outpath):
        print('Updated metadata file was detected and will not be overwritten.')
        return None
    else:
        print('Updated metadata file was not detected. Preparing...')
    if args.updated_metadata_dir == "inferred":
        updated_metadata_dir = os.path.realpath(os.path.join(args.out_dir, 'metadata/updated_metadata'))
    else:
        updated_metadata_dir = os.path.realpath(args.updated_metadata_dir)
    metadata = get_updated_metadata(metadata, updated_metadata_dir)
    quant_dir = os.path.join(args.out_dir, 'quant')
    metadata = get_mapping_rate(metadata, quant_dir)
    print('Writing updated metadata: {}'.format(outpath))
    metadata.df.to_csv(outpath, sep='\t', index=False)

def get_updated_metadata(metadata, updated_metadata_dir):
    if os.path.exists(updated_metadata_dir):
        print('Updated metadata directory found: {}'.format(updated_metadata_dir))
        files = os.listdir(updated_metadata_dir)
        updated_metadata_files = [ f for f in files if f.startswith('metadata')&f.endswith('.tsv') ]
        metadata_rows = list()
        for file in updated_metadata_files:
            file_path = os.path.join(updated_metadata_dir, file)
            tmp = pandas.read_csv(file_path, header=0, index_col=None, sep='\t')
            metadata_rows.append(tmp)
        tmp_concat = pandas.concat(metadata_rows, axis=0, ignore_index=True)
        cols = ['run',] + [ col for col in tmp_concat.columns if (col not in metadata.df.columns)&(col!='index') ]
        metadata.df = pandas.merge(metadata.df, tmp_concat.loc[:,cols], on='run', how='left')
    else:
        txt = 'Updated metadata directory not found. Some information may not appear in the plot: {}\n'
        sys.stderr.write(txt.format(updated_metadata_dir))
    return metadata

def get_mapping_rate(metadata, quant_dir):
    if os.path.exists(quant_dir):
        print('quant directory found: {}'.format(quant_dir))
        metadata.df.loc[:,'mapping_rate'] = numpy.nan
        sra_ids = metadata.df.loc[:,'run'].values
        sra_dirs = [ d for d in os.listdir(quant_dir) if d in sra_ids ]
        print('Number of quant sub-directories that matched to metadata: {:,}'.format(len(sra_dirs)))
        for sra_id in sra_dirs:
            run_info_path = os.path.join(quant_dir, sra_id, sra_id+'_run_info.json')
            if not os.path.exists(run_info_path):
                sys.stderr.write('run_info.json not found. Skipping {}.\n'.format(sra_id))
                continue
            is_sra = (metadata.df.loc[:,'run']==sra_id)
            with open(run_info_path) as f:
                run_info = json.load(f)
            metadata.df.loc[is_sra,'mapping_rate'] = run_info['p_pseudoaligned']
    else:
        txt = 'quant directory not found. Mapping rate cutoff will not be applied: {}\n'
        sys.stderr.write(txt.format(quant_dir))
    return metadata

def check_rscript():
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print(",<ERROR> Rscript is not installed.")
        sys.exit(1)

def curate_main(args):
    check_rscript()
    metadata = load_metadata(args)
    curate_dir = os.path.join(args.out_dir, 'curate')
    if not os.path.exists(curate_dir):
        os.mkdir(curate_dir)
    new_metadata_path = os.path.realpath(os.path.join(curate_dir, 'metadata.tsv'))
    write_updated_metadata(metadata, new_metadata_path, args)
    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    correlation_threshold = args.correlation_threshold
    intermediate = args.plot_intermediate
    curate_group = get_curate_group(args)
    curate_path = os.path.dirname(os.path.realpath(__file__))
    r_script_path = curate_path + '/transcriptome_curation.r'
    batch_script_path = curate_path + '/batch_curate.sh'
    quant_out = os.path.realpath(args.infile)
    eff_len_out = "NA"
    len_files = []
    count_files = []

    # check if cstmm output is used when --norm == tpm, because TPM undoes tmm normalization
    if args.norm == 'tpm':
        substring = "cstmm"
        try:
            quant_out.index(substring)
        except ValueError:
            warnings.warn("WARNING: TPM NORMALIZATION AND TMM NORMALIZATION ARE INCOMPATIBLE. IF INPUT DATA IS TMM NORMALIZED, PLEASE SWITCH --norm TO 'fpkm' INSTEAD.")
        else:
            raise ValueError("ERROR: AMALGKIT CSTMM NORMALIZED INPUT FILES DETECTED WHILE NORMALIZATION METHOD IS 'TPM'. TMM NORMALIZATION AND TPM NORMALIZATION ARE INCOMPATIBLE! PLEASE SWITCH --norm TO 'fpkm' INSTEAD." )


    # Input checks #
    # if single species mode active
    if args.batch is None:
        # Here, the user needs to provide single files for counts and effective lengths,
        # or just one file for abundances if the user already calculated fpkm/tpm
        # so there needs to be a warning, if there are insufficient user inputs
        if args.infile and args.eff_len_file:
            print("Single species mode")
            print("Both counts and effective length provided.")
            print("Calculating: ", args.norm)
            eff_len_out = os.path.realpath(args.eff_len_file)
        if args.infile and not args.eff_len_file:
            print("Single species mode")
            print("Only expression values provided. ")
            print("Assuming normalized expression values (like fpkm or tpm).")

        if not args.infile:
            print("No expression data provided. Please set Either --infile or --infile AND --eff_len_file")
            sys.exit(1)

        subprocess.call(['Rscript',
                         r_script_path,
                         quant_out,
                         new_metadata_path,
                         os.path.realpath(args.out_dir),
                         eff_len_out,
                         dist_method,
                         str(mr_cut),
                         '0',
                         str(intermediate),
                         curate_group,
                         str(args.norm),
                         str(args.one_outlier_per_iter),
                         str(correlation_threshold),
                         ])

    # if multiple species mode active
    if args.batch is not None:

        # Here, we still need arguments for counts and effective lengths,
        # But they need to be directories
        if args.infile_dir and args.eff_len_dir:
            print("Batch mode")
            print("Both counts and effective length provided. ")
            print("Transformation method:", args.norm)
            quant_dir = os.path.realpath(args.infile_dir)
            eff_len_dir = os.path.realpath(args.eff_len_dir)

            if not os.path.isdir(quant_dir):
                print(quant_dir, " is no directory")
                sys.exit(1)
            if not os.path.isdir(eff_len_dir):
                print(eff_len_dir, " is no directory")
                sys.exit(1)
            count_files = [f for f in os.listdir(quant_dir) if os.path.isfile(os.path.join(quant_dir, f))]
            len_files = [f for f in os.listdir(eff_len_dir) if os.path.isfile(os.path.join(eff_len_dir, f))]

        if args.infile_dir and not args.eff_len_dir:
            print("Batch mode")
            print("Only expression values provided.")
            print("Assuming normalized expression values (like log-fpkm or log-tpm).")
            quant_dir = os.path.realpath(args.infile_dir)

            if not os.path.isdir(quant_dir):
                print(quant_dir, " is no directory")
                sys.exit(1)

        if not args.infile_dir:
            print("No expression data provided. Please set Either --infile_dir or --infile_dir AND --eff_len_file_dir")
            sys.exit(1)

        # print("Trying to identify species in raw count directory ", quant_dir, " :")
        for f in count_files:
            split_fn = f.split("_")
            species = split_fn[0] + " " + split_fn[1]
            print("Species detected: ", species)
            sp = split_fn[0] + "_" + split_fn[1]
            len_file = [x for x in len_files if re.match(sp, x)]
            count_file = [x for x in count_files if re.match(sp, x)]
            print(count_file)
            if len(len_file) > 1:
                len_file = len_file[0]
            if len(count_file) > 1:
                count_file = len_file[0]
            path_metadata = os.path.realpath(args.metadata)
            export_string = " --export=QUANT=" + os.path.realpath(str(count_file)) + ",META=" + path_metadata + ",WORK=" + os.path.realpath(args.out_dir) + ",LEN=" + os.path.realpath(
                str(len_file)) + ",DIST=" + dist_method + ",CUT=" + str(mr_cut) + ",INTER=" + str(
                intermediate) + ",CURATE_GROUP=" + f'"{curate_group}"' + " " + batch_script_path + ",NORM=" + str(
                args.norm)
            # print(export_string)

            submit_command = ("sbatch " + "--job-name=" + sp + export_string)

            print(submit_command)  # Uncomment this line when done testing to use the submit command created
            # uncomment the following 3 lines when done testing to submit the jobs
            exit_status = subprocess.call(submit_command, shell=True)
            if exit_status != 0:  # Check to make sure the job submitted
                print("Job failed to submit".format(submit_command))
        print("Done submitting jobs!")
