from amalgkit.util import *

import re
import subprocess
import os
import sys

def get_tissues(args):
    if args.tissues is None:
        metadata = load_metadata(args)
        tissues = metadata.df.loc[:,'tissue'].dropna().unique()
    else:
        tissues = re.findall(r"[\w]+", args.tissues)
    tissues = '|'.join(tissues)
    print('Tissues to be included: {}'.format(', '.join(tissues)))
    return tissues

def curate_main(args):

    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print(",<ERROR> Rscript is not installed.")
        sys.exit(1)

    meta_out = os.path.join(args.work_dir, args.metadata)

    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    intermediate = args.cleanup
    tissues = get_tissues(args)
    curate_path = os.path.dirname(os.path.realpath(__file__))
    r_script_path = curate_path + '/transcriptome_curation.r'
    batch_script_path = curate_path + '/batch_curate.sh'

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
            quant_out = os.path.join(args.work_dir, args.infile)
            eff_len_out = os.path.join(args.work_dir, args.eff_len_file)
        if args.infile and not args.eff_len_file:
            print("Single species mode")
            print("Only expression values provided. ")
            print("Assuming normalized expression values (like fpkm or tpm).")
            quant_out = os.path.join(args.work_dir, args.infile)
            eff_len_out = "NA"
        if not args.infile:
            print("No expression data provided. Please set Either --infile or --infile AND --eff_len_file")
            sys.exit(1)

        subprocess.check_call(['Rscript',
                               r_script_path,
                               os.path.realpath(quant_out),
                               os.path.realpath(meta_out),
                               os.path.realpath(args.work_dir),
                               os.path.realpath(eff_len_out),
                               dist_method,
                               '0',
                               str(mr_cut),
                               str(intermediate),
                               tissues,
                               str(args.norm)])
    # if multiple species mode active
    if args.batch is not None:

        # Here, we still need arguments for counts and effective lengths,
        # But they need to be directories
        if args.infile_dir and args.eff_len_dir:
            print("Batch mode")
            print("Both counts and effective length provided. ")
            print("Calculating: ", args.norm)
            quant_dir = os.path.join(args.work_dir, args.infile_dir)
            eff_len_dir = os.path.join(args.work_dir, args.eff_len_dir)

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
            quant_dir = os.path.join(args.work_dir, args.infile_dir)
            eff_len_dir = "NA"
            if not os.path.isdir(quant_dir):
                print(quant_dir, " is no directory")
                sys.exit(1)

        if not args.infile_dir:
            print("No expression data provided. Please set Either --infile_dir or --infile_dir AND --eff_len_file_dir")
            sys.exit(1)


        #print("Trying to identify species in raw count directory ", quant_dir, " :")
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
            export_string=" --export=QUANT=" + os.path.realpath(str(count_file)) + ",META=" + os.path.realpath(meta_out) + ",WORK=" + os.path.realpath(args.work_dir) + ",LEN=" + os.path.realpath(str(len_file)) + ",DIST=" + dist_method + ",CUT=" + str(mr_cut) + ",INTER=" + str(intermediate) + ",TISSUES=" +f'"{tissues}"' + " " + batch_script_path + ",NORM=" +str(args.norm)
            #print(export_string)

            submit_command = ("sbatch " + "--job-name=" + sp +export_string)

            print(submit_command)  # Uncomment this line when done testing to use the submit command created
            # uncomment the following 3 lines when done testing to submit the jobs
            exit_status = subprocess.call(submit_command, shell=True)
            if exit_status != 0:  # Check to make sure the job submitted
                print("Job failed to submit".format(submit_command))
        print("Done submitting jobs!")
