import glob, subprocess,sys, numpy
from amalgkit.util import *

def check_quant_output(sra_id, output_dir, args):
    out_path = os.path.join(output_dir, sra_id+'_abundance.tsv')
    is_output_present = os.path.exists(out_path)
    if is_output_present:
        print('Output file detected: {}'.format(out_path))
        if (args.redo=='yes'):
            print('Continued. The output will be overwritten. Set "--redo no" to exit.')
            return None
        else:
            print('Exiting. Set "--redo yes" to overwrite.')
            sys.exit()
    else:
        print('Output file was not detected: {}'.format(out_path))
        return None

def call_kallisto(args,in_files,metadata,sra_id, output_dir, index):

    lib_layout = metadata.df.loc[(metadata.df.loc[:,'run'] == sra_id), 'lib_layout']
    lib_layout = lib_layout.iloc[0]
    print(lib_layout)
    if lib_layout == 'single':
        print("single end reads detected. Proceeding in single mode")
        if len(in_files) == 1:
            nominal_length = metadata.df.loc[:,'nominal_length'].values[0]
            # set nominal length to 200 if below 200 ...
            if nominal_length:
                if nominal_length < 200 or numpy.isnan(nominal_length):
                    nominal_length = 200
            # ...or if undefined.
            else:
                print("Could not find nominal length in metadata. Assuming fragment length.")
                nominal_length = 200
            print("fragment length set to: ", nominal_length)
            fragment_sd = nominal_length / 10
            print("fragment length standard deviation set to:", fragment_sd)

            kallist_cmd = ["kallisto", "quant", "--threads", str(args.threads), "--index", index, "-o", output_dir,
                           "--single", "-l", str(nominal_length), "-s", str(fragment_sd), in_files[0]]
            kallisto_out = subprocess.run(kallist_cmd)
        else:
            txt = "Library layout: {} and expected 1 input file. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
    elif lib_layout == 'paired':
        if len(in_files) == 2:
            print("paired end reads detected. Running in paired read mode.")
            kallisto_out = subprocess.run( ["kallisto", "quant", "--threads", str(args.threads), "-i", index, "-o",
                                            output_dir, in_files[0], in_files[1]])
        else:
            txt = "Library layout: {} and expected 2 input files. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))

    # TODO: Switch to try/except for error handling
    assert (kallisto_out.returncode == 0), "kallisto did not finish safely: {}".format(kallisto_out.stdout.decode('utf8'))

    # move output to results with unique name
    os.rename(os.path.join(output_dir, "run_info.json"), os.path.join(output_dir, sra_id + "_run_info.json"))
    os.rename(os.path.join(output_dir, "abundance.tsv"), os.path.join(output_dir, sra_id + "_abundance.tsv"))
    os.rename(os.path.join(output_dir, "abundance.h5"), os.path.join(output_dir, sra_id + "_abundance.h5"))

    return kallisto_out

def quant_main(args):

    #check kallisto dependency
    try:
        subprocess.run(['kallisto', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(",<ERROR> kallisto is not installed.")
        sys.exit(1)
    if args.id is not None:
        print('--id is specified.')
        sra_id = args.id
    if args.metadata is not None:
        print('--metadata is specified. Reading existing metadata table.')
        metadata = load_metadata(args) # loads single-row metadata according to --batch
        if args.batch is not None:
            sra_id = metadata.df.loc[:,'run'].values[0]
        else:
            sra_id = args.id
    print('SRA Run ID: {}'.format(sra_id))

    index = args.index

    # make results directory, if not already there
    output_dir = os.path.join(args.out_dir, 'quant', sra_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    check_quant_output(sra_id, output_dir, args)

    sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra=None)
    output_dir_getfastq = os.path.join(args.out_dir, 'getfastq', sra_id)
    ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir_getfastq)
    in_files = glob.glob(os.path.join(args.out_dir, 'getfastq', sra_id, sra_id + "*" + ext))

    # start quantification process.
    # throws exception, if in_files still empty.
    assert in_files, '{}: Fastq file not found. Check {}'.format(sra_id, output_dir)
    print('Input fastq detected:', ', '.join(in_files))

    call_kallisto(args, in_files, metadata,sra_id,output_dir,index)

    if (args.clean_fastq=='yes')&(os.path.exists(os.path.join(output_dir, sra_id + "_abundance.tsv"))):
        for in_file in in_files:
            print('Output file detected. Safely removing fastq:', in_file)
            os.remove(in_file)
            placeholder = open(in_file+'.safely_removed', "w")
            placeholder.write("This fastq file was safely removed after `amalgkit quant`.")
            placeholder.close()
