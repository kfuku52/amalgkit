import subprocess
from amalgkit.util import *

def check_quant_output(sra_id, output_dir, args):
    out_path = os.path.join(output_dir, sra_id + '_abundance.tsv')
    is_output_present = os.path.exists(out_path)
    if is_output_present:
        print('Output file detected: {}'.format(out_path))
        if args.redo == 'yes':
            print('Continued. The output will be overwritten. Set "--redo no" to not overwrite results.')
            return None
        else:
            print('Continued. The output will not be overwritten. If you want to overwrite the results, set "--redo yes".')
    else:
        print('Output file was not detected: {}'.format(out_path))
        return None


def call_kallisto(args, in_files, metadata, sra_stat, output_dir, index):
    sra_id = sra_stat['sra_id']
    lib_layout = sra_stat['layout']
    kallisto_cmd = ''
    if lib_layout == 'single':
        print("single end reads detected. Proceeding in single mode")
        if len(in_files) != 1:
            txt = "Library layout: {} and expected 1 input file. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        nominal_length = metadata.df.loc[:, 'nominal_length'].values[0]
        if nominal_length:
            print('Nominal length in metadata is unusually small ({}). Setting it to 200.'.format(nominal_length))
            if nominal_length < 200 or numpy.isnan(nominal_length):
                nominal_length = 200
        else:
            print("Could not find nominal length in metadata. Assuming fragment length.")
            nominal_length = 200
        print("Fragment length set to: {}".format(nominal_length))
        fragment_sd = nominal_length / 10
        print("Fragment length standard deviation set to: {}".format(fragment_sd))
        kallisto_cmd = ["kallisto", "quant", "--threads", str(args.threads), "--index", index, "-o", output_dir,
                        "--single", "-l", str(nominal_length), "-s", str(fragment_sd), in_files[0]]
    elif lib_layout == 'paired':
        if len(in_files) != 2:
            txt = "Library layout: {} and expected 2 input files. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        print("paired end reads detected. Running in paired read mode.")
        kallisto_cmd = ["kallisto", "quant", "--threads", str(args.threads), "-i", index, "-o",
                        output_dir, in_files[0], in_files[1]]

    print('Command: {}'.format(' '.join(kallisto_cmd)))
    kallisto_out = subprocess.run(kallisto_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print('kallisto stdout:')
    print(kallisto_out.stdout.decode('utf8'))
    print('kallisto stderr:')
    print(kallisto_out.stderr.decode('utf8'))
    if kallisto_out.returncode != 0:
        sys.stderr.write("kallisto did not finish safely.\n")
        if 'zero reads pseudoaligned' in kallisto_out.stderr.decode('utf8'):
            sys.stderr.write('No reads are mapped to the reference. This sample will be removed by `amalgkit curate`.')

    # move output to results with unique name
    try:
        os.rename(os.path.join(output_dir, "run_info.json"), os.path.join(output_dir, sra_id + "_run_info.json"))
        os.rename(os.path.join(output_dir, "abundance.tsv"), os.path.join(output_dir, sra_id + "_abundance.tsv"))
        os.rename(os.path.join(output_dir, "abundance.h5"), os.path.join(output_dir, sra_id + "_abundance.h5"))
    except FileNotFoundError:
        pass

    return kallisto_out


def check_kallisto_dependency():
    try:
        subprocess.run(['kallisto', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        raise Exception("kallisto is not installed.")


def check_layout_mismatch(sra_stat, output_dir):
    if sra_stat['layout'] == 'paired':
        fastq_files = glob.glob(os.path.join(output_dir, sra_stat['sra_id'] + '*.fastq*'))
        if len(fastq_files) == 1:
            sys.stderr.write('Single-end fastq was detected even though layout = {}\n'.format(sra_stat['layout']))
            sys.stderr.write('This sample will be treated as single-end sequencing.\n')
            sra_stat['layout'] = 'single'
    return sra_stat


def run_quant(args, metadata, sra_id, index):
    # make results directory, if not already there
    output_dir = os.path.join(args.out_dir, 'quant', sra_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    check_quant_output(sra_id, output_dir, args)

    output_dir_getfastq = os.path.join(args.out_dir, 'getfastq', sra_id)
    sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra=None)
    sra_stat = check_layout_mismatch(sra_stat, output_dir_getfastq)
    try:
        ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir_getfastq)
        if ext == '.safely_removed':
            print('These files have been deleted. If you wish to reobtain the .fastq file(s), run: getfastq -e email@adress.com --id ', sra_id, ' -w ', args.out_dir, '--redo yes --gcp yes --aws yes --ncbi yes')
            print('skipping.')
            return

    except FileNotFoundError:
        print('ERROR: could not find fastq file(s) in:', output_dir_getfastq)
        print(
            'If you wish to obtain the .fastq file(s), run: getfastq -e email@adress.com --id ',
            sra_id, ' -w ', args.out_dir, '--redo yes --gcp yes --aws yes --ncbi yes')
        print('skipping.')
        return
    in_files = glob.glob(os.path.join(args.out_dir, 'getfastq', sra_id, sra_id + "*" + ext))

    # start quantification process.
    # throws exception, if in_files still empty.
    assert in_files, '{}: Fastq file not found. Check {}'.format(sra_id, output_dir)
    print('Input fastq detected:', ', '.join(in_files))

    call_kallisto(args, in_files, metadata, sra_stat, output_dir, index)

    if (args.clean_fastq == 'yes') & (os.path.exists(os.path.join(output_dir, sra_id + "_abundance.tsv"))):
        for in_file in in_files:
            print('Output file detected. Safely removing fastq:', in_file)
            os.remove(in_file)
            placeholder = open(in_file + '.safely_removed', "w")
            placeholder.write("This fastq file was safely removed after `amalgkit quant`.")
            placeholder.close()


def get_index(args, sci_name):
    if args.index_dir is not None:
        index_dir = args.index_dir
    else:
        index_dir = os.path.join(args.out_dir, 'Index')
    if os.path.exists(index_dir):
        index = glob.glob(os.path.join(index_dir, sci_name + '*'))
        if len(index) > 1:
            raise ValueError(
                "found multiple index files for species. Please make sure there is only one index file for this species.")
        elif len(index) == 0:
            raise FileNotFoundError("Could not find Index file.")
        index = index[0]
    else:
        raise FileNotFoundError("could not find index folder")
    print("Index file found: {}".format(index))
    return index


def quant_main(args):
    check_kallisto_dependency()
    metadata = load_metadata(args)  # loads single-row metadata according to --batch
    if args.id is None:
        for i in metadata.df.index:
            print('')
            sra_id = metadata.df.at[i, 'run']
            sci_name = metadata.df.at[i, 'scientific_name']
            print('Species: {}'.format(sci_name))
            print('Run ID: {}'.format(sra_id))
            sci_name = sci_name.replace(" ", "_")
            print('looking for index folder in ', args.out_dir)
            index = get_index(args, sci_name)
            run_quant(args, metadata, sra_id, index)
    # if args.id is not specified, it will run the whole metadata sheet one by one
    else:
        sra_id = args.id
        sci_name = metadata.df.loc[metadata.df['run'] == args.id, 'scientific_name']
        print('Species: {}'.format(sci_name))
        print('Run ID: {}'.format(sra_id))
        sci_name = sci_name.replace(" ", "_")
        index = get_index(args, sci_name)
        run_quant(args, metadata, sra_id, index)
