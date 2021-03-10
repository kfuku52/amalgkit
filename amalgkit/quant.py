import re, glob, subprocess, os, sys
from Bio import Entrez
from amalgkit.getfastq import getfastq_getxml, getfastq_search_term
from amalgkit.util import *

def is_quant_output_present(sra_id, output_dir):
    out_path = os.path.join(output_dir, sra_id+'_abundance.tsv')
    is_output_present = os.path.exists(out_path)
    return is_output_present

def quant_main(args):

    #check kallisto dependency
    try:
        subprocess.run(['kallisto', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(",<ERROR> kallisto is not installed.")
        sys.exit(1)

    assert (args.id is None)!=(args.metadata is None), 'Either --id or --metadata should be specified.'
    if args.id is not None:
        print('--id is specified.')
        sra_id = args.id
    if args.metadata is not None:
        print('--metadata is specified. Reading existing metadata table.')
        assert (args.batch is not None), '--batch should be specified.'
        metadata = load_metadata(args)
        sra_id = metadata.df.loc[:,'run'].values[0]
    print('SRA ID:', sra_id)

    if args.index is None:
        index = args.ref + ".idx"
    else:
        index = args.index


    # build index via kallisto index

    if args.build_index == "yes":
        if not args.ref:
            raise ValueError("--build_index enabled, but no reference sequences given.")
        else:
            if os.path.exists(index):
                print('kallisto index was detected:', index)
            else:
                print('kallisto index was not detected. Creating:', index)
                kallisto_out = subprocess.run(["kallisto", "index", "--i", index, args.ref])
                # TODO: Switch to try/except for error handling
                assert (kallisto_out.returncode == 0), "kallisto did not finish safely: {}".format(kallisto_out.stdout.decode('utf8'))

    # prefer amalgkit processed files over others.

    in_files = glob.glob(os.path.join(args.work_dir, 'getfastq', sra_id, sra_id + "*.amalgkit.fastq.gz"))
    if not in_files:
        in_files = glob.glob(os.path.join(args.work_dir, sra_id) + "*.fastq*")

    # make results directory, if not already there
    output_dir = os.path.join(args.work_dir, 'quant', sra_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if (is_quant_output_present(sra_id, output_dir))&(args.redo=='no'):
        print('Output file(s) detected. Exiting. Set "--redo yes" for reanalysis.')
        sys.exit()

    # start quantification process.
    # throws exception, if in_files still empty.
    assert in_files, '{}: Fastq file not found. Check {}'.format(sra_id, output_dir)
    print('Input fastq detected:', ', '.join(in_files))
    # paired end read kallisto quant, if in_files == 2
    if len(in_files) == 2:
        print("paired end reads detected. Running in paired read mode.")
        kallisto_out = subprocess.run(["kallisto", "quant", "-i", index, "-o", output_dir, in_files[0], in_files[1]])
        # TODO: Switch to try/except for error handling
        assert (kallisto_out.returncode == 0), "kallisto did not finish safely: {}".format(
            kallisto_out.stdout.decode('utf8'))

    # throws exception, if more than 2 files in in_files. Could be expanded to handle more than 2 files.
    elif len(in_files) > 2:
        raise ValueError("more than 2 input files given. Refer to kallisto quant -h for more info")

    else:
        print("single end reads detected. Proceeding in --single mode")

        # kallisto needs fragment length and fragment length standard deviation.
        # if there is none supplied, check if ID matches SRA ID
        if args.fragment_length is None:

            print("No fragment length set.")
            SRA_ID_pattern= re.compile("[DES]R[RXP][0-9]{7}")
            # if it does match an SRA ID pattern, try to fetch fragment length from metadata.
            # this uses getfastq_getxml and getfastq_search_term from getfastq, as well as Metadata from metadata
            if SRA_ID_pattern.match(sra_id):
                print("SRA-ID detected. Trying to fetch fragment length from metadata.")
                Entrez.email = "test@test.com"
                search_term = getfastq_search_term(sra_id)
                print('Entrez search term:', search_term)
                xml_root = getfastq_getxml(search_term)
                metadata = Metadata.from_xml(xml_root)
                metadata.df = metadata.df.loc[(metadata.df['lib_layout'] == "single"), :]
                nominal_length = metadata.df['nominal_length'][0]
                # set nominal length to 200 if below 200 ...
                if nominal_length:
                    if nominal_length < 200:
                        nominal_length = 200
                # ...or if undefined.
                else:
                    nominal_length = 200
            # if sra_id isn't in SRA ID format, fragment length of 200 is assumed
            else:
                print("No SRA-ID detected. Assuming fragment length.")
                nominal_length = 200

            print("fragment length set to: ", nominal_length)
            fragment_sd = nominal_length/10
            print("fragment length standard deviation set to:", fragment_sd)
            kallisto_out = subprocess.run(["kallisto", "quant", "--index", index, "-o", output_dir, "--single", "-l", str(nominal_length), "-s", str(fragment_sd), in_files[0]])
            # TODO: Switch to try/except for error handling
            assert (kallisto_out.returncode == 0), "kallisto did not finish safely: {}".format(
                kallisto_out.stdout.decode('utf8'))
        # if fragment length is supplied by the user, kallisto quant can be run immediately
        else:
            print("fragment length set to: ", args.fragment_length)
            fragment_sd = args.fragment_length/10
            print("fragment length standard deviation set to:", fragment_sd)
            kallisto_out = subprocess.run(["kallisto", "quant", "--index", index,  "-o", output_dir, "--single", "-l", str(args.frament_length), "-s", str(fragment_sd), in_files[0]])
            # TODO: Switch to try/except for error handling
            assert (kallisto_out.returncode == 0), "kallisto did not finish safely: {}".format(
                kallisto_out.stdout.decode('utf8'))
    # move output to results with unique name
    os.rename(os.path.join(output_dir, "run_info.json"), os.path.join(output_dir, sra_id + "_run_info.json"))
    os.rename(os.path.join(output_dir, "abundance.tsv"), os.path.join(output_dir, sra_id + "_abundance.tsv"))
    os.rename(os.path.join(output_dir, "abundance.h5"), os.path.join(output_dir, sra_id + "_abundance.h5"))

    if (args.clean_fastq=='yes')&(os.path.exists(os.path.join(output_dir, sra_id + "_abundance.tsv"))):
        for in_file in in_files:
            print('Output file detected. Safely removing fastq:', in_file)
            os.remove(in_file)
            placeholder = open(in_file+'.safely_removed', "w")
            placeholder.write("This fastq file was safely removed after `amalgkit quant`.")
            placeholder.close()
