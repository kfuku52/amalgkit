import re, glob, subprocess, os, sys
from Bio import Entrez
from amalgkit.metadata import Metadata
from amalgkit.getfastq import getfastq_getxml, getfastq_search_term
from amalgkit.metadata import create_run_dir

def quant_main(args):

    #check kallisto depency
    try:
        subprocess.run(['kallisto', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(",<ERROR> kallisto is not installed.")
        sys.exit(1)

    if args.index is None:
        index = args.id + ".idx"
    else:
        index = args.index

    # build index via kallisto index

    if args.build_index == "yes":
        if not args.ref:
            raise ValueError("--build_index enabled, but no reference sequences given.")
        else:
            subprocess.run(["kallisto", "index", "--i", index, args.work_dir + args.ref])

    # prefer amalgkit processed files over others.

    in_files = glob.glob(os.path.join(args.work_dir, args.id + "*.amalgkit.fastq.gz"))

    if not in_files:
        in_files = glob.glob(os.path.join(args.work_dir, args.id) + "*.fastq*")


    # make results directory, if not already there
    output_dir = ''
    if args.auto_dir == 'yes':
        if not os.path.exists(os.path.join(args.out_dir, 'quant_output')):
            output_dir = create_run_dir(os.path.join(args.out_dir, 'quant_output'+args.id))
    else:
        if not os.path.exists(os.path.join(args.out_dir)):
            output_dir = os.path.join(args.out_dir)

    # start quantification process.
    # throws exception, if in_files still empty.

    if in_files:
        # paired end read kallisto quant, if in_files == 2
        if len(in_files) == 2:
            print("paired end reads detected. Running in paired read mode.")
            subprocess.run(["kallisto", "quant", "-i", index, "-o", output_dir, in_files[0], in_files[1]])

        # throws exception, if more than 2 files in in_files. Could be expanded to handle more than 2 files.
        elif len(in_files) > 2:
            raise ValueError("more than 2 input files given. Refer to kallisto quant -h for more info")

        else:
            print("single end reads detected. Proceeding in --single mode")

            # kallisto needs fragment length and fragment length standard deviation.
            # if there is none supplied, check if ID matches SRA ID
            if args.fragment_length is None:

                print("No fragment length set.")
                SRA_ID_pattern= re.compile("SR[RXP][0-9]{7}")
                # if it does match an SRA ID pattern, try to fetch fragment length from metadata.
                # this uses getfastq_getxml and getfastq_search_term from getfastq, as well as Metadata from metadata
                if SRA_ID_pattern.match(args.id):
                    print("SRA-ID detected. Trying to fetch fragment length from metadata.")
                    Entrez.email = "test@test.com"
                    search_term = getfastq_search_term(args.id)
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
                # if args.id isn't in SRA ID format, fragment length of 200 is assumed
                else:
                    print("No SRA-ID detected. Assuming fragment length.")
                    nominal_length = 200

                print("fragment length set to: ", nominal_length)
                fragment_sd = nominal_length/10
                print("fragment length standard deviation set to:", fragment_sd)
                subprocess.run(["kallisto", "quant", "--index", index, "-o", output_dir, "--single", "-l", str(nominal_length), "-s", str(fragment_sd), in_files[0]])

            # if fragment length is supplied by the user, kallisto quant can be run immediately
            else:
                print("fragment length set to: ", args.fragment_length)
                fragment_sd = args.fragment_length/10
                print("fragment length standard deviation set to:", fragment_sd)
                subprocess.run(["kallisto", "quant", "--index", index,  "-o", output_dir, "--single", "-l", str(args.frament_length), "-s", str(fragment_sd), in_files[0]])
    else:
        raise ValueError("ID ", args.id, "not found in working directory", args.work_dir)



    # move output to results with unique name
    os.rename(os.path.join(output_dir, "run_info.json"), os.path.join(output_dir, args.id + "_run_info.json"))
    os.rename(os.path.join(output_dir, "abundance.tsv"), os.path.join(output_dir, args.id + "_abundance.tsv"))
    os.rename(os.path.join(output_dir, "abundance.h5"), os.path.join(output_dir, args.id + "_abundance.h5"))
