import re
import subprocess
import warnings

from amalgkit.util import *


def get_curate_group(args, metadata):
    if args.curate_group is None:
        curate_group = metadata.df.loc[:, 'curate_group'].dropna().unique()
    else:
        curate_group = re.findall(r"[\w]+", args.curate_group)
    if (len(curate_group)==0):
        txt = 'The "curate_group" column in --metadata ({}) is not filled. Please manually edit the file.\n'
        txt += '`amalgkit curate` recognizes samples with the same string in this columns to belong the same group.\n'
        txt += 'Exiting.\n'
        sys.stderr.write(txt.format(args.metadata))
        sys.exit(1)
    print('Tissues to be included: {}'.format(', '.join(curate_group)))
    curate_group = '|'.join(curate_group)
    return curate_group

def check_rscript():
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print(",<ERROR> Rscript is not installed.")
        sys.exit(1)


def run_curate_r_script(args, new_metadata_path, metadata, sp):
    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    correlation_threshold = args.correlation_threshold
    intermediate = args.plot_intermediate
    curate_group = get_curate_group(args, metadata)
    curate_path = os.path.dirname(os.path.realpath(__file__))
    r_script_path = curate_path + '/transcriptome_curation.r'
    if args.input_dir is None:
        input_dir = os.path.join(args.out_dir, 'merge')
        print("no input_dir given. Assuming path is: ", input_dir)
    else:
        input_dir = args.input_dir

    # check if cstmm output is used when --norm == tpm, because TPM undoes tmm normalization
    if args.norm is not None and 'tpm' in args.norm:
        substring = "cstmm"
        try:
            input_dir.index(substring)
        except ValueError:
            warnings.warn(
                "WARNING: TPM NORMALIZATION AND TMM NORMALIZATION ARE INCOMPATIBLE. IF INPUT DATA IS TMM NORMALIZED, PLEASE SWITCH --norm TO ANY OF THE 'fpkm' NORMALISATIONS INSTEAD: (logn|log2|lognp1|log2p1|none)-(fpkm)")
        else:
            raise ValueError(
                "ERROR: AMALGKIT CSTMM NORMALIZED INPUT FILES DETECTED WHILE NORMALIZATION METHOD IS 'TPM'. TMM NORMALIZATION AND TPM NORMALIZATION ARE INCOMPATIBLE! PLEASE SWITCH --norm TO ANY OF THE 'fpkm' NORMALISATIONS INSTEAD: (logn|log2|lognp1|log2p1|none)-(fpkm)")

    len_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_eff_length.tsv')
    if 'cstmm' in input_dir:
        count_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_cstmm_counts.tsv')
    else:
        count_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_est_counts.tsv')

    if os.path.exists(count_file) and os.path.exists(len_file):
        print("Both counts and effective length files found.")
    else:
        sys.stderr.write("No expression data found. Please make sure `amalgkit merge` or `amalgkit cstmm` ran correctly and you provided the correct directory PATH.\n")
        if not os.path.exists(count_file):
            sys.stderr.write('Expected but undetected PATH of the count file: {}\n'.format(count_file))
        if not os.path.exists(len_file):
            sys.stderr.write('Expected but undetected PATH of the effective length file: {}\n'.format(len_file))
        sys.exit(1)
    print("Starting Rscript to obtain curated {} values.".format(args.norm))
    subprocess.call(['Rscript',
                     r_script_path,
                     count_file,
                     new_metadata_path,
                     os.path.realpath(args.out_dir),
                     len_file,
                     dist_method,
                     str(mr_cut),
                     '0',
                     str(intermediate),
                     curate_group,
                     str(args.norm),
                     str(args.one_outlier_per_iter),
                     str(correlation_threshold),
                     ])
    return


def curate_main(args):
    check_rscript()
    metadata = load_metadata(args)
    # figure out unique species
    spp = metadata.df.loc[:, 'scientific_name'].drop_duplicates().values
    # make directory
    curate_dir = os.path.join(args.out_dir, 'curate')
    if not os.path.exists(curate_dir):
        os.mkdir(curate_dir)
    new_metadata_path = os.path.realpath(os.path.join(curate_dir, 'metadata.tsv'))

    print('Found a total number of ', len(spp), ' species in this metadata table:')
    print('____________________________')
    for sp in spp:
        print(sp)
    print('____________________________')

    if args.batch is None:
        write_updated_metadata(metadata, new_metadata_path, args)
        # enter "normal" mode. Process all species, 1 at a time
        for sp in spp:
            sp = sp.replace(" ", "_")
            run_curate_r_script(args, new_metadata_path, metadata, sp)
    else:
        # enter Batch mode, only process 1 species
        print('Entering --batch mode. processing 1 species')
        txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table.'
        print(txt.format(args.batch, len(spp)))
        sp = spp[args.batch - 1]
        print('processing species number ', args.batch, ' : ', sp)
        metadata.df = metadata.df.loc[metadata.df['scientific_name'] == sp]
        sp = sp.replace(" ", "_")
        write_updated_metadata(metadata, new_metadata_path, args)
        run_curate_r_script(args, new_metadata_path, metadata, sp)
