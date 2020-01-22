import re
import subprocess
import os
import sys
from amalgkit.metadata import create_run_dir


def curate_main(args):

    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print(",<ERROR> Rscript is not installed.")
        sys.exit(1)

    quant_out = os.path.join(args.work_dir, args.infile)
    meta_out = os.path.join(args.work_dir, args.metafile)
    out_dir = ''
    if args.auto_dir == 'yes':
        if not os.path.exists(os.path.join(args.out_dir, 'curate_output')):
            out_dir = create_run_dir(os.path.join(args.out_dir, 'curate_output'))
    else:
        if not os.path.exists(os.path.join(args.out_dir)):
            out_dir = os.makedirs(os.path.join(args.out_dir))
        else:
            out_dir = os.path.join(args.out_dir)
    #out_dir = create_run_dir(os.path.join(args.out_dir, 'curate_output'))
    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    intermediate = args.cleanup
    tissues = re.findall(r"[\w]+", args.tissues)
    tissues = '|'.join(tissues)
    curate_path = os.path.dirname(os.path.realpath(__file__))
    r_script_path = curate_path+'/transcriptome_curation.r'
    subprocess.check_call(['Rscript', r_script_path, os.path.realpath(quant_out), os.path.realpath(meta_out), os.path.realpath(out_dir), dist_method, '0', str(mr_cut), str(intermediate), tissues])
