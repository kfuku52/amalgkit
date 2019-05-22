import re, glob, subprocess, os, sys
from amalgkit.metadata import create_run_dir

def curate_main(args):

    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(",<ERROR> Rscript is not installed.")
        sys.exit(1)

    quant_out = os.path.join(args.work_dir, args.infile)
    meta_out = os.path.join(args.work_dir,args.metafile)
    out_dir = create_run_dir(os.path.join(args.out_dir, 'getfastq_output'))
    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    intermediate = args.cleanup
    tissues = re.findall(r"[\w]+", args.tissues)
    tissues = '|'.join(tissues)
    subprocess.check_call(["Rscript", 'util/transcriptome_curation.r', quant_out, meta_out, out_dir, dist_method, '0', str(mr_cut), str(intermediate), tissues])

