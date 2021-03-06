import subprocess
import os
import sys
from os import listdir
from os.path import isfile, join


def cstmm_main(args):

    # Check if R is installed
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        print(",<ERROR> Rscript is not installed.")
        sys.exit(1)

    dir_count = os.path.join(args.work_dir, args.count)
    dir_ortho = os.path.join(args.work_dir, args.ortho)
    dir_work = args.work_dir
    count_files = [f for f in listdir(dir_count) if isfile(join(dir_count, f))]

    print("Trying to identify species in raw count directory ", dir_count, " :")
    for f in count_files:
        split_fn = f.split("_")
        species = split_fn[0]+" "+split_fn[1]
        print("Species detected: ", species)

    # Check if all files are there
    # TODO: do this more elegantly. :D
    if os.path.exists(os.path.join(dir_ortho, 'Orthogroups.tsv')):
        print("Orthogroups.tsv found")
    else:
        print("Could not find Orthogroups.tsv. Did you provide the correct OrthoFinder folder?")
        sys.exit()

    if os.path.exists(os.path.join(dir_ortho, 'Orthogroups_SingleCopyOrthologues.txt')):
        print("Orthogroups_SingleCopyOrthologues.txt found")
    else:
        print("Could not find Orthogroups_SingleCopyOrthologues.txt. Did you provide the correct OrthoFinder path?")
        sys.exit()

    cstmm_path = os.path.dirname(os.path.realpath(__file__))
    r_script_path = cstmm_path + '/cstmm.r'
    subprocess.check_call(['Rscript', r_script_path, os.path.realpath(dir_count), os.path.realpath(dir_ortho), os.path.realpath(dir_work)])
