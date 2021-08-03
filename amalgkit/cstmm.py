import subprocess
import os
import sys

def check_cstmm_dependency():
    try:
        subprocess.run(['Rscript', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as e:
        print(e)
        raise Exception("R (Rscript) is not installed.")

def get_count_files(dir_count):
    sciname_dirs = os.listdir(dir_count)
    count_files = list()
    for sciname_dir in sciname_dirs:
        sciname_path = os.path.join(dir_count, sciname_dir)
        if not os.path.isdir(sciname_path):
            continue
        files = os.listdir(sciname_path)
        for file in files:
            sciname_file = os.path.join(sciname_path, file)
            if not os.path.isfile(sciname_file):
                continue
            if not file.endswith('est_counts.tsv'):
                continue
            count_files.append(sciname_file)
        sciname_count_files = [ f for f in count_files if '/'+sciname_dir+'/' in f ]
        num_sciname_count_file = len(sciname_count_files)
        if (num_sciname_count_file==0):
            sys.stderr.write('No est_counts.tsv file found in: {}\n'.format(sciname_path))
        elif (num_sciname_count_file==1):
            print('One est_counts.tsv file found in {}'.format(sciname_path), flush=True)
        elif (num_sciname_count_file>=2):
            raise Exception('Multiple est_counts.tsv files found in: {}\n'.format(sciname_path))
    if (len(count_files)==0):
        raise Exception('No est_counts.tsv file was detected.')
    return count_files

def check_orthofinder_outputs(dir_ortho):
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

def cstmm_main(args):
    check_cstmm_dependency()

    # TODO still truncated PATHs are expected. Please update.
    dir_count = os.path.realpath(os.path.join(args.out_dir, args.count))
    dir_ortho = os.path.realpath(os.path.join(args.out_dir, args.ortho))
    dir_work = os.path.realpath(args.out_dir)
    count_files = get_count_files(dir_count)
    if (len(count_files)==1):
        txt = 'Only one species was detected. Standard TMM normalization will be applied.'
        print(txt, flush=True)
        mode_tmm = 'single_species'
    else:
        txt = 'Multiple species were detected. ' \
              'Cross-species TMM normalization will be applied with single-copy orthologs.'
        print(txt, flush=True)
        mode_tmm = 'multi_species'
        check_orthofinder_outputs(dir_ortho)
    cstmm_path = os.path.dirname(os.path.realpath(__file__))
    r_script_path = cstmm_path + '/cstmm.r'
    r_command = ['Rscript', r_script_path, dir_count, dir_ortho, dir_work, mode_tmm]
    print('')
    print('Starting R script: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)
