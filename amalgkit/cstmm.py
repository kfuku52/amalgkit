import subprocess
import os
import sys

from amalgkit.util import *

def resolve_species_count_file(sciname_path):
    sciname_count_files = []
    with os.scandir(sciname_path) as entries:
        for entry in entries:
            if (not entry.is_file()) or (not entry.name.endswith('est_counts.tsv')):
                continue
            sciname_count_files.append(entry.path)
    num_sciname_count_file = len(sciname_count_files)
    if num_sciname_count_file == 0:
        sys.stderr.write('No est_counts.tsv file found in: {}\n'.format(sciname_path))
        return None
    if num_sciname_count_file == 1:
        return sciname_count_files[0]
    raise ValueError('Multiple est_counts.tsv files found in: {}\n'.format(sciname_path))


def get_count_files(dir_count):
    if not os.path.exists(dir_count):
        raise FileNotFoundError('Count directory not found: {}'.format(dir_count))
    if not os.path.isdir(dir_count):
        raise NotADirectoryError('Count path exists but is not a directory: {}'.format(dir_count))
    count_files = []
    missing_species = []
    with os.scandir(dir_count) as species_entries:
        for species_entry in sorted(species_entries, key=lambda e: e.name):
            if (not species_entry.is_dir()) or species_entry.name.startswith('.') or species_entry.name.startswith('tmp.'):
                continue
            count_file = resolve_species_count_file(species_entry.path)
            if count_file is not None:
                count_files.append(count_file)
            else:
                missing_species.append(species_entry.name)
    if len(missing_species) > 0:
        raise FileNotFoundError(
            'No est_counts.tsv file found for species directory(ies): {} (under {}).'.format(
                ', '.join(missing_species),
                dir_count,
            )
        )
    if len(count_files) == 0:
        raise FileNotFoundError('No est_counts.tsv file was detected in: {}'.format(dir_count))
    return sorted(count_files)

def filepath2spp(file_paths):
    return [os.path.basename(path).split('_est_counts.tsv', 1)[0] for path in file_paths]

def cstmm_main(args):
    dir_out = os.path.realpath(args.out_dir)
    if os.path.exists(dir_out) and (not os.path.isdir(dir_out)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(dir_out))
    dir_cstmm = os.path.join(dir_out, 'cstmm')
    if os.path.exists(dir_cstmm) and (not os.path.isdir(dir_cstmm)):
        raise NotADirectoryError('cstmm path exists but is not a directory: {}'.format(dir_cstmm))
    check_rscript()
    os.makedirs(dir_cstmm, exist_ok=True)
    if args.dir_count=='inferred':
        dir_count = os.path.join(dir_out, 'merge')
    else:
        dir_count = os.path.realpath(args.dir_count)
    count_files = get_count_files(dir_count)
    if len(count_files) == 1:
        txt = 'Only one species was detected. Standard TMM normalization will be applied.'
        print(txt, flush=True)
        mode_tmm = 'single_species'
        file_orthogroup_table = ''
        file_genecount = ''
    else:
        txt = 'Multiple species were detected. ' \
              'Cross-species TMM normalization will be applied with single-copy orthologs.'
        print(txt, flush=True)
        mode_tmm = 'multi_species'
        check_ortholog_parameter_compatibility(args)
        if args.dir_busco is not None:
            file_orthogroup_table = os.path.join(dir_cstmm, 'cstmm_multispecies_busco_table.tsv')
            generate_multisp_busco_table(dir_busco=args.dir_busco, outfile=file_orthogroup_table)
        elif args.orthogroup_table is not None:
            file_orthogroup_table = os.path.realpath(args.orthogroup_table)
            if not os.path.exists(file_orthogroup_table):
                raise FileNotFoundError('Orthogroup table not found: {}'.format(file_orthogroup_table))
            if not os.path.isfile(file_orthogroup_table):
                raise IsADirectoryError(
                    'Orthogroup table path exists but is not a file: {}'.format(file_orthogroup_table)
                )
        file_genecount = os.path.join(dir_cstmm, 'cstmm_orthogroup_genecount.tsv')
        spp = filepath2spp(count_files)
        orthogroup2genecount(file_orthogroup=file_orthogroup_table, file_genecount=file_genecount, spp=spp)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    r_cstmm_path = os.path.join(dir_amalgkit_script, 'cstmm.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    r_command = ['Rscript', r_cstmm_path, dir_count, file_orthogroup_table, file_genecount, dir_cstmm, mode_tmm, r_util_path]
    print('')
    print('Starting R script: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)
    cleanup_tmp_amalgkit_files(work_dir='.')
