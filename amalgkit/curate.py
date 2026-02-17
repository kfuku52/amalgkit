import datetime
import os
import re
import shutil
import subprocess
import sys
import warnings

from amalgkit.util import *


def get_sample_group(args, metadata):
    if args.sample_group is None:
        sample_group = metadata.df.loc[:, 'sample_group'].dropna().unique()
    else:
        sample_group = re.findall(r"[\w]+", args.sample_group)
    if (len(sample_group)==0):
        txt = 'The "sample_group" column in --metadata ({}) is not filled. Please manually edit the file.\n'
        txt += '`amalgkit curate` recognizes samples with the same string in this columns to belong the same group.\n'
        txt += 'Exiting.\n'
        sys.stderr.write(txt.format(args.metadata))
        sys.exit(1)
    print('Tissues to be included: {}'.format(', '.join(sample_group)))
    sample_group = '|'.join(sample_group)
    return sample_group

def get_curate_completion_flag_path(curate_dir, sp):
    return os.path.join(curate_dir, sp, 'curate_completion_flag.txt')

def write_curate_completion_flag(curate_dir, sp):
    file_curate_completion_flag = get_curate_completion_flag_path(curate_dir, sp)
    with open(file_curate_completion_flag, 'w') as f:
        f.write('amalgkit curate completed at {}\n'.format(datetime.datetime.now()))

def register_curate_exit_status(curate_dir, sp, exit_status, failed_species):
    if exit_status == 0:
        write_curate_completion_flag(curate_dir=curate_dir, sp=sp)
    else:
        failed_species.add(sp)
        sys.stderr.write('amalgkit curate failed for species: {}\n'.format(sp))

def run_curate_r_script(args, metadata, sp, input_dir):
    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    correlation_threshold = args.correlation_threshold
    intermediate = args.plot_intermediate
    sample_group = get_sample_group(args, metadata)
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    r_script_path = os.path.join(dir_amalgkit_script, 'curate.r')
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    path_curate_input_metadata = os.path.join(input_dir, 'metadata.tsv')
    len_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_eff_length.tsv')
    if 'cstmm' in input_dir:
        count_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_cstmm_counts.tsv')
    else:
        count_file = os.path.join(os.path.abspath(input_dir), sp, sp + '_est_counts.tsv')
    if os.path.exists(count_file) and os.path.exists(len_file):
        print("Both counts and effective length files found: {}".format(sp), flush=True)
    else:
        if not os.path.exists(count_file):
            sys.stderr.write('Expected but undetected PATH of the count file: {}\n'.format(count_file))
        if not os.path.exists(len_file):
            sys.stderr.write('Expected but undetected PATH of the effective length file: {}\n'.format(len_file))
        sys.stderr.write('Skipping {}\n'.format(sp))
        print('Skipping {}'.format(sp), flush=True)
        return 1
    print("Starting Rscript to obtain curated {} values.".format(args.norm), flush=True)
    curate_r_result = subprocess.run([
            'Rscript',
            r_script_path,
            count_file,
            path_curate_input_metadata,
            os.path.realpath(args.out_dir),
            len_file,
            dist_method,
            str(mr_cut),
            '0',
            str(int(intermediate)),
            sample_group,
            args.sample_group_color,
            str(args.norm),
            str(int(args.one_outlier_per_iter)),
            str(correlation_threshold),
            str(args.batch_effect_alg),
            str(int(args.clip_negative)),
            str(int(args.maintain_zero)),
            os.path.realpath(r_util_path),
            str(int(args.skip_curation))
         ], check=False)
    return curate_r_result.returncode

def resolve_curate_input(args):
    if args.input_dir != 'inferred':
        print('Input_directory: {}'.format(args.input_dir))
        metadata = load_metadata(args, dir_subcommand=os.path.basename(args.input_dir))
        return metadata, args.input_dir

    dir_merge = os.path.realpath(os.path.join(args.out_dir, 'merge'))
    dir_cstmm = os.path.realpath(os.path.join(args.out_dir, 'cstmm'))
    if os.path.exists(dir_cstmm):
        print('Subdirectory for amalgkit cstmm will be used as input: {}'.format(dir_cstmm))
        metadata = load_metadata(args, dir_subcommand='cstmm')
        return metadata, dir_cstmm
    print('Subdirectory for amalgkit merge will be used as input: {}'.format(dir_merge))
    metadata = load_metadata(args, dir_subcommand='merge')
    return metadata, dir_merge


def list_selected_species(metadata):
    is_selected = (metadata.df['exclusion'] == 'no')
    selected_species = metadata.df.loc[is_selected, 'scientific_name'].drop_duplicates().values
    return [species.replace(' ', '_') for species in selected_species]


def collect_pending_species_for_curate(args, curate_dir, selected_species):
    pending_species = []
    for species in selected_species:
        flag_path = get_curate_completion_flag_path(curate_dir=curate_dir, sp=species)
        if not os.path.exists(flag_path):
            pending_species.append(species)
            continue
        if args.redo:
            print('Output file detected. Will be overwritten: {}'.format(species), flush=True)
            shutil.rmtree(os.path.join(curate_dir, species))
            pending_species.append(species)
            continue
        print('Skipping. Output file detected: {}'.format(species), flush=True)
    return pending_species


def run_curate_species_jobs(args, metadata, input_dir, curate_dir, pending_species, failed_species, species_jobs):
    if (species_jobs == 1) or (len(pending_species) <= 1):
        for species in pending_species:
            print('Starting: {}'.format(species), flush=True)
            exit_status = run_curate_r_script(args, metadata, species, input_dir)
            register_curate_exit_status(
                curate_dir=curate_dir,
                sp=species,
                exit_status=exit_status,
                failed_species=failed_species,
            )
        return

    max_workers = min(species_jobs, len(pending_species))
    print('Running curate for {:,} species with {:,} parallel jobs.'.format(len(pending_species), max_workers), flush=True)
    exit_status_by_species, failures = run_tasks_with_optional_threads(
        task_items=pending_species,
        task_fn=lambda species: run_curate_r_script(args, metadata, species, input_dir),
        max_workers=max_workers,
    )
    for species, exc in failures:
        failed_species.add(species)
        sys.stderr.write('amalgkit curate failed for species: {} ({})\n'.format(species, exc))
    for species, exit_status in exit_status_by_species.items():
        register_curate_exit_status(
            curate_dir=curate_dir,
            sp=species,
            exit_status=exit_status,
            failed_species=failed_species,
        )


def curate_main(args):
    species_jobs, _ = resolve_worker_allocation(
        requested_workers=getattr(args, 'species_jobs', 1),
        cpu_budget=getattr(args, 'cpu_budget', 0),
        worker_option_name='species_jobs',
        context='curate:',
    )
    check_rscript()
    metadata, input_dir = resolve_curate_input(args)
    if ('tpm' in args.norm) and ('cstmm' in input_dir):
            txt = ("TPM and TMM are incompatible. "
                   "If input data are CSTMM-normalized, "
                   "please switch --norm to any of the 'fpkm' normalization methods instead.")
            sys.stderr.write(txt)
    selected_species = list_selected_species(metadata)
    curate_dir = os.path.join(args.out_dir, 'curate')
    os.makedirs(curate_dir, exist_ok=True)
    failed_species = set()
    print(
        'Number of species in the selected metadata table ("exclusion"=="no"): {}'.format(len(selected_species)),
        flush=True,
    )
    pending_species = collect_pending_species_for_curate(args, curate_dir, selected_species)
    run_curate_species_jobs(
        args=args,
        metadata=metadata,
        input_dir=input_dir,
        curate_dir=curate_dir,
        pending_species=pending_species,
        failed_species=failed_species,
        species_jobs=species_jobs,
    )
    if len(failed_species) > 0:
        txt = 'amalgkit curate failed for {}/{} species: {}\n'
        sys.stderr.write(txt.format(len(failed_species), len(selected_species), ', '.join(sorted(failed_species))))
        sys.exit(1)
